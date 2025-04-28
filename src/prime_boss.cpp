#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cassert>
#include <omp.h>
#include <zlib.h>
#include <sstream>
#include <fstream>
#include "../include/BigInt.hpp"  // Assuming your BigInt class is here

const int SEGMENT_SIZE = 10'000'000;

struct CompressedBuffer {
    std::vector<unsigned char> data;
};

// Simple sieve for small integers (always up to sqrt(n) which fits in long long)
std::vector<int> simple_sieve(long long limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (long long i = 2; i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (long long j = i * i; j <= limit; j += i)
                is_prime[j] = false;
        }
    }
    std::vector<int> primes;
    for (long long i = 2; i <= limit; ++i)
        if (is_prime[i])
            primes.push_back(static_cast<int>(i));
    return primes;
}

// Delta-of-delta encode primes in [low, high)
template<typename T>
void sieve_segment(const T& low, const T& high, const std::vector<int>& base_primes,
                   std::vector<short>& encoded_primes, T& last_prime, T& last_delta) {
    T size = high - low;
    std::vector<bool> is_prime(static_cast<size_t>(size.to_long()), true);

    for (int p : base_primes) {
        T prime = static_cast<T>(p);
        if (prime * prime >= high) break;
        T start = prime * prime;
        if (start < low)
            start = ((low + prime - 1) / prime) * prime;
        for (T j = start; j < high; j = j + prime) {
            is_prime[static_cast<size_t>((j - low).to_long())] = false;
        }
    }

    for (size_t i = 0; i < static_cast<size_t>(size.to_long()); ++i) {
        if (is_prime[i]) {
            T prime = low + static_cast<T>(i);
            if (prime > 1) {
                if (last_prime == -1) {
                    encoded_primes.push_back(static_cast<short>(prime.to_long()));  // Assume prime fits in 16 bits here
                } else if (last_delta == -1) {
                    T delta = prime - last_prime;
                    encoded_primes.push_back(static_cast<short>(delta.to_long()));
                    last_delta = delta;
                } else {
                    T delta = prime - last_prime;
                    T delta_of_delta = delta - last_delta;
                    assert(delta_of_delta.to_long() >= -32768 && delta_of_delta.to_long() <= 32767);
                    encoded_primes.push_back(static_cast<short>(delta_of_delta.to_long()));
                    last_delta = delta;
                }
                last_prime = prime;
            }
        }
    }
}

// Compress a vector of shorts into a compressed memory buffer
CompressedBuffer compress_buffer(const std::vector<short>& deltas) {
    uLong source_size = deltas.size() * sizeof(short);
    uLong dest_size = compressBound(source_size);

    CompressedBuffer compressed;
    compressed.data.resize(dest_size);

    int res = compress2(compressed.data.data(), &dest_size,
                        reinterpret_cast<const Bytef*>(deltas.data()), source_size, Z_BEST_COMPRESSION);
    if (res != Z_OK) {
        std::cerr << "Compression failed!\n";
        exit(1);
    }

    compressed.data.resize(dest_size);
    return compressed;
}

// Lock-free multithreaded sieve and compression
template<typename T>
void multithreaded_sieve_and_compress_lockfree(const T& n, const std::string& output_filename, int num_threads = 4) {
    long long limit = static_cast<long long>(sqrt(n.to_long())) + 1;  // Base primes fit in long long
    auto base_primes = simple_sieve(limit);

    T two = 2;
    T segment_size = SEGMENT_SIZE;
    T total_segments = (n - two) / segment_size + 1;
    
    std::vector<CompressedBuffer> compressed_segments(static_cast<size_t>(total_segments.to_long()));

    omp_set_num_threads(num_threads);

    #pragma omp parallel for schedule(dynamic)
    for (long long segment_idx = 0; segment_idx < total_segments.to_long(); ++segment_idx) {
        T low = two + segment_size * segment_idx;
        T high = low + segment_size;
        if (high > n + 1) {
            high = n + 1;
        }

        std::vector<short> encoded;
        T last_prime = -1;
        T last_delta = -1;

        sieve_segment<T>(low, high, base_primes, encoded, last_prime, last_delta);

        compressed_segments[segment_idx] = compress_buffer(encoded);
    }

    // Merge all compressed segments sequentially
    std::ofstream out(output_filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return;
    }

    for (auto& segment : compressed_segments) {
        out.write(reinterpret_cast<const char*>(segment.data.data()), segment.data.size());
    }

    out.close();
}


// Read entire file into buffer
std::vector<unsigned char> read_entire_file(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for reading.");
    }
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<unsigned char> buffer(size);
    file.read(reinterpret_cast<char*>(buffer.data()), size);
    return buffer;
}

// Decompress buffer
std::vector<unsigned char> decompress_buffer(const std::vector<unsigned char>& compressed_data) {
    uLongf uncompressed_size = compressed_data.size() * 20;  // rough guess

    std::vector<unsigned char> uncompressed(uncompressed_size);

    int res = uncompress(uncompressed.data(), &uncompressed_size,
                         compressed_data.data(), compressed_data.size());

    if (res == Z_BUF_ERROR) {
        uncompressed_size *= 2;
        uncompressed.resize(uncompressed_size);
        res = uncompress(uncompressed.data(), &uncompressed_size,
                         compressed_data.data(), compressed_data.size());
    }

    if (res != Z_OK) {
        throw std::runtime_error("Failed to decompress buffer.");
    }

    uncompressed.resize(uncompressed_size);
    return uncompressed;
}


// Reconstruct primes from delta-of-delta encoded data
template<typename T>
std::vector<T> reconstruct_primes(const std::vector<unsigned char>& uncompressed_data) {
    const short* deltas = reinterpret_cast<const short*>(uncompressed_data.data());
    size_t num_deltas = uncompressed_data.size() / sizeof(short);

    std::vector<T> primes;
    primes.reserve(num_deltas);

    T running_sum(0);
    T last_delta(0);

    for (size_t i = 0; i < num_deltas; ++i) {
        short delta = deltas[i];
        if (i == 0) {
            running_sum = T(delta);
            primes.push_back(running_sum);
        } else if (i == 1) {
            running_sum = running_sum + T(delta);
            last_delta = T(delta);
            primes.push_back(running_sum);
        } else {
            last_delta = last_delta + T(delta);
            running_sum = running_sum + last_delta;
            primes.push_back(running_sum);
        }
    }

    return primes;
}

// Master function: load primes from compressed gz file
template<typename T>
std::vector<T> load_primes_from_compressed_file(const std::string& filename) {
    auto compressed_data = read_entire_file(filename);
    auto uncompressed_data = decompress_buffer(compressed_data);
    auto primes = reconstruct_primes<T>(uncompressed_data);
    return primes;
}




int main() {
    BigInt n("1000000");  
    std::string output_file = "primes_bigint_lockfree_compressed.gz";
    int num_threads = 8;

    std::cout << "Generating primes up to " << n.to_long() << " with " << num_threads << " threads (BigInt lock-free)...\n";
    multithreaded_sieve_and_compress_lockfree<BigInt>(n, output_file, num_threads);
    std::cout << "Done. Output saved to " << output_file << "\n";

    try {
        auto primes = load_primes_from_compressed_file<BigInt>(output_file);

        std::cout << "Loaded " << primes.size() << " primes.\n";
        for (size_t i = 0; i < std::min<size_t>(primes.size(), 20); ++i) {
            std::cout << primes[i] << " ";
        }
        std::cout << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
    }

    return 0;
}
