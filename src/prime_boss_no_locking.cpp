#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>         // memset
#include <cassert>
#include <omp.h>           // OpenMP
#include <zlib.h>          // zlib compression
#include <sstream>         // ostringstream
#include <fstream>

const int SEGMENT_SIZE = 10'000'000;

struct CompressedBuffer {
    std::vector<unsigned char> data;
};

// Simple sieve to find primes up to sqrt(n)
std::vector<int> simple_sieve(int limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (int j = i * i; j <= limit; j += i)
                is_prime[j] = false;
        }
    }
    std::vector<int> primes;
    for (int i = 2; i <= limit; ++i)
        if (is_prime[i])
            primes.push_back(i);
    return primes;
}

// Delta-of-delta encode primes in [low, high)
void sieve_segment(int low, int high, const std::vector<int>& base_primes,
                   std::vector<short>& encoded_primes, int& last_prime, int& last_delta) {
    int size = high - low;
    std::vector<bool> is_prime(size, true);

    for (int p : base_primes) {
        if (p * p >= high) break;
        int start = std::max(p * p, ((low + p - 1) / p) * p);
        for (int j = start; j < high; j += p)
            is_prime[j - low] = false;
    }

    for (int i = 0; i < size; ++i) {
        if (is_prime[i]) {
            int prime = low + i;
            if (prime > 1) {
                if (last_prime == -1) {
                    encoded_primes.push_back(static_cast<short>(prime));  // first prime
                } else if (last_delta == -1) {
                    int delta = prime - last_prime;
                    encoded_primes.push_back(static_cast<short>(delta));  // first delta
                    last_delta = delta;
                } else {
                    int delta = prime - last_prime;
                    int delta_of_delta = delta - last_delta;
                    assert(delta_of_delta >= -32768 && delta_of_delta <= 32767);
                    encoded_primes.push_back(static_cast<short>(delta_of_delta));
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

    int res = compress2(compressed.data.data(), &dest_size, reinterpret_cast<const Bytef*>(deltas.data()), source_size, Z_BEST_COMPRESSION);
    if (res != Z_OK) {
        std::cerr << "Compression failed!\n";
        exit(1);
    }

    compressed.data.resize(dest_size);
    return compressed;
}

// Lock-free multithreaded sieve + compression
void multithreaded_sieve_and_compress_lockfree(long long n, const std::string& output_filename, int num_threads = 4) {
    int limit = static_cast<int>(std::sqrt(n)) + 1;
    auto base_primes = simple_sieve(limit);

    int num_segments = (n - 1) / SEGMENT_SIZE + 1;
    std::vector<CompressedBuffer> compressed_segments(num_segments);

    omp_set_num_threads(num_threads);

    #pragma omp parallel for schedule(dynamic)
    for (int segment_idx = 0; segment_idx < num_segments; ++segment_idx) {
        long long low = 2LL + segment_idx * SEGMENT_SIZE;
        long long high = std::min(low + SEGMENT_SIZE, n + 1);

        std::vector<short> encoded;
        int last_prime = -1;
        int last_delta = -1;

        sieve_segment(static_cast<int>(low), static_cast<int>(high), base_primes, encoded, last_prime, last_delta);

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


// Load the entire compressed file into a memory buffer
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

// Decompress zlib compressed buffer
std::vector<unsigned char> decompress_buffer(const std::vector<unsigned char>& compressed_data) {
    uLongf uncompressed_size = compressed_data.size() * 20; // Rough estimate (primes are sparse)

    std::vector<unsigned char> uncompressed(uncompressed_size);

    int res = uncompress(uncompressed.data(), &uncompressed_size,
                         compressed_data.data(), compressed_data.size());

    if (res == Z_BUF_ERROR) {
        // Buffer too small; try a bigger buffer
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
std::vector<int> reconstruct_primes(const std::vector<unsigned char>& uncompressed_data) {
    const short* deltas = reinterpret_cast<const short*>(uncompressed_data.data());
    size_t num_deltas = uncompressed_data.size() / sizeof(short);

    std::vector<int> primes;
    primes.reserve(num_deltas);

    int running_sum = 0;
    int last_delta = 0;

    for (size_t i = 0; i < num_deltas; ++i) {
        short delta = deltas[i];
        if (i == 0) {
            running_sum = delta;
            primes.push_back(running_sum);
        } else if (i == 1) {
            running_sum += delta;
            last_delta = delta;
            primes.push_back(running_sum);
        } else {
            last_delta += delta;
            running_sum += last_delta;
            primes.push_back(running_sum);
        }
    }

    return primes;
}

// Master function: load and reconstruct primes from file
std::vector<int> load_primes_from_compressed_file(const std::string& filename) {
    auto compressed_data = read_entire_file(filename);
    auto uncompressed_data = decompress_buffer(compressed_data);
    auto primes = reconstruct_primes(uncompressed_data);
    return primes;
}


int main() {
    long long n = 1e6; // Try 1e9, 1e10, etc.
    std::string output_file = "primes_lockfree_compressed.gz";
    int num_threads = 4;

    std::cout << "Generating primes up to " << n << " using " << num_threads << " threads (lock-free)...\n";
    multithreaded_sieve_and_compress_lockfree(n, output_file, num_threads);
    std::cout << "Done. Output saved to " << output_file << "\n";

    try {
        auto primes = load_primes_from_compressed_file(output_file);
        std::cout << "Loaded " << primes.size() << " primes.\n";
        // Optional: print first few primes
        for (size_t i = 0; i < std::min<size_t>(primes.size(), 2000); ++i) {
            std::cout << primes[i] << " ";
        }
        std::cout << "\n";
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
    }

    return 0;
}
