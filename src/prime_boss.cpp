#include <iostream>
#include <vector>
#include <cmath>
#include <zlib.h>          // For gzip writing
#include <omp.h>           // For OpenMP multithreading
#include <cstring>         // For memset
#include <fstream>
#include <cassert>

const int SEGMENT_SIZE = 10'000'000;

// Write a vector of signed short deltas to gzip file
void write_deltas_gzip(const std::vector<short>& deltas, gzFile& out) {
    gzwrite(out, deltas.data(), deltas.size() * sizeof(short));
}

// Simple sieve of Eratosthenes up to limit
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

// Sieve one segment and encode primes into delta-of-delta format
void sieve_segment(int low, int high, const std::vector<int>& base_primes, 
                   std::vector<short>& output, int& last_prime, int& last_delta) {
    int size = high - low;
    std::vector<bool> is_prime(size, true);

    for (int p : base_primes) {
        if (p * p >= high) break;
        int start = std::max(p * p, ((low + p - 1) / p) * p);
        for (int j = start; j < high; j += p) {
            is_prime[j - low] = false;
        }
    }

    for (int i = 0; i < size; ++i) {
        if (is_prime[i]) {
            int prime = low + i;
            if (prime > 1) {
                if (last_prime == -1) {
                    output.push_back(static_cast<short>(prime));  // first prime
                } else if (last_delta == -1) {
                    int delta = prime - last_prime;
                    output.push_back(static_cast<short>(delta));  // first delta
                    last_delta = delta;
                } else {
                    int delta = prime - last_prime;
                    int delta_of_delta = delta - last_delta;
                    assert(delta_of_delta >= -32768 && delta_of_delta <= 32767);
                    output.push_back(static_cast<short>(delta_of_delta));
                    last_delta = delta;
                }
                last_prime = prime;
            }
        }
    }
}

// Full multithreaded segmented sieve and compression
void multithreaded_sieve_and_compress(long long n, const std::string& output_filename, int num_threads = 4) {
    int limit = static_cast<int>(std::sqrt(n)) + 1;
    auto base_primes = simple_sieve(limit);

    gzFile out = gzopen(output_filename.c_str(), "wb");
    if (!out) {
        std::cerr << "Failed to open output file.\n";
        return;
    }

    omp_set_num_threads(num_threads);

    //int last_prime = -1;
    //int last_delta = -1;

    // Synchronize across threads writing into gzip
    #pragma omp parallel
    {
        std::vector<short> local_buffer;
        #pragma omp for schedule(dynamic)
        for (long long low = 2; low <= n; low += SEGMENT_SIZE) {
            long long high = std::min(low + SEGMENT_SIZE, n + 1);
            std::vector<short> segment_output;
            int segment_last_prime = -1;
            int segment_last_delta = -1;

            sieve_segment(static_cast<int>(low), static_cast<int>(high), base_primes, segment_output, segment_last_prime, segment_last_delta);

            #pragma omp critical
            {
                write_deltas_gzip(segment_output, out);
            }
        }
    }

    gzclose(out);
}

int main() {
    long long n = 1000000;  // You can try 1e9, 1e12, etc.
    std::string output_file = "primes_multithreaded_delta_of_delta.gz";
    int num_threads = 4;  // Adjust for your CPU cores

    std::cout << "Generating primes up to " << n << " using " << num_threads << " threads...\n";
    multithreaded_sieve_and_compress(n, output_file, num_threads);
    std::cout << "Done. Output saved to " << output_file << "\n";

    return 0;
}
