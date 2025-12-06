#include "../include/fft.h"

// --- Precompute Roots of Unity ---
// Returns an array where roots[k] = G^k
std::vector<uint64_t> precompute_roots(int n) {
    std::vector<uint64_t> roots(n);
    // We need the N-th root of unity.
    // g_n = G^((MOD-1)/n)
    uint64_t g_n = power(G, (MOD - 1) / n);
    
    roots[0] = 1;
    for (int i = 1; i < n; i++) {
        roots[i] = (roots[i-1] * g_n) % MOD;
    }
    return roots;
}

// --- Bit Reversal ---
void bit_reverse_permute(std::vector<uint64_t>& a) {
    int n = a.size();
    int log_n = 0;
    while ((1 << log_n) < n) log_n++;

    // This is embarrassingly parallel
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int rev = 0;
        int temp = i;
        for (int j = 0; j < log_n; j++) {
            rev = (rev << 1) | (temp & 1);
            temp >>= 1;
        }
        if (i < rev) {
            std::swap(a[i], a[rev]);
        }
    }
}

// --- Parallel NTT ---
void ntt_openmp(std::vector<uint64_t>& a, bool invert, const std::vector<uint64_t>& all_roots) {
    int n = a.size();
    
    // 1. Bit Reversal
    bit_reverse_permute(a);

    // 2. The Butterfly Stages
    // len is the length of the current subarray being merged (2, 4, 8... N)
    for (int len = 2; len <= n; len <<= 1) {
        
        int half_len = len / 2;
        int stride = n / len; // Step size for accessing the precomputed roots

        // HEURISTIC: Switch loop strategies based on workload
        // If we have enough blocks to satisfy all threads, parallelize the blocks (i).
        // If blocks are few (late stages), parallelize the inner offsets (j).
        
        int num_blocks = n / len;
        
        // This threshold depends on your hardware threads. 
        // If blocks < threads, we are underutilizing the CPU if we parallelize 'i'.
        if (num_blocks >= omp_get_max_threads()) {
            
            // STRATEGY A: Many Blocks (Early Stages)
            // Parallelize the outer loop 'i'.
            // Each thread takes a full block (e.g., indices 0..len)
            
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; i += len) {
                // Inside a block, 'j' iterates sequentially
                for (int j = 0; j < half_len; j++) {
                    
                    // Direct lookup: O(1) access to w
                    // If inverting, we need w^-k = w^(N-k). 
                    // Or simpler: use roots[(N - k) % N].
                    int root_idx = j * stride;
                    if (invert && root_idx > 0) root_idx = n - root_idx;
                    
                    uint64_t w = all_roots[root_idx];
                    
                    uint64_t u = a[i + j];
                    uint64_t v = (a[i + j + half_len] * w) % MOD;
                    
                    a[i + j] = add(u, v);
                    a[i + j + half_len] = sub(u, v);
                }
            }
            
        } else {
            
            // STRATEGY B: Few Blocks / Large Blocks (Late Stages)
            // The 'i' loop is too small. We invert the loops.
            // Parallelize 'j' (indices *inside* the blocks).
            
            // Loop Order Inversion: 'j' is now the outer, parallel loop
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < half_len; j++) {
                
                // Precompute w for this 'j' once
                int root_idx = j * stride;
                if (invert && root_idx > 0) root_idx = n - root_idx;
                uint64_t w = all_roots[root_idx];
                
                // Apply this 'j' offset to every block 'i'
                // Since num_blocks is small, this inner loop is short.
                for (int i = 0; i < n; i += len) {
                    uint64_t u = a[i + j];
                    uint64_t v = (a[i + j + half_len] * w) % MOD;
                    
                    a[i + j] = add(u, v);
                    a[i + j + half_len] = sub(u, v);
                }
            }
        }
    }

    // 3. Normalization (if Inverse)
    if (invert) {
        uint64_t n_inv = modInverse(n);
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            a[i] = (a[i] * n_inv) % MOD;
        }
    }
}

int main() {
    // Determine size (Power of 2)
    int n = 1 << 20; // 1 Million elements
    
    // Precompute roots once (Sequential or Parallel)
    std::cout << "Precomputing roots..." << std::endl;
    std::vector<uint64_t> roots = precompute_roots(n);
    
    // Create data
    std::vector<uint64_t> data(n);
    for(int i=0; i<n; i++) data[i] = i;

    // Run NTT
    std::cout << "Running Parallel NTT on " << n << " elements..." << std::endl;
    double start = omp_get_wtime();
    
    ntt_openmp(data, false, roots);
    
    double end = omp_get_wtime();
    std::cout << "Done in " << (end - start) << " seconds." << std::endl;

    return 0;
}
