#include "../include/maths.h"


uint64_t power(uint64_t base, uint64_t exp) {
    uint64_t res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = mul(res, base);
        base = mul(base, base);
        exp /= 2;
    }
    return res;
}

uint64_t modInverse(uint64_t n) {
    return power(n, MOD - 2);
}

// --- NTT Core (Parallelized) ---

// Precompute roots of unity [1, w, w^2, ...]
std::vector<uint64_t> precompute_roots(int n) {
    std::vector<uint64_t> roots(n);
    uint64_t g_n = power(G, (MOD - 1) / n);
    
    roots[0] = 1;
    for (int i = 1; i < n; i++) {
        roots[i] = mul(roots[i-1], g_n);
    }
    return roots;
}

// In-place Bit Reversal
void bit_reverse_permute(std::vector<uint64_t>& a) {
    int n = a.size();
    int log_n = 0;
    while ((1 << log_n) < n) log_n++;

    #pragma omp parallel for schedule(static)
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

// The Transform Function
void ntt(std::vector<uint64_t>& a, bool invert, const std::vector<uint64_t>& all_roots) {
    int n = a.size();
    
    bit_reverse_permute(a);

    // Iterate through stages (length = 2, 4, 8 ...)
    for (int len = 2; len <= n; len <<= 1) {
        int half_len = len / 2;
        int stride = n / len; 
        int num_blocks = n / len;

        // Adaptive Parallelism Strategy
        if (num_blocks >= omp_get_max_threads()) {
            // Strategy A: Parallelize Outer Blocks
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; i += len) {
                for (int j = 0; j < half_len; j++) {
                    int root_idx = j * stride;
                    if (invert && root_idx > 0) root_idx = n - root_idx;
                    
                    uint64_t w = all_roots[root_idx];
                    uint64_t u = a[i + j];
                    uint64_t v = mul(a[i + j + half_len], w);
                    
                    a[i + j] = add(u, v);
                    a[i + j + half_len] = sub(u, v);
                }
            }
        } else {
            // Strategy B: Parallelize Inner Offsets
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < half_len; j++) {
                int root_idx = j * stride;
                if (invert && root_idx > 0) root_idx = n - root_idx;
                uint64_t w = all_roots[root_idx];
                
                for (int i = 0; i < n; i += len) {
                    uint64_t u = a[i + j];
                    uint64_t v = mul(a[i + j + half_len], w);
                    
                    a[i + j] = add(u, v);
                    a[i + j + half_len] = sub(u, v);
                }
            }
        }
    }

    // Normalization for Inverse Transform
    if (invert) {
        uint64_t n_inv = modInverse(n);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) {
            a[i] = mul(a[i], n_inv);
        }
    }
}

// --- High Level Multiplication ---

std::string multiply_large_numbers(const std::string& num1, const std::string& num2) {
    // 1. Determine size (Power of 2)
    // Result can have at most len1 + len2 digits
    size_t n = 1;
    while (n < num1.size() + num2.size()) n <<= 1;

    // 2. Precompute Roots (Expensive, do once per size)
    // Note: In a real library, you'd cache these based on 'n'
    std::vector<uint64_t> roots = precompute_roots(n);

    // 3. Convert Strings to Polynomials (Integer Vectors)
    // We process input in reverse order so index 0 is the 1s place
    std::vector<uint64_t> a(n, 0), b(n, 0);
    
    // Parallel Parse is tricky due to string indexing, keeping it serial or simple parallel
    #pragma omp parallel for
    for (size_t i = 0; i < num1.size(); i++) {
        a[i] = num1[num1.size() - 1 - i] - '0';
    }
    
    #pragma omp parallel for
    for (size_t i = 0; i < num2.size(); i++) {
        b[i] = num2[num2.size() - 1 - i] - '0';
    }

    // 4. Perform Forward NTT
    // We can run these two in parallel using sections
    #pragma omp parallel sections
    {
        #pragma omp section
        ntt(a, false, roots);
        
        #pragma omp section
        ntt(b, false, roots);
    }

    // 5. Pointwise Multiplication (Convolution Theorem)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        a[i] = mul(a[i], b[i]);
    }

    // 6. Perform Inverse NTT
    ntt(a, true, roots);

    // 7. Carry Propagation (The "Schoolbook" cleanup)
    // This step is inherently sequential because carry depends on previous result
    std::vector<int> result;
    result.reserve(n);
    
    uint64_t carry = 0;
    for (int i = 0; i < n; i++) {
        uint64_t val = a[i] + carry;
        result.push_back(val % 10);
        carry = val / 10;
    }
    
    // Handle remaining carry
    while (carry) {
        result.push_back(carry % 10);
        carry /= 10;
    }

    // 8. Format Output
    // Remove trailing zeros (which are leading zeros in the number)
    while (result.size() > 1 && result.back() == 0) {
        result.pop_back();
    }
    
    // Convert back to string (reverse logic)
    std::string res_str;
    res_str.resize(result.size());
    
    #pragma omp parallel for
    for(size_t i=0; i < result.size(); i++) {
        res_str[i] = result[result.size() - 1 - i] + '0';
    }

    return res_str;
}

#if 0

// --- Driver Code ---
int main() {
    // Fast I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    // Example: Square a large number
    std::string s1 = "123456789123456789123456789";
    std::string s2 = "987654321987654321987654321";
    
    // Or generate random massive numbers for testing
    // s1 = std::string(100000, '9'); // 100k digits
    // s2 = std::string(100000, '9');

    std::cout << "Multiplying strings of length " << s1.length() << " and " << s2.length() << "..." << std::endl;
    
    double start = omp_get_wtime();
    std::string product = multiply_large_numbers(s1, s2);
    double end = omp_get_wtime();

    std::cout << "Result length: " << product.length() << std::endl;
    std::cout << "Time taken: " << (end - start) << " seconds." << std::endl;
    
    // Print first/last few digits to verify
    if (product.length() > 20) {
        std::cout << "Result (truncated): " << product.substr(0, 10) 
                  << " ... " << product.substr(product.length()-10) << std::endl;
    } else {
        std::cout << "Result: " << product << std::endl;
    }

    return 0;
}


#endif  

