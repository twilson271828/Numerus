#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>
#include <cstdint>
#include <iomanip>

// --- Configuration ---
// The "Golden Prime" for NTT
const uint64_t MOD = 998244353;
const uint64_t G = 3;

// --- Modular Arithmetic Helpers ---
inline uint64_t add(uint64_t a, uint64_t b) {
    return (a + b >= MOD) ? (a + b - MOD) : (a + b);
}

inline uint64_t sub(uint64_t a, uint64_t b) {
    return (a >= b) ? (a - b) : (a - b + MOD);
}

inline uint64_t mul(uint64_t a, uint64_t b) {
    return (a * b) % MOD;
}


uint64_t power(uint64_t base, uint64_t exp);
uint64_t modInverse(uint64_t n);
std::vector<uint64_t> precompute_roots(int n);
void bit_reverse_permute(std::vector<uint64_t>& a);
void ntt_openmp(std::vector<uint64_t>& a, bool invert, const std::vector<uint64_t>& all_roots);
void ntt(std::vector<uint64_t>& a, bool invert, const std::vector<uint64_t>& all_roots);
std::string multiply_large_numbers(const std::string& num1, const std::string& num2)