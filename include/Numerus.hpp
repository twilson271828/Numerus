#include <bitset>
#include <complex>
#include <vector>   
#include <iostream>


namespace Numerus {
std::bitset<4> convertToBinary(uint8_t &n); 

size_t convertToDecimal(std::bitset<4> const &b);

std::complex<double> exponentiate(size_t k, size_t n, size_t N); 

std::vector<std::complex<double>> n_roots_of_unity(int N);

std::complex<double> dft_coef(std::vector<std::complex<double>> &input,size_t n);

std::complex<double> dift(std::vector<std::complex<double>> &input,size_t n); 

std::complex<double> dft(std::vector<std::complex<double>> &input,size_t n);

} // namespace Numerus
