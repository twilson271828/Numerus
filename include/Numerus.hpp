#include <bitset>
#include <complex>
#include <vector>   
#include <iostream>
#include "BigInt.hpp"


std::bitset<4> convertToBinary(uint8_t &n); 

size_t convertToDecimal(std::bitset<4> const &b);

std::complex<double> exponentiate(size_t k, size_t n, size_t N); 

std::complex<double> dift(std::vector<std::complex<double>> &input,size_t n); 

std::complex<double> dft(std::vector<std::complex<double>> &input,size_t n);

  