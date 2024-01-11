#include "../include/Numerus.hpp"

std::bitset<4> Numerus::convertToBinary(uint8_t &n) {
  std::bitset<4> b;
  int i = 0;

  try {
    while (n > 0) {
      int r = n % 2;
      n /= 2;
      b.set(i, r);
      i += 1;
    }
  } catch (std::out_of_range &e) {
    std::cout << "std::out_of_range exception caught: " << e.what() << "\n";
  }

  return b;
}

size_t Numerus::convertToDecimal(std::bitset<4> const &b) {

  size_t n = (size_t)b.to_ulong();
  return n;
}


std::complex<double> Numerus::exponentiate(size_t k, size_t n, size_t N) {

  double x = 0.0;
  double y = 0.0;
  double theta = ((2 * M_PI * k * n) / N);

  x = std::cos(theta);
  y = std::sin(theta);
  std::complex<double> omega(x, y);
  return omega;
}

std::complex<double> Numerus::dift(std::vector<std::complex<double>> &input,
                                  size_t n) {

  size_t N = input.size();
  std::complex<double> coeff(0.0, 0.0);

  for (int k = 0; k < N; k++) {
    coeff += input[k] * exponentiate(k, n, N);
  }

  if (abs(coeff.real()) < 1e-6) {
    coeff.real(0.0);
  }

  if (abs(coeff.imag()) < 1e-6) {
    coeff.imag(0.0);
  }

  coeff /= N;
  return coeff;
}

std::complex<double> Numerus::dft(std::vector<std::complex<double>> &input,
                                 size_t n) {

  size_t N = input.size();
  std::complex<double> coeff(0.0, 0.0);

  for (int k = 0; k < N; k++) {
    coeff += input[k] / exponentiate(k, n, N);
  }

  if (abs(coeff.real()) < 1e-6) {
    coeff.real(0.0);
  }

  if (abs(coeff.imag()) < 1e-6) {
    coeff.imag(0.0);
  }

  return coeff;
}
std::complex<double> Numerus::dft_coef(std::vector<std::complex<double>> &input,
                              size_t n) {

  size_t N = input.size();
  std::complex<double> coeff(0.0, 0.0);

  for (int k = 0; k < N; k++) {
    coeff += input[k] / exponentiate(k, n, N);
  }

  if (abs(coeff.real()) < 1e-6) {
    coeff.real(0.0);
  }

  if (abs(coeff.imag()) < 1e-6) {
    coeff.imag(0.0);
  }

  return coeff;
}



std::vector<std::complex<double>> Numerus::n_roots_of_unity(int N) {

  std::vector<std::complex<double>> nroots;

  for (int k = 0; k < N; k++) {
    std::complex<double> root = exponentiate(k, 1, N);
    nroots.push_back(root);
  }
  return nroots;
}