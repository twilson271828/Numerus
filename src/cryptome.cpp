#include "../include/BigInt.hpp"
#include <algorithm>
#include <bitset>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

template <typename T> void printVector(std::vector<T> &x) {
  for (auto &i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}



void fft(std::vector<cd>& a, bool invert) {
  size_t n = a.size();
  if (n == 1) return;

  std::vector<cd> a0(n / 2), a1(n / 2);
  for (size_t i = 0; 2 * i < n; ++i) {
      a0[i] = a[i*2];
      a1[i] = a[i*2+1];
  }

  fft(a0, invert);
  fft(a1, invert);

  double ang = 2 * PI / n * (invert ? -1 : 1);
  cd w(1), wn(cos(ang), sin(ang));
  for (size_t i = 0; 2 * i < n; ++i) {
      a[i] = a0[i] + w * a1[i];
      a[i + n/2] = a0[i] - w * a1[i];
      if (invert) {
          a[i] /= 2;
          a[i + n/2] /= 2;
      }
      w *= wn;
  }
}


std::vector<int> multiply_fft(const std::vector<int>& a, const std::vector<int>& b) {
  std::vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
  size_t n = 1;
  while (n < a.size() + b.size()) n <<= 1;
  fa.resize(n); fb.resize(n);

  fft(fa, false);
  fft(fb, false);
  for (size_t i = 0; i < n; ++i)
      fa[i] *= fb[i];
  fft(fa, true);

  std::vector<int> result(n);
  long long carry = 0;
  for (size_t i = 0; i < n; ++i) {
      long long val = static_cast<long long>(std::round(fa[i].real())) + carry;
      result[i] = val % 10;
      carry = val / 10;
  }

  while (carry) {
      result.push_back(carry % 10);
      carry /= 10;
  }

  // Remove leading zeros
  while (result.size() > 1 && result.back() == 0)
      result.pop_back();

  return result;
}


BigInt Schonhage_Strassen(BigInt& x, BigInt& y) {
  
  if (x == 0 || y == 0) {
    return BigInt(0);
  }
  if (x.get_sign() == NEG && y.get_sign() == NEG) {
    x.set_sign(POS);
    y.set_sign(POS);
  } else if (x.get_sign()== NEG || y.get_sign() == NEG) {
    x.set_sign(POS);
    y.set_sign(POS);
  }
  if (x.get_sign() == UNDEFINED || y.get_sign() == UNDEFINED) {
    return BigInt(0);
  }
  if (x.get_sign() == _NULL || y.get_sign() == _NULL) {
    return BigInt(0);
  }
  if (x == 0 || y == 0) {
    return BigInt(0);
  }
  std::vector<uint8_t> x_numerus = x.getNumerus();
  std::vector<uint8_t> y_numerus = y.getNumerus();
  std::vector<int> a(x_numerus.begin(), x_numerus.end());
  std::vector<int> b(y_numerus.begin(), y_numerus.end());

  auto result = multiply_fft(a, b);

  BigInt product;
  product.setNumerus(std::vector<uint8_t>(result.begin(), result.end()));

  if (x.get_sign() == y.get_sign()) {
    product.set_sign(POS);
  } else {
    product.set_sign(NEG);
  }
  

  return product.trim_zeros();
}

std::vector<BigInt> split_number(BigInt x, int m) {

  std::vector<uint8_t> numerus = x.getNumerus();
  std::vector<BigInt> result;

  for (int i = numerus.size(); i >= 0; i -= m) {
    BigInt temp = x.slice(i - m, i - 1);
    // result.push_back(temp);
    result.insert(result.begin(), temp);
  }
  return result;
}

divmod10 div(BigInt &x, BigInt &y) {

  divmod10 d;
  long ylong = y.to_long();
  int m = 4;
  std::vector<BigInt> x_parts = split_number(x, m);
  printVector(x_parts);
  std::vector<BigInt> y_parts = split_number(y, m);
  BigInt remainder = 0;
  long quotient_part;
  std::vector<BigInt> quotient;
  std::vector<int> quotient_list;

  for (auto &part : x_parts) {
    std::cout << "remainder = " << remainder << "\n";
    std::cout << "part = " << part << "\n";

    remainder = remainder * std::pow(10, m) + part;
    // std::cout << "result = " << remainder << "\n";
    quotient_part = remainder.to_long() / ylong;
    quotient_list = BigInt(quotient_part).to_list();
    int nzeros = m - quotient_list.size();

    if (nzeros > 0) {
      for (int i = 0; i < nzeros; i++) {
        quotient_list.insert(quotient_list.begin(), 0);
      }
    }

    std::cout << "quotient_list = ";
    printVector(quotient_list);

    remainder = remainder - BigInt(quotient_part) * y;

    quotient.push_back(quotient_list);

    std::cout << "************************************************\n";
  }

  printVector(quotient);
  d.quotient = BigInt(quotient);
  d.remainder = remainder;

  return d;
}

void printVector(std::vector<int> &x) {
  for (auto &i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}
bool isLeadingZeroPresent(const std::string &c) {
  if (c.empty()) {
    return false;
  }
  if (c[0] == '0') {
    return true;
  }
  return false;
}

int main() {

BigInt a("1234567890");
BigInt b("987654321098");
BigInt result = Schonhage_Strassen(a, b);
std::cout << result << std::endl;
std::cout << "1219326311247340343220" << std::endl;


  return 0;
}
