#include <bitset>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <typeinfo>
#include <vector>

enum SIGN { POS, NEG, UNDEFINED, _NULL };

struct split;
struct divmod10;

class BigInt {

private:
  std::vector<uint8_t> numerus;

  SIGN sign;

  BigInt vsub(BigInt &x, BigInt &y) const;
  BigInt vadd(BigInt &x, BigInt &y) const;
  BigInt vmult(BigInt &x, BigInt &y) const;
  // bool is_digit_char(char c) const;

  // BigInt Schonhage_Strassen(BigInt &x, BigInt &y) const;
  // BigInt Toom3(BigInt &x, BigInt &y) const;
  // std::complex<double> exponentiate(size_t k, size_t n, size_t N);
  // std::complex<double> dift(std::vector<std::complex<double>> &input, size_t
  // n); std::complex<double> dft(std::vector<std::complex<double>> &input,
  // size_t n); size_t bitrev(size_t n);

public:
  BigInt();

  /// @brief
  /// @param c
  BigInt(const std::string c);

  BigInt(const BigInt &num);

  BigInt(const long &num);

  BigInt(const std::vector<int> &v, SIGN s = POS);
  BigInt(const std::vector<BigInt> &v, SIGN s = POS);
  void setNumerus(const std::vector<uint8_t> &source);
  BigInt(const std::vector<uint8_t> &num, SIGN s = POS);
  BigInt karatsuba(BigInt &x, BigInt &y) const;
  BigInt lshift(const int n) const;
  BigInt shift_n(const int n, bool add_to_front = false) const;
  divmod10 divmod(const long n) const;
  BigInt slice(int i, int j) const;
  uint8_t operator[](const int i) const;
  size_t size() const;
  void print_numerus() const;
  BigInt abs() const;
  SIGN get_sign() const;
  void set_sign(SIGN x);
  std::unique_ptr<std::vector<uint8_t>> numerus_ptr();
  std::vector<uint8_t> get_numerus();
  bool is_digit(char c) const;
  long to_long() const;
  std::vector<int> to_list() const;

  /// @brief
  void insert(const uint8_t &val, const int &ix);

  BigInt operator*(const BigInt &num);
  BigInt operator/(const BigInt &num) const;

  BigInt operator+(const BigInt &num) const;

  BigInt operator-(const BigInt &num) const;

  bool operator<(const BigInt &num) const;

  bool operator<=(const BigInt &num) const;

  bool operator>=(const BigInt &num) const;

  bool operator>(const BigInt &num) const;

  bool operator!=(const BigInt &num) const;

  bool operator==(const BigInt &num) const;

  void operator++();

  void operator--();

  split split_it(size_t m) const;

  friend std::ostream &operator<<(std::ostream &out, const BigInt &num);
};

struct split {
  BigInt xleft;
  BigInt xright;
  size_t m;
};

struct divmod10 {
  BigInt quotient;
  BigInt remainder;
};

inline std::ostream &operator<<(std::ostream &os, const divmod10 &d) {
  os << "(" << d.quotient << "," << d.remainder << ")";

  return os;
}
