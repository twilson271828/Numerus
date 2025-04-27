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
  divmod10 burnikel_ziegler(BigInt const &x, const BigInt &y) const;
  std::vector<BigInt> split_number(const BigInt x, const int m) const;
  BigInt vsub(BigInt &x, BigInt &y) const;
  BigInt vadd(BigInt &x, BigInt &y) const;
  BigInt vmult(BigInt &x, BigInt &y) const;

  // BigInt Schonhage_Strassen(BigInt &x, BigInt &y) const;
  // BigInt Toom3(BigInt &x, BigInt &y) const;
  // std::complex<double> exponentiate(size_t k, size_t n, size_t N);
  // std::complex<double> dift(std::vector<std::complex<double>> &input, size_t
  // n); std::complex<double> dft(std::vector<std::complex<double>> &input,
  // size_t n); size_t bitrev(size_t n);

public:
  BigInt();
  BigInt trim_zeros() const;
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
  BigInt sqrt_bigint(const BigInt& n);
  BigInt lshift(const int n) const;
  BigInt shift_n(const int n, bool add_to_front = false) const;
  divmod10 divmod(const long n) const;
  BigInt slice(int i, int j) const;
  int operator[](const int i) const;
  int size() const;
  void print_numerus() const;
  BigInt abs() const;
  SIGN get_sign() const;
  void set_sign(SIGN x);
  std::unique_ptr<std::vector<uint8_t>> numerus_ptr();
  std::vector<uint8_t> getNumerus() const;
  bool is_digit(char c) const;
  long to_long() const;
  std::vector<int> to_list() const;

  /// @brief
  void insert(const uint8_t &val, const int &ix);

  BigInt operator*(const BigInt &num);
  BigInt operator/(const long n) const;
  

  divmod10 div(const BigInt &num) const;

  BigInt operator+(const BigInt &num) const;

  BigInt operator-(const BigInt &num) const;

  bool operator<(const BigInt &num) const;

  bool operator<=(const BigInt &num) const;

  bool operator>=(const BigInt &num) const;

  bool operator>(const BigInt &num) const;

  bool operator!=(const BigInt &num) const;

  BigInt& operator=(const BigInt& other);
  
  bool operator==(const BigInt &num) const;

  void operator++();

  void operator--();

  split split_it(int m) const;

  friend std::ostream &operator<<(std::ostream &out, const BigInt &num);
};

BigInt operator/(const BigInt& a, const BigInt& b);

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
