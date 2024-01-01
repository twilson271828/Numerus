#include <bitset>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
//#include <armadillo>

enum SIGN { POS, NEG, UNDEFINED };

struct split;

class BigInt {

private:
  
  std::vector<std::bitset<4> > numerus;

  SIGN sign;
 
  BigInt vsub(BigInt &x, BigInt &y) const;
  BigInt vadd(BigInt &x, BigInt &y) const;
  BigInt vmult(BigInt &x, BigInt &y) const;
  //BigInt karatsuba(BigInt &x, BigInt &y) const;
  //BigInt Schonhage_Strassen(BigInt &x, BigInt &y) const;
  //BigInt Toom3(BigInt &x, BigInt &y) const;
  //std::complex<double> exponentiate(size_t k, size_t n, size_t N);
  //std::complex<double> dift(std::vector<std::complex<double>> &input, size_t n);
  //std::complex<double> dft(std::vector<std::complex<double>> &input, size_t n);
  //size_t bitrev(size_t n);

public:
  BigInt();

  /// @brief
  /// @param c
  BigInt(const std::string c);

  BigInt(const BigInt &num);

  BigInt(const long &num);

  BigInt(const std::vector<std::bitset<4>> &num);

  void numerus_ix(const int &ix, const ulong &val);
  
  BigInt m16(const int m, bool add_to_front = false) const;
  BigInt slice(int i, int j) const;
  std::bitset<4> operator[](const int i) const;
  size_t size() const;

  int get_sign() const;
  void set_numerus(const std::vector<std::bitset<4>> &num);
  std::vector<int> get_numerus();
  

  /// @brief
  //void negative();
  void insert(const std::bitset<4> &val,const int &ix);
  void insert(const int &val, const int &ix);

  BigInt operator*(const BigInt &num);
  BigInt operator/(const BigInt &num) const;

  BigInt operator+(const BigInt &num) const;

  BigInt operator-(const BigInt &num) const;

  bool operator<(const BigInt &num) const;

  //bool operator<=(const BigInt &num) const;

  //bool operator>=(const BigInt &num) const;

  bool operator>(const BigInt &num) const;

  //bool operator!=(const BigInt &num) const;

  bool operator==(const BigInt &num) const;

  //BigInt operator++();

  //BigInt operator--();

  //BigInt operator!() const;

  split split_it(size_t m) const;

  friend std::ostream &operator<<(std::ostream &out, const BigInt &num);
};

std::bitset<4> convertToBinary(uint8_t &n);
size_t convertToDecimal( std::bitset<4> const &n);

struct bitset_add {
  std::bitset<4> sum;
  bool carry;
};

struct split {
  BigInt xleft;
  BigInt xright;
  size_t m;
};
#if 0
inline void printSplit(const split &p) {
  std::cout << "high_p = " << p.xleft;
  std::cout << "low_p = " << p.xright;
  std::cout << "m = " << p.m << "\n";
}
#endif  