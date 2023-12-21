#include "../include/BigInt.hpp"

//std::vector<uint8_t> BigInt::get_numerus() { return this->numerus; }


std::bitset<4> convertToBinary(const uint8_t &n) {
  std::bitset<4> b;  
  int i = 0;
  size_t n1 = n;
  try {
    while (n > 0) {
      int r = n % 2;
      n1 /= 2;
      b.set(i, r);
      i += 1;
    }
  } catch (std::out_of_range &e) {
    std::cout << "std::out_of_range exception caught: " << e.what() << "\n";
  }

  return b;
}


size_t convertToDecimal(std::bitset<4> const &b) {

  size_t n = (size_t)b.to_ulong();
  return n;
}


#if 0

std::vector<std::bitset<4>> BigInt::get_binary_numerus() {

  std::vector<std::bitset<4>> binary_numerus;

  for (uint8_t n : this->numerus) {
    std::bitset<4> b = this->convertToBinary(n);
    binary_numerus.push_back(b);
  }

  return binary_numerus;
}

size_t BigInt::bitrev(size_t n) {

  size_t rev = 0;

  // traversing bits of 'n' from the right
  while (n > 0) {
    // bitwise left shift
    // 'rev' by 1
    rev <<= 1;

    // if current bit is '1'
    if (n & 1 == 1)
      rev ^= 1;

    // bitwise right shift
    // 'n' by 1
    n >>= 1;
  }

  // required number
  return rev;
}

std::complex<double> BigInt::exponentiate(size_t k, size_t n, size_t N) {

  double x = 0.0;
  double y = 0.0;
  double theta = ((2 * M_PI * k * n) / N);

  x = std::cos(theta);
  y = std::sin(theta);
  std::complex<double> omega(x, y);
  return omega;
}

std::complex<double> BigInt::dft(std::vector<std::complex<double>> &input,
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

std::complex<double> BigInt::dift(std::vector<std::complex<double>> &input,
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

BigInt BigInt::vsub(BigInt &x, BigInt &y) const {

  BigInt z;

  int n = x.size();
  int m = y.size();

  if (n != m) {
    if (n > m) {
      int d = n - m;
      y = y.m10(d, true);
    } else {
      int d = m - n;
      x = x.m10(d, true);
    }
  }

  for (int i = n - 1; i >= 0; i--) {
    int tot = x[i] - y[i];
    if (tot < 0) {
      int val = 0;
      if (i > 0) {
        val = x[i - 1] - 1;
        x.insert(val, i - 1);
      }
      // x[i-1]=x[i-1]-1;
      z.insert(10 + tot, 0);
    } else if (tot > 0) {
      z.insert(tot, 0);
    } else {
      if (i != 0) {
        z.insert(0, 0);
      }
    }
  }
  return z;
}

BigInt BigInt::Schonhage_Strassen(BigInt &x, BigInt &y) const { return x; }

BigInt BigInt::slice(int i, int j) const {
  BigInt z;

  if ((j - i) + i > this->size() - 1) {
    std::cout << "The slice range [i,j] = [" << i << "," << j
              << "] is greater than the length of " << *this << "\n";
    return BigInt("NAN");
  }

  if (i > j) {
    std::cout << "[i,j] = "
              << "[" << i << " , " << j << "]\n";
    std::cout << " The starting index for BigInt::slice must be less than the "
                 "ending index.\n";
    return BigInt("NAN");
  }

  if ((i < 0) || (j < 0)) {
    std::cout << "The starting and ending indices for the BigInt::slice "
                 "routine must be greater than or equal to 0\n";
    return BigInt("NAN");
  }

  // Starting and Ending iterators
  auto start = this->numerus.begin() + i;
  auto end = this->numerus.begin() + i + (j - i) + 1;

  // To store the sliced vector
  // std::vector<uint8_t> result(j - i + 1);
  std::vector<uint8_t> result(j - i + 1);
  z.numerus = result;

  std::copy(start, end, z.numerus.begin());

  return z;
}

split BigInt::split_it(size_t m) const {
  size_t n = this->size();

  BigInt r;
  BigInt c;

  split z;

  z.xright = this->slice(n - m, n - 1);
  z.xleft = this->slice(0, n - m - 1);
  z.m = m;
  return z;
}

BigInt BigInt::karatsuba(BigInt &x, BigInt &y) const {

  size_t n = x.size();
  size_t m = y.size();

  if (n > m) {

    y = y.m10(n - m, true);
  }
  if (n < m) {
    x = x.m10(m - n, true);
  }

  if (n < 2 || m < 2) {
    return x * y;
  }

  size_t k = std::max(n, m);
  size_t k2 = std::floor(k / 2);

  split split_x = x.split_it(k2);
  split split_y = y.split_it(k2);

  BigInt low_x = split_x.xright;
  BigInt low_y = split_y.xright;
  BigInt high_x = split_x.xleft;
  BigInt high_y = split_y.xleft;

  BigInt z2 = karatsuba(high_x, high_y);
  BigInt z0 = karatsuba(low_x, low_y);
  BigInt w1 = high_x + low_x;
  BigInt w2 = high_y + low_y;

  BigInt z1 = karatsuba(w1, w2);

  BigInt W = z1 - z2 - z0;

  BigInt P = z2.m10(k2 * 2, false) + W.m10(k2, false) + z0;

  return P;
}

BigInt BigInt::Toom3(BigInt &x, BigInt &y) const { return x; }

BigInt BigInt::vmult(BigInt &x, BigInt &y) const {

  int n = x.size();
  int m = y.size();

  if (n > 50 && m > 50) {
    return karatsuba(x, y);
  }

  int shift = 0;
  int carry = 0;
  int base = 10;
  int t = 0;
  std::vector<BigInt> vecs;

  for (int i = m - 1; i >= 0; i--) {
    BigInt z;
    carry = 0;
    for (int j = n - 1; j >= 0; j--) {
      t = x[j] * y[i] + carry;
      carry = 0;
      if (t >= 10) {
        auto dv = std::div(t, base);
        carry = (int)dv.quot;
        if (j == 0) {
          z.insert((int)dv.rem, 0);
          z.insert(carry, 0);
        } else {
          z.insert((int)dv.rem, 0);
        }
      } else {
        z.insert(t, 0);
      }
    }

    BigInt z1 = z.m10(shift);
    shift += 1;
    std::vector<BigInt>::iterator ix = vecs.begin();
    vecs.insert(ix, z1);
  }

  BigInt a;
  for (int i = 0; i < vecs.size(); i++) {
    a = vadd(a, vecs[i]);
  }
  return a;
}

BigInt BigInt::vadd(BigInt &x, BigInt &y) const {
  BigInt z;

  int n = x.size();
  int m = y.size();
  int k = std::max(n, m);
  if (n != m) {
    if (n > m) {
      int d = n - m;
      y = y.m10(d, true);
    } else {
      int d = m - n;
      x = x.m10(d, true);
    }
  }

  int c = 0;
  int tot = 0;
  for (int i = k - 1; i >= 0; i--) {
    tot = x[i] + y[i] + c;
    if (tot >= 10) {
      c = 1;
      z.insert(tot % 10, 0);
      if (k == 1) {
        z.insert(1, 0);
      }
      if (i == 0) {
        z.insert(c, 0);
      }
    } else {
      c = 0;
      z.insert(tot, 0);
    }
  }

  return z;
}
#endif 
BigInt::BigInt() {
  sign = POS;
  return;
}

BigInt::BigInt(const std::vector<std::bitset<4>> &num) {
  numerus =num;
  sign = POS;
  return;
}

BigInt::BigInt(const BigInt &num) {
  numerus = num.numerus;
  sign = num.sign;
}

BigInt::BigInt(const long &num) {
  BigInt z;

  long x = num;
  if (x < 0) {
    x *= -1;
    z.sign = NEG;
  }

  while (x > 0) {
    z.insert(x % 2, 0);
    x /= 2;
  }

  *this = z;
  
}

BigInt::BigInt(const std::string c) {
  sign = POS;
  if (c == "NAN") {
    sign = UNDEFINED;
    return;
  }
  try {

    int n = c.size();
    for (int i = 0; i < n; i++) {

      char ch = c[i];

      if (int(ch) < int('0') || int(ch) > int('9')) {
        if (i == 0 && int(ch) == 45) {
          sign = NEG;
        } else if (int(ch) == 43) {
          sign = POS;
        }

        else {
          throw std::runtime_error("Encountered a non-numeric char character "
                                   "when attempting to construct a BigInt");
        }
      } else {

        uint8_t x = int(ch) - int('0');
        
        numerus.push_back(convertToBinary(x));
      }
      
    }
  } catch (std::exception &e) {
    std::cout << "Caught Exception:" << e.what() << "\n";
    std::exit(0);
  }
}
#if 0
BigInt BigInt::m10(int m, bool add_to_front) const {
  BigInt z = *this;
  for (int i = 0; i < m; i++) {
    if (add_to_front) {
      z.numerus.insert(z.numerus.begin(), 0);
     
    } else {
      z.numerus.push_back(0);
    }
  }
  return z;
}


#endif 

size_t BigInt::operator[](const int i) const { 
  
  std::bitset<4> b = numerus[i];
  return convertToDecimal(b);
   }

 
int BigInt::get_sign() const {

  if (sign == NEG)
    return -1;
  else if (sign == POS)
    return 1;
  else
    return 0;
}


size_t BigInt::size() const { return numerus.size(); }
  


std::ostream &operator<<(std::ostream &out, const BigInt &num) {

  size_t n = num.size();
  if (num.sign == NEG) {
    out << "-";
  }

  for (auto x : num.numerus) {
      out << x;
    }

  int i = 0;
  
  if (num[i] == 0) {
  
    while (num[i] == 0 && i < n) {
      ++i;
    }
    while (i < n) {
      out << convertToDecimal(num[i]);
      i++;
    }
  
  } else {
    for (auto x : num.numerus) {
      out << convertToDecimal(x);
    }
  }   
  out << "\n";

  return out;

}


#if 0

void BigInt::insert(const int &val, const int &ix) {
  
  numerus.insert(numerus.begin() + ix, (uint8_t)val);
}
BigInt BigInt::operator/(const BigInt &num) const { return num; }

BigInt BigInt::operator-(const BigInt &num) const {
  BigInt x = *this;
  BigInt y = num;
  BigInt z;
  if (x == y) {
    return BigInt("0");
  }
  int nx = x.size();
  int ny = y.size();

  int n = std::max(nx, ny);
  int m = std::min(nx, ny);

  if (nx < ny) {
    BigInt zx = x.m10(n - m, true);
    x = zx;
  } else if (nx > ny) {
    BigInt zy = y.m10(n - m, true);
    y = zy;
  }

  if (x < y && x.sign == POS && y.sign == POS) {
    // tested
    // swap x and y

    BigInt temp = x;
    x = y;
    y = temp;
    z = vsub(x, y);
    z.sign = NEG;
    return z;
  }

  if (x < y && x.sign == NEG && y.sign == POS) {
    // cannot happen
  }

  if (x < y && x.sign == POS && y.sign == NEG) {
    // cannot happen
  }

  if (x < y && x.sign == NEG && y.sign == NEG) {
    // tested
    z = vsub(x, y);
    z.sign = NEG;
    return z;
  }

  if (x > y && x.sign == POS && y.sign == POS) {
    // tested
    z = vsub(x, y);
    z.sign = POS;
    return z;
  }

  if (x > y && x.sign == POS && y.sign == NEG) {
    // tested
    z = vadd(x, y);
    z.sign = POS;
    return z;
  }

  if (x > y && x.sign == NEG && y.sign == POS) {
    // cannot happen
  }

  if (x > y && x.sign == NEG && y.sign == NEG) {
    // tested
    // swap x and y
    BigInt temp = x;
    x = y;
    y = temp;

    z = vsub(x, y);
    z.sign = POS;
    return z;
  }

  return z;
}

BigInt BigInt::operator+(const BigInt &num) const {
  BigInt z;
  BigInt x = *this;
  BigInt y = num;
  int nx = x.size();
  int ny = y.size();
  int n = std::max(nx, ny);
  int m = std::min(nx, ny);
  if (nx < ny) {
    BigInt zx = x.m10(n - m, true);
    x = zx;
  } else if (nx > ny) {
    BigInt zy = y.m10(n - m, true);
    y = zy;
  }

  int xsign = x.get_sign();
  int ysign = y.get_sign();

  if (xsign < 0 && ysign > 0) {
    x.sign = POS;
    return y - x;
  }
  if (xsign > 0 && ysign < 0) {
    y.sign = POS;
    return x - y;
  }
  if (xsign < 0 && ysign < 0) {
    z.sign = NEG;
  }
  if (xsign > 0 && ysign > 0) {
    z.sign = POS;
  }

  int c = 0;
  int tot = 0;
  for (int i = n - 1; i >= 0; i--) {
    tot = x[i] + y[i] + c;
    if (tot >= 10) {
      c = 1;
      z.insert(tot % 10, 0);
      if (n == 1) {
        z.insert(1, 0);
      }
    } else {
      c = 0;
      z.insert(tot, 0);
    }
  }
  return z;
}

BigInt BigInt::operator*(const BigInt &num) {
  BigInt x = *this;
  BigInt y = num;
  BigInt z;

  if (y.size() > x.size()) {
    z = vmult(y, x);
  } else {
    z = vmult(x, y);
  }

  return z;
}

BigInt BigInt::operator!() const {
  BigInt z = *this;
  if (z.sign == POS) {
    z.sign = NEG;
  } else {
    z.sign = POS;
  }
  return z;
}
#endif




#if 0
bool BigInt::operator==(const BigInt &y) const {
  BigInt temp = y;
  int xsign = get_sign();
  int ysign = temp.get_sign();
  if (xsign != ysign)
    return false;

  int m = size();
  int n = temp.size();
  if (m != n)
    return false;
  for (int i = 0; i < n; i++) {
    if (numerus[i] != temp[i])
      return false;
  }

  return true;
}

bool BigInt::operator<(const BigInt &y) const {

  // -1 == True
  // 1 == False

  BigInt temp = y;
  int xsign = get_sign();
  int ysign = temp.get_sign();

  int n = size();
  int m = temp.size();

  if (xsign < 0 and ysign > 0)
    return true;
  else if (xsign > 0 and ysign < 0)
    return false;

  if (n > m) {
    if (xsign == -1 and ysign == -1)
      return true;
    else if (xsign == 1 and ysign == 1)
      return false;
  }

  if (n < m) {
    if (xsign == -1 and ysign == -1)
      return false;
    else if (xsign == 1 and ysign == 1) {

      return false;
    }
  }

  if (n == m) {
    for (int i = 0; i < n; i++) {
      if ((numerus[i] > temp[i]) && (xsign == -1)) {

        return true;
      } else if ((numerus[i] > temp[i]) && (xsign == 1))
        return false;
      else if ((temp[i] > numerus[i]) && (xsign == -1))
        return false;
      else if ((temp[i] > numerus[i]) && (xsign == 1)) {

        return true;
      }
    } // end for
  }

  return false;
}

bool BigInt::operator>(const BigInt &y) const {
  // 1 == True
  // -1 == False
  BigInt temp = y;

  int xsign = get_sign();
  int ysign = temp.get_sign();

  int n = size();
  int m = temp.size();

  if (xsign < 0 and ysign > 0)
    return false;
  else if (xsign > 0 and ysign < 0)
    return true;

  if (n > m) {
    if (xsign == -1 and ysign == -1)
      return false;
    else if (xsign == 1 and ysign == 1)
      return true;
  }

  if (n < m) {
    if (xsign == -1 and ysign == -1)
      return true;
    else if (xsign == 1 and ysign == 1)
      return false;
  }

  if (n == m) {
    for (int i = 0; i < n; i++) {
      if ((numerus[i] > temp[i]) && (xsign == -1))
        return false;
      else if ((numerus[i] > temp[i]) && (xsign == 1))
        return true;
      else if ((temp[i] > numerus[i]) && (xsign == -1))
        return true;
      else if ((temp[i] > numerus[i]) && (xsign == 1))
        return false;
    }
  }

  return false;
}
#endif

