#include "../include/BigInt.hpp"

#include <stdio.h>

std::vector<uint8_t> BigInt::get_numerus() {
  std::vector<uint8_t> v = numerus;
  return v;
}

void BigInt::setNumerus(const std::vector<uint8_t>& source) {
    numerus = source; // or numerus.assign(source.begin(), source.end());
}

BigInt BigInt::slice(int i, int j) const {
  BigInt z;
  z.set_sign(this->get_sign());

  int d = (j - 1) + 1;
  int ds = this->size() - 1;

  if (d > ds) {
    return BigInt("NaN");
  }

  if (i > j) {

    return BigInt("NaN");
  }

  if ((i < 0) || (j < 0)) {

    return BigInt("NaN");
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

  
  int n = x.size();
  int m = y.size();

  if (x.abs() < 10 && y.abs() < 10) {

    int result = (int)x[0] * (int)y[0];

    return BigInt(result);
  }

  int k = std::max(n, m);

  int k2 = std::floor(k / 2);

  divmod10 dx = x.divmod(k2);

  divmod10 dy = y.divmod(k2);
  BigInt x_high = dx.quotient;
  BigInt x_low = dx.remainder;
  BigInt y_high = dy.quotient;
  BigInt y_low = dy.remainder;

  BigInt z0 = karatsuba(x_low, y_low);
  BigInt c1 = x_low + x_high;
  BigInt c2 = y_low + y_high;

  BigInt z1 = karatsuba(c1, c2);

  BigInt z2 = karatsuba(x_high, y_high);
  BigInt z3 = z1 - z2 - z0;

  BigInt result = z2.lshift(2 * k2) + z3.lshift(k2) + z0;

  return result;
}

#if 0


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



BigInt BigInt::Schonhage_Strassen(BigInt &x, BigInt &y) const { return x; }

BigInt BigInt::Toom3(BigInt &x, BigInt &y) const { return x; }

#endif




BigInt::BigInt() {
  sign = _NULL;
  numerus = std::vector<uint8_t>(0);
  
}

BigInt::BigInt(const std::vector<uint8_t> &num) {
  numerus = num;
  sign = POS;
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
  } else {
    z.sign = POS;
  }

  if (x == 0) {
    z.insert(0, 0);
  }

  while (x > 0) {
    z.insert(x % 10, 0);
    x /= 10;
  }

  *this = z;
}

bool BigInt::is_digit(char ch) const {

  return (int(ch) >= int('0') && int(ch) <= int('9'));

}

BigInt::BigInt(const std::string c) {
  sign=POS;
  char ch;
  int n = c.size();
  if (n==0 || c == "NaN"){
    sign = UNDEFINED;
    return;
  }

  ch=c[0];  

  if (n==1){    
    if ( ! is_digit(ch)) {
      sign = UNDEFINED;
      return;
    }
    else {
      uint8_t x = int(ch) - int('0');
      numerus.push_back(x);
      return;
    }
  }

  if (int(ch) != 45 && int(ch) != 43 && !is_digit(ch)) {
    sign = UNDEFINED;
    
    return;
  } 
  if (int(ch) ==45){
      sign = NEG;
    }
  if (int(ch) == 43){
      sign = POS;
    }
  if (is_digit(ch) && int(ch) != 45 && int(ch) != 43) {
     uint8_t x = int(ch) - int('0');
    numerus.push_back(x);
  }
  
    for (int i = 1; i < n; i++) {     
      ch = c[i];                   
      if (is_digit(ch)) {        
        uint8_t x = int(ch) - int('0');        
        numerus.push_back(x);
      }
      else {
        sign = UNDEFINED;
        return;
      }  
  }

}

#if 0
BigInt::BigInt(const std::string c) {
  sign = POS;
  char ch;

  int n = c.size();
  if (n==0 || c == "Nan"){
    sign = UNDEFINED;
    return;
  }
  if (n==1){
    ch=c[0];  
    if ( ! is_digit_char(ch)) {
      sign = UNDEFINED;
      return;
    }
    
  }
    for (int i = 0; i < n; i++) {
      ch = c[i];

      if (is_digit_char(ch)) {
        if (i == 0 && int(ch) == 45) {
          sign = NEG;
        } else if (i==0 && int(ch) == 43) {
          sign = POS;
        }

        else {
          sign = UNDEFINED;
          return;
        }

      } else {
        uint8_t x = int(ch) - int('0');
        numerus.push_back(x);
      }
  }
}

#endif
divmod10 BigInt::divmod(const long n) const {
  BigInt z = *this;

  divmod10 result;
  int m = z.size();
  if (m <= n) {

    return divmod10{BigInt("0"), BigInt(z)};
  }

  int d = m - n;

  std::vector<uint8_t> quotient(d);
  std::vector<uint8_t> remainder(n);
  for (int i = 0; i <= d; i++) {
    quotient[i] = z[i];
  }

  for (int i = d; i < m; i++) {
    int j = i - d;
    remainder[j] = z[i];
  }

  if (remainder.size() > 1 && remainder[0] == 0) {
    remainder.erase(remainder.begin());
  }
  BigInt q(quotient);
  BigInt r(remainder);

  result.quotient = q;
  result.remainder = r;

  return result;
}


// add_to_front = true
BigInt BigInt::lshift(const int n) const {
  BigInt z = *this;
  if (z == 0) {
    return z;
  }
  for (int i = 0; i < n; i++) {
    uint8_t b = 0;
    z.numerus.push_back(b);
  }
  return z;
}

std::unique_ptr<std::vector<uint8_t>> BigInt::numerus_ptr() {
  return std::make_unique<std::vector<uint8_t>>(numerus);
}

uint8_t BigInt::operator[](const int i) const { return numerus[i]; }

void BigInt::set_sign(SIGN x) { sign = x; }
#if 0
void BigInt::operator!() {
  if (sign == NEG) {
    sign = POS;
  }
  else if (sign == POS) {
    sign = NEG;
  }
  else if (sign == UNDEFINED) {
    sign = UNDEFINED;
  }
 
}
#endif

BigInt BigInt::abs() const {
  BigInt z = *this;
  if (z.get_sign() == NEG) {
    z.set_sign(POS);
  }
  return z;
}

SIGN BigInt::get_sign() const { return sign; }

size_t BigInt::size() const { return numerus.size(); }

std::ostream &operator<<(std::ostream &out, const BigInt &num) {

  if (num.sign == UNDEFINED) {
    out << "NaN";
    return out;
  }

   if (num.sign == _NULL) {
    out << "_NULL";
    return out;
  }

  if (num.sign == NEG) {
    out << "-";
  }
  size_t n = num.size();
  
  int i = 0;

  if (num[i] == 0) {

    while (num[i] == 0 && i < n) {
      ++i;
    }
    while (i < n) {
      out << static_cast<int>(num[i]);
      i++;
    }
  } else {
    for (auto x : num.numerus) {
      out << (int)(x);
    }
  }

  return out;
}

void BigInt::insert(const uint8_t &val, const int &ix) {
  numerus.insert(numerus.begin() + ix, val);
}

#if 0
void BigInt::insert(const int &val, const int &ix) {

  uint8_t x = (uint8_t)val;
  numerus.insert(numerus.begin() + ix, x);
}
#endif

BigInt BigInt::vadd(BigInt &x, BigInt &y) const {
  BigInt z;

  int n = x.size();
  int m = y.size();
  int k = std::max(n, m);
  if (n != m) {
    if (n > m) {
      int d = n - m;
      y = y.shift_n(d, true);
    } else {
      int d = m - n;
      x = x.shift_n(d, true);
    }
  }

  std::vector<uint8_t> x_numerus = x.get_numerus();
  std::vector<uint8_t> y_numerus = y.get_numerus();

  int c = 0;
  int tot = 0;

  for (int i = k - 1; i >= 0; i--) {
    tot = x_numerus[i] + y_numerus[i] + c;

    if (tot >= 10) {
      c = 1;

      z.insert(tot % 10, 0);
      if (k == 1) {
        z.insert(c, 0);
      }
      if (i == 0) {
        z.insert(c, 0);
      }
    } else {
      c = 0;
      z.insert(tot, 0);
    }
  }

  int i = 0;
  while (z[i] == 0) {
    i++;
  }

  if (i > 0) {

    std::unique_ptr<std::vector<uint8_t>> ptr = z.numerus_ptr();
    ptr->erase(ptr->begin(), ptr->begin() + i);
  }

  return z;
}

void BigInt::print_numerus() const {
  for (auto x : this->numerus) {
    std::cout << (int)(x);
  }
  std::cout << "\n";
}
BigInt BigInt::vsub(BigInt &x, BigInt &y) const {

  int n = x.size();
  int m = y.size();
  int k = std::max(n, m);

  std::vector<uint8_t> result(k);
  std::fill(result.begin(), result.end(), uint8_t(0));
  std::vector<uint8_t> x_numerus = x.get_numerus();
  std::vector<uint8_t> y_numerus = y.get_numerus();

  int carry = 0;
  for (int i = k - 1; i >= 0; i--) {

    int diff = (int)x_numerus[i] - (int)y_numerus[i] - carry;

    if (diff < 0) {
      carry = 1;
      int j = i - 1;
      while ((int)x_numerus[j] == 0) {
        x_numerus[j] = 10 - carry;
        j--;
      }
      x_numerus[j] -= carry;
      diff += 10;
      carry = 0;
    } else {
      carry = 0;
    }

    result[i] = diff;
  }

  int i = 0;
  // Remove any leading zeros
  while (result[i] == 0) {
    i++;
  }
  if (i > 0) {
    result.erase(result.begin(), result.begin() + i);
  }

  BigInt z(result);
  return z;
}

BigInt BigInt::vmult(BigInt &x, BigInt &y) const {

  int n = x.size();
  int m = y.size();

  std::vector<uint8_t> x_numerus = x.get_numerus();
  std::vector<uint8_t> y_numerus = y.get_numerus();

  if (x_numerus[0] == 0 || y_numerus[0] == 0) {

    return BigInt(0);
  }

  if (n > 10 || m > 10) {
    
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
      t = x_numerus[j] * y_numerus[i] + carry;
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

    BigInt z1 = z.shift_n(shift);
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

BigInt BigInt::operator*(const BigInt &num) {
  BigInt x = *this;
  BigInt y = num;
  BigInt z;
  if (this->sign == UNDEFINED || num.sign == UNDEFINED || this->sign == _NULL || num.sign == _NULL) {
    return BigInt("NaN");
  }
  
  if (x.size() == 0 || y.size() == 0) {
    return BigInt("NaN");
  }
  if (x == 0 || y == 0) {
    return BigInt("0");
  }

  if (y.size() > x.size()) {
    z = vmult(y, x);
  } else {

    z = vmult(x, y);
  }

  SIGN xsign = x.get_sign();
  SIGN ysign = y.get_sign();

  if ((xsign == NEG and ysign == POS) or (xsign == NEG and ysign == POS)) {

    z.set_sign(NEG);
  }

  if (xsign == NEG and ysign == NEG) {
    z.set_sign(POS);
  }

  if (xsign == POS and ysign == POS) {
    z.set_sign(POS);
  }

  return z;
}

bool BigInt::operator!=(const BigInt &num) const {
  if (*this == num)
    return false;

  return true;
}

void BigInt::operator++() {
  BigInt z = *this;
  BigInt one("1");
  z = z + one;
  *this = z;
}

void BigInt::operator--() {
  BigInt z = *this;
  BigInt one("1");
  z = z - one;
  *this = z;
}
// void BigInt::operator!(){
//   BigInt z = *this;
//   z = -z;
// }

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
    x = x.shift_n(n - m, true);
  } else if (nx > ny) {
    y = y.shift_n(n - m, true);
  }

  if (x < y && x.sign == POS && y.sign == POS) {
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
    z = vsub(x, y);
    z.sign = NEG;

    return z;
  }

  if (x > y && x.sign == POS && y.sign == POS) {
    z = vsub(x, y);
    z.sign = POS;

    return z;
  }

  if (x > y && x.sign == POS && y.sign == NEG) {

    z = vadd(x, y);
    z.sign = POS;

    return z;
  }

  if (x > y && x.sign == NEG && y.sign == POS) {

    // cannot happen
  }

  if (x > y && x.sign == NEG && y.sign == NEG) {
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

BigInt BigInt::shift_n(int m, bool add_to_front) const {
  BigInt z = *this;
  for (int i = 0; i < m; i++) {
    uint8_t b = 0;
    if (add_to_front) {

      z.numerus.insert(z.numerus.begin(), b);
    } else {
      z.numerus.push_back(b);
    }
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
    BigInt zx = x.shift_n(n - m, true);
    x = zx;
  } else if (nx > ny) {
    BigInt zy = y.shift_n(n - m, true);
    y = zy;
  }

  SIGN xsign = x.get_sign();
  SIGN ysign = y.get_sign();

  if (xsign == NEG && ysign == POS) {
    x.sign = POS;
    return y - x;
  }
  if (xsign == POS && ysign == NEG) {
    y.sign = POS;
    return x - y;
  }
  if (xsign == NEG && ysign == NEG) {
    z.sign = NEG;
  }
  if (xsign == POS && ysign == POS) {
    z.sign = POS;
  }

  int c = 0;
  int tot = 0;

  for (int i = n - 1; i >= 0; i--) {
    int xi = x[i];
    int yi = y[i];
    tot = xi + yi + c;
    if (tot >= 10) {
      c = 1;
      z.insert(tot % 10, 0);
      if (n == 1 || i == 0) {
        z.insert(1, 0);
      }
    } else {
      c = 0;
      z.insert(tot, 0);
    }
  }
  return z;
}

#if 0
 void BigInt::operator!() {
  BigInt z = *this;
  if (z.sign == POS) {
    z.sign = NEG;
  } else {
    z.sign = POS;
  }
  
}
#endif
bool BigInt::operator==(const BigInt &y) const {
  BigInt temp = y;
  SIGN xsign = get_sign();
  SIGN ysign = temp.get_sign();
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

  SIGN xsign = this->get_sign();
  SIGN ysign = y.get_sign();

  int n = this->size();
  int m = y.size();

  if (xsign == POS and ysign == NEG) {
    return false;
  }
  if (xsign == NEG and ysign == POS) {
    return true;
  }
  if (n > m) {

    if (xsign == POS and ysign == POS) {
      return false;
    } else if (xsign == NEG and ysign == NEG)
      return true;
  }

  if (n < m) {

    if (xsign == NEG and ysign == NEG) {
      return false;
    } else if (xsign == POS and ysign == POS) {
      return true;
    }
  }

  if (n == m) {

    for (int i = 0; i < n; i++) {
      size_t numerus_i = numerus[i];
      size_t temp_i = y[i];
      if ((numerus_i > temp_i) && (xsign == NEG)) {
        return true;
      } else if ((numerus_i > temp_i) && (xsign == POS))
        return false;
      else if ((temp_i > numerus_i) && (xsign == NEG))
        return false;
      else if ((temp_i > numerus_i) && (xsign == POS)) {

        return true;
      }
    } // end for
  }

  return false;
}

bool BigInt::operator<=(const BigInt &y) const {
    if (*this == y) {
        return true;
    }

    SIGN xsign = get_sign();
    SIGN ysign = y.get_sign();

    if (xsign != ysign) {
        if (xsign == POS) {
            return false;
        } else {
            return true;
        }
      
    }

    int n = size();
    int m = y.size();

    if (n != m) {
        return (n > m) == (xsign == NEG);
      
    }

    for (int i = 0; i < n; ++i) {
        if (numerus[i] != y.numerus[i]) {
            return (numerus[i] > y.numerus[i]) == (xsign == NEG);
        }
    }
    std::cout <<"here" << std::endl;
    return false;
}

bool BigInt::operator>(const BigInt &y) const {
    if (*this == y) {
        return false;
    }

    SIGN xsign = get_sign();
    SIGN ysign = y.get_sign();

    if (xsign != ysign) {
        if (xsign == POS) {
            return true;
        } else {
            return false;
        }        
    }

    int n = size();
    int m = y.size();

    if (n != m) {
        return (n > m) == (xsign == POS);
      
    }

    for (int i = 0; i < n; ++i) {
        if (numerus[i] != y.numerus[i]) {
            return (numerus[i] > y.numerus[i]) == (xsign == POS);
        }
    }

    return false;

}

bool BigInt::operator>=(const BigInt &y) const {
    if (*this == y) {
        return true;
    }

    SIGN xsign = get_sign();
    SIGN ysign = y.get_sign();

    if (xsign != ysign) {
        if (xsign == POS) {
            return true;
        } else {
            return false;
        }
        
    }

    int n = size();
    int m = y.size();

    if (n != m) {
        return (n > m) == (xsign == POS);
      
    }

    for (int i = 0; i < n; ++i) {
        if (numerus[i] != y.numerus[i]) {
            return (numerus[i] > y.numerus[i]) == (xsign == POS);
        }
    }

    return false;
}

