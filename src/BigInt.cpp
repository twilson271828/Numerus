#include "../include/BigInt.hpp"
#include "../include/maths.h"
#include <stdio.h>



BigInt BigInt::Schonhage_Strassen(const std::string& num1, const std::string& num2) const {
 // 1. Determine size (Power of 2)
    // Result can have at most len1 + len2 digits
    size_t n = 1;
    while (n < num1.size() + num2.size()) n <<= 1;

    // 2. Precompute Roots (Expensive, do once per size)
    // Note: In a real library, you'd cache these based on 'n'
      std::vector<uint64_t> roots = precompute_roots(n);

    // 3. Convert Strings to Polynomials (Integer Vectors)
    // We process input in reverse order so index 0 is the 1s place
      std::vector<uint64_t> a(n, 0), b(n, 0);
    
    // Parallel Parse is tricky due to string indexing, keeping it serial or simple parallel
    #pragma omp parallel for
    for (size_t i = 0; i < num1.size(); i++) {
        a[i] = num1[num1.size() - 1 - i] - '0';
    }
    
    #pragma omp parallel for
    for (size_t i = 0; i < num2.size(); i++) {
        b[i] = num2[num2.size() - 1 - i] - '0';
    }

    // 4. Perform Forward NTT
    // We can run these two in parallel using sections
    #pragma omp parallel sections
    {
        #pragma omp section
        ntt(a, false, roots);
        
        #pragma omp section
        ntt(b, false, roots);
    }

    // 5. Pointwise Multiplication (Convolution Theorem)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        a[i] = mul(a[i], b[i]);
    }

    // 6. Perform Inverse NTT
    ntt(a, true, roots);

    // 7. Carry Propagation (The "Schoolbook" cleanup)
    // This step is inherently sequential because carry depends on previous result
    std::vector<int> result;
    result.reserve(n);
    
    uint64_t carry = 0;
    for (int i = 0; i < n; i++) {
        uint64_t val = a[i] + carry;
        result.push_back(val % 10);
        carry = val / 10;
    }
    
    // Handle remaining carry
    while (carry) {
        result.push_back(carry % 10);
        carry /= 10;
    }

    // 8. Format Output
    // Remove trailing zeros (which are leading zeros in the number)
    while (result.size() > 1 && result.back() == 0) {
        result.pop_back();
    }
    
    // Convert back to string (reverse logic)
    std::string res_str;
    res_str.resize(result.size());
    
    #pragma omp parallel for
    for(size_t i=0; i < result.size(); i++) {
        res_str[i] = result[result.size() - 1 - i] + '0';
    }
    BigInt final_result(res_str);
    return final_result;
}




BigInt operator/(const BigInt& a, const BigInt& b) {
  return a.divmod(b.to_long()).quotient;
}



void BigInt::setNumerus(const std::vector<uint8_t> &source) {
  numerus = source;
}

BigInt BigInt::sqrt_bigint(const BigInt& n) {
  if (n == BigInt(0) || n == BigInt(1))
      return n;

  BigInt low(0);
  BigInt high = n;
  BigInt ans(0);

  while (low <= high) {
      BigInt mid = (low + high) / 2;
      BigInt midsq = mid * mid;

      if (midsq == n)
          return mid;
      else if (midsq < n) {
          low = mid + BigInt(1);
          ans = mid;
      } else {
          high = mid - BigInt(1);
      }
  }

  return ans;
}


std::vector<uint8_t> BigInt::getNumerus() const {
  std::vector<uint8_t> v = numerus;
  return v;
}

std::vector<int> BigInt::to_list() const {
  std::vector<int> result;
  for (auto x : numerus) {
    result.push_back((int)x);
  }
  return result;
}

long BigInt::to_long() const {
  long result = 0;
  for (int i = 0; i < numerus.size(); i++) {
    result = result * 10 + numerus[i];
  }
  return result;
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
    if (i < 0) {
      i = 0;
    }
    if (j < 0) {
      return BigInt("NaN");
    }
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

split BigInt::split_it(int m) const {
  int n = this->size();

  BigInt r;
  BigInt c;

  split z;
  if (n < m) {
    z.xleft = BigInt("NaN");
    z.xright = BigInt("NaN");
    z.m = m;
    return z;
  }
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

BigInt::BigInt(const std::vector<BigInt> &v, SIGN s) {

  std::vector<uint8_t> result;
  std::vector<uint8_t> temp;

  for (auto x : v) {
    temp = x.getNumerus();
    for (auto y : temp) {
      result.push_back(y);
    }
  }
  BigInt z(result, s);

  *this = z;
}

BigInt::BigInt(const std::vector<int> &v, SIGN s) {

  std::vector<uint8_t> result;
  for (auto x : v) {
    result.push_back((uint8_t)x);
  }
  numerus = result;
  sign = s;
}

BigInt::BigInt() {
  sign = _NULL;
  std::vector<uint8_t> v;
  numerus = v;
}

BigInt::BigInt(const std::vector<uint8_t> &num, SIGN s) {
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

BigInt BigInt::trim_zeros() const {
  int n = numerus.size();
  BigInt z = *this;
  std::vector<uint8_t> numerus = z.getNumerus();

  if (n == 0 || n == 1) {
    return *this;
  }
  int i = 0;
  while (numerus[i] == 0) {
    i++;
  }
  if (i > 0) {

    numerus.erase(numerus.begin(), numerus.begin() + i);
    z.setNumerus(numerus);
  }
  return z;
}

bool BigInt::is_digit(char ch) const {

  return (int(ch) >= int('0') && int(ch) <= int('9'));
}

BigInt::BigInt(const std::string c) {
  sign = POS;
  char ch;
  int n = c.size();
  if (n == 0 || c == "NaN") {
    sign = UNDEFINED;
    return;
  }

  ch = c[0];

  if (n == 1) {
    if (!is_digit(ch)) {
      sign = UNDEFINED;
      return;
    } else {
      uint8_t x = int(ch) - int('0');
      numerus.push_back(x);
      return;
    }
  }

  if (int(ch) != 45 && int(ch) != 43 && !is_digit(ch)) {
    sign = UNDEFINED;

    return;
  }
  if (int(ch) == 45) {
    sign = NEG;
  }
  if (int(ch) == 43) {
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
    } else {
      sign = UNDEFINED;
      return;
    }
  }
}

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

int BigInt::operator[](const int i) const { return (int)numerus[i]; }

void BigInt::set_sign(SIGN x) { sign = x; }

BigInt BigInt::abs() const {
  BigInt z = *this;
  if (z.get_sign() == NEG) {
    z.set_sign(POS);
  }
  return z;
}

SIGN BigInt::get_sign() const { return sign; }

int BigInt::size() const { return numerus.size(); }

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
  if (n == 1 && num[0] == 0) {
    out << static_cast<int>(num[0]);
    return out;
  }

  if (num[i] == 0) {

    while (num[i] == 0 && i < n) {
      // out << static_cast<int>(num[i]);
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

  std::vector<uint8_t> x_numerus = x.getNumerus();
  std::vector<uint8_t> y_numerus = y.getNumerus();

  int c = 0;
  int tot = 0;

  for (int i = k - 1; i >= 0; i--) {
    tot = x_numerus[i] + y_numerus[i] + c;

    if (tot >= 10) {
      c = 1;

      z.insert(tot % 10, 0);

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

std::string BigInt::to_string() const {
  std::string result;
  for (auto x : this->numerus) {
    result += std::to_string(x);
  }
  return result;
}

void BigInt::print() const {
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
  std::vector<uint8_t> x_numerus = x.getNumerus();
  std::vector<uint8_t> y_numerus = y.getNumerus();

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
  const int gmp_threshold1 = 2000;
  const int gmp_threshold2 = 100000;
  int n = x.size();
  int m = y.size();
  long order = n+m;

  std::vector<uint8_t> x_numerus = x.getNumerus();
  std::vector<uint8_t> y_numerus = y.getNumerus();

  if (order > gmp_threshold1 && order < gmp_threshold2) {
    return karatsuba(x, y);
  }

  if (order >= gmp_threshold2) {
    std::string x_str = x.to_string();  
    std::string y_str = y.to_string();
    return Schonhage_Strassen(x_str, y_str);
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
  if (this->sign == UNDEFINED || num.sign == UNDEFINED || this->sign == _NULL ||
      num.sign == _NULL) {
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

std::vector<BigInt> BigInt::split_number(const BigInt x, const int m) const {

  std::vector<uint8_t> numerus = x.getNumerus();
  std::vector<BigInt> result;

  for (int i = numerus.size(); i >= 0; i -= m) {
    BigInt temp = x.slice(i - m, i - 1);
    result.insert(result.begin(), temp);
  }
  return result;
}

divmod10 BigInt::burnikel_ziegler(const BigInt &x, const BigInt &y) const {

  divmod10 d;
  long ylong = y.to_long();
  const int m =
      4; // change this to some optimized value determined the the inputs.
  std::vector<BigInt> x_parts = split_number(x, m);

  std::vector<BigInt> y_parts = split_number(y, m);
  BigInt remainder = 0;
  long quotient_part;
  std::vector<BigInt> quotient;
  std::vector<int> quotient_list;

  for (auto &part : x_parts) {
    remainder = remainder * std::pow(10, m) + part;
    quotient_part = remainder.to_long() / ylong;
    quotient_list = BigInt(quotient_part).to_list();
    int nzeros = m - quotient_list.size();

    if (nzeros > 0) {
      for (int i = 0; i < nzeros; i++) {
        quotient_list.insert(quotient_list.begin(), 0);
      }
    }
    remainder = remainder - BigInt(quotient_part) * y;
    quotient.push_back(quotient_list);
  }

  d.quotient = BigInt(quotient).trim_zeros();
  d.remainder = remainder.trim_zeros();

  return d;
}


BigInt BigInt::operator/(const long n) const{
  BigInt x = *this;
  divmod10 d = x.divmod(n);
  return d.quotient;
}


divmod10 BigInt::div(const BigInt &num) const {
  BigInt x = *this;
  BigInt y = num;
  divmod10 d = burnikel_ziegler(x, y);

  return d;
}


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



BigInt& BigInt::operator=(const BigInt& other) {
  if (this != &other) { // Protect against self-assignment
      numerus = other.numerus;
      sign = other.sign;
  }
  return *this;
}

bool BigInt::operator==(const BigInt &y) const {

  BigInt x = *this;
  BigInt temp = y;
  SIGN xsign = get_sign();
  SIGN ysign = temp.get_sign();
  if (xsign != ysign) {
    return false;
  }
  // BigInt x1 = x.trim_zeros();
  // BigInt temp1 = temp.trim_zeros();
  int m = x.size();
  int n = temp.size();

  if (m != n) {

    return false;
  }

  for (int i = 0; i < n; i++) {
    if (x[i] != temp[i]) {
      return false;
    }
  }

  return true;
}

bool BigInt::operator<(const BigInt &y) const {

  if (*this == y) {
    return false;
  }

  SIGN xsign = get_sign();
  SIGN ysign = y.get_sign();
  if (xsign == _NULL || ysign == _NULL) {
    return false;
  }

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

  return false;
}

bool BigInt::operator<=(const BigInt &y) const {
  if (*this == y) {
    return true;
  }

  SIGN xsign = get_sign();
  SIGN ysign = y.get_sign();
  if (xsign == _NULL || ysign == _NULL) {
    return false;
  }

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

  return false;
}

bool BigInt::operator>(const BigInt &y) const {
  if (*this == y) {
    return false;
  }

  SIGN xsign = get_sign();
  SIGN ysign = y.get_sign();
  if (xsign == _NULL || ysign == _NULL) {
    return false;
  }

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

  if (xsign == _NULL || ysign == _NULL) {

    return false;
  }

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
