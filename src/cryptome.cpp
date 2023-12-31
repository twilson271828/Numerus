#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <limits>
#include <iostream>


#if 0

void extend(BigInt &x, BigInt &y){

  unsigned long long n = x.size();
  unsigned long long m = y.size();

  unsigned long long p = 1;

  while(p < (n + m - 1)  ){
    p *= 2;
  }
 
  std::cout << "p-n = " << p-n << "\n";
  std::cout << "p-m = " << p-m << "\n";
  
  x=x.m10(p-n,true);
  y=y.m10(p-m,true);

 
}


template <size_t N>
bool bitset_less(const std::bitset<N>& lhs, const std::bitset<N>& rhs) {
  for (size_t i = N; i-- > 0; ) {
    if (lhs[i] && !rhs[i]) return false;
    if (!lhs[i] && rhs[i]) return true;
  }
  return false;
}

template <size_t N>
bool bitset_greater(const std::bitset<N>& lhs, const std::bitset<N>& rhs) {
  for (size_t i = N; i-- > 0; ) {
    if (lhs[i] && !rhs[i]) return  true;
    if (!lhs[i] && rhs[i]) return false;
  }
  return false;
}



std::bitset<4> convertToBinary(uint8_t &n) {

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

  // std::cout << "b = " << b << "\n";
  return b;
}

std::complex<double> exponentiate(size_t k, size_t n, size_t N) {

  double x = 0.0;
  double y = 0.0;
  double theta = ((2 * M_PI * k * n) / N);

  x = std::cos(theta);
  y = std::sin(theta);
  std::complex<double> omega(x, y);
  return omega;
}

std::complex<double> dift_coef(std::vector<std::complex<double>> &input,
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

std::vector<std::complex<double>> dift(std::vector<std::complex<double>> &X) {

  std::vector<std::complex<double>> y;
  size_t N = X.size();
  for (int ix = 0; ix < N; ix++) {
    std::complex<double> coef = dift_coef(X, ix);
    y.push_back(coef);
  }

  return y;
}

std::vector<std::complex<double>> n_roots_of_unity(int N) {

  std::vector<std::complex<double>> nroots;

  for (int k = 0; k < N; k++) {
    std::complex<double> root = exponentiate(k, 1, N);
    nroots.push_back(root);
  }
  return nroots;
}

std::complex<double> dft_coef(std::vector<std::complex<double>> &input,
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

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> &X) {

  std::vector<std::complex<double>> y;
  size_t N = X.size();
  for (int ix = 0; ix < N; ix++) {
    std::complex<double> coef = dft_coef(X, ix);
    y.push_back(coef);
  }

  return y;
}

template <typename T>
std::vector<T> filter(std::vector<T> &x, std::vector<int> &ix) {

  std::vector<T> y;
  for (auto &i : ix) {
    y.push_back(x[i]);
  }
  return y;
}

std::vector<std::complex<double>> initialize_y(int N) {

  std::vector<std::complex<double>> y;
  for (int i = 0; i < N; i++) {
    std::complex<double> v(0.0, 0.0);
    y.push_back(v);
  }
  return y;
}


std::vector<std::complex<double>> complexify_numerus(BigInt &x) {

  std::vector<std::complex<double>> complex_numerus;
  int N = x.size();
  for (int i = 0; i < N; i++) {
    std::complex<double> val(x[i], 0.0);
    complex_numerus.push_back(val);
  }
  return complex_numerus;
}


uint8_t convertToDec(std::bitset<4> x){

  return static_cast<uint8_t>(x.to_ulong());
}

BigInt find_nearest_power_of_2(BigInt &x) {

  std::vector<std::bitset<4>> binary_numerus = x.get_binary_numerus();
  BigInt y1;
  std::vector<std::bitset<4>> y = binary_numerus;
  size_t n = binary_numerus.size();
  for (std::bitset<4> x: y) {

    x.reset();

  }
  y[n-1].set(3);
  

  
  

return y1;
  
}


int find_nearest_power_of_2(uint64_t x) {

  int n = std::ceil(std::log(x) / std::log(2));
  return n;
}

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> &x,
                                      std::complex<double> omega) {

  int N = x.size();
  if (N == 1) {
    return x;
  } else {
    std::vector<int> even;
    std::vector<int> odd;

    for (int i = 0; i < N; i++) {
      if ((i % 2) == 0) {
        even.push_back(i);
      } else {
        odd.push_back(i);
      }
    }

    std::vector<std::complex<double>> a_even = filter(x, even);
    std::vector<std::complex<double>> a_odd = filter(x, odd);

    std::vector<std::complex<double>> y_even = fft(a_even, omega * omega);
    std::vector<std::complex<double>> y_odd = fft(a_odd, omega * omega);
    std::vector<std::complex<double>> y = initialize_y(N);

    int n2 = std::floor(N / 2);
    std::complex<double> x(1, 0);
    for (int i = 0; i < n2; i++) {
      y[i] = y_even[i] + x * y_odd[i];
      y[i + n2] = y_even[i] - x * y_odd[i];
      x = x * omega;
    }

    return y;
  }
}

BigInt Cooley_Tukey(BigInt &x, BigInt &y) { return x; }

BigInt Schonhage_Strassen(BigInt &x, BigInt &y) { return x; }

std::vector<std::complex<double>> mult(std::vector<std::complex<double>> X1,
                                       std::vector<std::complex<double>> X2) {

  std::vector<std::complex<double>> c;
  int n1 = X1.size();
  int n2 = X2.size();
  if (n1 != n2) {
    throw std::runtime_error("Attempted to element-wise multiply two vectors "
                             "with unequal dimensions.");
  }

  for (int i = 0; i < n1; i++) {
    c.push_back(X1[i] * X2[i]);
  }

  return c;
}

std::vector<std::complex<double>>
convolution(std::vector<std::complex<double>> X1,
            std::vector<std::complex<double>> X2) {

  std::vector<std::complex<double>> Z1 = mult(dft(X1), dft(X2));
  std::vector<std::complex<double>> Z2 = dift(Z1);
  return Z1;
}



#endif 


BigInt vadd(BigInt &x, BigInt &y) {
BigInt z;

  int n = x.size();
  int m = y.size();
  int k = std::max(n, m);
  if (n != m) {
    if (n > m) {
      int d = n - m;
      y = y.m16(d, true);
    } else {
      int d = m - n;
      x = x.m16(d, true);
    }
  }
  

  std::vector<int> x_numerus = x.get_numerus(); 
  std::vector<int> y_numerus = y.get_numerus();

  
  int c = 0;
  int tot = 0;

  for(int i = k-1;i>=0;i--) {
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

  return z;  
  }




  void printVector(std::vector<int> &x) {
    for (auto &i : x) {
      std::cout << i << " ";
    }
    std::cout << "\n";
  }

  BigInt vsub(BigInt &x, BigInt &y)  {
  
  int n = x.size();
  int m = y.size();
  int k = std::max(n, m);


  std::vector<std::bitset<4>> result(k);
  std::fill(result.begin(), result.end(), std::bitset<4>(0));

  if (n != m) {
    if (n > m) {
      int d = n - m;
      y = y.m16(d, true);
    } else {
      int d = m - n;
      x = x.m16(d, true);
    }
  }

  
std::vector<int> x_numerus = x.get_numerus(); 
std::vector<int> y_numerus = y.get_numerus();

  

  for (int i = k - 1; i >= 0; i--) {
    if (x_numerus[i] < y_numerus[i]) {                 
      int val = x_numerus[i-1] -1;           
      x_numerus[i-1] = val;      
       result[i] = std::bitset<4>(x_numerus[i] +10 - y_numerus[i]);      
      }
      else{
        result[i] = std::bitset<4>(x_numerus[i] - y_numerus[i]);        
      }        
    }

    int i =0;
    //Remove any leading zeros
    while(result[i].to_ulong() == 0){
      i++;
    }
    if (i > 0) {
    result.erase(result.begin(),result.begin()+i);

    }

   
  BigInt z(result);
  return z;
  }




int main() {

  BigInt x("200");
  BigInt y("2");

  

  
  return 0;
}