#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <limits>


std::complex<double> exponentiate(size_t k, size_t n, size_t N) {

  double x = 0.0;
  double y = 0.0;
  double theta = ((2 * M_PI * k * n) / N);

  x = std::cos(theta);
  y = std::sin(theta);
  std::complex<double> omega(x, y);
  return omega;
}

std::complex<double> dift_coef(std::vector<std::complex<double>> &input, size_t n) {
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


std::vector<std::complex<double>> dift(std::vector<std::complex<double>> &X){

  std::vector<std::complex<double> > y;
  size_t N = X.size();
  for (int ix =0;ix < N;ix++) {
    std::complex<double> coef = dift_coef(X,ix);
    y.push_back(coef);

  }

  return y;

}


std::vector<std::complex<double> > n_roots_of_unity(int N){

  std::vector<std::complex<double> > nroots;

  for (int k = 0; k < N; k++) {
    std::complex<double> root = exponentiate(k,1,N);
    nroots.push_back(root);
  }
  return nroots;

}

std::complex<double> dft_coef(std::vector<std::complex<double>> &input, size_t n) {

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


std::vector<std::complex<double>> dft(std::vector<std::complex<double>> &X){

std::vector<std::complex<double> > y;
  size_t N = X.size();
  for (int ix =0;ix < N;ix++) {
    std::complex<double> coef = dft_coef(X,ix);
    y.push_back(coef);

  }

  return y;

}


template <typename T> std::vector<T>  filter(std::vector<T>  &x,std::vector<int> &ix) {

    std::vector<T> y;
    for(auto & i:ix){
      y.push_back(x[i]);
    }
  return y;
}

std::vector<std::complex<double> > initialize_y(int N) {

  std::vector<std::complex<double> > y;
  for (int i = 0;i < N; i++){
    std::complex<double> v(0.0,0.0);
    y.push_back(v);
    }
    return y;

}


std::vector<std::complex<double> > complexify_numerus(BigInt &x) {

  std::vector<std::complex<double> > complex_numerus;
  int N = x.size();
  for (int i = 0;i < N;i++){
    std::complex<double> val(x[i],0.0);
    complex_numerus.push_back(val);

  }
  return complex_numerus;
}

std::vector<std::complex<double> >  fft(std::vector<std::complex<double> > &x, std::complex<double> omega) {


   
   int N = x.size();
   if (N == 1) {
    return x;
   }
   else {
   std::vector<int> even;
   std::vector<int> odd;

   for (int i =0;i < N; i++){
    if ( (i % 2) == 0) {
    even.push_back(i);
    }
    else {
      odd.push_back(i);
    }
  } 

  std::vector<std::complex<double>>  a_even = filter(x,even);
  std::vector<std::complex<double> >  a_odd = filter(x,odd);

  std::vector<std::complex<double> > y_even = fft(a_even,omega*omega);
  std::vector<std::complex<double> > y_odd = fft(a_odd,omega*omega);
  std::vector<std::complex<double> > y = initialize_y(N);

  int n2 = std::floor(N/2);
  std::complex<double> x(1,0);
  for (int i =0; i < n2;i++ ) {
    y[i]=y_even[i]+x*y_odd[i];
    y[i+n2]=y_even[i]- x*y_odd[i];
    x = x*omega;
  }
   
   return y;
   }

}

BigInt Cooley_Tukey(BigInt &x, BigInt &y) {
  return x;
}

BigInt Schonhage_Strassen(BigInt &x, BigInt &y) { return x; }

std::vector<std::complex<double> > mult(std::vector<std::complex<double> >X1,std::vector<std::complex<double> >X2) {

  std::vector<std::complex<double> > c;
  int n1 = X1.size();
  int n2 = X2.size();
  if (n1 != n2){
    throw std::runtime_error("Attempted to element-wise multiply two vectors with unequal dimensions.");
  }

  for (int i = 0; i < n1; i++) {
    c.push_back(X1[i]*X2[i]);
  }

  return c;

}

std::vector<std::complex<double> > convolution(std::vector<std::complex<double>> X1,std::vector<std::complex<double>> X2) {

    std::vector<std::complex<double>> Z1 = mult(dft(X1),dft(X2));
    std::vector<std::complex<double>> Z2 = dift(Z1);
    return Z1;    


}




int main() {

std::vector< std::complex<double> > X1 ={{1,0},{2,-1 },{0,-1},{-1,2}};
std::vector< std::complex<double> > X2 ={{2,0},{4,-2 },{0,-2},{-2,4}};

std::vector<std::complex<double>> Z1 = convolution(X1,X2);

for (auto x : Z1) {
      std::cout << x;
    }
    std::cout << "\n\n";




  std::vector< std::complex<double> > X ={{1,0},{2,-1 },{0,-1},{-1,2}};
  std::cout << "X[0] = " << X[0]<< "\n";
  std::cout << "X[1] = " << X[1]<< "\n";
  std::cout << "X[2] = " << X[2]<< "\n";
  std::cout << "X[3] = " << X[3]<< "\n";

  std::vector<std::complex<double>> Y = dft(X);

  std::cout << "Y[0] = " << Y[0] << "\n";
  std::cout << "Y[1] = " << Y[1] << "\n";
  std::cout << "Y[2] = " << Y[2] << "\n";
  std::cout << "Y[3] = " << Y[3] << "\n";
 

  std::vector<std::complex<double>> Z = dift(Y);

  std::cout << "Z[0] = " << Z[0] << "\n";
  std::cout << "Z[1] = " << Z[1] << "\n";
  std::cout << "Z[2] = " << Z[2] << "\n";
  std::cout << "Z[3] = " << Z[3] << "\n";

  #if 0
  int N = 10;
  std::vector<std::complex<double> > nroots = n_roots_of_unity(N);
  
  std::vector<int> evens;
  std::vector<int> odds;

  for (int i =0;i < N; i++){
    if ( (i % 2) == 0) {
    evens.push_back(i);
    }
    else {
      odds.push_back(i);
    }
  } 

#endif
 return 0;
 
}