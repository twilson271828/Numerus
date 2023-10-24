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


std::vector<std::complex<double> >  fft(BigInt &x, std::complex<double> omega) {


   std::vector<std::complex<double>> fft_numerus;
   std::vector<uint8_t> numerus = x.get_numerus();

   int N = numerus.size();
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

  std::vector<std::complex<double> > a_even = filter(numerus,even);

   return fft_numerus;

    

}


BigInt Schonhage_Strassen(BigInt &x, BigInt &y) { return x; }



int main() {

  #if 0
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

  #endif
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

 

  for (auto & ix: odds) {
    std::cout << ix << "\n";
  }

  for (auto & ix: evens) {
    std::cout << ix << "\n";
  }

   std::vector<std::complex<double> > even_nroots = filter(nroots,evens);

  int ix = 0;
  for (auto & root: nroots){
    std::cout << ix << " : " << root << "\n";
    ix+=1;
  }
std::cout << "********************\n";
ix = 0;
   for(auto & root: even_nroots){
    std::cout << ix << " : " << root << "\n";
    ix+=1;
   }

  
 return 0;
 
}