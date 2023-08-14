#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <limits>



size_t bitrev(size_t n) {

    //const size_t K = sizeof(n)*8;
    const size_t K = 4;
    std::cout << "K = " << K << "\n";
    std::bitset<K> b(n);

    std::cout << n << " in binary is " << b << "\n";

    return 0;
}

std::complex<double> exponentiate(size_t k,size_t n,size_t N) {

    double x = 0.0;
    double y = 0.0;
    double theta = ((2* M_PI *k * n)/N);
   
    x = std::cos(theta);
    y = std::sin(theta);
    std::complex<double> omega(x,y);
    return omega;    
}


BigInt Schonhage_Strassen(BigInt &x,BigInt&y) {

    return x;
}



std::complex<double> dift(std::vector<std::complex<double>>& input, size_t n) {

    size_t N = input.size();
    std::complex<double> coeff(0.0,0.0);
    
    for (int k = 0; k < N; k++) {   
            coeff += input[k]*exponentiate(k,n,N);       
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

std::complex<double> dft(std::vector<std::complex<double>>& input,size_t n) {

    size_t N = input.size();
    std::complex<double> coeff(0.0,0.0);
    
    for (int k = 0; k < N; k++) {   
            coeff += input[k]/exponentiate(k,n,N);       
    }

    if (abs(coeff.real()) < 1e-6) {
        coeff.real(0.0);
    }

    if (abs(coeff.imag()) < 1e-6) {
        coeff.imag(0.0);    
    }

    return coeff;
} 

int main() {

  BigInt z("6789353555355553535353553535");
  std::cout << "z is " << z.size() << " digits long. \n";
  BigInt zslice1("67893535553");
  BigInt z1 = z.slice(-1,-1);
  BigInt z2 = z.slice(5,5);
  BigInt z3 = z.slice(0,1000);
  BigInt z4 = z.slice(27,34);
  BigInt z5 = z.slice(27,28);

  std::cout << "z2 = " << z2 <<"\n";
  std::cout << "z3 = " << z3 <<"\n";
  std::cout << "z4 = " << z4 <<"\n";
  std::cout << "z5 = " << z5 <<"\n";


  




  
}