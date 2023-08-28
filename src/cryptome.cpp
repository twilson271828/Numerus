#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <limits>


BigInt karatsuba(BigInt &x, BigInt &y) {

    size_t n = x.size(); 
    size_t m = y.size();
    size_t k = std::max(n,m);
    size_t k2 = std::floor(k/2);

    if (n == 1) {
        return x*y;
    }

    split split_x = x.split_it(k2);
    split split_y = y.split_it(k2);
 
    BigInt low_x = split_x.xright;
    BigInt low_y = split_y.xright;
    BigInt high_x = split_x.xleft;
    BigInt high_y = split_y.xleft;
    
    BigInt xsum = low_x+high_x;
    BigInt ysum = low_y+high_y;


    BigInt U = karatsuba(high_x,high_y);
    BigInt V = karatsuba(low_x,low_y);
    BigInt w1 = high_x - low_x;
    BigInt w2 = high_y - low_y;

    BigInt W = karatsuba(w1,w2);
    
   
    BigInt Z = U+V-W;

    BigInt P = U.m10(k2*2,false) + Z.m10(k2,false)+V;

    return P;

}



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

  BigInt z1("12345976");
  BigInt z2("67897787");

  size_t n = z1.size(); 
  size_t m = z2.size();
  size_t k = std::max(n,m);
 
  size_t k2 = 3;

  split z3 = z1.split_it(k2);
  std::cout << "z3.xright = " << z3.xright << "\n";
  std::cout << "z3.xleft = " << z3.xleft << "\n";
  std::cout << "z3.m = " << z3.m << "\n";

  split z4 = z2.split_it(k2);
  std::cout << "z4.xleft = " << z4.xleft << "\n";
  std::cout << "z4.xright = " << z4.xright << "\n";
  std::cout << "z4.m = " << z4.m << "\n";

  BigInt z5 = karatsuba(z1,z2);
  std::cout << "z5 = " << z5 <<"\n";

  BigInt z6 = z1*z2;

  std::cout << "z6 = " << z6 << "\n";


}