#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <limits>


BigInt karatsuba(BigInt &x, BigInt &y) {

    size_t n = x.size(); 
    size_t m = y.size();
    size_t k = std::max(n,m);
    size_t k2 = std::floor(k/2);

    if (n < 2 || m < 2) {

        return x*y;
    }

   

    split split_x = x.split_it(k2);
    split split_y = y.split_it(k2);
    //BigInt z2 = split_x.x1*split_y.x1;
    //BigInt z1 = split_x.x1*split_y.x0+split_x.x0*split_y.x1;
    //BigInt z0 =split_x.x0*split_y.x0;
    BigInt low_x = split_x.x0;
    BigInt low_y = split_y.x0;
    BigInt high_x = split_x.x1;
    BigInt high_y = split_y.x1;
    
    BigInt xsum = low_x+high_x;
    BigInt ysum = low_y+high_y;

    BigInt z0 = karatsuba(low_x,low_y);
    BigInt z1 = karatsuba(xsum,ysum);
    BigInt z2 = karatsuba(high_x,high_y);

    BigInt first = z2.m10(k2*2,false);
    BigInt second = z1-z2-z0;
    second = second.m10(k2,false);
    BigInt third = z0;

    return first + second + third;

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

  BigInt z1("12345");
  BigInt z2("6789");
  BigInt z3 = karatsuba(z1,z2);
  std::cout << "z3 = " << z3 << "\n";


}