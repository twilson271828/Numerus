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

BigInt karatsuba(BigInt &x, BigInt &y)  {

    size_t n = x.size(); 
    size_t m = y.size();

    if (n < 10 && y < 10) {

        return x*y;
    }

    size_t k = std::max(n,m);
    size_t k2 = std::floor(k/2);

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

    
    return x;
}


int main() {

  size_t k = 10;
  BigInt z("6789424643665457123213125523442134324234242352342380724234242");

  split split10 = z.split_it(k);

  std:: cout << "split10.x0 = "<< split10.x0 << "\n";
  std:: cout << "split10.x1 = "<< split10.x1 << "\n";

  size_t m = 0;
  split split0=z.split_it(m);
  std:: cout << "split0.x0 = "<< split0.x0 << "\n";
  std::cout << "split0.x0 sign = " << split0.x0.get_sign() << "\n";
  std:: cout << "split0.x1 = "<< split0.x1 << "\n";
  std::cout << "z.size() = " << z.size() << "\n";


  size_t n = z.size();
  split split_n = z.split_it(n);
  std:: cout << "split_n.x0 = "<< split_n.x0 << "\n";
  std:: cout << "split_n.x1 = "<< split_n.x1 << "\n";

}