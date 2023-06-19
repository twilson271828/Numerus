#include "../include/BigInt.hpp"
#include <bitset>
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

    std::complex<double> c0(1,0);
    std::complex<double> c1(2,-1);
    std::complex<double> c2(0,-1);
    std::complex<double> c3(-1,2);

    std::complex<double> c4(2,0);
    std::complex<double> c5(-2,-2);
    std::complex<double> c6(0,-2);
    std::complex<double> c7(4,4);

std::vector<std::complex<double>> v = {c0,c1,c2,c3};
std::vector<std::complex<double>> V = {c4,c5,c6,c7};

std::complex<double> x0 = dft(v,0);
std::complex<double> x1 = dft(v,1);
std::complex<double> x2 = dft(v,2);
std::complex<double> x3 = dft(v,3);

std::cout << "x0 = " << x0 << "\n";
std::cout << "x1 = " << x1 << "\n";
std::cout << "x2 = " << x2 << "\n";
std::cout << "x3 = " << x3 << "\n";
   
std::complex<double> x4 = dift(V,0);
std::complex<double> x5 = dift(V,1);
std::complex<double> x6 = dift(V,2);
std::complex<double> x7 = dift(V,3);
std::cout << "\n";
std::cout << "x4 = " << x4 << "\n";
std::cout << "x5 = " << x5 << "\n";
std::cout << "x6 = " << x6 << "\n";
std::cout << "x7 = " << x7 << "\n";
   
  
}