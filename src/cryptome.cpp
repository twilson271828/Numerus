#include "../include/BigInt.hpp"
#include <bitset>



//bitrev(j, 8) gives 0, 4, 2, 6, 1, 5, 3, 7 for j = 0, . . . , 7


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

std::complex<double> dift(std::vector<std::complex<double>>& input, size_t n) {

    size_t N = input.size();
    std::complex<double> coeff(0.0,0.0);
    
    for (int k = 0; k < N; k++) {   
            coeff += input[k]*exponentiate(k,n,N);       
    }
    //coeff /= N;
    return coeff;
} 

std::complex<double> dft(std::vector<std::complex<double>>& input,size_t n) {

    size_t N = input.size();
    std::complex<double> coeff(0.0,0.0);
    
    for (int k = 0; k < N; k++) {   
            coeff += input[k]*exponentiate(-k,n,N);       
    }
    //coeff /= N;
    return coeff;
} 

int main() {

    std::complex<double> c0(1,0);
    std::complex<double> c1(2,-1);
    std::complex<double> c2(0,-1);
    std::complex<double> c3(-1,2);

std::vector<std::complex<double>> v = {c0,c1,c2,c3};

std::complex<double> x0 = dft(v,0);
std::complex<double> x1 = dft(v,1);
std::complex<double> x2 = dft(v,2);
std::complex<double> x3 = dft(v,3);

std::cout << "x0 = " << x0 << "\n";
std::cout << "x1 = " << x1 << "\n";
std::cout << "x2 = " << x2 << "\n";
std::cout << "x3 = " << x3 << "\n";
   
  
}