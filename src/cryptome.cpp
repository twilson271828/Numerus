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

int main() {

    std::complex<double> c0(1,0);
    std::complex<double> c1(1,-2);
    std::complex<double> c2(0,-1);
    std::complex<double> c3(-1,2);

std::vector<std::complex<double>> v = {c0,c1,c2,c3};
   

 
  
}