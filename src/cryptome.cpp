#include "../include/BigInt.hpp"
#include <bitset>

size_t bitrev(size_t n,size_t k) {

 size_t rev = 0;
 
    // traversing bits of 'n' from the right
    for (size_t i =0;i < k;i++) {
      ;
    }
 
    // required number
    return rev;
}
int bitwise_operators(){

 // a = 5(00000101), b = 9(00001001)
    int a = 5, b = 9;
 
    // The result is 00000001
    std::cout<<"a = " << a <<","<< " b = " << b <<std::endl;
    std::cout << "a & b = " << (a & b) << std::endl;
 
    // The result is 00001101
    std::cout << "a | b = " << (a | b) << std::endl;
 
    // The result is 00001100
    std::cout << "a ^ b = " << (a ^ b) << std::endl;
 
    // The result is 11111010
    std::cout << "~a = " << (~a) << std::endl;
 
    // The result is 00010010
    std::cout<<"b << 1" <<" = "<< (b << 1) <<std::endl;
 
    // The result is 00000100

    return 0;
}

int main() {

    bitwise_operators();
    //for (size_t i = 0;i < 8;i++){
    //std::cout << bitrev(i,64) <<"\n";
    //}
    //return 0;
}