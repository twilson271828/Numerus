#include "../include/BigInt.hpp"
#include <bitset>

//bitrev(j, 8) gives 0, 4, 2, 6, 1, 5, 3, 7 for j = 0, . . . , 7

size_t bitrev(size_t n) {

    unsigned int NO_OF_BITS = sizeof(n) * 8;
    size_t reverse_num = 0;
    int i;
    for (i = 0; i < NO_OF_BITS; i++) {
        if ((n & (1 << i)))
            reverse_num |= 1 << ((NO_OF_BITS - 1) - i);
    }
    return reverse_num;

 
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

    //bitwise_operators();
    for (int i = 0; i < 8;i++) {
      std::cout <<bitrev(i) <<std::endl;

    }
}