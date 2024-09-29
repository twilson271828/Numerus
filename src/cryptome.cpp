#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

divmod10 div(BigInt &x,BigInt &y){

  divmod10 d;
  d.quotient=x;
  d.remainder=y;
  return d;
}

void printVector(std::vector<int> &x) {
  for (auto &i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}
bool isLeadingZeroPresent(const std::string &c) {
  if (c.empty()) {
    return false;
  }
  if (c[0] == '0') {
    return true;
  }
  return false;
}

std::vector<BigInt> split_number(BigInt x,int m) {

  std::vector<uint8_t> numerus=x.get_numerus();
  std::vector<BigInt> result;
  
   for (int i = numerus.size() - 1; i >= 0; i -= m) {
    std::cout << "[i-m,i] = [" << i-m << "," << i  << "]\n";
    BigInt temp = x.slice(i- m,i);
    std::cout << "temp = " <<temp << "\n";
      
    }

  
  
  
  return result;



}


void print_v(std::vector<BigInt> &x) {
  for (auto &i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}


int main() {

  BigInt x("7294372378472835723758");
  int n = x.size();
  std::cout << "n = " << n << "\n";
  int m = 3;
  std::vector<BigInt> result = split_number(x, m);
  //BigInt v = x.slice(n-3,n-1);
 // std::cout << "v = " << v << "\n";
  return 0;
}
