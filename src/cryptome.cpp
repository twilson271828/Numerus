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

std::vector<BigInt> split_number(BigInt x,size_t m) {

  std::vector<uint8_t> numerus=x.get_numerus();
  std::vector<BigInt> result;

   for (int i = numerus.size() - 1; i >= 0; --i) {
        std::cout << (int)numerus[i] << " ";
    }
  std::cout << "\n";
  
  return result;



}



int main() {

  BigInt x("7294372378472835723758");
  int m = 3;
  std::vector<BigInt> result = split_number(x, m);
 
  return 0;
}
