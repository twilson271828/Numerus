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



int main() {

  BigInt x(0);
  std::cout << "x[0] = "<< (int)x[0] << "\n";
 
  return 0;
}
