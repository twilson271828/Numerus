#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

void printVector(std::vector<BigInt> &x) {
  for (auto &i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}

std::vector<BigInt> split_number(BigInt x,int m) {

  std::vector<uint8_t> numerus=x.get_numerus();
  std::vector<BigInt> result;
  
   for (int i = numerus.size(); i >= 0; i -= m) {    
    BigInt temp = x.slice(i- m,i-1);
    //result.push_back(temp);
    result.insert(result.begin(),temp);      
    }
  return result;

}


divmod10 div(BigInt &x,BigInt &y){

  divmod10 d;
  int m = 4;
  std::vector<BigInt> x_parts = split_number(x, m);
  std::vector<BigInt> y_parts = split_number(y, m);
  BigInt remainder = 0;
  std::vector<BigInt> quotient_part;

  for (auto &part : x_parts) {
    remainder = remainder * pow(10, m) + part;
    std::cout << "remainder: " << remainder << "\n";
  }


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




void print_v(std::vector<BigInt> &x) {
  for (auto &i : x) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}


int main() {

  BigInt x("7294372378472835723758");
  BigInt y("2568");
  int m = 4;
  divmod10 z = div(x,y);
  
  
  
  return 0;
}
