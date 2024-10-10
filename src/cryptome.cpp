#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>

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
  long ylong = y.to_long();
  int m = 4;
  std::vector<BigInt> x_parts = split_number(x, m);
  printVector(x_parts);
  std::vector<BigInt> y_parts = split_number(y, m);
  BigInt remainder = 0;
  long quotient_part;
  std::vector<BigInt> quotient_list;

  for (auto &part : x_parts) {
    std::cout << "remainder = " << remainder << "\n";
    std::cout << "part = " << part << "\n";
    //std::cout << "10**m = " << std::pow(10, m) << "\n";
    //std::cout << "part = " << part << "\n";
    remainder = remainder *std::pow(10,m) + part;
    //std::cout << "result = " << remainder << "\n";
    quotient_part = remainder.to_long() % ylong;
    remainder = quotient_part;
    //std::cout << "quotient_part = " << quotient_part << "\n";
    quotient_list.insert(quotient_list.begin(),BigInt(quotient_part));


    std::cout << "************************************************\n";
  }

  int nzeros = m-quotient_list.size();

  if (nzeros>0){
            for (int i = 0; i < nzeros; i++) {
                quotient_list.insert(quotient_list.begin(),BigInt(0));
            }
  }

  printVector(quotient_list);

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
