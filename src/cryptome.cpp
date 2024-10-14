#include "../include/BigInt.hpp"
#include <bitset>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <cmath>


template <typename T>
void printVector(std::vector<T> &x) {
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
  std::vector<BigInt> quotient;
  std::vector<int> quotient_list;

  for (auto &part : x_parts) {
    std::cout << "remainder = " << remainder << "\n";
    std::cout << "part = " << part << "\n";
    
    remainder = remainder *std::pow(10,m) + part;
    //std::cout << "result = " << remainder << "\n";
    quotient_part = remainder.to_long() / ylong;
    quotient_list = BigInt(quotient_part).to_list();
    int nzeros = m-quotient_list.size();

  if (nzeros>0){
            for (int i = 0; i < nzeros; i++) {
                quotient_list.insert(quotient_list.begin(),0);
            }
  }
    
    std::cout << "quotient_list = ";
    printVector(quotient_list);
    
    remainder = remainder - BigInt(quotient_part) * y;

    quotient.push_back(quotient_list);
   
    
    std::cout << "************************************************\n";
  }

  

  printVector(quotient);
  d.quotient=BigInt(quotient);
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

  BigInt x("7294372378472835723758");
  BigInt y("2568");
  
  int m = 4;

  divmod10 z = div(x,y);
  std::cout << "quotient = " << z.quotient << "\n";
  std::cout << "remainder = " << z.remainder << "\n";



  
  
  
  return 0;
}
