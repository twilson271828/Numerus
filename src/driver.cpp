#include "../include/BigInt.hpp"

int main() {

  BigInt c1("314159");
  BigInt c2("271828");

  BigInt z1 = c1 * c2;
  // BigInt z2 = c2*c1;

  // std::cout << "c1 = " << c1 << "\n";
  // std::cout << "c2 = " << c2 << "\n";
  std::cout << "z1 = " << z1 << "\n";
}