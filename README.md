# pchmatexport
A small header-only script to extract system matrices from a Nastran punch file (pch) and store on the memory

Example use (C++11):
```
#include <iostream>
#include <string>
#include "pchmatexport.h"

int main() {
  std::string pchAddress("C:/myPCHfile.pch");
  
  auto sysmat = PCH::read_matrix(pchAddress, PCH::MATRIX::STIF); // to extract the stiffness matrix
  auto sysmat = PCH::read_matrix(pchAddress, PCH::MATRIX::MASS); // to extract the mass matrix
  
  std::cout << sysmat << std::endl; // possible since sysmat is an Eigen::MatrixXd
  
  return 0;
}

```
The only dependency is Eigen 3.3.7 or newer, which is a header-only library. 
