# pchmatexport
A small header-only script to extract system matrices from a Nastran punch file (pch) and store on the memory

Example use (C++11):
```
#include <iostream>
#include <string>
#include "pchmatexport.h"

int main() {
  std::string pchAddress("C:/myPCHfile.pch");
  
  auto stifmat = PCH::read_matrix(pchAddress, PCH::MATRIX::STIF); // to extract the stiffness matrix
  auto massmat = PCH::read_matrix(pchAddress, PCH::MATRIX::MASS); // to extract the mass matrix
  
  std::cout << stifmat << std::endl; // possible since sysmat is an Eigen::MatrixXd
  
  return 0;
}

```
The only dependency is Eigen 3.3.7 or newer, which is a header-only library. 

**Important Note:**

Some PCH files Nastran generates may have outputs where the `*` lines that come after the `DMIG*` line have the last two columns contiguous as below: 

```
DMIG*   KAAX                           2               1   
*                      1               1-1.866666647D+07  <-- 1 and -1.866666647D+07 are contiguous
```
In this case, those two columns need to be separated as below before running the `read_matrix()` function:

```
DMIG*   KAAX                           2               1   
*                      1               1     -1.866666647D+07  <-- 1 and -1.866666647D+07 are not contiguous anymore
```

I will modify the script later to automatically account for this.
