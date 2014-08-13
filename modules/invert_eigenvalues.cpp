#include <array>
#include <cstdlib>
#include <iostream>

#include "read_write.h"
#include "recover_spec.h"
int main() {
  
  //Set up data structure
  std::array<double, 120> evals;
  std::array<double, 120> test;
  read_eigenvalues_from_binary("eigenvalues_bin",600,0,test);

  //Spectrum Recovery
  recover_spectrum("eigenvalues_ape",600,0,evals);
  
  //Test output
  for (int k = 0; k < 120; ++k) std::cout << test.at(k) << std::endl;
  return 0;

}


