#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include "statistics.h"

int main() {
  int ev_num = 120;
  std::vector<double> average (ev_num,0);
  average = ev_average("/hiskp2/helmes/A60_0600_L120/eigensystems/hyp_05_07_03/eigenvalues", ev_num);
  for(auto& v : average) std::cout << v << std::endl;
  return 0;

}
