#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <Eigen/Eigen>
#include "variables.h"
#include "read_write.h"

int main() {
  const int dim_row = MAT_ENTRIES;
  const int dim_col = 120;// number of eigenvectors
  //Set up data structure
  Eigen::MatrixXcd V;
  Eigen::MatrixXcd V_dagger;
  Eigen::MatrixXcd S;
  V = Eigen::MatrixXcd::Zero( dim_row, dim_col);
  V_dagger = Eigen::MatrixXcd::Zero( dim_row, dim_col);
  S = Eigen::MatrixXcd::Zero( dim_row, dim_col);
  read_eigenvectors_from_binary_ts("eigenvectors", 600, 0, V);
  //read_eigenvectors_from_binary_ts("eigenvectors", 600, 0);
  std::cout << "Calculating v v_dagger" << std::endl;
  V_dagger = V.adjoint();
 
  for(int i = 0; i < dim_row; ++i ){
    for(int j = 0; j < dim_col; ++j) {
      //std::complex<double> entry =  V.row(i) * V_dagger.col(j);
      if (std::abs(V(i,j)) > 1e-15)  
      std::cout << i << " " << j << " " << V(i,j) << std::endl;
    }
  }
}

