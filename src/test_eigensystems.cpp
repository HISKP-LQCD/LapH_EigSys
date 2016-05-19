#include <cstdlib>
#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "eigensystem.h"
#include "navigation.h"
#include "par_io.h"
#include "read_write.h"

int main(int argc, char **argv) {


  //get number of configuration from last argument to main
  int cfg = atoi( argv[ (argc-1) ] );
  --argc;
  IO* para = IO::getInstance();
  para -> IO::set_values("parameters.txt");
  int num_evecs = para -> get_int("NEV");
  //lookup tables
  //Set up navigation
  Nav* lookup = Nav::getInstance();
  lookup -> init();
  Tslice* slice = Tslice::getInstance();
  slice -> Tslice::init();
  
  //build all combinations  
  std::complex<double> tr_uv_dag, sum_uv_dag, sum_uv_dag_abs, sum_uv_dag_abs2;
  std::complex<double> tr_uu_dag, sum_uu_dag;
  std::complex<double> tr_vu_dag, sum_vu_dag;
  std::complex<double> tr_vv_dag, sum_vv_dag;
for(int ts = 0; ts < 64; ++ts){ 
  Eigen::MatrixXcd M (para -> get_int("MAT_ENTRIES"),num_evecs);
  Eigen::MatrixXcd N (num_evecs, num_evecs);

  std::cout << "\n\n///////////////////////////////// ts = " << ts << "/////////////////////////////\n" << std::endl;
 //judge eigensystem
  Eigen::MatrixXcd U = M;
  Eigen::MatrixXcd U_fix = M;
  std::vector<double> phase_u;
  std::cout << "set up U: ";

  std::string in_1 = "cheb_static8/eigenvectors";
  read_evectors_bin_ts(in_1.c_str(), cfg, ts, num_evecs, U);
  //fix_phase(U, U_fix, phase_u);

  
  //transformed eigensystem
  Eigen::MatrixXcd V = M;
  Eigen::MatrixXcd V_fix = M;
  std::vector<double> phase_v;
  std::cout << "set up V: ";

  std::string in_2 = "cheb_dyn8/eigenvectors";
  read_evectors_bin_ts(in_2.c_str(),cfg, ts, num_evecs, V);
  std::cout << M.block(0,0,6,6) << std::endl;
  std::cout << " traces and sums of combinations: \n";
 
  //U_dag V
  N = ( U.adjoint() ) * V;
  tr_uv_dag = N.trace();
  sum_uv_dag = N.sum();
  sum_uv_dag_abs = (N.cwiseAbs()).sum();
  sum_uv_dag_abs2 = (N.cwiseAbs2()).sum();

  std::cout << "c-wise tr( U_dag * V) = " << tr_uv_dag << "\n" << 
    " sum( U_dag * V ) = " << sum_uv_dag << "\n" <<
    " sum.cwiseAbs( U_dag * V ) = " << sum_uv_dag_abs << "\n" <<
    " sum.cwiseAbs2( U_dag * V ) = " << sum_uv_dag_abs2 << std::endl; 
  
  std::cout << " (U.real().cwiseAbs()-V.real().cwiseAbs()).sum = " << 
          (U.real().cwiseAbs() - V.real().cwiseAbs()).sum() << std::endl;
  std::cout << " (U.imag().cwiseAbs()-V.imag().cwiseAbs()).sum = " <<
          (U.imag().cwiseAbs() - V.imag().cwiseAbs()).sum() << std::endl;

  //U_dag * U
  N = ( U.adjoint() ) * U;
  tr_uu_dag = N.trace();
  sum_uu_dag = N.sum();

  std::cout << "tr( U_dag * U) = " << tr_uu_dag 
    << " sum( U_dag * U ) = " << sum_uu_dag << std::endl; 

  //V_dag * V
  N = ( V.adjoint() ) * V;
  tr_vv_dag = N.trace();
  sum_vv_dag = N.sum();

  std::cout << "tr( V_dag * V) = " << tr_vv_dag 
    << " sum( V_dag * V ) = " << sum_vv_dag << std::endl; 
  
}
  return 0;
}
