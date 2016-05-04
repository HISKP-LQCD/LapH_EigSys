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

  std::cout << "\n\n///////////////////////////////// ts = " << ts << "/////////////////////////////\n" << std::endl;
 //judge eigensystem
  Eigen::MatrixXcd U = M;
  Eigen::MatrixXcd U_fix = M;
  std::vector<double> phase_u;
  std::cout << "set up U: ";

  std::string in_1 = para -> get_path("ev_1");
  in_1 += "/eigenvectors";
  read_evectors_bin_ts(in_1.c_str(), cfg, ts, num_evecs, U);
  //fix_phase(U, U_fix, phase_u);

  
  //transformed eigensystem
  Eigen::MatrixXcd V = M;
  Eigen::MatrixXcd V_fix = M;
  std::vector<double> phase_v;
  std::cout << "set up V: ";

  std::string in_2 = para -> get_path("ev_2");
  in_2 += "/eigenvectors";
  read_evectors_bin_ts(in_2.c_str(),cfg, ts, num_evecs, U);
  //fix_phase(V, V_fix, phase_v);
  //M = U_fix-V_fix;
  //for(int c = 0; c < M.cols(); ++c){
  //  for(int r = 0; r < M.rows(); ++r){
  //    if (fabs(M(r,c)) > 10e-6){
  //      std::cout << r << " " << c << " " << M(r,c) << std::endl;
  //    }
  //  }
  //}
  std::cout << M.block(0,0,6,6) << std::endl;
  std::cout << " traces and sums of combinations: \n";
 
  //U_dag V
  M = ( U.adjoint() ) * V;
  tr_uv_dag = M.trace();
  sum_uv_dag = M.sum();
  sum_uv_dag_abs = (M.cwiseAbs()).sum();
  sum_uv_dag_abs2 = (M.cwiseAbs2()).sum();

  std::cout << "c-wise tr( U_dag * V) = " << tr_uv_dag << "\n" << 
    " sum( U_dag * V ) = " << sum_uv_dag << "\n" <<
    " sum.cwiseAbs( U_dag * V ) = " << sum_uv_dag_abs << "\n" <<
    " sum.cwiseAbs2( U_dag * V ) = " << sum_uv_dag_abs2 << std::endl; 
  
  std::cout << (U.cwiseAbs2() - V.cwiseAbs2()).sum() << "\n";

  //U_dag * U
  M = ( U.adjoint() ) * U;
  tr_uu_dag = M.trace();
  sum_uu_dag = M.sum();

  std::cout << "tr( U_dag * U) = " << tr_uu_dag 
    << " sum( U_dag * U ) = " << sum_uu_dag << std::endl; 

  //V_dag * V
  M = ( V.adjoint() ) * V;
  tr_vv_dag = M.trace();
  sum_vv_dag = M.sum();

  std::cout << "tr( V_dag * V) = " << tr_vv_dag 
    << " sum( V_dag * V ) = " << sum_vv_dag << std::endl; 
  
}
  return 0;
}
