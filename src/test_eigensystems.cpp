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
  std::string inarch = "/work/hbn28/hbn282/eigensystems/D30arch/";
  std::string intest = "/work/hbn28/hbn282/eigensystems/D30newsingle/";
  std::string intest1 = "/work/hbn28/hbn282/eigensystems/D30new/";
  std::string intest2 = "/work/hbn28/hbn282/eigensystems/D30old/";
  std::string intest3 = "/work/hbn28/hbn282/eigensystems/D30newlc/";
  
  //build all combinations  
  std::complex<double> tr_uv_dag, sum_uv_dag, sum_uv_dag_abs, sum_uv_dag_abs2;
  std::complex<double> tr_uu_dag, sum_uu_dag;
  std::complex<double> tr_vu_dag, sum_vu_dag;
  std::complex<double> tr_vv_dag, sum_vv_dag;
for(int ts = 0; ts < 1; ++ts){ 
  Eigen::MatrixXcd M (para -> get_int("MAT_ENTRIES"),num_evecs);

  std::cout << "\n\n///////////////////////////////// ts = " << ts << "/////////////////////////////\n" << std::endl;
 //judge eigensystem
  Eigen::MatrixXcd U = M;
  Eigen::MatrixXcd U_fix = M;
  std::vector<double> phase_u;
  std::cout << "set up U: ";

  para -> set_inputpath(inarch);
  read_evectors_bin_ts("eigenvectors", cfg, ts, num_evecs, U);
  //fix_phase(U, U_fix, phase_u);

  
  //transformed eigensystem
  Eigen::MatrixXcd V = M;
  Eigen::MatrixXcd V_fix = M;
  std::vector<double> phase_v;
  std::cout << "set up V: ";

  para -> set_inputpath(intest3);
  read_evectors_bin_ts("eigenvectors",cfg, ts, num_evecs, V);
  //fix_phase(V, V_fix, phase_v);
  //M = U_fix-V_fix;
  //for(int c = 0; c < M.cols(); ++c){
  //  for(int r = 0; r < M.rows(); ++r){
  //    if (fabs(M(r,c)) > 10e-6){
  //      std::cout << r << " " << c << " " << M(r,c) << std::endl;
  //    }
  //  }
  //}
  Eigen::MatrixXcd W = M;
  Eigen::MatrixXcd W_fix = M;
  std::vector<double> phase_w;
  std::cout << "set up W: ";

  para -> set_inputpath(intest);
  read_evectors_bin_ts("eigenvectors",cfg, ts, num_evecs, W);

  std::cout << M.block(0,0,6,6) << std::endl;
  std::cout << " traces and sums of combinations: \n";
 
  //U_dag V
  Eigen::MatrixXcd N (num_evecs,num_evecs);
  std::cout << "U vs V" << std::endl;
  N = ( U.adjoint() ) * V;
  tr_uv_dag = N.trace();
  sum_uv_dag = N.sum();
  sum_uv_dag_abs = (N.cwiseAbs()).sum();
  sum_uv_dag_abs2 = (N.cwiseAbs2()).sum();

  std::cout << "c-wise tr( U_dag * V) = " << tr_uv_dag << "\n" << 
    " sum( U_dag * V ) = " << sum_uv_dag << "\n" <<
    " sum.cwiseAbs( U_dag * V ) = " << sum_uv_dag_abs << "\n" <<
    " sum.cwiseAbs2( U_dag * V ) = " << sum_uv_dag_abs2 << std::endl; 
  
  M = U.real().cwiseAbs() - V.real().cwiseAbs();
  std::cout << "Re((U.cwiseAbs() - V.cwiseAbs()).sum()) = " << M.sum() << std::endl;
  //std::cout << "colwise sum : " << M.colwise().sum() << std::endl;
  M = U.imag().cwiseAbs() - V.imag().cwiseAbs();
  std::cout << "Im((U.cwiseAbs() - V.cwiseAbs()).sum()) = " << M.sum() << std::endl;
  //std::cout << "colwise sum : " << M.colwise().sum() << std::endl;

  //U_dag W
  std::cout << "U vs W" << std::endl;
  N = ( U.adjoint() ) * W;
  tr_uv_dag = N.trace();
  sum_uv_dag = N.sum();
  sum_uv_dag_abs = (N.cwiseAbs()).sum();
  sum_uv_dag_abs2 = (N.cwiseAbs2()).sum();

  std::cout << "c-wise tr( U_dag * W) = " << tr_uv_dag << "\n" << 
    " sum( U_dag * W ) = " << sum_uv_dag << "\n" <<
    " sum.cwiseAbs( U_dag * W ) = " << sum_uv_dag_abs << "\n" <<
    " sum.cwiseAbs2( U_dag * W ) = " << sum_uv_dag_abs2 << std::endl; 
  
  M = U.real().cwiseAbs() - W.real().cwiseAbs();
  std::cout << "Re((U.cwiseAbs() - W.cwiseAbs()).sum()) = " << M.sum() << std::endl;
  //std::cout << "colwise sum : " << M.colwise().sum() << std::endl;
  M = U.imag().cwiseAbs() - W.imag().cwiseAbs();
  std::cout << "Im((U.cwiseAbs() - W.cwiseAbs()).sum()) = " << M.sum() << std::endl;
  //std::cout << "colwise sum : " << M.colwise().sum() << std::endl;
  //V_dag W
  std::cout << "V vs W" << std::endl;
  N = ( V.adjoint() ) * W;
  tr_uv_dag = N.trace();
  sum_uv_dag = N.sum();
  sum_uv_dag_abs = (N.cwiseAbs()).sum();
  sum_uv_dag_abs2 = (N.cwiseAbs2()).sum();

  std::cout << "c-wise tr( V_dag * W) = " << tr_uv_dag << "\n" << 
    " sum( V_dag * W ) = " << sum_uv_dag << "\n" <<
    " sum.cwiseAbs( V_dag * W ) = " << sum_uv_dag_abs << "\n" <<
    " sum.cwiseAbs2( V_dag * W ) = " << sum_uv_dag_abs2 << std::endl; 
  
  M = V.real().cwiseAbs() - W.real().cwiseAbs();
  //auto Mcol = M.colwise().sum();
  std::cout << "Re((V.cwiseAbs() - W.cwiseAbs()).sum()) = " << M.sum() << std::endl;
  //std::cout << "colwise sum : " << Mcol << std::endl;
  M = V.imag().cwiseAbs() - W.imag().cwiseAbs();
  //Mcol = M.colwise().sum();
  std::cout << "Im((V.cwiseAbs() - W.cwiseAbs()).sum()) = " << M.sum() << std::endl;
  //std::cout << "colwise sum : " << Mcol << std::endl;

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

  //W_dag * W
  N = ( W.adjoint() ) * W;
  tr_vv_dag = N.trace();
  sum_vv_dag = N.sum();

  std::cout << "tr( W_dag * W) = " << tr_vv_dag 
    << " sum( W_dag * W ) = " << sum_vv_dag << std::endl; 
  
}
  return 0;
}
