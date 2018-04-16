// Reset the phase of a newly calculated eigensystem, such that it resembles a
// deleted eigensystem
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
  const int num_evecs = para -> get_int("NEV");
  const int T = 35;
  //Set up needed objects
  Eigen::MatrixXcd V_old (para -> get_int("MAT_ENTRIES"),num_evecs);
  Eigen::MatrixXcd V (para -> get_int("MAT_ENTRIES"),num_evecs);
  Eigen::MatrixXcd V_fix (para -> get_int("MAT_ENTRIES"),num_evecs);
  Eigen::MatrixXcd V_reset (para -> get_int("MAT_ENTRIES"),num_evecs);
  std::vector<double> phase_old;
  std::vector<double> phase_new;

  // Read in data
  read_evectors_bin_ts("eigenvectors",cfg,T,num_evecs,V);
  read_evectors_bin_ts(para -> get_path("ev1").c_str(),"eigenvectors",cfg,T,num_evecs,V_old);
  // build filename as prefix
  std::string path_old_phases("/hiskp4/hiskp2/eigensystems/A80.24_L24_T48_beta190_mul0080_musig150_mudel190_kappa1632600/hyp_062_062_3/nev_120/phases/phases"); 
  read_eigenvalues_bin(path_old_phases.c_str(),cfg,T,num_evecs,phase_old);
  //for (auto po : phase_old) std::cout << po << std::endl;
  // Reset and fix the phases of the new system
  reset_phase(V,V_reset,phase_old);
  fix_phase(V_reset,V_fix,phase_new);
  write_eig_sys_bin("eigenvectors_reset",cfg,T,num_evecs,V_fix);
  write_eigenvalues_bin("eigenvectors_old_phase",cfg,T,num_evecs,phase_old);
  std::cout << "\nFirst 10 entries of Eigenvector 0:" << std::endl;
  std::cout << "\nold eigenvector:" << std::endl;
  for (size_t i = 0; i < 100; ++i) {
    std::cout << V_old(i,0) <<std::endl;
  }
  //std::cout << "\nnew eigenvectors:" << std::endl;
  //for (size_t i = 0; i < 10; ++i) {
  //  std::cout << V(i,0) <<std::endl;
  //}

  std::cout << "\nnew fixed eigenvectors:" << std::endl;
  for (size_t i = 0; i < 100; ++i) {
    std::cout << V_reset(i,0) <<std::endl;
  }
  return 0;
}
  
