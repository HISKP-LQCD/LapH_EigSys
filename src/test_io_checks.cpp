#include <cstdlib>
//measure time
#include <chrono>
#include <limits>
#include <string>
#include <Eigen/Core>
#include "config_utils.h"
#include "eigensystem.h"
#include "io.h"
#include "navigation.h"
#include "par_io.h"
#include "timeslice.h"
#include "recover_spec.h"
#include "read_write.h"

int main(int argc, char **argv) {
  //--------------------------------------------------------------------------//
  //                             Local variables                              //
  //--------------------------------------------------------------------------//
  std::cout << "Maximal streamsize: " << std::numeric_limits<std::streamsize>::max() << std::endl; 
  std::cout << "Minimal streamsize: " << std::numeric_limits<std::streamsize>::min() << std::endl; 
  std::cout << "Maximal long uint: " << std::numeric_limits<unsigned long int>::max() << std::endl; 
  std::cout << "Minimal long uint: " << std::numeric_limits<unsigned long int>::min() << std::endl; 
  //Handling infile
  IO* pars = IO::getInstance();
  pars -> set_values("parameters.txt");
  pars -> print_summary();
  //in and outpaths
  std::string GAUGE_FIELDS = pars -> get_path("conf");
  //lattice layout from infile
  const int L0 = pars -> get_int("LT");
  const int L1 = pars -> get_int("LX");
  const int L2 = pars -> get_int("LY");
  const int L3 = pars -> get_int("LZ");
  const int V3 = pars -> get_int("V3");
  const int V_TS = pars -> get_int("V_TS");
  //calculation parameters from infile
  const int NEV = pars -> get_int("NEV");
  const int V_4_LIME = pars -> get_int("V4_LIME");
  const int MAT_ENTRIES = pars -> get_int("MAT_ENTRIES");
  //chebyshev parameters
  const double LAM_L = pars -> get_float("lambda_l");
  const double LAM_C = pars -> get_float("lambda_c"); 
  //hyp-smearing parameters
  const double ALPHA_1 = pars -> get_float("alpha_1");
  const double ALPHA_2 = pars -> get_float("alpha_2");
  const int ITER = pars -> get_int("iter");

  //get number of configuration from last argument to main
  int config;
  config = atoi( argv[ (argc-1) ] );
  --argc;
  int ts = 62;
  int nev = NEV;
  Eigen::MatrixXcd eigensystem(5, 4);

  //std::chrono::high_resolution_clock::time_point now_rd = std::chrono::system_clock::now();
  //read_evectors_bin_ts("eigenvectors", config, ts, nev, eigensystem);
  //std::chrono::high_resolution_clock::time_point then_rd = std::chrono::system_clock::now();
  //auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (then_rd-now_rd).count();
  //std::cout << "Read takes: " << duration << " ms" << std::endl;

  std::chrono::high_resolution_clock::time_point now = std::chrono::system_clock::now();
  write_eig_sys_bin("eigenvectors", config, ts, nev, eigensystem);
  std::chrono::high_resolution_clock::time_point then = std::chrono::system_clock::now();
  auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds> (then-now).count();
  std::cout << "Write takes: "<< duration1 << " ms" << std::endl;
  std::cout << "Write complete" << std::endl;

  return 0;
}
