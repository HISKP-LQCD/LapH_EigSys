#include "eigensystem.h"
static IO* const pars = IO::getInstance();
//fix phase of eigensystem and store phase of first entry of each eigenvector
void fix_phase(Eigen::MatrixXcd& V, Eigen::MatrixXcd& V_fix, std::vector<double>& phase) {
   const int V3 = pars -> get_int("V3");
  //helper variables:
  //Number of eigenvectors
  int n_ev;
  //negative imaginary
  std::complex<double> i_neg (0,-1);
  //tmp factor and phase
  std::complex<double> fac (1.,1.);
  double tmp_phase = 0;
  //get sizes right, resize if necessary
  n_ev = V.cols();
  if (phase.size() != n_ev) phase.resize(n_ev);
  if (V_fix.cols() != n_ev) V_fix.resize(3*V3,n_ev);
  //loop over all eigenvectors of system
  for (int n = 0; n < n_ev; ++n) {

    tmp_phase = std::arg(V(0,n));
    phase.at(n) = tmp_phase;
    fac = std::exp(i_neg*tmp_phase);
    //Fix phase of eigenvector with negative polar angle of first entry
    V_fix.col(n) = fac * V.col(n); 

  }
}
