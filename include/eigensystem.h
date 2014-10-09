#ifndef _EIGENSYSTEM_H_
#define _EIGENSYSTEM_H_
#include <cstdlib>
#include <Eigen/Eigen>
#include <complex>
#include "par_io.h"
#include "read_write.h"
//fix phase of eigensystem and store phase of first entry of each eigenvector
void fix_phase(Eigen::MatrixXcd& V, Eigen::MatrixXcd& V_fix, std::vector<double>& phase);

#endif
