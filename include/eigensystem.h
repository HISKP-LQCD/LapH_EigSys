#ifndef _EIGENSYSTEM_H_
#define _EIGENSYSTEM_H_
#include <cstdlib>
#include <Eigen/Eigen>
#include <complex>
#include "par_io.h"
#include "read_write.h"
//fix phase of eigensystem and store phase of first entry of each eigenvector
void fix_phase(Eigen::MatrixXcd& V, Eigen::MatrixXcd& V_fix,
               std::vector<double>& phase);
//reset the phase of an eigensystem to the given one
//reset_phase calculates the phases of V, changes them to phase_want and stores
//the result in V_reset
void reset_phase(Eigen::MatrixXcd& V, Eigen::MatrixXcd& V_reset,
                 std::vector<double>& phase_want);
#endif
