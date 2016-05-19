#ifndef _RECOVER_SPEC_H_
#define _RECOVER_SPEC_H_
#include <vector>
#include <Eigen/Eigen>

//Recover original eigenvalues from Chebyshev and map
void recover_spectrum(const int nb_ev, const std::vector<double>& evals_in, std::vector<double>& evals_out);
void recover_spectrum_large(const int nb_ev, const std::vector<double>& evals_in, std::vector<double>& evals_out);

//Invert 8th Chebyshev-Polynom of first kind
double invert_T8(const double value_in);

//Invert map
double invert_B(const double value_in);
double invert_Bprime(const double value_in);

#endif //_RECOVER_SPEC_H_
