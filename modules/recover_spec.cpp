#include <array>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "par_io.h"
#include "recover_spec.h"
#include "read_write.h"

static IO* const pars = IO::getInstance();

//Invert map
double invert_B(const double value_in) {

  //Parameters for projection
  double lambda_l = pars -> get_float("lambda_l");
  double lambda_c = pars -> get_float("lambda_c");
  return( (lambda_l-lambda_c)*(value_in-1)*0.5 - lambda_c );

}

//TODO: Implement variable recovery
//Invert nth Chebyshev-Polynomial of first kind
double invert_Tn(const double value_in){
  const int n = pars -> get_int("DEG");
  //const int n = 8;
  double inv = 0;
  //double tmp = fabs(fmod(value_in,M_PI)-1);
  inv = cosh(1./n*acosh(value_in));
  return inv;
}
//Invert 8th Chebyshev-Polynomial of first kind
double invert_T8(const double value_in) {

  double value_out;
  value_out = sqrt( value_in + 1 );
  value_out *= sqrt(2.);
  value_out = sqrt( 2 + sqrt( 2. + value_out ) );
  value_out *= 0.5;
  return value_out;

}

//Recover original eigenvalues from Chebyshev and map
void recover_spectrum(const int nb_ev, const std::vector<double>& evals_in, std::vector<double>& evals_out) {
  if (evals_out.size() != nb_ev) evals_out.resize(nb_ev, 0);
  for (unsigned int i = 0; i < evals_in.size(); ++i) {

    evals_out.at(i) = invert_Tn( evals_in.at(i) );
    evals_out.at(i) = invert_B(evals_out.at(i));

  }

}
