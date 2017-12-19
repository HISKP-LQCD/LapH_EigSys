#ifndef SOURCESHAPE_H_
#define SOURCESHAPE_H

#include <complex>
#include <array>
#include <fstream>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>
#include "Eigen/Eigen"

#include "navigation.h"
#include "structs.h"
#include "par_io.h"
//#include "variables.h"


//calculate source shape
//std::complex<double> source_shape(Eigen::MatrixXcd& V, int r, int dir);
void source_shape_complete(const Eigen::MatrixXcd& V, const int nb_ev, std::vector<shp>& result);
//void hopping3d(int iup[][3], int idown[][3]);

//Map sourceshape to correct coordinates
void afterburner( std::vector<shp>& results );

//get averages of psi over same radii
void avg_radii(const std::vector<shp>& results, std::vector<std::pair<double, double> >& avg_shp);

//Weight eigenvectors with square root of eigenvalues
void weight_eigenvectors(const std::vector<double>& evalues, const int nb_ev, Eigen::MatrixXcd& V);

#endif /* SOURCESHAPE_H_ */
