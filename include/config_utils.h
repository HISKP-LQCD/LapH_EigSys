/*
 * configs.h
 * Declares: - lookup-tables in 3d (hopping3d)
 *					 - mapping from ildg to Eigen Array (map_timeslice_to_eigen)
 *					 - building of Laplacian in E4- and C-space (BuildLaplacian)
 *					 - gauge transform (transform_ts)
 *					 - check for gauge invariance (check_gauge)
 * Created on: Aug 26, 2013
 * Author: christopher helmes
 */

#ifndef CONFIGS_H_
#define CONFIGS_H_
#include <complex>
#include <iomanip>
#include <iostream>
#include <array>
#include <slepceps.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>
#include <petsctime.h>

#include "read_write.h"
#include "variables.h"


//navigation through timeslice
void hopping3d(int iup[][3], int idown[][3]);

//write timeslice of ildg-message to array of eigen-matrices
void map_timeslice_to_eigen(Eigen::Matrix3cd **eigen, double *timeslice);

//Write Laplacian in euclidean and colour space (ts) to PETSc-Matrix Lap 
PetscErrorCode BuildLaplacian(Mat Lap, Eigen::Matrix3cd **ts, int iup[][3], int idown[][3]);
void smearing_stout(Eigen::Matrix3cd **eigen_timeslice, double rho, int iter);
void smearing_ape(Eigen::Matrix3cd **eigen_timeslice, double alpha, int iter);
void smearing_hyp(Eigen::Matrix3cd **eigen_timeslice, double alpha_1, double alpha_2, int iter);
std::array< int, 2 > get_dirs(int mu);

//Right displacement of one eigensystem
void right_displacement_one_dir(Eigen::Matrix3cd** config, const int iup[][3],
    const int idown[][3], const int dir, Eigen::MatrixXcd& V, Eigen::MatrixXcd& W );

//----------------------------------------------------------------------------//
//                                 Debugging                                  //
//----------------------------------------------------------------------------//
//Plaquette of timeslice
//calculate plaquettes at given point i
double get_plaque_nx_ny(int x, int y, int &i, int iup[][3], int idown[][3],
                        Eigen::Matrix3cd **U );
double plaque_timeslice(Eigen::Matrix3cd **ts, int iup[][3], int idown[][3]);

//checking gauge invariance of read in configuration
void check_gauge(Eigen::Matrix3cd **eigen_timeslice, Eigen::Matrix3cd **eigen_ts_gauge);
//transform configuration with SU(3)-Matrices
void transform_ts( Eigen::Matrix3cd **config,
    int up[][3], Eigen::Matrix3cd* Omega, Eigen::Matrix3cd **gauge_config);

//build gauge array
void build_gauge_array(Eigen::Matrix3cd* gauge);

//transform matrix of eigenvectors with gauge array
void transform_ev(const int nb_ev, Eigen::Matrix3cd* Omega, Eigen::MatrixXcd& V);

//helper functions for transformation
void construct_update_matrices( Eigen::Matrix3cd& x, double eps );

//Build a SU(3)-matrix from a random complex matrix and a back-projection
void construct_random_su3(Eigen::Matrix3cd& x);
void su_2_matrix( double eps, std::vector< std::complex< double > >& vec_temp_4d );

#endif /* CONFIGS_H_ */
