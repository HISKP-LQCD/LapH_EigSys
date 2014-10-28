/*
 * shell_matop.h
 *
 *  Created on: Aug 26, 2013
 *      Author: helmes
 */

#ifndef SHELL_MATOP_H_
#define SHELL_MATOP_H_

#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <slepceps.h>
#include <petscblaslapack.h>
#include <Eigen/Eigen>
//#include <petsctime.h>
#include "config_utils.h"
#include "navigation.h"
#include "timeslice.h"
#include "par_io.h"
PetscErrorCode MatMult_Laplacian2D( Mat A,Vec x,Vec y); //replaces MATOP_MULT of PETSc
PetscErrorCode MatGetDiagonal_Laplacian2D(Mat A,Vec diag);//replaces MATOP_DIAF of PETSc

#endif /* SHELL_MATOP_H_ */
