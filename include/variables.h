  /*
   * variables.h
   * Holds needed lattice info : 
   * 	- L0...L3: euclidean dimensions
 * 	- V3: 3d-Volume-size
 *	- NDIR: # directions
 *	- NCOL: # colors
 *	- V_TS: # entries in ildg-timeslice
 *	- MAT_ENTRIES: # rows/cols for PETSc-Matrix
 * Created on: Aug 26, 2013
 * Author: christopher helmes
 */

#ifndef VARIABLES_H_
#define VARIABLES_H_
//Geometry (set up for equal spatial extent)
const int L0 = 4; //nt
const int L1 = 4; //nx
const int L2 = L1; //ny
const int L3 = L1; //nz
const int V3 = L1*L2*L3;
const int NDIR = 4;
const int NCOL = 3;
const int V_TS = L1*L2*L3*NDIR*NCOL*NCOL*2; //2 is for complex
const int V_4_LIME = V_TS * L0;
const int MAT_ENTRIES = NCOL*V3;
//desired number of eigenvalues
const int NEV = 120;

//IO-Paths
//path to gauge fields
const std::string GAUGE_FIELDS = "/hiskp2/gauges/test4x4x4x4";
//Chebyshev-acceleration
const double LAM_L = 11.8;
const double LAM_C = 6.0;

//Hyp-Smearing Parameters
const double ALPHA_1 = 0.62;
const double ALPHA_2 = 0.62;
const int ITER = 3;

//random seed
const unsigned int RND_SEED = 1227;

//Global Declarations
extern int up_3d[V3][3], down_3d[V3][3];
//Time Slice of Configuration
extern Eigen::Matrix3cd **eigen_timeslice;

#endif /* VARIABLES_H_ */
