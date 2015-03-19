//Quick and dirty handling of timeslice for global handling
//Declare as singleton for write protection and global access
#ifndef _TIMESLICE_H_
#define _TIMESLICE_H_
#include <array>
#include <complex>
#include "Eigen/Eigen"
#include <unsupported/Eigen/MatrixFunctions> 
#include "navigation.h"
#include "par_io.h"
class Tslice {
  private:
    Tslice (){};
    ~Tslice ();
    //Eigen Array
    Eigen::Matrix3cd** eigen_timeslice;
  public:

    static Tslice* getInstance();
    //initialize timeslice with identity matrices
    void init(); 
    //Get SU(3)-Matrices from timeslice and sort them into Eigen Array
    void map_timeslice_to_eigen( double* timeslice);
    void transform_ts(Eigen::Matrix3cd* Omega);
    //obtain one entry of gaugefield timeslice
    void set_gauge(const int spat, const int dir, Eigen::Matrix3cd el);
    Eigen::Matrix3cd get_gauge(const int spat, const int dir);
    //Hyp-smear the timeslice
    void smearing_hyp(const double a_1, const double a_2, const int it);
    void smearing_stout(double rho, int iter);
    void smearing_ape( double alpha, int iter);
};
#endif//_TIMESLICE_H_
