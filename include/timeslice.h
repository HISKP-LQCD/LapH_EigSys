//Quick and dirty handling of timeslice for global handling
//Declare as singleton for write protection and global access
#ifndef _TIMESLICE_H_
#define _TIMESLICE_H_
#include "Eigen/Eigen"
class Tslice {
  private:
    Tslice () {};
    ~Tslice () {};
  //Eigen Array
  Eigen::Matrix3cd **eigen_timeslice;
  public:

  static Tslice* getInstance();
  //Fill Eigenarray with one timeslice of eigenvalues
  
  //Get SU(3)-Matrices from timeslice and sort them into Eigen Array
  void map_timeslice_to_eigen(Eigen::Matrix3cd **eigen, double *timeslice);
  //obtain one entry of gaugefield timeslice
  Eigen::Matrix3cd get_gauge(const int spat, const int dir); 
};
#endif//_TIMESLICE_H_