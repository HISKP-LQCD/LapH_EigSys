#include "par_io.h"
static IO* const pars=IO::getInstance();

Tslice* Tslice::getInstance(){
  static Tslice theInstance;
  return & theInstance;
}
//Constructor
Tslice::Tslice() {
  pars -> get_int("V3");
  Eigen::Matrix3cd **eigen_timeslice = new Eigen::Matrix3cd *[V3];
  //Allocate Eigen Array to hold timeslice
  for ( auto i = 0; i < V3; ++i ) {
    eigen_timeslice[i] = new Eigen::Matrix3cd[3];
    for (auto dir = 0; dir < 3; ++dir) {
      eigen_timeslice[i][dir] = Eigen::Matrix3cd::Identity();
    }
  }
}

//mapping from gauge config to Eigen 3x3 complex matrix arrays
void Tslice::map_timeslice_to_eigen( double *timeslice) {
  int L1 = pars -> get_int("LX");
  int L2 = pars -> get_int("LY");
  int L3 = pars -> get_int("LZ");
  
  int NDIR = pars -> get_int("NDIR");
  int NCOL = pars -> get_int("NCOL");
  int V3 = pars -> get_int("V3");

  int V_TS = pars -> get_int("V_TS");

  //read in elements
  int el_input = 0;
  for (int z = 0; z < L3; ++z) {//spatial loops
    for (int y = 0; y < L2; ++y) {
      for (int x = 0; x < L1; ++x) {
        for (int mu = 1; mu < 4; ++mu) {//direction loop
          std::complex< double > array[9];
          for (int a = 0; a < 3; ++a) {//colour loops
            for (int b = 0; b < 3; ++b) {
              //timeslice index of real part
              int ind_r = z*V_TS/L3+y*V_TS/(L3*L2)+x*V_TS/(V3)+
                mu*V_TS/(V3*NDIR)+a*V_TS/(V3*NDIR*NCOL)
                +b*V_TS/(V3*NDIR*NCOL*NCOL)+0;
              //timeslice index of imaginary part
              int ind_i = z*V_TS/L3+y*V_TS/(L3*L2)+x*V_TS/(V3)+
                mu*V_TS/(V3*NDIR)+a*V_TS/(V3*NDIR*NCOL)
                +b*V_TS/(V3*NDIR*NCOL*NCOL)+1;
              std::complex<double> pair(timeslice[ind_r], timeslice[ind_i]);
              //array to be mapped to Eigen Array
              array[3*b+a] = pair;
              ++el_input;
            }
          }
          Eigen::Map<Eigen::Matrix3cd> dummy(array);
          //spatial index
          int ind = z*L2*L1+y*L1+x;
          eigen_timeslice[ind][mu-1] = dummy;
        }
      }
    }
  }
  //std::cout << el_input << " doubles read in from ildg timeslice " << std::endl;
}

Eigen::Matrix3cd Tslice::get_gauge(const int spat, const int dir) {
  return eigen_timeslice[spat][dir];
}
