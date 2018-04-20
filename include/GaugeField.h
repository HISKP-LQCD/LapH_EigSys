/*! \file Gaugefield.h
    \brief Contains the class for gaugefield handling

    The Gauge fields of a range of timeslices are read into a vector of
    Multiarrays. Functions: - 
 */

#ifndef GAUGEFIELD_H_
#define GAUGEFIELD_H_
//!
//! \class GaugeField
//! \brief Class for handling the Gaugefield, displacement routines and Smearing
//! methods.
//!                                                                             
//! This Class reads in timeslices from Gaugefields in lime format into a
//! Eigen::Matrix Structure. Eigenvectors or Matrices can be displaced in a
//! symmetrized and an unsymmetrized way. Furthermore Gaugetransformations can
//! be conducted either of the Gaugefield itself or passed Eigensystems.
//!                                                                             
//! Smearings implemented so far are Hyp-Smearing, Ape-Smearing and Stout
//! smearing, acting directly on the read in gaugefield
//!
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

//#include "Eigen/Dense"
//#include "Eigen/Core"
//#include "Eigen/SparseCore"
#include <unsupported/Eigen/MatrixFunctions>
#include "boost/multi_array.hpp" 

//#include "lime.h"
//#include "config_utils.h"
#include "EigenVector.h"
#include "io_utils.h"
//#include "propagator_io.h"
//#include "quark.h"
#include "RandomVector.h"
#include "typedefs.h"


//! \typedef 2dim array for lookup tables
typedef boost::multi_array<int, 2> look;
//! \typedef 2dim array for gauge field matrices of one timeslice.
//! 
typedef boost::multi_array<Eigen::Matrix3cd, 2> array_3cd_d2_eigen;
class GaugeField {
  public:
    //! \brief Constructor

    //!t0 is first, tf is last timeslice to be read in
    GaugeField(const int Lt, const int Lx, const int Ly, const int Lz, 
               const std::string config_path, const size_t t0, 
               const size_t tf, const size_t ndir);
    ~GaugeField() {};
    //! \brief overloaded access operator
    //!
    //! returns the complex 3x3 link variable of timeslice t at volumeindex v in
    //! direction dir
    //!
    inline const Eigen::Matrix3cd& operator()(const size_t t, const size_t v,
                                              const size_t dir) const {
        return (tslices.at(t))[v][dir];
    };
    //! \brief read in gauge field
    //! 
    //! slice_i is the initial, slice_f the final timeindex of the lattice.
    //! The read in lime gaugefield is internally splitted into slice_f -
    //! slice_i 3d timeslices and mapped to the eigenformat.
    //!
    void read_gauge_field(const size_t config_i, const size_t slice_i,
                          const size_t slice_f);

    //! \brief Hypercubic Blocking in 3 dimensions
    //! 
    //! t is the timeslice index, alpha_1 the inner weight, alpha_2 the outer
    //! weight and iter the number of iterations the smearing is applied
    void smearing_hyp(const size_t t, const double alpha_1,
                      const double alpha_2, const size_t iter);
    //! \brief Stout Smearing of Gaugelinks
    //!
    //! t is the timeslice index, rho is the staple weight and iter the number
    //! of applications of the smearing
    void smearing_stout(const size_t t, const double rho, const size_t iter);

    //! \brief APE Smearing of one timeslice of Gaugelinks
    //! 
    //! t is the timeslice index, alpha_1 is the staple weight, iter is the
    //! number of applications
    void smearing_ape(const size_t t, const double alpha_1, const size_t iter);

    //! \brief Initialize lookup tables for 3d navigation through timeslices 
    //! 
    //! LX: Lattice extent in X direction, LY: Lattice extent in Y direction,
    //! LZ: Lattice extent in Z direction 
    void init(const size_t LX, const size_t LY, const size_t LZ);

    //Test functions for navigation
    int get_up(const int pos, const int dir);
    int get_dn(const int pos, const int dir);

    //! Summ displacement vectors for phasefactor
    //!
    //!
    Eigen::Vector3f summed_displacement(const DisplacementDirection);

    //! \brief Forward displacement acting to the right once in direction dir
    //! 
    //! v is the eigensystem at timeslice t
    //! t is the timeslice index
    //! dir is the direction of the derivative
    Eigen::MatrixXcd forward_uv(const Eigen::MatrixXcd& v, const size_t t,
                                                  const char dir,
                                                  const size_t verbose);

    //! \brief Backward displacement acting to the right once in direction dir
    //! 
    //! v is the eigensystem at timeslice t
    //! t is the timeslice index
    //! dir is the direction of the derivative
    Eigen::MatrixXcd backward_uv(const Eigen::MatrixXcd& v, const size_t t,
                                                  const char dir,
                                                  const size_t verbose);

    //! \brief Displace Matrix of eigenvectors by looping over displacement
    //! vectors
    //!
    //! v is the Eigensystem at a specific timeslice
    //! t is the timeslice to take into account
    //! disp is a vector of displacement pairs
    //! verbose is an integer controlling the debug level
    //! The eigenvectors of timeslice t are displaced using the vector of
    //! displacement pairs. The First entry of each pair specifies the direction
    //! of the right acting derivative, ">" meaning forward and "<" meaning backward
    Eigen::MatrixXcd displace_eigenvectors(const Eigen::MatrixXcd& v,
                                           const size_t t,
                                           const DisplacementDirection disp,
                                           const size_t verbose);
  
    //! \brief Returns displaced vector or matrix
    //! 
    //! v is the address of the Object to be displaced, t is the timeslice
    //! index, dir is one of 0,1 or 2 meaning x-,y- or z-direction respectively
    //! sym inidcates whether a symmetrized derivative should be used.
    //Eigen::MatrixXcd disp(const Eigen::MatrixXcd& v, const size_t t,
    //                      const size_t dir, bool sym);
    //! brief Shift eigenvectors about one step up or down in one direction
    //!
    //! v is the address of the Object to be shifted (one timeslice of an
    //! Eigensystem). step chooses if we want to shift up or down (+1 for up,-1
    //! for down), dir chooses the x,y or z direction 0,1 or 2, respectively.
    //Eigen::MatrixXcd shift(const Eigen::MatrixXcd& v, const size_t step,
    //                       const size_t dir);

    //! brief Returns displaced vector or matrix
    //
    //! v is the address of the Object to be displaced, t is the timeslice
    //! index, dir is one of 0,1 or 2 meaning x-,y- or z-direction respectively
    //! The derivative is taken in one operator $\bar{\Psi_{f}}
    //! \overleftrightarrow{D}\Psi_{f'}$.
    //Eigen::MatrixXcd symmetric_derivative(const Eigen::MatrixXcd& v,
    //                                      const size_t t, const size_t dir);
    //! brief subfunction for symmetric_derivative
    //!
    //! calculates product U_\mu(x) * V(x)
    //Eigen::MatrixXcd Umu_times_V(const Eigen::MatrixXcd& v,
    //                                              const size_t t,
    //                                              const size_t dir,
    //                                              const size_t verbose);
    //! brief subfunction for symmetric_derivative
    //!
    //! calculates product U_\mu(x) * V(x+\mu)
    //Eigen::MatrixXcd Umu_times_shiftedV(const Eigen::MatrixXcd& v,       
    //                                              const size_t t,
    //                                              const size_t dir,
    //                                              const size_t verbose);

    //! \brief Returns symmetric 2 times displaced vector or matrix
    //! 
    //! v is the address of the Object to be displaced, t is the timeslice
    //! index, dir is one of 0,1 or 2 meaning x-,y- or z-direction respectively
    //Eigen::MatrixXcd disp_2(const Eigen::MatrixXcd& v, const size_t t,
    //                      const size_t dir);
    //! \brief Gauge Transformation of timeslices
    //!
    //! For generating the transformation fields indices of the initial
    //! timeslice t0 and the final timeslice tf need to be passed
    void trafo(const size_t t0, const size_t tf);

    //! \brief Returns a gaugetransformed transformed Eigenvector or LapH-Matrix
    //!
    //! If no transformation fields are generated, one is generated, otherwise
    //! first one is used
    Eigen::MatrixXcd trafo_ev(const Eigen::MatrixXcd& eig_sys);

    //! \brief Returns plaquette for one timeslice
    //!
    //! t is the timeslice index of the gaugefield
    double plaque_ts(const size_t t);

  private:
  
    //! \brief project one timeslice of lime Gaugefield to 3d timeslice
    //!
    //! t is the timeslice to be mapped timeslice points to the array of the
    //! lime Gaugefield.
    void map_timeslice_to_eigen(const size_t t, const double* timeslice);

    //! \brief tmlQCD function to read in gaugefields to array of doubles
    //!
    //! gaugefield is a pointer to the storage for the gaugefield
    //! filename indicates the path to the configuration
    //! slice_i is the initial, slice_f the final timeslice to be read
    void read_lime_gauge_field_doubleprec_timeslices(double* gaugefield,
                                                     const char* filename,
                                                     const size_t slice_i,
                                                     const size_t slice_f);

    //! \brief Build gauge array
    //!
    //! Constructs trange gaugefields, stored internally
    void build_gauge_array(const size_t trange);

    //! \brief Returns palquette at one space-time point
    //! 
    //! mu and nu are plaquette directions, vol is the volume index, t is the
    //! timeslice index
    //! Used internally by plaque_ts
    double plaque_pnt(const size_t mu, const size_t nu, const size_t vol, const size_t t);

  //! \brief Parameters specifying size of the Lattice
  const int Lt, Lx, Ly, Lz, V3, dim_row, V_TS, V_for_lime;
  //! \brief Path to gauge configurations
  const std::string config_path;
  //! \brief Vector holding boost multi_arrays for range of timeslices
  std::vector<array_3cd_d2_eigen> tslices;
  //! \brief One 3d timeslice of Gaugelinks as 2d boost multi_array of 3x3
  //! Matrices
  array_3cd_d2_eigen omega;
  //! \brief 2d boost_multiarray for indices in up direction
  look iup;
  
  //! \brief 2d boost_multiarray for indices in down direction
  look idown;
};

#endif //GAUGEFIELD_H_
