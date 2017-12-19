//calculates the source shape psi(r)

#include "sourceshape_funcs.h"
#include "config_utils.h"
//#include "variables.h"
static Nav* const lookup = Nav::getInstance();
static IO* const pars = IO::getInstance();

// get index from distance dst to point pnt in direction dir (x = 0, y = 1, z = 2)
static int get_dst_ind(const int pnt, const int dst, const int dir) {
  int coord = pnt;
  //for ( int cnt = 0; cnt < dst; ++cnt ) coord = up_3d[coord][dir];
  for ( int cnt = 0; cnt < dst; ++cnt ) coord = lookup -> get_up(coord, dir);
  return coord;

}

// get index for vector case of r (r1,r2,r3) Taking one spatial index and going
// ri steps in each direction
static int get_vectorial_dst_ind(const int pnt, const Eigen::Vector3i r) {

  int dst_tmp = get_dst_ind(pnt, r(0), 0);
  dst_tmp = get_dst_ind(dst_tmp, r(1), 1);
  dst_tmp = get_dst_ind(dst_tmp, r(2), 2);
  return dst_tmp;

}

//construct one element of Distillation operator, i.e. a 3x3 complex matrix in
//one direction for r 
static std::complex<double> dll_op(const Eigen::MatrixXcd& V,
    const int spatial_index,const int nb_ev, const Eigen::RowVector3i r) {

  //get spatial indices in every direction (x,y,z) at distance
  //Factors of three are for N_col
  int dst_ind = 3 * get_vectorial_dst_ind(spatial_index, r);
  int spatial_ind = 3 * spatial_index;

  //Every V.block() comprises a matrix of size  N_c x N_evectors either at spatial point
  //x or distant point x+r. Using cyclic invariance the computation can be sped
  //up significantly. Spatial and Distant index multiplied by three because of
  //Factor N_c in rows.
  Eigen::Matrix3cd tmp_1 = ( V.block(dst_ind,0,3,nb_ev) )* 
              ( ( V.block(spatial_ind,0,3,nb_ev) ).adjoint() );

 /* Eigen::Matrix3cd tmp_2 = ( V.block(spatial_ind,0,3,nb_ev) ) *
              ( ( V.block(dst_ind,0,3,nb_ev) ).adjoint() );
 */
 // std::cout << ( tmp_1.adjoint() ) - tmp_2 << "\n\n";
  return sqrt( (tmp_1 * ( tmp_1.adjoint() ) ).trace() );

}

//Calculates the source shape in dependance of an vector r
void  source_shape_complete(const Eigen::MatrixXcd& V, const int nb_ev,
    std::vector<shp>& results) {

  for (int r1 = 0; r1 < pars -> get_int("L1"); r1+=2) {
    for (int r2 = 0; r2 < pars -> get_int("L2"); r2+=2) {
      //keep in x-y-plane
      //int r3 = 0;
      for (int r3 = 0; r3 < pars -> get_int("L3"); r3+=2 ) {
        std::complex<double> psi (0.0,0.0);
        shp element;
      Eigen::RowVector3i distance (r1,r2,r3);
        for (int space = 0; space < pars -> get_int("V3"); ++space) {
          // creating the source shape at every spatial lattice point
          psi += dll_op(V, space, nb_ev, distance);
          //std::cout << "psi_tmp(r) =  " << psi  << std::endl;
        }
        //write results
        element.r = distance;
//        element.length = sqrt(distance.squaredNorm());
        element.shape = psi.real();
        results.push_back(element);
//        std::cout << element.shape << std::endl;
      }//end r3
    }//end r2
  }//end r1  

}

static double flip_sign(double value) {
  if(value < 0.) value *= (-1.);
  return value;
}

//Weight eigenvectors with square root of eigenvalues
void weight_eigenvectors(const std::vector<double>& evalues, const int nb_ev, Eigen::MatrixXcd& V){

  //Multiply every eigenvector with the squareroot of its real eigenvalue
  for (int idx = 0; idx < nb_ev; ++idx) {
    double ev = flip_sign(evalues.at(idx));
    V.col(idx) *= sqrt(ev);
  }
}


//get averages of psi over same radii
void avg_radii(const std::vector<shp>& results, std::vector<std::pair<double, double> >& avg_shp) {

  //construct vector of tuples. Each entry holds a tuple with cnt, r and convoluted psi
  std::vector<std::tuple<int,double,double> > tmp;
  //search whole vector of results
  for(auto i=results.begin(); i != results.end(); ++i) {
    shp res = *i;
    bool copied = false;
    //get radius from each element shp
    double radius = res.length;
    for(auto t=tmp.begin(); t != tmp.end(); ++t) {
      if (std::abs(radius - std::get<1>(*t)) < 10e-12) {
        ++std::get<0>(*t);
        std::get<2>(*t) += res.shape; 
        copied = true;
//        std::cout << std::get<0>(*t) << " " << std::get<1>(*t) << " " << std::get<2>(*t) << std::endl;
//        std::cout << tmp.size() << std::endl;
      }
    }
    
    if (!copied) {
      std::tuple<int,double,double> new_tmp;
      new_tmp = std::make_tuple(1,res.length, res.shape);
      tmp.push_back(new_tmp);
//      std::cout << tmp.size() << std::endl;
    }
  }
  
  //copy averages to avg_shp
  for(auto s:tmp) {
    //std::cout << std::get<0>(s) << " " << std::get<1>(s) << " " << std::get<2>(s) << std::endl; 
    double r_avg = std::get<1>(s);
    double psi_avg = std::get<2>(s)/std::get<0>(s);
    avg_shp.push_back(std::make_pair(r_avg,psi_avg));
  }
  
}

//Massage result to get source shape in all directions, proper radii
//If ri is between 12 and 23 subtract Li
void afterburner( std::vector<shp>& results ){

  //loop over results
  for (auto& element:results){
    //check ri 
    for (int indx = 0 ; indx < 3; ++indx) {
      //Here L1 is placeholder for any spatial extension
      if ( (element.r)(indx) >= (pars->get_int("L1"))/2 ){
        (element.r)(indx) -= pars->get_int("L1");
      }
    }
    element.length = sqrt((element.r).squaredNorm());
  }
}
