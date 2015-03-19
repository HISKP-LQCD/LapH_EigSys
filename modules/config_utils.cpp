/*
 * configs.cpp
 * Source code to config navigation and config mappings
 * Defines: - lookup-tables in 3d (hopping3d)
 *					- mapping from ildg to Eigen Array (map_timeslice_to_eigen)
 *					- building of Laplacian in E4- and C-space (BuildLaplacian)
 *					- gauge transform (transform_ts)
 *					- check for gauge invariance (check_gauge)
 * Created on: Aug 26, 2013
 * Author: christopher helmes
 */
#include "config_utils.h"

static IO* const pars = IO::getInstance();
static Nav* const lookup = Nav::getInstance();

//Maps SU(3)-Matrices of ts to sparse Laplacian matrix Lap
PetscErrorCode BuildLaplacian(Mat Lap, Eigen::Matrix3cd **ts) {
  int V3 = pars -> get_int("V3");
  PetscErrorCode error;
  //# of transferred elements
  size_t n_el = 0;
  PetscScalar el_diag = 6;
  for (int i = 0; i < V3; ++i) {
    //k corresponds to index of SU(3)-matrix
    register const size_t k = 3*i;
    //diagonal values
    error = MatSetValue(Lap, k+0, k+0, el_diag, INSERT_VALUES);
    CHKERRQ(error);
    error = MatSetValue(Lap, k+1, k+1, el_diag, INSERT_VALUES);
    CHKERRQ(error);
    error = MatSetValue(Lap, k+2, k+2, el_diag, INSERT_VALUES);
    CHKERRQ(error);
    n_el += 3;
    for ( int mu = 0; mu < 3; ++mu ) {//direction of matrix
      register const int up = lookup -> get_up(i,mu);
      register const int down = lookup -> get_dn(i,mu);
      for ( int a = 0; a < 3; ++a ) {//colour a
        register const size_t k_prime = k+a;
        for ( int b = 0; b < 3; ++b ) {//colour b
          //if ( (3*up+b) > k_prime ) {
          error = MatSetValue(Lap, k_prime, 3*up+b, -(ts[i][mu])(a,b),
              INSERT_VALUES);
          CHKERRQ(error);
          //}
          //if ( (3*down+b) > k_prime ) {
          error = MatSetValue(Lap, k_prime, 3*down+b,
              -(ts[down][mu].adjoint())(a,b), INSERT_VALUES);
          CHKERRQ(error);
          //}
          n_el += 2;
        }
      }	
    }
  }
  std::cout << n_el << std::endl;
  return(error);
}

static Eigen::Matrix3cd proj_to_su3_imp(Eigen::Matrix3cd& in){
  //avoid possible aliasing issues:
  Eigen::Matrix3cd out = ( (in.adjoint() * in).sqrt() ).inverse();
  in = in*out;
  std::complex<double> det = 1./pow(in.determinant(),1./3.);
  in *= det;
  return in;
  
  /*
  in *= ( ( (in.adjoint()) * in ).sqrt() ).inverse();
  std::complex<double> det = 1./pow(in.determinant(),1./3.);
  in *= det;
  return in;
  */
}

//std::array< int, 2 > get_dirs(int mu) {
//  std::array<int,2> dirs;  
//  if (mu == 0) {
//    dirs.at(0) = 1;
//    dirs.at(1) = 2;
//  }
//  else if (mu == 1) {
//    dirs.at(0) = 0;
//    dirs.at(1) = 2; 
//  }
//  else if (mu == 2) {
//    dirs.at(0) = 0;
//    dirs.at(1) = 1; 
//  }
//  
//  return(dirs);
//}
//
////storage map for dec_timeslice
//static int decor (int dir, int smear_plane) {
//  int ret;
//  if (dir == 0) {
//    if (smear_plane == 1) ret = 0;
//    else ret = 3; 
//  }
//  else if (dir == 1) {
//    if (smear_plane == 0) ret = 1;
//    else ret = 4;
//  }
//  else {
//    if (smear_plane == 0) ret = 2;
//    else ret = 5;
//  }
//  return ret;
//}
//
////Smearing schemes
////Stout-Smearing
//
//void smearing_stout(const int up_3d[][3], const int down_3d[][3], Eigen::Matrix3cd **eigen_timeslice, double rho, int iter) {
//
//  int V3 = pars -> get_int("V3");
//  std::complex<double> im_half(0,0.5);
//  Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
//  for ( auto i = 0; i < V3; ++i ) {
//    eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
//  }
//  for (int j = 0; j < iter; ++j) {
//    for (int i = 0; i < V3; ++i) {
//      for (int dir = 0; dir < 3; ++dir) {
//        Eigen::Matrix3cd staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link
//        //filling each element of smearer using 3 links
//        //For each link calculate staples summing them up
//        for (int not_dir = 0; not_dir < 3; ++not_dir) {
//          if (dir != not_dir) {
//            int mu = up_3d[i][not_dir];
//            int nu = up_3d[mu][dir];
//            int eta = down_3d[nu][not_dir];
//            //Staples in positive direction
//            staple += eigen_timeslice[i][not_dir]*
//              (eigen_timeslice[mu][dir]*(eigen_timeslice[eta][not_dir].adjoint()));
//
//            mu = down_3d[i][not_dir];
//            nu = up_3d[mu][dir];
//            //Staples in negative direction
//            staple += (eigen_timeslice[mu][not_dir].adjoint())*
//              (eigen_timeslice[mu][dir]*eigen_timeslice[nu][dir]);
//          }
//        }
//        Eigen::Matrix3cd omega = ( staple * ( rho/4. ) ) * ( eigen_timeslice[i][dir].adjoint() );
//        Eigen::Matrix3cd q = ( ( omega.adjoint() - omega ) * 0.5 )
//          - ( ( ( ( omega.adjoint() - omega ).trace() ) * Eigen::Matrix3cd::Identity() ) * (1./6.) );
//        eigen_timeslice_ts[i][dir] = ( ( q*(-1) ).exp() ) * eigen_timeslice[i][dir];
//      }
//    }
//    for ( auto i = 0; i < V3; ++i ) {
//      for ( auto mu = 0; mu < 3; ++mu) {
//        eigen_timeslice[i][mu] = eigen_timeslice_ts[i][mu];
//      }
//    }
//  }
//  //clean up
//  for (int k = 0; k < V3; ++k) {
//    delete[] eigen_timeslice_ts[k];
//  }
//  delete eigen_timeslice_ts;
//
//}
//
////Ape-Smearing
//void smearing_ape(const int up_3d[][3], const int down_3d[][3], Eigen::Matrix3cd **eigen_timeslice, double alpha, int iter){
//
//  int V3 = pars -> get_int("V3");
//  Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
//  for ( auto i = 0; i < V3; ++i ) {
//    eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
//  }
//  for (int j = 0; j < iter; ++j) {
//    for (int i = 0; i < V3; ++i) {
//      for (int dir = 0; dir < 3; ++dir) {
//        Eigen::Matrix3cd staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link
//        //filling each element of smearer using 3 links
//        //Position indices mu nu eta
//        for (int not_dir = 0; not_dir < 3; ++not_dir) {
//          if (dir != not_dir) {
//            int mu = up_3d[i][not_dir];
//            int nu = up_3d[mu][dir];
//            int eta = down_3d[nu][not_dir];
//            //Staples in positive direction
//            staple += eigen_timeslice[i][not_dir]*
//              (eigen_timeslice[mu][dir]*(eigen_timeslice[eta][not_dir].adjoint()));
//
//            mu = down_3d[i][not_dir];
//            nu = up_3d[mu][dir];
//            //Staples in negative direction
//            staple += (eigen_timeslice[mu][not_dir].adjoint())*
//              (eigen_timeslice[mu][dir]*eigen_timeslice[nu][dir]);
//          }
//        }
//        eigen_timeslice_ts[i][dir] = (eigen_timeslice[i][dir] * (1.-alpha)) + (staple * alpha/4.);
//      }
//    }
//    for ( auto i = 0; i < V3; ++i ) {
//      for ( auto mu = 0; mu < 3; ++mu) {
//        eigen_timeslice[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
//      }
//    }
//  }
//  //clean up
//  for (int k = 0; k < V3; ++k) {
//    delete[] eigen_timeslice_ts[k];
//  }
//  delete eigen_timeslice_ts;
//}
//
////HYP-Smearing
//
//void smearing_hyp(const int up_3d[][3], const int down_3d[][3], Eigen::Matrix3cd **eigen_timeslice, double alpha_1, double alpha_2, int iter) {
//  
//  int V3 = pars -> get_int("V3");
//  //temporal timeslice twice the size for decorated links 
//  Eigen::Matrix3cd **dec_timeslice = new Eigen::Matrix3cd *[V3];
//  for (auto vol = 0; vol < V3; ++vol) {
//    dec_timeslice[vol] = new Eigen::Matrix3cd[6];
//  }
//
//  //temporal timeslice from decorated links
//  Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
//  for ( auto i = 0; i < V3; ++i ) {
//    eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
//  }
//  
//  
//  //temporal integers holding directions for decorated smearing
//  int mu, nu, eta;
//  for (int run = 0; run < iter; ++run) {
//    //calculate inner staple from original timeslice, store in dec_timeslice each link can get smeared in two planes
//    for (auto vol = 0; vol < V3; ++vol) {
//      for (auto dir = 0; dir < 3; ++dir) {
//        //inner staple
//        Eigen::Matrix3cd inner_staple = Eigen::Matrix3cd::Zero();
//        std::array< Eigen::Matrix3cd, 2 > tmp_staples;
//        std::array<int, 2> perpendics = get_dirs( dir );
//        for (auto it_perp_dir = perpendics.begin(); it_perp_dir != perpendics.end(); ++it_perp_dir ) {
//          int perp_dir = *it_perp_dir;
//          //up-type smearing_indices
//          mu = up_3d[vol][perp_dir];
//          nu = up_3d[mu][dir];
//          eta = up_3d[vol][dir];
//
//          //product of up matrices
//          inner_staple = eigen_timeslice[vol][perp_dir] * ( eigen_timeslice[mu][dir] *
//                        ( eigen_timeslice[eta][perp_dir].adjoint() ) ); 
//          //down-type smearing indices
//          mu = down_3d[vol][perp_dir];
//          nu = up_3d[mu][dir];
//
//          //eta is same endpoint no adjoint necessary here
//          //product of down matrices
//          inner_staple += ( eigen_timeslice[mu][perp_dir].adjoint() ) *
//                          ( eigen_timeslice[mu][dir] * eigen_timeslice[nu][perp_dir] );
//
//          //Careful placement of decorated links in dec_timeslices:
//          //dir=0 has placement in dec_dir = 0 (smeared in 1 plane)
//          //                       dec_dir = 3 (smeared in 2 plane)
//          //dir=1 has placement in dec_dir = 1 (smeared in 0 plane)
//          //                       dec_dir = 4 (smeared in 2 plane)
//          //dir=2 has placement in dec_dir = 2 (smeared in 0 plane)
//          //                       dec_dir = 5 (smeared in 1 plane)
//          Eigen::Matrix3cd stac = ( eigen_timeslice[vol][dir] *
//                              (1-alpha_2) ) +  ( inner_staple * alpha_2/2.);  
//          int n_el = it_perp_dir - perpendics.begin();
//          tmp_staples.at(n_el) = proj_to_su3_imp(stac);
//          //tmp_staples.at(n_el) = stac;//without SU(3) projection
//        }
//
//        //staple link in direction dir in non participating and negative directions
//        dec_timeslice[vol][dir] = tmp_staples.at(0);
//        dec_timeslice[vol][dir+3] = tmp_staples.at(1);
//      }
//    }
//    //just debugging
//    //if(run == 0) std::cout << "Fat link a x,y,z (11,19,29): \n" << dec_timeslice[11*L2*L1+19*L1+29][1]<< "\n\n";
//    //calculate outer staple from dec_timeslice as modified ape-smearing
//
//      for (int i = 0; i < V3; ++i) {
//        for (int dir = 0; dir < 3; ++dir) {
//          Eigen::Matrix3cd outer_staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link
//
//          //filling each element of smearer using 3 links
//          //debugging link
//          Eigen::Matrix3cd outer_staple_test;
//          for (int not_dir = 0; not_dir < 3; ++not_dir) {
//            if (dir != not_dir) {
//
//              //calculate plane in which was smeared
//              int plane = ( (dir+1) ^ (not_dir+1) ) - 1;
//              mu = up_3d[i][not_dir];
//              nu = up_3d[mu][dir];
//              eta = down_3d[nu][not_dir];
//
//              //Staples in positive direction
//              //replace directions by appropriate decor 
//              int a,b;
//              a = decor(not_dir,plane);
//              b = decor(dir,plane);
//              outer_staple += dec_timeslice[i][a]*
//                (dec_timeslice[mu][b]*(dec_timeslice[eta][a].adjoint()));
///*           if(i == (11*L2*L3+19*L3+29)&& run == 0)   std::cout << "dir, i: " << dir << " " << i << "\n" << dec_timeslice[i][a]*
//                (dec_timeslice[mu][b]*(dec_timeslice[eta][a].adjoint())) << "\n\n";
//                */
//              mu = down_3d[i][not_dir];
//              nu = up_3d[mu][dir];
//            
//              //Staples in negative direction
//              outer_staple += (dec_timeslice[mu][a].adjoint())*
//                (dec_timeslice[mu][b]*dec_timeslice[nu][a]); //structure has to be a_dag, b, a
//  /*         if(i == (11*L2*L3+19*L3+29)&& run == 0)   std::cout << "dir, i: " << dir << " " << i << "\n" <<(dec_timeslice[mu][a].adjoint())*
//                (dec_timeslice[mu][b]*dec_timeslice[nu][a]) << "\n\n";
//                */
// 
//            }
//          }
//          eigen_timeslice_ts[i][dir] = (eigen_timeslice[i][dir] * (1.-alpha_1)) + (outer_staple * alpha_1/4.);
//        }
//      }
//      for ( auto i = 0; i < V3; ++i ) {
//        for ( auto mu = 0; mu < 3; ++mu) {
//          eigen_timeslice[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
//          //eigen_timeslice[i][mu] = eigen_timeslice_ts[i][mu];//without SU(3)-projection
//        }
//      }
//  }
//  //clean up
//  for (int k = 0; k < 3; ++k) {
//    delete[] dec_timeslice[k];
//    delete[] dec_timeslice[k+3];
//    delete[] eigen_timeslice_ts[k];
//  } 
//  delete dec_timeslice;
//  delete eigen_timeslice_ts;
//  //printf("Timeslice successfully HYP-smeared with (a_1, a_2, iterations): %f, %f, %d \n",
//      //alpha_1, alpha_2, iter);
//}

//Displacements
//TODO: replace lookup tables
//displacement in one direction i acting to the right
void right_displacement_one_dir(Eigen::Matrix3cd** config, const int iup[][3],
    const int idown[][3], const int dir, Eigen::MatrixXcd& V, Eigen::MatrixXcd& W ) {

  int V3 = pars -> get_int("V3");
  //Information on Matrix size
  const int num_cols = V.cols();
  const int num_rows = V.rows();

  //Loop over all eigenvectors in 
  for (int i = 0; i < num_cols; ++i ) {
//    std::cout << "eigenvector: " << i << std::endl;
    //storing eigenvector
    Eigen::VectorXcd in(num_rows);
    Eigen::VectorXcd out(num_rows);

    in = V.col(i);

    //Displace eigenvector
    for (int spatial_ind = 0; spatial_ind < V3; ++spatial_ind) {
      //std::cout << "x: " << spatial_ind << std::endl;
      Eigen::Vector3cd tmp;
      Eigen::Vector3cd quark_up;
      Eigen::Vector3cd quark_down;

      //determine needed indices from lookup tables;
      int up_ind = iup[spatial_ind][dir];
      int down_ind = idown[spatial_ind][dir];

      quark_up = in.segment(3*up_ind,3);
      quark_down = in.segment(3*down_ind,3);
      tmp = 0.5*( (config[spatial_ind][dir] * quark_up) - ((config[down_ind][dir].adjoint()) * quark_down) ); 
      out.segment(3*spatial_ind,3) = tmp;

    }//end spatial loop
    //write displaced eigenvector to W
        W.col(i) = out;
  }//end eigenvector loop

}

//debugging functions

//Plaquette of timeslice
//calculate plaquettes at given point i
double get_plaque_nx_ny(int x, int y, int &i, int iup[][3], int idown[][3],
                        Eigen::Matrix3cd **U ) { 

        Eigen::Matrix3cd plaque;
	unsigned j,k,l;
	double P;
	j=iup[i][x];  //convention(up): 0:= n0+1, 1:= n1+1, 2:= n2+1, 3:= n3+1
	k=iup[j][y]; //convention(down): 0:= n0-1, 1:= n1-1, 2:= n2-1, 3:= n3-1
	l=idown[k][x];
	Eigen::Matrix3cd a,b,c,d;
	a=U[i][x]; //convention links: 0:= n0+1, 1:= n1+1, 2:= n2+1, 3:= n3+1
	b=U[j][y];
	c=U[l][x].adjoint();
	d=U[i][y].adjoint();
	plaque = a*(b*(c*d));
	P=1./3.*real(plaque.trace());

	return P;
}

double plaque_timeslice(Eigen::Matrix3cd **ts, int iup[][3] , int idown[][3]){

  int V3 = pars -> get_int("V3");

  double plaquette = 0;
  double cnt = 0;
  for (int i = 0; i < V3; ++i) {
    for (int dir_1 = 0; dir_1 < 3; ++dir_1) {
      for (int dir_2 = 0; dir_2 < 3; ++dir_2) {
        //make sure not to run in same direction
        if (dir_1 != dir_2) {
          plaquette += get_plaque_nx_ny(dir_1, dir_2, i, iup, idown, ts);
          ++cnt;
        }
      }//loop dir_2
    }//loop_dir1
  }//loop V3
  return (plaquette/cnt);
}

void check_gauge(Eigen::Matrix3cd **eigen_timeslice,
    Eigen::Matrix3cd **eigen_ts_gauge) {
  int V3 = pars -> get_int("V3");
  std::complex< double > tr_ts, tr_gauge;
  for ( int j = 0; j < V3; ++j ) {
    for (int dir = 0; dir < 3; ++dir ) {
      tr_ts = (eigen_timeslice[j][dir]*eigen_timeslice[j][dir].adjoint()).sum();
      tr_gauge = (eigen_ts_gauge[j][dir]*eigen_ts_gauge[j][dir].adjoint()).sum();
      if (abs(tr_ts - tr_gauge) > 10e-13 ) {
        std::cout << "Error in Gauge Transformation at "<< j << " " << dir << std::endl;
        std::cout << "tr_ts = " << tr_ts << " tr_gauge = " << tr_gauge << std::endl;
        exit(1);
      }
      /*std::cout << "ts-matrix:\n" << eigen_timeslice[j][dir] << "\n\ngauged ts-matrix:\n" << eigen_ts_gauge[j][dir] << "\n\n";*/
    }
  }
}

//Gauge-transform config and store transformed in gauge_config
void transform_ts( Eigen::Matrix3cd **config,
    int up[][3], Eigen::Matrix3cd* Omega, Eigen::Matrix3cd **gauge_config) {

  int V3 = pars -> get_int("V3");
  /*
  //SU(3)-gauge field of same size as config
  Eigen::Matrix3cd* Omega = new Eigen::Matrix3cd [V3];
  build_gauge_array(Omega);
  write_gauge_matrices("ts_trafo_log.bin", Omega);
  */
  for ( int i = 0; i < V3; ++i ) {
    for ( int mu = 0; mu < 3; ++mu ) {
      int j = up[i][mu];
      gauge_config[i][mu] = (Omega[i].adjoint())*(config[i][mu]*Omega[j]);
    //std::cout << "gauge config: \n" << gauge_config[i][mu] << "\n\n";
     // std::cout << "config: \n" << config[i][mu] << "\n\n\n";
    }
  }
//  delete[] Omega;
  /*for ( auto i = 0; i < V3; ++i ) {
    for ( auto mu = 0; mu < 3; ++mu ) {
      config[i][mu] = gauge_config[i][mu];
    }
  }*/

}

//transform matrix of eigenvectors with gauge array
void transform_ev(const int nb_ev, Eigen::Matrix3cd* Omega, Eigen::MatrixXcd& V) {
  int V3 = pars -> get_int("V3");
  /*
  //declare gauge trafo
  Eigen::Matrix3cd* Omega = new Eigen::Matrix3cd[V3];
  //initialize gauge trafo
  build_gauge_array(Omega);
  */
  //write_gauge_matrices("ev_trafo_log.bin",Omega);
  for (int nev = 0; nev < nb_ev ; ++nev) {
    for (int i = 0; i < V3; ++i) {
      int ind_c = 3 * i;
      Eigen::Vector3cd tmp = Omega[i].adjoint() * (V.col(nev)).segment(ind_c,3); 
      (V.col(nev)).segment(ind_c,3) = tmp; 
    }//end V3
  }//end nb_ev
  //delete[] Omega;
}  

//build gauge array
void build_gauge_array(Eigen::Matrix3cd* gauge) {
  srand(1227);
  int V3 = pars -> get_int("V3");
  //for (int vol = 0; vol < V3; ++vol) gauge[vol] = Eigen::Matrix3cd::Identity();
  for (int vol = 0; vol < V3; ++vol) construct_random_su3(gauge[vol]);

}

//checks if su is SU(3) matrix by tr(x*x^dagger) and sum(x*x^dagger)
static void check_su3(const Eigen::Matrix3cd& su) {

  std::complex<double> test1 =  (su * su.adjoint() ).trace();
  std::complex<double> test2 =  (su * su.adjoint() ).sum();
  
  if ( ( (abs(test1) - 3.) > 10e-11 ) || ( (abs(test2) - 3.) > 10e-11 ) ) {
    std::cout << std::setprecision(13) << " det:  " << su.determinant() << "\n tr: " << test1 
      << "\n sum: " << test2 << "\n" << std::endl;
  }
}

//Build a SU(3)-matrix from a random complex matrix and a back-projection
void construct_random_su3(Eigen::Matrix3cd& x) {
  
  std::array<std::complex<double>, 9 > random;
  //9 complex random numbers
  for (int c = 0; c < random.size(); ++c) {
    random.at(c) = std::complex<double>( (double(rand()) / RAND_MAX ) - 0.5,(double(rand()) / RAND_MAX ) - 0.5 ); 
        
  }
  x << random.at(0), random.at(1), random.at(2),
      random.at(3), random.at(4), random.at(5),
      random.at(6), random.at(7), random.at(8);
  x = proj_to_su3_imp(x);
  check_su3(x);

}


//Build a SU(3)-matrix from random SU(2)-matrix
void construct_update_matrices( Eigen::Matrix3cd& x, double eps ) {
  Eigen::Matrix3cd r;
  Eigen::Matrix3cd s;
  Eigen::Matrix3cd t;
  std::vector <std::complex < double > > vector;
  //embedding SU(2) in SU(3)
  su_2_matrix( eps, vector );//vector characterising SU(2) matrices
  r << vector.at(0), vector.at(1), 0, vector.at(2), vector.at(3), 0, 0, 0, 1;

  su_2_matrix( eps, vector );//vector characterising SU(2) matrices
  s << vector.at(0), 0, vector.at(1), 0, 1, 0, vector.at(2), 0, vector.at(3);

  su_2_matrix( eps, vector );//vector characterising SU(2) matrices
  t << 1, 0, 0, 0, vector.at(0), vector.at(1), 0, vector.at(2), vector.at(3);

  //Building Update matrix X
  x = r*(s*t);
  if ( (abs(x.trace()) > 10e-12) ) std::cout << x << "\n" << x.determinant() << "\n" << x.trace() << std::endl;

}

//Build random SU(2)-matrix
void su_2_matrix( double eps, std::vector< std::complex< double > >& vec_temp_4d ) {
  vec_temp_4d.clear();
  unsigned i;
  double x_0 = sqrt( 1 - eps*eps );
  double random[3];
  for ( i = 0; i < 3; ++i ) {
    random[i] = (double(rand()) / RAND_MAX ) - 0.5;
  }
  double mod_r;
  mod_r = sqrt(random[0]*random[0]+random[1]*random[1]+random[2]*random[2]);
  std::complex< double > z_0(x_0, (eps/mod_r) * random[2]);
  vec_temp_4d.push_back(z_0);
  std::complex< double > z_1((eps/mod_r) * random[1], (eps/mod_r) * random[0]);
  vec_temp_4d.push_back(z_1);
  std::complex< double > z_2((-1)*(eps/mod_r) * random[1], (eps/mod_r) * random[0]);
  vec_temp_4d.push_back(z_2);
  vec_temp_4d.push_back(conj(z_0));
}


