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
//smearing and associated stuff
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

std::array< int, 2 > get_dirs(int mu) {
  std::array<int,2> dirs;  
  if (mu == 0) {
    dirs.at(0) = 1;
    dirs.at(1) = 2;
  }
  else if (mu == 1) {
    dirs.at(0) = 0;
    dirs.at(1) = 2; 
  }
  else if (mu == 2) {
    dirs.at(0) = 0;
    dirs.at(1) = 1; 
  }
  
  return(dirs);
}

//storage map for dec_timeslice
static int decor (int dir, int smear_plane) {
  int ret;
  if (dir == 0) {
    if (smear_plane == 1) ret = 0;
    else ret = 3; 
  }
  else if (dir == 1) {
    if (smear_plane == 0) ret = 1;
    else ret = 4;
  }
  else {
    if (smear_plane == 0) ret = 2;
    else ret = 5;
  }
  return ret;
}

//Smearing schemes
//Stout-Smearing

void Tslice::smearing_stout(const int up_3d[][3], const int down_3d[][3], double rho, int iter) {

  int V3 = pars -> get_int("V3");
  std::complex<double> im_half(0,0.5);
  Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
  for ( auto i = 0; i < V3; ++i ) {
    eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
  }
  for (int j = 0; j < iter; ++j) {
    for (int i = 0; i < V3; ++i) {
      for (int dir = 0; dir < 3; ++dir) {
        Eigen::Matrix3cd staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link
        //filling each element of smearer using 3 links
        //For each link calculate staples summing them up
        for (int not_dir = 0; not_dir < 3; ++not_dir) {
          if (dir != not_dir) {
            int mu = up_3d[i][not_dir];
            int nu = up_3d[mu][dir];
            int eta = down_3d[nu][not_dir];
            //Staples in positive direction
            staple += eigen_timeslice[i][not_dir]*
              (eigen_timeslice[mu][dir]*(eigen_timeslice[eta][not_dir].adjoint()));

            mu = down_3d[i][not_dir];
            nu = up_3d[mu][dir];
            //Staples in negative direction
            staple += (eigen_timeslice[mu][not_dir].adjoint())*
              (eigen_timeslice[mu][dir]*eigen_timeslice[nu][dir]);
          }
        }
        Eigen::Matrix3cd omega = ( staple * ( rho/4. ) ) * ( eigen_timeslice[i][dir].adjoint() );
        Eigen::Matrix3cd q = ( ( omega.adjoint() - omega ) * 0.5 )
          - ( ( ( ( omega.adjoint() - omega ).trace() ) * Eigen::Matrix3cd::Identity() ) * (1./6.) );
        eigen_timeslice_ts[i][dir] = ( ( q*(-1) ).exp() ) * eigen_timeslice[i][dir];
      }
    }
    for ( auto i = 0; i < V3; ++i ) {
      for ( auto mu = 0; mu < 3; ++mu) {
        eigen_timeslice[i][mu] = eigen_timeslice_ts[i][mu];
      }
    }
  }
  //clean up
  for (int k = 0; k < V3; ++k) {
    delete[] eigen_timeslice_ts[k];
  }
  delete eigen_timeslice_ts;

}

//Ape-Smearing
void Tslice::smearing_ape(const int up_3d[][3], const int down_3d[][3], double alpha, int iter){

  int V3 = pars -> get_int("V3");
  Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
  for ( auto i = 0; i < V3; ++i ) {
    eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
  }
  for (int j = 0; j < iter; ++j) {
    for (int i = 0; i < V3; ++i) {
      for (int dir = 0; dir < 3; ++dir) {
        Eigen::Matrix3cd staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link
        //filling each element of smearer using 3 links
        //Position indices mu nu eta
        for (int not_dir = 0; not_dir < 3; ++not_dir) {
          if (dir != not_dir) {
            int mu = up_3d[i][not_dir];
            int nu = up_3d[mu][dir];
            int eta = down_3d[nu][not_dir];
            //Staples in positive direction
            staple += eigen_timeslice[i][not_dir]*
              (eigen_timeslice[mu][dir]*(eigen_timeslice[eta][not_dir].adjoint()));

            mu = down_3d[i][not_dir];
            nu = up_3d[mu][dir];
            //Staples in negative direction
            staple += (eigen_timeslice[mu][not_dir].adjoint())*
              (eigen_timeslice[mu][dir]*eigen_timeslice[nu][dir]);
          }
        }
        eigen_timeslice_ts[i][dir] = (eigen_timeslice[i][dir] * (1.-alpha)) + (staple * alpha/4.);
      }
    }
    for ( auto i = 0; i < V3; ++i ) {
      for ( auto mu = 0; mu < 3; ++mu) {
        eigen_timeslice[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
      }
    }
  }
  //clean up
  for (int k = 0; k < V3; ++k) {
    delete[] eigen_timeslice_ts[k];
  }
  delete eigen_timeslice_ts;
}

//HYP-Smearing

void Tslice::smearing_hyp(const int up_3d[][3], const int down_3d[][3], double alpha_1, double alpha_2, int iter) {
  
  int V3 = pars -> get_int("V3");
  //temporal timeslice twice the size for decorated links 
  Eigen::Matrix3cd **dec_timeslice = new Eigen::Matrix3cd *[V3];
  for (auto vol = 0; vol < V3; ++vol) {
    dec_timeslice[vol] = new Eigen::Matrix3cd[6];
  }

  //temporal timeslice from decorated links
  Eigen::Matrix3cd **eigen_timeslice_ts = new Eigen::Matrix3cd *[V3]; 
  for ( auto i = 0; i < V3; ++i ) {
    eigen_timeslice_ts[i] = new Eigen::Matrix3cd[3];
  }
  
  
  //temporal integers holding directions for decorated smearing
  int mu, nu, eta;
  for (int run = 0; run < iter; ++run) {
    //calculate inner staple from original timeslice, store in dec_timeslice each link can get smeared in two planes
    for (auto vol = 0; vol < V3; ++vol) {
      for (auto dir = 0; dir < 3; ++dir) {
        //inner staple
        Eigen::Matrix3cd inner_staple = Eigen::Matrix3cd::Zero();
        std::array< Eigen::Matrix3cd, 2 > tmp_staples;
        std::array<int, 2> perpendics = get_dirs( dir );
        for (auto it_perp_dir = perpendics.begin(); it_perp_dir != perpendics.end(); ++it_perp_dir ) {
          int perp_dir = *it_perp_dir;
          //up-type smearing_indices
          mu = up_3d[vol][perp_dir];
          nu = up_3d[mu][dir];
          eta = up_3d[vol][dir];

          //product of up matrices
          inner_staple = eigen_timeslice[vol][perp_dir] * ( eigen_timeslice[mu][dir] *
                        ( eigen_timeslice[eta][perp_dir].adjoint() ) ); 
          //down-type smearing indices
          mu = down_3d[vol][perp_dir];
          nu = up_3d[mu][dir];

          //eta is same endpoint no adjoint necessary here
          //product of down matrices
          inner_staple += ( eigen_timeslice[mu][perp_dir].adjoint() ) *
                          ( eigen_timeslice[mu][dir] * eigen_timeslice[nu][perp_dir] );

          //Careful placement of decorated links in dec_timeslices:
          //dir=0 has placement in dec_dir = 0 (smeared in 1 plane)
          //                       dec_dir = 3 (smeared in 2 plane)
          //dir=1 has placement in dec_dir = 1 (smeared in 0 plane)
          //                       dec_dir = 4 (smeared in 2 plane)
          //dir=2 has placement in dec_dir = 2 (smeared in 0 plane)
          //                       dec_dir = 5 (smeared in 1 plane)
          Eigen::Matrix3cd stac = ( eigen_timeslice[vol][dir] *
                              (1-alpha_2) ) +  ( inner_staple * alpha_2/2.);  
          int n_el = it_perp_dir - perpendics.begin();
          tmp_staples.at(n_el) = proj_to_su3_imp(stac);
          //tmp_staples.at(n_el) = stac;//without SU(3) projection
        }

        //staple link in direction dir in non participating and negative directions
        dec_timeslice[vol][dir] = tmp_staples.at(0);
        dec_timeslice[vol][dir+3] = tmp_staples.at(1);
      }
    }
    //just debugging
    //if(run == 0) std::cout << "Fat link a x,y,z (11,19,29): \n" << dec_timeslice[11*L2*L1+19*L1+29][1]<< "\n\n";
    //calculate outer staple from dec_timeslice as modified ape-smearing

      for (int i = 0; i < V3; ++i) {
        for (int dir = 0; dir < 3; ++dir) {
          Eigen::Matrix3cd outer_staple = Eigen::Matrix3cd::Zero(); //Holding all smearing matrices for one link

          //filling each element of smearer using 3 links
          //debugging link
          Eigen::Matrix3cd outer_staple_test;
          for (int not_dir = 0; not_dir < 3; ++not_dir) {
            if (dir != not_dir) {

              //calculate plane in which was smeared
              int plane = ( (dir+1) ^ (not_dir+1) ) - 1;
              mu = up_3d[i][not_dir];
              nu = up_3d[mu][dir];
              eta = down_3d[nu][not_dir];

              //Staples in positive direction
              //replace directions by appropriate decor 
              int a,b;
              a = decor(not_dir,plane);
              b = decor(dir,plane);
              outer_staple += dec_timeslice[i][a]*
                (dec_timeslice[mu][b]*(dec_timeslice[eta][a].adjoint()));
/*           if(i == (11*L2*L3+19*L3+29)&& run == 0)   std::cout << "dir, i: " << dir << " " << i << "\n" << dec_timeslice[i][a]*
                (dec_timeslice[mu][b]*(dec_timeslice[eta][a].adjoint())) << "\n\n";
                */
              mu = down_3d[i][not_dir];
              nu = up_3d[mu][dir];
            
              //Staples in negative direction
              outer_staple += (dec_timeslice[mu][a].adjoint())*
                (dec_timeslice[mu][b]*dec_timeslice[nu][a]); //structure has to be a_dag, b, a
  /*         if(i == (11*L2*L3+19*L3+29)&& run == 0)   std::cout << "dir, i: " << dir << " " << i << "\n" <<(dec_timeslice[mu][a].adjoint())*
                (dec_timeslice[mu][b]*dec_timeslice[nu][a]) << "\n\n";
                */
 
            }
          }
          eigen_timeslice_ts[i][dir] = (eigen_timeslice[i][dir] * (1.-alpha_1)) + (outer_staple * alpha_1/4.);
        }
      }
      for ( auto i = 0; i < V3; ++i ) {
        for ( auto mu = 0; mu < 3; ++mu) {
          eigen_timeslice[i][mu] = proj_to_su3_imp(eigen_timeslice_ts[i][mu]);
          //eigen_timeslice[i][mu] = eigen_timeslice_ts[i][mu];//without SU(3)-projection
        }
      }
  }
  //clean up
  for (int k = 0; k < 3; ++k) {
    delete[] dec_timeslice[k];
    delete[] dec_timeslice[k+3];
    delete[] eigen_timeslice_ts[k];
  } 
  delete dec_timeslice;
  delete eigen_timeslice_ts;
  //printf("Timeslice successfully HYP-smeared with (a_1, a_2, iterations): %f, %f, %d \n",
      //alpha_1, alpha_2, iter);
}
