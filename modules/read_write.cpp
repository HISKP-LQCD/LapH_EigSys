//everything to read and write from/to files

#include "read_write.h"

static IO* const pars = IO::getInstance();
/********************************Input from files*****************************/

//Reads in Eigenvectors from one Timeslice in binary format to V
void read_evectors_bin_ts(const char* prefix, const int config_i, const int t,
    const int nb_ev, Eigen::MatrixXcd& V) {
  int V3 = pars -> get_int("V3");
  const int dim_row = 3 * V3;
  //buffer for read in
  std::complex<double>* eigen_vec = new std::complex<double>[dim_row];
  //setting up file
  char filename[200];
  sprintf(filename, "%s.%04d.%03d", prefix, config_i, t);
  std::cout << "Reading file: " << filename << std::endl;
  std::ifstream infile(filename, std::ifstream::binary);
  for (int nev = 0; nev < nb_ev; ++nev) {
    infile.read( reinterpret_cast<char*> (eigen_vec), 2*dim_row*sizeof(double));
    V.col(nev) = Eigen::Map<Eigen::VectorXcd, 0 >(eigen_vec, dim_row);
    //    std::cout << eigen_vec[5] << std::endl;
  }
  //clean up
  delete[] eigen_vec;
  infile.close();
}

//Reads in eigenvalues from ascii file to std::vector
//Beware of arguments ordering in sprintf
void read_eigenvalues_ascii(const char* prefix, const int config_i, const int t,
    const int nb_ev, std::vector<double>& ev) {

  //setting filename 
  char filename[200];
  sprintf(filename, "%s.%03d.%04d", prefix, t, config_i);
  std::ifstream infile(filename);
  if (infile) {
    ev.resize(nb_ev);
    for (int nev = 0; nev < nb_ev; ++nev) {
      infile >> ev.at(nev);
      //std::cout << ev.at(nev) << std::endl;
    }
    infile.close();
  }
}
//Reads in Array of gauge-trafo matrices from binary file to Array of
//Eigen::3cd matrices
void read_gauge_matrices (const char* prefix, Eigen::Matrix3cd* G) {
  int V3 = pars -> get_int("V3");
  const int entries = 9;
  std::ifstream infile(prefix, std::ifstream::binary);
  std::complex<double>* matrix_su3 = new std::complex<double>[entries];
  //read from binary file named "prefix"
  for (auto vol = 0; vol < V3; ++vol ) {
    infile.read( reinterpret_cast<char*> (matrix_su3), 2*entries*sizeof(double) );
    G[vol] = Eigen::Map<Eigen::Matrix3cd> (matrix_su3);
  }//end vol

  //clean up
  delete[] matrix_su3;
  infile.close();

}

//Reads in eigenvalues from binary to std::vector
void read_eigenvalues_bin( const char* prefix, const int config_i, const int t,
    const int nb_ev, std::vector<double>& ev) {

  //Build filename
  char filename[200];
  sprintf(filename, "%s.%04d.%03d", prefix, config_i, t);
  std::ifstream infile(filename, std::ifstream::binary);
  if (infile) {
    ev.resize(nb_ev);
    infile.read (reinterpret_cast<char*>(&ev[0]), ev.size()*sizeof(double));
    //if(infile) std::cout << filename << " read successfully" << std::endl;
    //else std::cout << "error: only" << infile.gcount() << "could be read" << std::endl;
    infile.close();
  }  

}
/*
//Reads in sourceshape from binary
//Layout is x,y,z,r,psi(r)
void read_sourceshape_bin(const char* filename, std::vector<shp> sourceshape) {
std::ifstream source_in(filename, std::ifstream::binary);
if (source_in){
source_in.seekg (0, source_in.end);
int shp_lngth = source_in.tellg();
source_in.seekg (0, is.beg);
char* buffer = new char [shp_length];
source_in.read(buffer,shp_length);
int elmnt_length = 3*sizeof(int)+2*sizeof(double);
std::cout << shp_length/elmnt_length << std::endl;
sourceshape.resize(shp_length/elmnt_length);
for(auto i = sourceshape.begin(); i<=sourceshape.end(); ++i) {
sourceshape.at(i) = buffer[i+elmnt_length];
}
delete[] buffer;
}
}
*/
/****************************Output to files**********************************/

//Writes eigenvectors from one Timeslice to file in binary format

void write_eig_sys_bin(const char* prefix, const int config_i, const int t, const int nb_ev, Eigen::MatrixXcd& V) {
  const int V3 = pars -> get_int("V3");
  //set up filename
  char file [200];
  sprintf(file, "%s.%04d.%03d", prefix, config_i, t);
  std::cout << "Writing to file:" << file << std::endl;
  std::ofstream outfile(file, std::ofstream::binary);
  outfile.write(reinterpret_cast<char*> (V.data()), 2*3*V3*nb_ev*sizeof(double));
  outfile.close();

}


//Write Results for source shape to ASCII-file
void write_sourceshape_ascii(const char* prefix, const int config,
    const int tslice, const int nb_ev, const std::vector<std::pair<double,double> >& results) {

  char filename[200];
  sprintf(filename, "%s_nev%d.%04d.%03d.txt", prefix, nb_ev, config, tslice);
  std::ofstream write_file(filename);
  for (auto& element:results) {
    write_file << std::setprecision(12) << std::get<0>(element) << " " << std::get<1>(element) << std::endl;
  }
  write_file.close();

}

//Write Results for source shape to binary file
void write_sourceshape_bin(const char* prefix, const int config,
    const int tslice, const int nb_ev, const std::vector<std::pair<double,double> >& results) {

  //build filename
  char filename[200];
  sprintf(filename, "%s_nev%d.%04d.%03d.bin", prefix, nb_ev, config, tslice);
  std::ofstream eigenvalues(filename, std::ofstream::binary);
  eigenvalues.write(reinterpret_cast<const char*>(&results[0]),
      results.size()*sizeof(std::pair<double,double>));
  eigenvalues.close();
}

//Write eigenvalues from std::vector to binary
void write_eigenvalues_bin( const char* prefix, const int config_i, const int t,
    const int nb_ev, std::vector<double>& ev) {

  //Build filename
  char filename[200];
  sprintf(filename, "%s.%04d.%03d", prefix, config_i, t);

  std::ofstream outfile(filename, std::ofstream::binary);
  if(outfile) {
    outfile.write(reinterpret_cast<char*>(&ev[0]), ev.size()*sizeof(double));
    outfile.close();
    std::cout << "recovery written to: " << filename << std::endl;
  }
}

//write gauge trafo matrices to binary file
void write_gauge_matrices(const char* prefix, Eigen::Matrix3cd* G) {
  int V3 = pars -> get_int("V3");

  std::ofstream outfile (prefix, std::ofstream::binary);
  for (auto vol = 0; vol < V3; ++vol) {
    outfile.write(reinterpret_cast<const char*>(&(G[vol])), 2*9*sizeof(double));
  }//end vol
  outfile.close();
}
//write gauge link matrices of one timeslice to binary file
void write_link_matrices_ts(const char* prefix, Eigen::Matrix3cd** links) {

  int V3 = pars -> get_int("V3");
  std::ofstream outfile (prefix, std::ofstream::binary);
  for (auto vol = 0; vol < V3; ++vol) {
    for (auto mu = 0; mu < 3; ++mu) {
      outfile.write(reinterpret_cast<const char*>(&(links[vol][mu])), 2*9*sizeof(double));
    }//end dir
  }//end vol
  outfile.close();
}

