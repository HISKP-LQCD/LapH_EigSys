//everything to read and write from/to files

#include "read_write.h"
static IO* const pars = IO::getInstance();
static Tslice* const eigen_timeslice = Tslice::getInstance();

/******************************** Helper functions ****************************/
// check if all eigenvectors are read at the end of the file. Programm exits if
// not
static void eof_check(const int t, const int nev, const int nb_ev,
                      const bool file_end){
  if (file_end){
    std::cout << "Timeslice: " << t << ": " << nev 
              << " eigenvectors read from file" << std::endl;
    if (nev != nb_ev){
      std::cout << "Error: Wrong number of eigenvectors read, exiting" 
                << std::endl;
      exit(1);
    }
  }
}

static bool check_trace(const Eigen::MatrixXcd& V, const int nb_ev){
  bool read_state = true;
  Eigen::MatrixXcd VdV(nb_ev,nb_ev);
  VdV = V.adjoint() * V;
  double eps = 10e-10;
  std::complex<double> trace (0.,0.), sum(0.,0.);
  trace = VdV.trace();
  sum = VdV.sum();
  std::cout << trace.real() << std::endl;
  if ( fabs( trace.real()) - nb_ev > eps ||
       fabs(trace.imag()) > eps ||
       fabs(sum.real()) - nb_ev > eps ||
       fabs(sum.imag()) > eps)
    read_state = false;
  return read_state; 
}

/********************************Input from files*****************************/

//Reads in Eigenvectors from one Timeslice in binary format to V
void read_evectors_bin_ts(const std::string& prefix, const int config_i,
                          const int t, const int nb_ev, Eigen::MatrixXcd& V){
  int V3 = pars -> get_int("V3");
  //bool thorough = pars -> get_int("strict");
  const int dim_row = 3 * V3;
  //TODO: Change path getting to something keyword independent
  std::string path = pars -> get_path("in_path");
  //buffer for read in
  std::complex<double>* eigen_vec = new std::complex<double>[dim_row];
  //setting up file
  std::string filename = (boost::format("%s/%s.%04d.%03d")% path % prefix 
                          % config_i % t).str();
  std::cout << "Reading file: " << filename << std::endl;
  std::ifstream infile(filename, std::ifstream::binary);
  for (int nev = 0; nev < nb_ev; ++nev) {
    infile.read( reinterpret_cast<char*> (eigen_vec), 2*dim_row*sizeof(double));
    V.col(nev) = Eigen::Map<Eigen::VectorXcd, 0 >(eigen_vec, dim_row);
    eof_check(t,nev,nb_ev,infile.eof());
  }
  if(check_trace(V, nb_ev) != true){
    std::cout << "Timeslice: " << t 
              << ": Eigenvectors damaged, exiting" << std::endl;
    exit(0);
  }

  //clean up
  delete[] eigen_vec;
  infile.close();
}

//Reads in Eigenvectors from one Timeslice in binary format to V
void read_evectors_bin_ts(const std::string& path, const std::string& prefix,
                          const int config_i, const int t, const int nb_ev,
                          Eigen::MatrixXcd& V){
  int V3 = pars -> get_int("V3");
  //bool thorough = pars -> get_int("strict");
  const int dim_row = 3 * V3;
  //TODO: Change path getting to something keyword independent
  //buffer for read in
  std::complex<double>* eigen_vec = new std::complex<double>[dim_row];
  //setting up file
  std::string filename = (boost::format("%s/%s.%04d.%03d")% path % prefix
                          % config_i % t).str();
  std::cout << "Reading file: " << filename << std::endl;
  std::ifstream infile(filename, std::ifstream::binary);
  for (int nev = 0; nev < nb_ev; ++nev) {
    infile.read( reinterpret_cast<char*> (eigen_vec), 2*dim_row*sizeof(double));
    V.col(nev) = Eigen::Map<Eigen::VectorXcd, 0 >(eigen_vec, dim_row);
    eof_check(t,nev,nb_ev,infile.eof());
  }
  if(check_trace(V, nb_ev) != true){
    std::cout << "Timeslice: " << t 
              << ": Eigenvectors damaged, exiting" << std::endl;
    exit(0);
  }

  //clean up
  delete[] eigen_vec;
  infile.close();
}
//Reads in eigenvalues from ascii file to std::vector
//Beware of arguments ordering in sprintf
void read_eigenvalues_ascii(const std::string& path, const std::string& prefix,
                            const int config_i, const int t, const int nb_ev,
                            std::vector<double>& ev){

  //setting filename 
  std::string filename = (boost::format("%s/%s.%04d.%03d")% path % prefix
                          % config_i % t).str();
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
void read_gauge_matrices (const std::string& prefix, Eigen::Matrix3cd* G) {
  int V3 = pars -> get_int("V3");
  const int entries = 9;
  std::ifstream infile(prefix, std::ifstream::binary);
  std::complex<double>* matrix_su3 = new std::complex<double>[entries];
  //read from binary file named "prefix"
  for (int vol = 0; vol < V3; ++vol ) {
    infile.read( reinterpret_cast<char*> (matrix_su3), 2*entries*sizeof(double) );
    G[vol] = Eigen::Map<Eigen::Matrix3cd> (matrix_su3);
  }//end vol

  //clean up
  delete[] matrix_su3;
  infile.close();

}

//Reads in eigenvalues from binary to std::vector
void read_eigenvalues_bin(const std::string& path, const std::string& prefix,
                          const int config_i, const int t, const int nb_ev,
                          std::vector<double>& ev){
  //Build filename
  std::string filename = (boost::format("%s/%s.%04d.%03d")% path % prefix 
                          % config_i % t).str();
  std::cout<< filename <<std::endl;
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
void read_sourceshape_bin(const std::string& filename, std::vector<shp> sourceshape) {
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

void write_eig_sys_bin(const std::string& path, const std::string& prefix,
                       const int config_i, const int t, const int nb_ev,
                       Eigen::MatrixXcd& V){
  const int V3 = pars -> get_int("V3");
  //set up filename
  std::string filename = (boost::format("%s/%s.%04d.%03d")% path % prefix 
                          % config_i % t).str();
  if(check_trace(V, nb_ev) != true){
    std::cout << "Timeslice: " << t 
              << ": Eigenvectors damaged, abort writing" << std::endl;
    exit(1);
  }
  std::cout << "Writing to file:" << filename << std::endl;
  std::ofstream outfile(filename, std::ofstream::binary);
  std::streamsize  begin = outfile.tellp();
  std::streamsize eigsys_bytes =2*3*V3*nb_ev*sizeof(double); 
  outfile.write(reinterpret_cast<char*> (V.data()), eigsys_bytes);
  std::streamsize end = outfile.tellp();
  if ( (end - begin)/eigsys_bytes != 1 ){
    std::cout << "Timeslice:  " << t 
              << " Error: write incomplete, exiting" << std::endl;
    std::cout << (end-begin) << " bytes instead of expected "
              << eigsys_bytes << " bytes" << std::endl;
    exit(1);
  } 
  //std::cout << end - begin << " bytes written" << std::endl;
  outfile.close();

}


//Write Results for source shape to ASCII-file
void write_sourceshape_ascii(const std::string& path, const std::string& prefix,
                             const int config, const int tslice,
                             const int nb_ev,
                             const std::vector<std::pair<double,double> >& results){

  std::string filename = (boost::format("%s/%s_nev%d.%04d.%03d")% path % prefix
                          % nb_ev %config % tslice).str();
  std::ofstream write_file(filename);
  // A header is a brilliant idea
  write_file << "#r\tpsi(r)"<< std::endl;
  for (auto& element:results) {
    write_file << std::setprecision(12) << std::get<0>(element) 
               << " " << std::get<1>(element) << std::endl;
  }
  write_file.close();

}

//Write Results for source shape to binary file
void write_sourceshape_bin(const std::string& path,const std::string& prefix,
                           const int config, const int tslice, const int nb_ev,
                           const std::vector<std::pair<double,double> >& results){
  //build filename
  std::string filename = (boost::format("%s/%s_nev%d.%04d.%03d")% path % prefix
                          % nb_ev %config % tslice).str();
  std::ofstream eigenvalues(filename, std::ofstream::binary);
  //eigenvalues.write(reinterpret_cast<const std::string&>(&results[0]),
  //    results.size()*sizeof(std::pair<double,double>));
  for (auto el:results){ 
    eigenvalues << el.first << "\t" << el.second << "\n";
  }
  eigenvalues.close();
}

//Write eigenvalues from std::vector to binary
void write_eigenvalues_bin(const std::string& path,const std::string& prefix,
                           const int config_i, const int t, const int nb_ev,
                           std::vector<double>& ev){
  //Build filename
  std::string filename = (boost::format("%s/%s.%04d.%03d")% path % prefix 
                          % config_i % t).str();

  std::ofstream outfile(filename, std::ofstream::binary);
  if(outfile) {
    outfile.write(reinterpret_cast<char*>(&ev[0]), ev.size()*sizeof(double));
    outfile.close();
    std::cout << "recovery written to: " << filename << std::endl;
  }
}

//write gauge trafo matrices to binary file
void write_gauge_matrices(const std::string& prefix, Eigen::Matrix3cd* G) {
  int V3 = pars -> get_int("V3");

  std::ofstream outfile (prefix, std::ofstream::binary);
  for (int vol = 0; vol < V3; ++vol) {
    //outfile.write(reinterpret_cast<const std::string&>(&(G[vol])), 2*9*sizeof(double));
    outfile << G << "\n";
  }//end vol
  outfile.close();
}
//write gauge link matrices of one timeslice to binary file
void write_link_matrices_ts(const std::string& prefix) {

  int V3 = pars -> get_int("V3");
  std::ofstream outfile (prefix, std::ofstream::binary);
  for (int vol = 0; vol < V3; ++vol) {
    for (int mu = 0; mu < 3; ++mu) {
      Eigen::Matrix3cd tmp = eigen_timeslice -> get_gauge(vol,mu); 
      //outfile.write(reinterpret_cast<const std::string&>(&tmp), 2*9*sizeof(double));
      outfile << tmp << "\n";
    }//end dir
  }//end vol
  outfile.close();
}

