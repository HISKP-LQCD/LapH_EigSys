#include <tuple>
#include <utility>
//#include "config_utils.h"
#include "sourceshape_funcs.h"
#include "structs.h"
#include "par_io.h"
//#include "variables.h"
#include "read_write.h"

//int up_3d[V3][3];
//int down_3d[V3][3];

int main(int argc, char* argv[]) {

  // handle input file
  IO* pars = IO::getInstance();
  pars->set_values("parameters.txt");
  if (argc != 4) std::cout << "No parameters entered. Arguments are: configuration timeslice dim(V)" << std::endl;
  //Parameters from command line and variables
  int config = atoi(argv[1]);
  int time = atoi(argv[2]);
  int dim_col = atoi(argv[3]);
  std::cout << "Entered parameters: " << "config: " << config << " ts: " << time << " nev: " << dim_col << std::endl;
  const int dim_row = 3 * pars -> get_int("V3");
  //hopping3d(up_3d, down_3d);

  //initialize memory for timeslice and results
  Eigen::MatrixXcd V (dim_row, dim_col);

  std::vector<double> ev;
  ev.reserve(dim_col);

  std::vector<shp> result;
  result.reserve(pars-> get_int("V3"));
  int size1 = result.size();

  //temporary vector for averaging
  std::vector<std::pair<double,double> > shp_avg;

  //Read in eigenvectors of one timeslice
  read_evectors_bin_ts("eigenvectors",
      config, time, dim_col, V);
  std::cout << V(0,0) << std::endl;
  
  //Read in eigenvalues of one timeslice
  //read_eigenvalues_bin("/hiskp2/helmes/A60_0600_L120/eigensystems/hyp_05_07_03/eigenvalues",
  //    config, time, dim_col, ev);

  //for(auto& value : ev) std::cout << value << std::endl;
  
  std::complex<double> reference = (V.col(2))(dim_row-1);
  //weight_eigenvectors(ev, dim_col, V);
  
 //calculate sourceshape
  std::cout << "calculating, please wait..." << std::endl;
  source_shape_complete(V, dim_col, result);
  afterburner(result);
  
  //average psi with same r afterburner has to complete before
  avg_radii(result, shp_avg);
  std::cout << "avg complete" << std::endl; 
  for (auto& s:shp_avg) std::cout << std::get<0>(s) << " "<< std::get<1>(s) << std::endl;

  //save sourceshape
  write_sourceshape_ascii("sourceshape_A40_L20",config,time,dim_col,shp_avg);
  write_sourceshape_bin("sourceshape_A40_L20",config,time,dim_col,shp_avg);

  //Clean-up
  return 0;

}
