#include "par_io.h"

IO* IO::getInstance(){
  static IO theInstance;
return &theInstance;
}

void IO::set_inputpath(std::string infile) {
  input_path = infile;
}

void IO::set_values(const char* infile) {
  //one line of config file
  std::string config_line;
  //key and value to divide string to
  std::string keyword;
  std::string value;
  
  //three map class instances for int double and const char
  std::map <std::string, std::string> names;
  std::map <std::string, int> integers;
  std::map <std::string, double> floatings;

  //stream for config file
  std::fstream config;
  config.open(infile, std::fstream::in);

  //get each line of configuration into config_line
  while (!std::getline(config, config_line).eof()){
    if(config_line[0] != '#') {
        //find first occurence of '=' in each line and divide strings
        unsigned pos = config_line.find_first_of("=");
        keyword = config_line.substr(0,(pos-1));
        value = config_line.substr( (pos+2), std::string::npos);

        //String valued map variables
        if (keyword == "config_path") {
          names["config_path"] = value;
        }
        if (keyword == "result_path") {
          names["result_path"] = value;
        }
        if (keyword == "in_path") {
          names["in_path"] = value;
        }

        //Integer valued map variables
        if (keyword == "LT" ) {
          integers["LT"] = atoi(value.c_str());
        }
        if (keyword == "LX" ) {
          integers["LX"] = atoi(value.c_str());
        }
        if (keyword == "NEV") {
          integers["NEV"] = atoi(value.c_str());
        }
        if (keyword == "iter") {
          integers["iter"] = atoi(value.c_str());
        }
        if (keyword == "DEG") {
          integers["DEG"] = atoi(value.c_str());
        }

        //floating point valued map variables
        if (keyword == "alpha_1") {
          floatings["alpha_1"] = atof(value.c_str());
        }
        if (keyword == "alpha_2") {
          floatings["alpha_2"] = atof(value.c_str());
        }
        if (keyword == "LAM_L") {
          floatings["LAM_L"] = atof(value.c_str());
        }
        if (keyword == "LAM_C") {
          floatings["LAM_C"] = atof(value.c_str());
        }
    }
  }
  config.close();

  //set the values
  config_path = names["config_path"];
  result_path = names["result_path"];
  input_path = names["in_path"];

  LT = integers["LT"];
  LX = integers["LX"];
  LY = LX;
  LZ = LX;
  
  NEV = integers["NEV"];
  NCOL = 3;
  NDIR = 4;
  V3 = LX*LY*LZ;
  MAT_ENTRIES=V3*NCOL;

  V_TS = LX*LY*LZ*NDIR*NCOL*NCOL*2; //2 is for complex
  V4_LIME = V_TS * LT;
  
  iter = integers["iter"];
  alpha_1 = floatings["alpha_1"];
  alpha_2 = floatings["alpha_2"];

  LAM_L = floatings["LAM_L"];
  LAM_C = floatings["LAM_C"];
  DEG = integers["DEG"];

}

void IO::print_summary() {
  std::cout << "Configuration Summary for the eigensystem calculation:" << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "Lattice size is:                                   " << LT << " x " << LX << " x " << LY << " x " << LZ << std::endl;
  std::cout << "Number of eigenvectors is:                         " << NEV << std::endl;
  std::cout << "Hyp-smearing parameters (alpha_1, alpha_2, iter):  " << alpha_1 << ", " << alpha_2 << ", " << iter << std::endl;
  std::cout << "Using Chebyshev parameters:                        Lambda_l = " << LAM_L << " Lambda_c = " << LAM_C << std::endl;
  std::cout << "Degree of Chebyshev polynomial:                    " << DEG << std::endl;
  std::cout << "Configurations are read in from:                   " << config_path << std::endl;
  std::cout << "Results are stored in:                             " << result_path << std::endl;
  std::cout << "Entries in Matrix:				   " << MAT_ENTRIES << std::endl;
  std::cout << "number of colours:				   " << NCOL << std::endl;
  std::cout << "number of directions:				   " << NDIR << std::endl;
  std::cout << "Volume in gauge matrices:			   " << V3 << std::endl;
  std::cout << "lime config volume in doubles:			   " << V4_LIME << std::endl;
}

int IO::get_int(std::string spec) {
  if (spec == "LT") {
    return LT;
  }
  //so far only isotropic volumes considered
  if (spec == "LX" || spec == "LY" || spec == "LZ") {
    return LX;
  }

  if (spec == "NEV") {
    return NEV;
  }

  if (spec == "NDIR") {
    return NDIR;
  }

  if (spec == "NCOL") {
    return NCOL;
  }

  if (spec == "V3") {
    return V3;
  }

  if (spec == "V4_LIME") {
    return V4_LIME;
  }

  if (spec == "V_TS") {
    return V_TS;
  }

  if (spec == "iter") {
    return iter;
  }
  if (spec == "MAT_ENTRIES") {
    return MAT_ENTRIES;
  }
  if (spec == "DEG") {
    return DEG;
  }
  else return -1;

}

double IO::get_float(std::string spec) {
  if (spec == "alpha_1") {
    return alpha_1;
  }

  if (spec == "alpha_2") {
    return alpha_2;
  }
  
  if (spec == "lambda_l") {
    return LAM_L;
  }

  if (spec == "lambda_c") {
    return LAM_C;
  }
  else return -1;
}

std::string IO::get_path(std::string spec) {
  if (spec == "conf") {
    return config_path;
  }

  if (spec == "res") {
    return result_path;
  }

  if (spec == "input") {
    return input_path;
  }

  else return "not found";
}



