
#include "statistics.h"

static double flip_sign(double value) {
  if(value < 0.) value *= (-1.);
  return value;
}

//calculate average of all eigenvalues 
std::vector<double> ev_average(const char* filename, const int nb_ev) { 
  //vector to hold result fixed size and set to zero
  std::vector<double> avg (nb_ev, 0);
  //counter for files (1/n in average)
  int cnt = 0;

  for (int config = 600; config < 792; config += 8) {//loop over configurations
    for (int tslice = 0; tslice < 48; ++tslice) {//loop over timeslices

      //temporary vector to hold eigenvalues of one file
      std::vector<double> tmp;
      read_eigenvalues_bin(filename, config, tslice, nb_ev, tmp);
      //ensure size
      if (tmp.size() == nb_ev) { 

        //incrementally add to avg
        for (int element = 0; element < nb_ev; ++element) {
          avg.at(element) += flip_sign(tmp.at(element));
        }
        ++cnt;

      }
    }
  }

  //divide by files read in
  double n = cnt;
  for (auto& average : avg) average /= n;
  return avg;

}

//calculate standard deviation of average over all eigenvalues
std::vector<double> std_dev (const char* filename, const int nb_ev, const std::vector<double>& avg_ev) {
  //vector holding standard deviations fixed size 0
  std::vector<double> stdev (nb_ev, 0);
  //counter for files (1/n in average)
  int cnt = 0;

  for (int config = 600; config < 792; config += 8) {//loop over configurations
    for (int tslice = 0; tslice < 48; ++tslice) {//loop over timeslices

      //temporary vector to hold eigenvalues of one file
      std::vector<double> tmp;
      read_eigenvalues_bin(filename, config, tslice, nb_ev, tmp);
      //ensure size
      if (tmp.size() == nb_ev) { 

        //incrementally add to avg
        for (int element = 0; element < nb_ev; ++element) {
          stdev.at(element) += (flip_sign(tmp.at(element))-avg_ev.at(element));
        }
        ++cnt;

      }
    }
  }

  //divide by files read in
  double n = cnt - 1.;
  for (auto& average : stdev) average /= n;
return stdev;
}
