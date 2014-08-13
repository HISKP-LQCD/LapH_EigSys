#ifndef _STATISTICS_H_
#define _STATISTICS_H_
#include <vector>
#include "read_write.h"

//Calculates the average of all eigenvalues in the ensemble
//Inputs: (in order of appearance) the prefix of the filename, Number of Eigenvalues
//in each file
//Returns: std::vector of averages of eigenvalues
std::vector<double> ev_average(const char* filename, const int nb_ev);

//Calculates the standard deviation of all eigenvalues if average is supplied
//Inputs: (in order of appearance) prefix of filenmae, number of eigenvalues, reference to vector of averages
//Returns: std::vector of standard deviations of eigenvalues
std::vector<double> std_dev (const char* filename, const int nb_ev, const std::vector<double>& avg_ev);
#endif //_STATISTICS_H_
