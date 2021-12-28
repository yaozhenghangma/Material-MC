#ifndef RANDOM_FUNCTION
#define RANDOM_FUNCTION

#include <random>
#include <vector>

std::vector<int> RandomSite(int n_x, int n_y, int n_z, int base_n);
std::vector<double> RandomSpin(double scaling);
double RandomFloat();

#endif