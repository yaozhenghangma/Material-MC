#ifndef RESULT_OUT
#define RESULT_OUT

#include <fmt/core.h>
#include <fmt/os.h>

#include "MC_structure.h"

int WriteOutput(MonteCarlo & monte_carlo, std::vector<double> energy, std::vector<double> Cv, \
std::vector<double> moment, std::vector<double> chi, \
std::vector<double> moment_x, std::vector<double> chi_x, \
bool field, std::string output_file);

#endif