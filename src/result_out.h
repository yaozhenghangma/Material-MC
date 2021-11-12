#ifndef RESULT_OUT
#define RESULT_OUT

#include <fmt/core.h>
#include <fmt/os.h>

#include "MC_structure.h"

int WriteOutput(MonteCarlo & monte_carlo, std::vector<double> energy, std::vector<double> Cv, \
std::vector<double> moment, std::vector<double> Ki, \
std::vector<double> moment_x, std::vector<double> moment_y, std::vector<double> moment_z, \
std::vector<double> Ki_x, std::vector<double> Ki_y, std::vector<double> Ki_z, \
std::string output_file);

#endif