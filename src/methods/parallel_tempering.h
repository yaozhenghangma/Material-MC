#ifndef PARALLEL_TEMPERING
#define PARALLEL_TEMPERING

#include <mpi.h>

#include "../constants.h"
#include "../MC_structure.h"
#include "../random_function.h"
#include "../spin_out.h"
#include "local_update.h"

int ParallelTemperingMonteCarlo(MPI_Comm world,
MonteCarlo & monte_carlo, Supercell & supercell, std::string & spin_structure_file_prefix,
std::vector<double> & energy, std::vector<double> & Cv,
std::vector<double> & moment, std::vector<double> & chi,
std::vector<double> & moment_projection, std::vector<double> & chi_projection);

#endif