#ifndef CLASSICAL
#define CLASSICAL

#include <boost/mpi.hpp>

#include "../constants.h"
#include "../MC_structure.h"
#include "../random_function.h"
#include "../spin_out.h"
#include "local_update.h"

int ClassicalMonteCarlo(boost::mpi::environment & env, boost::mpi::communicator & world, 
MonteCarlo & monte_carlo, Supercell & supercell, std::string & spin_structure_file_prefix, 
std::vector<double> & energy, std::vector<double> & Cv,
std::vector<double> & moment, std::vector<double> & chi,
std::vector<double> & moment_projection, std::vector<double> & chi_projection);

#endif