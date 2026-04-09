#ifndef PARALLEL_TEMPERING
#define PARALLEL_TEMPERING

#include <mpi.h>

#include <string>
#include <vector>

#include "../constants.h"
#include "../MC_structure.h"
#include "../random_function.h"
#include "../spin_out.h"
#include "local_update.h"

/**
 * @brief Runs parallel tempering Monte Carlo sampling with MPI replica exchange.
 *
 * @param world MPI communicator.
 * @param monte_carlo Monte Carlo schedule and exchange interval settings.
 * @param supercell Simulation supercell.
 * @param spin_structure_file_prefix Spin output prefix (currently unused).
 * @param energy Output average energy per spin.
 * @param Cv Output heat capacity per spin.
 * @param moment Output average magnetic moment per spin.
 * @param chi Output magnetic susceptibility per spin.
 * @param moment_projection Output projected moment (when field enabled).
 * @param chi_projection Output projected susceptibility (when field enabled).
 * @return int Returns 0 on completion.
 */
int ParallelTemperingMonteCarlo(MPI_Comm world,
MonteCarlo & monte_carlo, Supercell & supercell, std::string & spin_structure_file_prefix,
std::vector<double> & energy, std::vector<double> & Cv,
std::vector<double> & moment, std::vector<double> & chi,
std::vector<double> & moment_projection, std::vector<double> & chi_projection);

#endif