#ifndef LOCAL_UPDATE
#define LOCAL_UPDATE

#include "../MC_structure.h"
#include "../random_function.h"

/**
 * @brief Performs one local Metropolis update for a Heisenberg spin.
 *
 * @param supercell Simulation supercell.
 * @param one_site Site selected for update.
 * @param T Temperature.
 * @return int Returns 1 after one trial update.
 */
int LocalUpdateHeisenberg(Supercell & supercell, Site & one_site, double T);

/**
 * @brief Performs one local Metropolis update for an Ising spin.
 *
 * @param supercell Simulation supercell.
 * @param one_site Site selected for update.
 * @param T Temperature.
 * @return int Returns 1 after one trial update.
 */
int LocalUpdateIsing(Supercell & supercell, Site & one_site, double T);

#endif