#include "local_update.h"

/**
 * @file local_update.cpp
 * @brief Local Metropolis spin-update kernels for Heisenberg and Ising models.
 */

/**
 * @brief Performs one Metropolis trial update for a Heisenberg spin.
 *
 * A new random spin direction is proposed, local energy difference is computed,
 * and the proposal is accepted with Boltzmann probability exp(-dE / (k_B T)).
 *
 * @param supercell Simulation supercell.
 * @param one_site Site to update.
 * @param T Temperature.
 * @return int Returns 1 after one trial update.
 */
int LocalUpdateHeisenberg(Supercell & supercell, Site & one_site, double T) {
    // Flip one spin for Heisenberg model.
    // Energy and spin before flip.
    double energy_old = supercell.Hamiltonian(supercell.base_site, one_site);
    std::vector<double> old_spin = one_site.spin;

    // Energy and spin after proposing a random reorientation.
    one_site.spin = RandomSpin(*one_site.spin_scaling);
    double energy_new = supercell.Hamiltonian(supercell.base_site, one_site);
    double de = energy_new - energy_old;

    // Metropolis accept/reject step.
    if (RandomFloat() > exp(-de/(KB*T))) {
        one_site.spin = old_spin;
    } else {
        supercell.lattice.total_energy += de;
    }

    return 1;
}

/**
 * @brief Reverses a spin vector in-place.
 *
 * @param spin Spin vector {sx, sy, sz}.
 * @return int Returns 0 on completion.
 */
int FlipSpin(std::vector<double> & spin) {
    spin[0] = -spin[0];
    spin[1] = -spin[1];
    spin[2] = -spin[2];
    return 0;
}

/**
 * @brief Performs one Metropolis trial update for an Ising spin.
 *
 * The proposal is a full spin inversion. Acceptance follows the same
 * Boltzmann criterion as the Heisenberg local update.
 *
 * @param supercell Simulation supercell.
 * @param one_site Site to update.
 * @param T Temperature.
 * @return int Returns 1 after one trial update.
 */
int LocalUpdateIsing(Supercell & supercell, Site & one_site, double T) {
    // Flip one spin for Ising model.
    // Energy and spin before flip.
    double energy_old = supercell.Hamiltonian(supercell.base_site, one_site);
    std::vector<double> old_spin = one_site.spin;

    // Energy and spin after proposing an Ising inversion.
    FlipSpin(one_site.spin);
    double energy_new = supercell.Hamiltonian(supercell.base_site, one_site);
    double de = energy_new - energy_old;

    // Metropolis accept/reject step.
    if (RandomFloat() > exp(-de/(KB*T))) {
        one_site.spin = old_spin;
    } else {
        supercell.lattice.total_energy += de;
    }

    return 1;
}