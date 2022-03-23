#include "local_update.h"

int LocalUpdateHeisenberg(Supercell & supercell, Site & one_site, double T) {
    // Flip one spin for Heisenberg model.
    // Energy and spin before flip.
    double energy_old = supercell.Hamiltonian(supercell.base_site, one_site);
    std::vector<double> old_spin = one_site.spin;

    // Energy and spin after flip.
    one_site.spin = RandomSpin(*one_site.spin_scaling);
    double energy_new = supercell.Hamiltonian(supercell.base_site, one_site);
    double de = energy_new - energy_old;

    // Judge whether to flip.
    if (RandomFloat() > exp(-de/(KB*T))) {
        one_site.spin = old_spin;
    } else {
        supercell.lattice.total_energy += de;
    }

    return 1;
}

int FlipSpin(std::vector<double> & spin) {
    spin[0] = -spin[0];
    spin[1] = -spin[1];
    spin[2] = -spin[2];
    return 0;
}

int LocalUpdateIsing(Supercell & supercell, Site & one_site, double T) {
    // Flip one spin for Ising model.
    // Energy and spin before flip.
    double energy_old = supercell.Hamiltonian(supercell.base_site, one_site);
    std::vector<double> old_spin = one_site.spin;

    // Energy and spin after flip.
    FlipSpin(one_site.spin);
    double energy_new = supercell.Hamiltonian(supercell.base_site, one_site);
    double de = energy_new - energy_old;

    // Judge whether to flip.
    if (RandomFloat() > exp(-de/(KB*T))) {
        one_site.spin = old_spin;
    } else {
        supercell.lattice.total_energy += de;
    }

    return 1;
}