#include "local_update.h"

int LocalUpdate(Supercell & supercell, Site & one_site, double T) {
    // Flip one spin.
    // Energy and spin before flip.
    double energy_old = supercell.Hamiltonian(supercell.base_site, one_site);
    std::vector<double> old_spin = one_site.spin;

    // Energy and spin after flip.
    one_site.spin = RandomSpin(*one_site.spin_scaling);
    double energy_new = supercell.Hamiltonian(supercell.base_site, one_site);
    double de = energy_new - energy_old;

    // Judge whether to flip.
    double crition = RandomFloat();
    if (crition > exp(-de/(KB*T))) {
        one_site.spin = old_spin;
    } else {
        supercell.lattice.total_energy += de;
    }

    return 1;
}