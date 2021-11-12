#include "Hamiltonion.h"

double Heisenberg(BaseSite & base_site, Site & site) {
    double energy = 0;
    std::vector<double> spin_sum;
    for(int i=0; i<*site.neighbor_number; i++) {
        spin_sum = {0, 0, 0};
        for(int j=0; j<site.neighbor[i].size(); j++) {
            spin_sum[0] += (*site.neighbor[i][j]).spin[0];
            spin_sum[1] += (*site.neighbor[i][j]).spin[1];
            spin_sum[2] += (*site.neighbor[i][j]).spin[2];
        }
        energy += (*site.super_exchange_parameter)[i] * (site.spin[0]*spin_sum[0] + site.spin[1]*spin_sum[1] + site.spin[2]*spin_sum[2]);
    }

    energy += 2 * base_site.anisotropic_factor_D * (*site.anisotropic_factor)*site.spin[2]*site.spin[2];
    return energy;
}

double Heisenberg_3_anisotropy(BaseSite & base_site, Site & site) {
    double energy = 0;
    std::vector<double> spin_sum;
    for(int i=0; i<*site.neighbor_number; i++) {
        spin_sum = {0, 0, 0};
        for(int j=0; j<site.neighbor[i].size(); j++) {
            spin_sum[0] += (*site.neighbor[i][j]).spin[0];
            spin_sum[1] += (*site.neighbor[i][j]).spin[1];
            spin_sum[2] += (*site.neighbor[i][j]).spin[2];
        }
        energy += (*site.super_exchange_parameter)[i] * (site.spin[0]*spin_sum[0] + site.spin[1]*spin_sum[1] + site.spin[2]*spin_sum[2]);
    }

    energy += 2 * base_site.anisotropic_factor_D * (*site.anisotropic_factor)*site.spin[2]*site.spin[2];
    energy += 2 * base_site.anisotropic_factor_En * site.spin[0]*site.spin[0];
    energy -= 2 * base_site.anisotropic_factor_En * site.spin[1]*site.spin[1];
    return energy;
}

double Heisenberg_external_field(BaseSite & base_site, Site & site) {
    double energy = 0;
    std::vector<double> spin_sum;

    // Exchange energy
    for(int i=0; i<*site.neighbor_number; i++) {
        spin_sum = {0, 0, 0};
        for(int j=0; j<site.neighbor[i].size(); j++) {
            spin_sum[0] += (*site.neighbor[i][j]).spin[0];
            spin_sum[1] += (*site.neighbor[i][j]).spin[1];
            spin_sum[2] += (*site.neighbor[i][j]).spin[2];
        }
        energy += (*site.super_exchange_parameter)[i] * (site.spin[0]*spin_sum[0] + site.spin[1]*spin_sum[1] + site.spin[2]*spin_sum[2]);
    }

    // Magnetic anisotropy energy
    energy += 2 * base_site.anisotropic_factor_D * (*site.anisotropic_factor)*site.spin[2]*site.spin[2];
    energy += 2 * base_site.anisotropic_factor_En * site.spin[0]*site.spin[0];
    energy -= 2 * base_site.anisotropic_factor_En * site.spin[1]*site.spin[1];

    // External magnetic field
    energy -= 4 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}