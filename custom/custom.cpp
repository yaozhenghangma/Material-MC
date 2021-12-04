#include "custom.h"

double Hamiltonian_custom_base(BaseSite & base_site, Site & site) {
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

    return energy * 0.5;
}

double Hamiltonian_custom(BaseSite & base_site, Site & site) {
    double energy = 0;
    static double K1 = -43.2;
    static double K2 = 48.6;
    std::vector<double> spin_sum;
    double spin_x_square = site.spin[0]*site.spin[0];
    double spin_y_square = site.spin[1]*site.spin[1];
    double spin_z_square = site.spin[2]*site.spin[2];
    for(int i=0; i<*site.neighbor_number; i++) {
        spin_sum = {0, 0, 0};
        for(int j=0; j<site.neighbor[i].size(); j++) {
            spin_sum[0] += (*site.neighbor[i][j]).spin[0];
            spin_sum[1] += (*site.neighbor[i][j]).spin[1];
            spin_sum[2] += (*site.neighbor[i][j]).spin[2];
        }
        energy += (*site.super_exchange_parameter)[i] * (site.spin[0]*spin_sum[0] + site.spin[1]*spin_sum[1] + site.spin[2]*spin_sum[2]);
    }
    energy += K1 * (spin_x_square * spin_y_square + spin_y_square * spin_z_square + spin_z_square * spin_x_square);
    energy += K2 * (spin_x_square * spin_y_square * spin_z_square);

    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}
