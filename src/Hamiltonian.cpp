#include "Hamiltonian.h"

#include <cstdlib>
#include <iostream>

namespace {

void ResolveKhCyclicComponents(char direction_label,
                               int & gamma,
                               int & alpha,
                               int & beta,
                               int shell,
                               int entry) {
    switch (direction_label) {
        case 'x':
            gamma = 0;
            alpha = 1;
            beta = 2;
            break;
        case 'y':
            gamma = 1;
            alpha = 2;
            beta = 0;
            break;
        case 'z':
            gamma = 2;
            alpha = 0;
            beta = 1;
            break;
        default:
            std::cerr << "Invalid KH bond direction label in Hamiltonian evaluation: "
                      << "shell " << shell
                      << ", entry " << entry
                      << ", label \"" << direction_label << "\".\n";
            exit(-1);
    }
}

} // namespace

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

    return energy;
}

double Heisenberg_xyz_anisotropy(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    return energy;
}

double Heisenberg_with_field(BaseSite & base_site, Site & site) {
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

    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Heisenberg_xyz_anisotropy_with_field(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Heisenberg_x_anisotropy(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    
    return energy;
}

double Heisenberg_y_anisotropy(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    
    return energy;
}

double Heisenberg_z_anisotropy(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    return energy;
}

double Heisenberg_xy_anisotropy(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    
    return energy;
}

double Heisenberg_yz_anisotropy(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    return energy;
}

double Heisenberg_zx_anisotropy(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    return energy;
}

double Heisenberg_x_anisotropy_with_field(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    
    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Heisenberg_y_anisotropy_with_field(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    
    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Heisenberg_z_anisotropy_with_field(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Heisenberg_xy_anisotropy_with_field(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    
    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Heisenberg_yz_anisotropy_with_field(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[1] * (*site.anisotropic_ratio)[1] * site.spin[1] * site.spin[1];
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Heisenberg_zx_anisotropy_with_field(BaseSite & base_site, Site & site) {
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
    energy += base_site.anisotropy[0] * (*site.anisotropic_ratio)[0] * site.spin[0] * site.spin[0];
    energy += base_site.anisotropy[2] * (*site.anisotropic_ratio)[2] * site.spin[2] * site.spin[2];

    energy -= 2 * (base_site.B[0]*site.spin[0] + base_site.B[1]*site.spin[1] + base_site.B[2]*site.spin[2]);

    return energy;
}

double Kitaev_Heisenberg(BaseSite & base_site, Site & site) {
    double energy = 0;
    if(site.neighbor.size() != site.neighbor_direction.size()) {
        std::cerr << "Inconsistent KH neighbor shell metadata during Hamiltonian evaluation: "
                  << site.neighbor.size() << " neighbor shells but "
                  << site.neighbor_direction.size() << " direction-label shells.\n";
        exit(-1);
    }

    for(int i=0; i<*site.neighbor_number; i++) {
        if(site.neighbor[i].size() != site.neighbor_direction[i].size()) {
            std::cerr << "Inconsistent KH neighbor metadata during Hamiltonian evaluation: "
                      << "shell " << i
                      << " has " << site.neighbor[i].size() << " neighbors but "
                      << site.neighbor_direction[i].size() << " direction labels.\n";
            exit(-1);
        }

        for(int j=0; j<site.neighbor[i].size(); j++) {
            int gamma = 0;
            int alpha = 0;
            int beta = 0;
            ResolveKhCyclicComponents(site.neighbor_direction[i][j], gamma, alpha, beta, i, j);

            Site & neighbor = *site.neighbor[i][j];
            energy += base_site.kh_j * (site.spin[0]*neighbor.spin[0] + site.spin[1]*neighbor.spin[1] + site.spin[2]*neighbor.spin[2]);
            energy += base_site.kh_k * site.spin[gamma] * neighbor.spin[gamma];
            energy += base_site.kh_g * (site.spin[alpha]*neighbor.spin[beta] + site.spin[beta]*neighbor.spin[alpha]);
            energy += base_site.kh_gp * (
                site.spin[alpha]*neighbor.spin[gamma]
                + site.spin[gamma]*neighbor.spin[alpha]
                + site.spin[beta]*neighbor.spin[gamma]
                + site.spin[gamma]*neighbor.spin[beta]);
        }
    }

    return energy;
}

double Kitaev_Heisenberg_base(BaseSite & base_site, Site & site) {
    return Kitaev_Heisenberg(base_site, site) * 0.5;
}

double Heisenberg_base(BaseSite & base_site, Site & site) {
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