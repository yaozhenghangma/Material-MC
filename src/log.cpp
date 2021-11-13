#include "log.h"

int WriteHamiltonion(Supercell & supercell, std::shared_ptr<spdlog::logger> logger) {
    switch (supercell.lattice.hamiltonion_type) {
        case HamiltonionType::Heisenberg :
            logger->info("Hamiltonion: Heisenberg model.");
            logger->info("Hamiltonion function: H = J*S_i*S_j");
            break;
        case HamiltonionType::Heisenberg_x_anisotropy :
            logger->info("Hamiltonion: Heisenberg model with 1 axis anisotropy.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_x*(S_i*S_i)_x");
            break;
        case HamiltonionType::Heisenberg_y_anisotropy :
            logger->info("Hamiltonion: Heisenberg model with 1 axis anisotropy.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_y*(S_i*S_i)_y");
            break;
        case HamiltonionType::Heisenberg_z_anisotropy :
            logger->info("Hamiltonion: Heisenberg model with 1 axis anisotropy.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_z*(S_i*S_i)_z");
            break;
        case HamiltonionType::Heisenberg_xy_anisotropy :
            logger->info("Hamiltonion: Heisenberg model with 2 axis anisotropy.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y");
            break;
        case HamiltonionType::Heisenberg_yz_anisotropy :
            logger->info("Hamiltonion: Heisenberg model with 2 axis anisotropy.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z");
            break;
        case HamiltonionType::Heisenberg_zx_anisotropy :
            logger->info("Hamiltonion: Heisenberg model with 2 axis anisotropy.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x");
            break;
        case HamiltonionType::Heisenberg_xyz_anisotropy :
            logger->info("Hamiltonion: Heisenberg model with r axis anisotropy.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z");
            break;
        case HamiltonionType::Heisenberg_with_field :
            logger->info("Hamiltonion: Heisenberg model in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + B*S_i");
            break;
        case HamiltonionType::Heisenberg_x_anisotropy_with_field :
            logger->info("Hamiltonion: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + B*S_i");
            break;
        case HamiltonionType::Heisenberg_y_anisotropy_with_field :
            logger->info("Hamiltonion: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + B*S_i");
            break;
        case HamiltonionType::Heisenberg_z_anisotropy_with_field :
            logger->info("Hamiltonion: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonionType::Heisenberg_xy_anisotropy_with_field :
            logger->info("Hamiltonion: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + B*S_i");
            break;
        case HamiltonionType::Heisenberg_yz_anisotropy_with_field :
            logger->info("Hamiltonion: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonionType::Heisenberg_zx_anisotropy_with_field :
            logger->info("Hamiltonion: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x + B*S_i");
            break;
        case HamiltonionType::Heisenberg_xyz_anisotropy_with_field :
            logger->info("Hamiltonion: Heisenberg model with r axis anisotropy in external magnetic field.");
            logger->info("Hamiltonion function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z + B*S_i");
            break;
        default:
            break;
    }
    return 0;
}

int WriteLog(Supercell & supercell, MonteCarlo & monte_carlo, std::shared_ptr<spdlog::logger> logger) {
    logger->info("Lattice constants:");
    logger->info("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.lattice.a[0], supercell.lattice.a[1], supercell.lattice.a[2]);
    logger->info("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.lattice.b[0], supercell.lattice.b[1], supercell.lattice.b[2]);
    logger->info("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.lattice.c[0], supercell.lattice.c[1], supercell.lattice.c[2]);
    WriteHamiltonion(supercell, logger);
    logger->info("Anisotropy:");
    logger->info("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.base_site.anisotropic_factor[0], supercell.base_site.anisotropic_factor[1], supercell.base_site.anisotropic_factor[2]);
    logger->info("Magnetic field (meV/spin):");
    logger->info("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.base_site.B[0], supercell.base_site.B[1], supercell.base_site.B[2]);
    logger->info("Super-exchange parameters:");
    std::string tmp_string;
    for(int i=0; i<supercell.base_site.number; i++) {
        tmp_string = fmt::format("Element {}, ", supercell.base_site.elements[i]);
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{:12.5f} ", supercell.base_site.super_exchange_parameter[i][j]);
        }
        logger->info(tmp_string);
    }
    logger->info("Coordination number:");
    for(int i=0; i<supercell.base_site.number; i++) {
        tmp_string = fmt::format("Element {}, ", supercell.base_site.elements[i]);
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{} ", supercell.site[0][0][0][i].neighbor[j].size());
        }
        logger->info(tmp_string);
    }
    logger->info("Initialization:");
    logger->info("");
    return 0;
}