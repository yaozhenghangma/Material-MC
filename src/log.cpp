#include "log.h"

int WriteHamiltonian(Supercell & supercell, std::shared_ptr<spdlog::logger> logger) {
    switch (supercell.lattice.hamiltonian_type) {
        case HamiltonianType::Heisenberg :
            logger->info("Hamiltonian: Heisenberg model.");
            logger->info("Hamiltonian function: H = J*S_i*S_j");
            break;
        case HamiltonianType::Heisenberg_x_anisotropy :
            logger->info("Hamiltonian: Heisenberg model with 1 axis anisotropy.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x");
            break;
        case HamiltonianType::Heisenberg_y_anisotropy :
            logger->info("Hamiltonian: Heisenberg model with 1 axis anisotropy.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y");
            break;
        case HamiltonianType::Heisenberg_z_anisotropy :
            logger->info("Hamiltonian: Heisenberg model with 1 axis anisotropy.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z");
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy :
            logger->info("Hamiltonian: Heisenberg model with 2 axis anisotropy.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y");
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy :
            logger->info("Hamiltonian: Heisenberg model with 2 axis anisotropy.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z");
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy :
            logger->info("Hamiltonian: Heisenberg model with 2 axis anisotropy.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x");
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy :
            logger->info("Hamiltonian: Heisenberg model with r axis anisotropy.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z");
            break;
        case HamiltonianType::Heisenberg_with_field :
            logger->info("Hamiltonian: Heisenberg model in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + B*S_i");
            break;
        case HamiltonianType::Heisenberg_x_anisotropy_with_field :
            logger->info("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + B*S_i");
            break;
        case HamiltonianType::Heisenberg_y_anisotropy_with_field :
            logger->info("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + B*S_i");
            break;
        case HamiltonianType::Heisenberg_z_anisotropy_with_field :
            logger->info("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy_with_field :
            logger->info("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + B*S_i");
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy_with_field :
            logger->info("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy_with_field :
            logger->info("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x + B*S_i");
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy_with_field :
            logger->info("Hamiltonian: Heisenberg model with r axis anisotropy in external magnetic field.");
            logger->info("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonianType::Heisenberg_custom :
            logger->info("Hamiltonian: Custom Hamiltonian")
            break;
        default:
            logger->info("Known Hamiltonian.")
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
    WriteHamiltonian(supercell, logger);
    logger->info("Anisotropy:");
    logger->info("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.base_site.anisotropy[0], supercell.base_site.anisotropy[1], supercell.base_site.anisotropy[2]);
    logger->info("Magnetic field (meV/spin):");
    logger->info("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.base_site.B[0], supercell.base_site.B[1], supercell.base_site.B[2]);
    logger->info("Cell information:");
    std::string tmp_string;
    for(int i=0; i<supercell.base_site.number; i++) {
        logger->info("Element {}, spin {}:", supercell.base_site.elements[i], supercell.base_site.spin_scaling[i]);
        logger->info("Anisotropic ratio: {:12.5f} {:12.5f} {:12.5f}", \
        supercell.base_site.anisotropic_ratio[i][0], \
        supercell.base_site.anisotropic_ratio[i][1], \
        supercell.base_site.anisotropic_ratio[i][2]);
        logger->info("Super-exchange parameters:");
        tmp_string = fmt::format("Super-exchange parameters: ");
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{:12.5f} ", supercell.base_site.super_exchange_parameter[i][j]);
        }
        logger->info(tmp_string);
        tmp_string = fmt::format("Coordination number: ");
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{} ", supercell.site[0][0][0][i].neighbor[j].size());
        }
        logger->info(tmp_string);
    }
    logger->info("Initialization:");
    logger->info("");
    return 0;
}