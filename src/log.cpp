#include "log.h"

int WriteHamiltonian(Supercell & supercell, fmt::v8::ostream & logger) {
    switch (supercell.lattice.hamiltonian_type) {
        case HamiltonianType::Heisenberg :
            logger.print("Hamiltonian: Heisenberg model.");
            logger.print("Hamiltonian function: H = J*S_i*S_j");
            break;
        case HamiltonianType::Heisenberg_x_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x");
            break;
        case HamiltonianType::Heisenberg_y_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y");
            break;
        case HamiltonianType::Heisenberg_z_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z");
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y");
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z");
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x");
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with r axis anisotropy.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z");
            break;
        case HamiltonianType::Heisenberg_with_field :
            logger.print("Hamiltonian: Heisenberg model in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + B*S_i");
            break;
        case HamiltonianType::Heisenberg_x_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + B*S_i");
            break;
        case HamiltonianType::Heisenberg_y_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + B*S_i");
            break;
        case HamiltonianType::Heisenberg_z_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + B*S_i");
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x + B*S_i");
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with r axis anisotropy in external magnetic field.");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z + B*S_i");
            break;
        case HamiltonianType::Heisenberg_custom :
            logger.print("Hamiltonian: Custom Hamiltonian");
            break;
        default:
            logger.print("Known Hamiltonian.");
            break;
    }
    return 0;
}

int WriteLog(Supercell & supercell, MonteCarlo & monte_carlo, fmt::v8::ostream & logger) {
    logger.print("Lattice constants:");
    logger.print("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.lattice.a[0], supercell.lattice.a[1], supercell.lattice.a[2]);
    logger.print("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.lattice.b[0], supercell.lattice.b[1], supercell.lattice.b[2]);
    logger.print("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.lattice.c[0], supercell.lattice.c[1], supercell.lattice.c[2]);
    WriteHamiltonian(supercell, logger);
    logger.print("Anisotropy:");
    logger.print("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.base_site.anisotropy[0], supercell.base_site.anisotropy[1], supercell.base_site.anisotropy[2]);
    logger.print("Magnetic field (meV/spin):");
    logger.print("{:12.5f} {:12.5f} {:12.5f}", \
    supercell.base_site.B[0], supercell.base_site.B[1], supercell.base_site.B[2]);
    logger.print("Cell information:");
    std::string tmp_string;
    for(int i=0; i<supercell.base_site.number; i++) {
        logger.print("Element {}, spin {}:", supercell.base_site.elements[i], supercell.base_site.spin_scaling[i]);
        logger.print("Anisotropic ratio: {:12.5f} {:12.5f} {:12.5f}", \
        supercell.base_site.anisotropic_ratio[i][0], \
        supercell.base_site.anisotropic_ratio[i][1], \
        supercell.base_site.anisotropic_ratio[i][2]);
        logger.print("Super-exchange parameters:");
        tmp_string = fmt::format("Super-exchange parameters: ");
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{:12.5f} ", supercell.base_site.super_exchange_parameter[i][j]);
        }
        logger.print(tmp_string);
        tmp_string = fmt::format("Coordination number: ");
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{} ", supercell.site[0][0][0][i].neighbor[j].size());
        }
        logger.print(tmp_string);
    }
    logger.print("Initialization:");
    logger.print("");
    return 0;
}