#include "log.h"

/**
 * @brief Write a human-readable Hamiltonian summary to the log stream.
 *
 * The selected runtime Hamiltonian type is mapped to a text description and
 * a symbolic formula string.
 *
 * @param supercell Runtime container that stores selected hamiltonian_type.
 * @param logger Open fmt output stream for log.txt.
 * @return Always 0.
 */
int WriteHamiltonian(Supercell & supercell, fmt::v8::ostream & logger) {
    // KH energy is dispatched by model_type rather than hamiltonian_type.
    if(supercell.lattice.model_type == ModelType::Kitaev_Heisenberg) {
        logger.print("Hamiltonian: Kitaev-Heisenberg model (bond-dependent J/K/G/Gp terms).\n");
        return 0;
    }

    // Print Hamiltonian family and the corresponding symbolic expression.
    switch (supercell.lattice.hamiltonian_type) {
        case HamiltonianType::Heisenberg :
            logger.print("Hamiltonian: Heisenberg model.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j\n");
            break;
        case HamiltonianType::Heisenberg_x_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x\n");
            break;
        case HamiltonianType::Heisenberg_y_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y\n");
            break;
        case HamiltonianType::Heisenberg_z_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z\n");
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y\n");
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z\n");
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x\n");
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy :
            logger.print("Hamiltonian: Heisenberg model with r axis anisotropy.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z\n");
            break;
        case HamiltonianType::Heisenberg_with_field :
            logger.print("Hamiltonian: Heisenberg model in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_x_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_y_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_z_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 1 axis anisotropy in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_y*(S_i*S_i)_y + D_z*(S_i*S_i)_z + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with 2 axis anisotropy in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_z*(S_i*S_i)_z + D_x*(S_i*S_i)_x + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy_with_field :
            logger.print("Hamiltonian: Heisenberg model with r axis anisotropy in external magnetic field.\n");
            logger.print("Hamiltonian function: H = J*S_i*S_j + D_x*(S_i*S_i)_x + D_y*(S_i*S_i)_y + + D_z*(S_i*S_i)_z + B*S_i\n");
            break;
        case HamiltonianType::Heisenberg_custom :
            logger.print("Hamiltonian: Custom Hamiltonian\n");
            break;
        default:
            logger.print("Known Hamiltonian.\n");
            break;
    }
    return 0;
}

/**
 * @brief Write simulation metadata summary to the log stream.
 *
 * This function logs lattice vectors, algorithm/model selection, Hamiltonian
 * description, anisotropy/field values, and per-element interaction summaries.
 *
 * @param supercell Runtime structure containing lattice/base-site/site data.
 * @param monte_carlo Runtime Monte Carlo method selection.
 * @param logger Open fmt output stream for log.txt.
 * @return Always 0.
 */
int WriteLog(Supercell & supercell, MonteCarlo & monte_carlo, fmt::v8::ostream & logger) {
    // Lattice vectors from parsed POSCAR.
    logger.print("Lattice constants:\n");
    logger.print("{:12.5f} {:12.5f} {:12.5f}\n", \
    supercell.lattice.a[0], supercell.lattice.a[1], supercell.lattice.a[2]);
    logger.print("{:12.5f} {:12.5f} {:12.5f}\n", \
    supercell.lattice.b[0], supercell.lattice.b[1], supercell.lattice.b[2]);
    logger.print("{:12.5f} {:12.5f} {:12.5f}\n", \
    supercell.lattice.c[0], supercell.lattice.c[1], supercell.lattice.c[2]);
    // Runtime Monte Carlo algorithm dispatch.
    logger.print("Monte Carlo Method: ");
    switch (monte_carlo.methods) {
        case Methods::classical:
            logger.print("Classical Monte Carlo.\n");
            break;
        case Methods::parallel_tempering:
            logger.print("Parallel Tempering Monte Carlo.\n");
            break;
        default:
            logger.print("Unknown method.");
            break;
    }
    // Spin model family used by local updates.
    logger.print("Model Type: ");
    switch (supercell.lattice.model_type) {
        case ModelType::Heisenberg:
            logger.print("Heisenberg model.\n");
            break;
        case ModelType::Ising:
            logger.print("Ising model.\n");
            break;
        case ModelType::Kitaev_Heisenberg:
            logger.print("Kitaev-Heisenberg model.\n");
            break;
        default:
            logger.print("Unknown model.\n");
            break;
    }
    if(supercell.lattice.model_type == ModelType::Kitaev_Heisenberg) {
        logger.print("KH global couplings (J, K, G, Gp): {:12.5f} {:12.5f} {:12.5f} {:12.5f}\n", \
        supercell.base_site.kh_j, supercell.base_site.kh_k, supercell.base_site.kh_g, supercell.base_site.kh_gp);
        logger.print("KH convention: H_ij = J*Si*Sj + K*Si_gamma*Sj_gamma + G*(Si_alpha*Sj_beta + Si_beta*Sj_alpha) + Gp*(Si_alpha*Sj_gamma + Si_gamma*Sj_alpha + Si_beta*Sj_gamma + Si_gamma*Sj_beta).\n");
        logger.print("KH cyclic mapping: gamma=x => (alpha,beta)=(y,z), gamma=y => (alpha,beta)=(z,x), gamma=z => (alpha,beta)=(x,y).\n");
        if(supercell.base_site.kh_bond_type_direction.size() == 3) {
            logger.print("KH bond type mapping (type1, type2, type3): {} {} {}\n",
            supercell.base_site.kh_bond_type_direction[0],
            supercell.base_site.kh_bond_type_direction[1],
            supercell.base_site.kh_bond_type_direction[2]);
        }
    }
    // Hamiltonian summary and formula text.
    WriteHamiltonian(supercell, logger);
    logger.print("Anisotropy: ");
    logger.print("{:12.5f} {:12.5f} {:12.5f}\n", \
    supercell.base_site.anisotropy[0], supercell.base_site.anisotropy[1], supercell.base_site.anisotropy[2]);
    logger.print("Magnetic field (meV/spin): ");
    logger.print("{:12.5f} {:12.5f} {:12.5f}\n", \
    supercell.base_site.B[0], supercell.base_site.B[1], supercell.base_site.B[2]);
    // Per-element interaction summary in the primitive cell.
    logger.print("Cell information:\n");
    std::string tmp_string;
    for(int i=0; i<supercell.base_site.number; i++) {
        logger.print("Element {}, spin {}:\n", supercell.base_site.elements[i], supercell.base_site.spin_scaling[i]);
        logger.print("Anisotropic ratio: {:12.5f} {:12.5f} {:12.5f}\n", \
        supercell.base_site.anisotropic_ratio[i][0], \
        supercell.base_site.anisotropic_ratio[i][1], \
        supercell.base_site.anisotropic_ratio[i][2]);
        // List configured exchange parameters by neighbor shell.
        tmp_string = fmt::format("Super-exchange parameters: ");
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{:12.5f} ", supercell.base_site.super_exchange_parameter[i][j]);
        }
        tmp_string += "\n";
        logger.print(tmp_string);
        // Coordination numbers are taken from initialized neighbor links at (0,0,0).
        tmp_string = fmt::format("Coordination number: ");
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            tmp_string += fmt::format("{} ", supercell.site[0][0][0][i].neighbor[j].size());
        }
        tmp_string += "\n";
        logger.print(tmp_string);
    }
    logger.print("Initialization:\n");
    logger.print("\n");
    return 0;
}