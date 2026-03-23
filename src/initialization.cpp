#include "initialization.h"

#include <iostream>

/**
 * @brief Parse command-line options for input/output file paths.
 *
 * Supported options:
 * - -c <file>: structure file (POSCAR)
 * - -i <file>: input TOML file
 * - -o <file>: output table file
 * - -s <prefix>: spin-structure output prefix
 * - -h: print usage and exit
 *
 * @param argc Argument count from main.
 * @param argv Argument vector from main.
 * @param cell_structure_file Output: path of structure file.
 * @param input_file Output: path of input TOML.
 * @param output_file Output: path of output data file.
 * @param spin_structure_file_prefix Output: prefix for spin output files.
 * @return Always 0.
 */
int ReadOptions(int argc, char** argv, std::string & cell_structure_file, std::string & input_file, std::string & output_file, std::string & spin_structure_file_prefix) {
    // Process options from command line.
    char option;
    while ((option = getopt(argc, argv, "c:i:o:s:h")) != -1) {
        switch (option) {
            case 'c':
                cell_structure_file = optarg; break;
            case 'i':
                input_file = optarg; break;
            case 'o':
                output_file = optarg; break;
            case 's':
                spin_structure_file_prefix = optarg; break;
            case 'h':
                usage(argv[0]);
        }
    }

    return 0;
}

/**
 * @brief Print command-line usage and terminate.
 *
 * @param program Program name from argv[0].
 */
void usage(char* program) {
    std::cout
        << "Material-MC: Monte Carlo simulation for magnetic materials\n\n"
        << "Usage:\n"
        << "  " << program << " [options]\n\n"
        << "Options:\n"
        << "  -c <file>   Path to structure file (POSCAR format).\n"
        << "              Default: POSCAR\n"
        << "  -i <file>   Path to simulation parameter file (TOML).\n"
        << "              Default: input.toml\n"
        << "  -o <file>   Path to thermodynamic output file.\n"
        << "              Default: output.txt\n"
        << "  -s <prefix> Prefix for spin-structure output files.\n"
        << "              Default: spin\n"
        << "  -h          Print this help message and exit.\n\n"
        << "Examples:\n"
        << "  mpirun -np 8 " << program << "\n"
        << "  mpirun -np 8 " << program
        << " -c example/Ba2NiIrO6/POSCAR -i example/Ba2NiIrO6/input.toml"
        << " -o Ba2NiIrO6_output.txt -s Ba2NiIrO6_spin\n\n"
        << "Notes:\n"
        << "  - Run with MPI (mpirun/mpiexec).\n"
        << "  - input.toml controls Monte Carlo, lattice, Hamiltonian, and initialization settings.\n"
        << std::endl;
    exit(0);
}

/**
 * @brief Build the supercell sites and initialize spin vectors.
 *
 * Initialization policy:
 * - (0,0,0) cell uses base-site spin_initialization scaled by spin_scaling.
 * - +a/+b/+c directions propagate spins by Euler rotations angleA/angleB/angleC.
 *
 * Per-site pointers to species-level parameters are attached after spin setup.
 *
 * @param supercell Runtime lattice/site container updated in-place.
 * @return Always 0.
 */
int EnlargeCell(Supercell & supercell) {
    // Enlarge the system with given number and initialize the spin.
    std::vector<std::vector<std::vector<Site>>> site1;
    std::vector<std::vector<Site>> site2;
    std::vector<Site> site3;
    Site site4;

    for(int i=0; i<supercell.lattice.n_x; i++) {
        supercell.site.push_back(site1);
        for(int j=0; j<supercell.lattice.n_y; j++) {
            supercell.site[i].push_back(site2);
            for(int k=0; k<supercell.lattice.n_z; k++) {
                supercell.site[i][j].push_back(site3);
                for(int l=0; l<supercell.base_site.number; l++) {
                    supercell.site[i][j][k].push_back(site4);
                    // Spin initialization is seeded at origin cell and propagated by axis rotations.
                    if(k == 0) {
                        if(j == 0) {
                            if(i == 0) {
                                supercell.site[i][j][k][l].spin = supercell.base_site.spin_initialization[l];
                                supercell.site[i][j][k][l].spin[0] *= supercell.base_site.spin_scaling[l];
                                supercell.site[i][j][k][l].spin[1] *= supercell.base_site.spin_scaling[l];
                                supercell.site[i][j][k][l].spin[2] *= supercell.base_site.spin_scaling[l];
                            } else {
                                supercell.site[i][j][k][l].spin = Rotation(supercell.initialization.angleA, supercell.site[i-1][j][k][l].spin);
                            }
                        } else {
                            supercell.site[i][j][k][l].spin = Rotation(supercell.initialization.angleB, supercell.site[i][j-1][k][l].spin);
                        }
                    } else {
                        supercell.site[i][j][k][l].spin = Rotation(supercell.initialization.angleC, supercell.site[i][j][k-1][l].spin);
                    }
                    // Bind per-site references to species-level immutable parameters.
                    supercell.site[i][j][k][l].spin_scaling = & supercell.base_site.spin_scaling[l];
                    supercell.site[i][j][k][l].anisotropic_ratio = & supercell.base_site.anisotropic_ratio[l];
                    supercell.site[i][j][k][l].super_exchange_parameter = & supercell.base_site.super_exchange_parameter[l];
                    supercell.site[i][j][k][l].neighbor_number = & supercell.base_site.neighbor_number[l];
                }
            }
        }
    }

    return 0;
}

/**
 * @brief Compute squared Cartesian distance between two basis sites across a lattice translation.
 *
 * Formula uses fractional offset (index + base_site2 - base_site1), then maps it
 * to Cartesian coordinates by lattice vectors a/b/c, and returns squared norm.
 *
 * @return Squared distance in Cartesian units.
 */
double Distance(std::vector<double> lattice_constant1, std::vector<double> lattice_constant2, std::vector<double> lattice_constant3, \
std::vector<int> index, std::vector<double> base_site1, std::vector<double> base_site2) {
    // Return squared distance of index-base_site1+base_site2.
    std::vector<double> total_index = {index[0]+base_site2[0]-base_site1[0], index[1]+base_site2[1]-base_site1[1], index[2]+base_site2[2]-base_site1[2]};
    double result = 0;
    for(int i=0; i<3; i++) {
        result += (total_index[0]*lattice_constant1[i] + total_index[1]*lattice_constant2[i] + total_index[2]*lattice_constant3[i]) \
        * (total_index[0]*lattice_constant1[i] + total_index[1]*lattice_constant2[i] + total_index[2]*lattice_constant3[i]);
    }
    return result;
}

/**
 * @brief Insert a candidate shell distance into sorted distance_list if not duplicate.
 *
 * Distances within min(d1, d2) * tolerance_percentage are treated as the same shell.
 * The list is maintained in ascending order (zeros are treated as empty slots).
 *
 * @param distance Candidate squared distance.
 * @param distance_list In/out sorted shell distances.
 * @param tolerance_percentage Relative tolerance for duplicate detection.
 * @return Always 0.
 */
int AddDistance(double distance, std::vector<double> & distance_list, double tolerance_percentage) {
    int s = distance_list.size()-1;
    // Check similar distance
    for(int i=0; i<s+1; i++) {
        if(distance_list[i] == 0) {
            break;
        } else if (std::abs(distance_list[i]-distance) < std::min(distance, distance_list[i]) * tolerance_percentage){
            return 0;
        }
    }

    // Add new distance and order
    if(distance_list[s] == 0 || distance_list[s] > distance) {
        distance_list[s] = distance;
        for(int i=s-1; i>=0; i--) {
            if(distance_list[i] == 0 || distance_list[i] > distance) {
                distance_list[i+1] = distance_list[i];
                distance_list[i] = distance;
            }
        }
    }
    return 0;
}

/**
 * @brief Finalize runtime function bindings and neighbor links for the supercell.
 *
 * Main steps:
 * 1) Bind Hamiltonian function pointer from selected HamiltonianType.
 * 2) Bind local-update kernel from selected ModelType.
 * 3) Build neighbor index templates around each base site.
 * 4) Materialize neighbor pointers for each site with periodic boundary conditions.
 * 5) Compute initial total energy.
 *
 * @param supercell Runtime container updated in-place.
 * @return Always 0.
 */
int InitializeSupercell(Supercell & supercell) {
    // Initialize Hamiltonian
    supercell.HamiltonianBase = Heisenberg_base;
    switch (supercell.lattice.hamiltonian_type) {
        case HamiltonianType::Heisenberg :
            supercell.Hamiltonian = Heisenberg;
            break;
        case HamiltonianType::Heisenberg_with_field :
            supercell.Hamiltonian = Heisenberg_with_field;
            break;
        case HamiltonianType::Heisenberg_x_anisotropy :
            supercell.Hamiltonian = Heisenberg_x_anisotropy;
            break;
        case HamiltonianType::Heisenberg_x_anisotropy_with_field :
            supercell.Hamiltonian = Heisenberg_x_anisotropy_with_field;
            break;
        case HamiltonianType::Heisenberg_y_anisotropy :
            supercell.Hamiltonian = Heisenberg_y_anisotropy;
            break;
        case HamiltonianType::Heisenberg_y_anisotropy_with_field :
            supercell.Hamiltonian = Heisenberg_y_anisotropy_with_field;
            break;
        case HamiltonianType::Heisenberg_z_anisotropy :
            supercell.Hamiltonian = Heisenberg_z_anisotropy;
            break;
        case HamiltonianType::Heisenberg_z_anisotropy_with_field :
            supercell.Hamiltonian = Heisenberg_z_anisotropy_with_field;
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy :
            supercell.Hamiltonian = Heisenberg_xy_anisotropy;
            break;
        case HamiltonianType::Heisenberg_xy_anisotropy_with_field :
            supercell.Hamiltonian = Heisenberg_xy_anisotropy_with_field;
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy :
            supercell.Hamiltonian = Heisenberg_yz_anisotropy;
            break;
        case HamiltonianType::Heisenberg_yz_anisotropy_with_field :
            supercell.Hamiltonian = Heisenberg_yz_anisotropy_with_field;
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy :
            supercell.Hamiltonian = Heisenberg_zx_anisotropy;
            break;
        case HamiltonianType::Heisenberg_zx_anisotropy_with_field :
            supercell.Hamiltonian = Heisenberg_zx_anisotropy_with_field;
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy :
            supercell.Hamiltonian = Heisenberg_xyz_anisotropy;
            break;
        case HamiltonianType::Heisenberg_xyz_anisotropy_with_field :
            supercell.Hamiltonian = Heisenberg_xyz_anisotropy_with_field;
            break;
        case HamiltonianType::Heisenberg_custom :
            supercell.Hamiltonian = Hamiltonian_custom;
            supercell.HamiltonianBase = Hamiltonian_custom_base;
        default:
            break;
    }

    // Initialize local spin-update kernel by model type.
    switch (supercell.lattice.model_type) {
        case ModelType::Heisenberg :
            supercell.Update = LocalUpdateHeisenberg;
            break;
        case ModelType::Ising :
            supercell.Update = LocalUpdateIsing;
            break;
        default:
            supercell.Update = LocalUpdateHeisenberg;
            break;
    }

    // Build neighbor templates in base-cell coordinates first.
    // neighbors_index[site][shell] -> list of {dx, dy, dz, base_site_id} entries.
    std::vector<std::vector<std::vector<std::vector<int>>>> neighbors_index;

    for(int i=0; i<supercell.base_site.number; i++) {
        std::vector<std::vector<std::vector<int>>> index1;
        neighbors_index.push_back(index1);
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            std::vector<std::vector<int>> index2;
            neighbors_index[i].push_back(index2);
        }
    }

    for(int i=0; i<supercell.base_site.number; i++) {
        double distance_square = 0;

        // Enumerate translation window and collect matched neighbors by (element, distance).
        for(int j=-supercell.base_site.neighbor_number[i]; j<supercell.base_site.neighbor_number[i]+1; j++) {
            for(int k=-supercell.base_site.neighbor_number[i]; k<supercell.base_site.neighbor_number[i]+1; k++) {
                for(int l=-supercell.base_site.neighbor_number[i]; l<supercell.base_site.neighbor_number[i]+1; l++) {
                    for(int m=0; m<supercell.base_site.number; m++) {
                        distance_square = Distance(supercell.lattice.a, supercell.lattice.b, supercell.lattice.c, {j, k, l}, \
                        supercell.base_site.coordinate[i], supercell.base_site.coordinate[m]);

                        if(distance_square == 0.0) {
                            continue;
                        } else {
                            for(int n=0; n<supercell.base_site.neighbor_number[i]; n++) {
                                if(supercell.base_site.elements[m] == supercell.base_site.neighbor_elements[i][n] \
                                && abs(distance_square - supercell.base_site.neighbor_distance_square[i][n]) \
                                < supercell.lattice.tolerance_percentage*std::min(distance_square, supercell.base_site.neighbor_distance_square[i][n])) {
                                    // Matched one configured shell for base site i.
                                    std::vector<int> ind = {j, k, l, m};
                                    neighbors_index[i][n].emplace_back(ind);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Materialize neighbor pointers for every real site using periodic wrapping.
    for(int i=0; i<supercell.lattice.n_x; i++) {
        for(int j=0; j<supercell.lattice.n_y; j++) {
            for(int k=0; k<supercell.lattice.n_z; k++) {
                for(int l=0; l<supercell.base_site.number; l++) {
                    std::vector<Site*> temp = {};
                    for(int m=0; m<supercell.base_site.neighbor_number[l]; m++) {
                        supercell.site[i][j][k][l].neighbor.push_back(temp);
                        for(int n=0; n<neighbors_index[l][m].size(); n++) {
                            supercell.site[i][j][k][l].neighbor[m].push_back( 
                            & supercell.site[(i+neighbors_index[l][m][n][0]%supercell.lattice.n_x+supercell.lattice.n_x) % supercell.lattice.n_x] \
                            [(j+neighbors_index[l][m][n][1]%supercell.lattice.n_y+supercell.lattice.n_y) % supercell.lattice.n_y] \
                            [(k+neighbors_index[l][m][n][2]%supercell.lattice.n_z+supercell.lattice.n_z) % supercell.lattice.n_z] \
                            [neighbors_index[l][m][n][3]]);
                        }
                    }
                }
            }
        }
    }

    // Cache initial total energy after topology/function bindings are ready.
    supercell.lattice.total_energy = supercell.energy();
    return 0;
}