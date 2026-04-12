#include "initialization.h"

namespace {

constexpr double kKhDirectionSortTolerance = 1e-9;

std::array<double, 3> BuildBondVectorCartesian(const Supercell& supercell,
                                               const std::vector<int>& neighbor_template,
                                               int source_base_site) {
    const std::vector<double>& source = supercell.base_site.coordinate[source_base_site];
    const std::vector<double>& target = supercell.base_site.coordinate[neighbor_template[3]];

    const std::vector<double> fractional = {
        static_cast<double>(neighbor_template[0]) + target[0] - source[0],
        static_cast<double>(neighbor_template[1]) + target[1] - source[1],
        static_cast<double>(neighbor_template[2]) + target[2] - source[2]
    };

    return {
        fractional[0]*supercell.lattice.a[0] + fractional[1]*supercell.lattice.b[0] + fractional[2]*supercell.lattice.c[0],
        fractional[0]*supercell.lattice.a[1] + fractional[1]*supercell.lattice.b[1] + fractional[2]*supercell.lattice.c[1],
        fractional[0]*supercell.lattice.a[2] + fractional[1]*supercell.lattice.b[2] + fractional[2]*supercell.lattice.c[2]
    };
}

double DotProduct(const std::array<double, 3>& left, const std::array<double, 3>& right) {
    return left[0]*right[0] + left[1]*right[1] + left[2]*right[2];
}

std::array<double, 3> NormalizeAndCanonicalizeDirection(const std::array<double, 3>& vector,
                                                         double direction_epsilon,
                                                         int base_site,
                                                         int shell,
                                                         int entry) {
    const double norm = std::sqrt(DotProduct(vector, vector));
    if(norm <= direction_epsilon) {
        std::cerr << "Invalid KH bond vector during geometry classification: "
                  << "zero-length vector at base site " << base_site
                  << ", shell " << shell
                  << ", entry " << entry << ".\n";
        exit(-1);
    }

    std::array<double, 3> direction = {
        vector[0] / norm,
        vector[1] / norm,
        vector[2] / norm
    };

    int sign = 1;
    for(int axis=0; axis<3; axis++) {
        if(std::abs(direction[axis]) > direction_epsilon) {
            if(direction[axis] < 0.0) {
                sign = -1;
            }
            break;
        }
    }

    if(sign < 0) {
        direction[0] = -direction[0];
        direction[1] = -direction[1];
        direction[2] = -direction[2];
    }

    return direction;
}

std::vector<char> ClassifyKhShellDirections(const Supercell& supercell,
                                            const std::vector<std::vector<int>>& shell_templates,
                                            int base_site,
                                            int shell) {
    const double direction_epsilon = supercell.base_site.kh_direction_epsilon;
    const double direction_tolerance = supercell.base_site.kh_direction_tolerance;

    if(direction_epsilon <= 0.0) {
        std::cerr << "Invalid KH direction epsilon in runtime data: expected > 0, got "
                  << direction_epsilon << ".\n";
        exit(-1);
    }

    if(direction_tolerance <= 0.0 || direction_tolerance >= 1.0) {
        std::cerr << "Invalid KH direction tolerance in runtime data: expected 0 < value < 1, got "
                  << direction_tolerance << ".\n";
        exit(-1);
    }

    if(shell_templates.empty()) {
        std::cerr << "Missing KH bonds for classification: "
                  << "base site " << base_site
                  << ", shell " << shell << " has no neighbors.\n";
        exit(-1);
    }

    if(supercell.base_site.kh_bond_type_direction.size() != 3) {
        std::cerr << "Invalid KH bond-type mapping size in runtime data: expected 3 labels, got "
                  << supercell.base_site.kh_bond_type_direction.size() << ".\n";
        exit(-1);
    }

    struct DirectionGroup {
        std::array<double, 3> axis;
        std::vector<int> entries;
    };

    std::vector<DirectionGroup> groups;
    std::vector<std::array<double, 3>> directions;
    directions.reserve(shell_templates.size());

    for(int entry=0; entry<shell_templates.size(); entry++) {
        std::array<double, 3> bond_vector = BuildBondVectorCartesian(supercell, shell_templates[entry], base_site);
        directions.push_back(NormalizeAndCanonicalizeDirection(
            bond_vector,
            direction_epsilon,
            base_site,
            shell,
            entry));
    }

    for(int entry=0; entry<directions.size(); entry++) {
        int matched_group = -1;
        for(int group=0; group<groups.size(); group++) {
            const double parallel_score = std::abs(DotProduct(directions[entry], groups[group].axis));
            if(parallel_score >= 1.0 - direction_tolerance) {
                if(matched_group != -1) {
                    std::cerr << "Ambiguous KH bond classification: "
                              << "base site " << base_site
                              << ", shell " << shell
                              << ", entry " << entry
                              << " matches multiple direction groups.\n";
                    exit(-1);
                }
                matched_group = group;
            }
        }

        if(matched_group == -1) {
            DirectionGroup new_group;
            new_group.axis = directions[entry];
            new_group.entries.push_back(entry);
            groups.push_back(new_group);
        } else {
            groups[matched_group].entries.push_back(entry);
        }
    }

    if(groups.size() != 3) {
        std::cerr << "Invalid KH bond classification for honeycomb scope: "
                  << "base site " << base_site
                  << ", shell " << shell
                  << " produced " << groups.size()
                  << " geometric bond types, expected exactly 3.\n";
        exit(-1);
    }

    std::vector<int> sorted_group_index(groups.size());
    std::iota(sorted_group_index.begin(), sorted_group_index.end(), 0);
    std::sort(sorted_group_index.begin(), sorted_group_index.end(),
              [&groups](int left, int right) {
                  for(int axis=0; axis<3; axis++) {
                      if(std::abs(groups[left].axis[axis] - groups[right].axis[axis]) > kKhDirectionSortTolerance) {
                          return groups[left].axis[axis] < groups[right].axis[axis];
                      }
                  }
                  return left < right;
              });

    std::vector<char> shell_labels(shell_templates.size(), '\0');
    for(int type_index=0; type_index<sorted_group_index.size(); type_index++) {
        const int group_index = sorted_group_index[type_index];
        const char mapped_label = supercell.base_site.kh_bond_type_direction[type_index];
        for(int entry : groups[group_index].entries) {
            shell_labels[entry] = mapped_label;
        }
    }

    for(int entry=0; entry<shell_labels.size(); entry++) {
        if(shell_labels[entry] == '\0') {
            std::cerr << "Incomplete KH bond classification: "
                      << "base site " << base_site
                      << ", shell " << shell
                      << ", entry " << entry
                      << " has no direction label.\n";
            exit(-1);
        }
    }

    return shell_labels;
}

std::vector<std::vector<std::vector<char>>> BuildKhNeighborDirectionTemplates(
    const Supercell& supercell,
    const std::vector<std::vector<std::vector<std::vector<int>>>>& neighbors_index) {
    std::vector<std::vector<std::vector<char>>> labels;
    labels.resize(neighbors_index.size());

    for(int base_site=0; base_site<neighbors_index.size(); base_site++) {
        labels[base_site].resize(neighbors_index[base_site].size());
        for(int shell=0; shell<neighbors_index[base_site].size(); shell++) {
            labels[base_site][shell] = ClassifyKhShellDirections(
                supercell,
                neighbors_index[base_site][shell],
                base_site,
                shell);
        }
    }

    return labels;
}

} // namespace

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
    if(supercell.lattice.model_type == ModelType::Kitaev_Heisenberg) {
        supercell.Hamiltonian = Kitaev_Heisenberg;
        supercell.HamiltonianBase = Kitaev_Heisenberg_base;
    } else {
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
    }

    // Initialize local spin-update kernel by model type.
    switch (supercell.lattice.model_type) {
        case ModelType::Heisenberg :
            supercell.Update = LocalUpdateHeisenberg;
            break;
        case ModelType::Ising :
            supercell.Update = LocalUpdateIsing;
            break;
        case ModelType::Kitaev_Heisenberg :
            // KH uses classical vector spins in this stage; keep Heisenberg local update.
            supercell.Update = LocalUpdateHeisenberg;
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
                                && std::abs(distance_square - supercell.base_site.neighbor_distance_square[i][n]) \
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

    // Pre-compute KH direction labels for each neighbor template entry.
    std::vector<std::vector<std::vector<char>>> kh_neighbor_direction_templates;
    if(supercell.lattice.model_type == ModelType::Kitaev_Heisenberg) {
        kh_neighbor_direction_templates = BuildKhNeighborDirectionTemplates(supercell, neighbors_index);
    }

    // Materialize neighbor pointers for every real site using periodic wrapping.
    for(int i=0; i<supercell.lattice.n_x; i++) {
        for(int j=0; j<supercell.lattice.n_y; j++) {
            for(int k=0; k<supercell.lattice.n_z; k++) {
                for(int l=0; l<supercell.base_site.number; l++) {
                    supercell.site[i][j][k][l].neighbor.clear();
                    supercell.site[i][j][k][l].neighbor_direction.clear();
                    std::vector<Site*> temp = {};
                    std::vector<char> temp_direction = {};
                    for(int m=0; m<supercell.base_site.neighbor_number[l]; m++) {
                        supercell.site[i][j][k][l].neighbor.push_back(temp);
                        supercell.site[i][j][k][l].neighbor_direction.push_back(temp_direction);
                        for(int n=0; n<neighbors_index[l][m].size(); n++) {
                            supercell.site[i][j][k][l].neighbor[m].push_back(
                            & supercell.site[(i+neighbors_index[l][m][n][0]%supercell.lattice.n_x+supercell.lattice.n_x) % supercell.lattice.n_x] \
                            [(j+neighbors_index[l][m][n][1]%supercell.lattice.n_y+supercell.lattice.n_y) % supercell.lattice.n_y] \
                            [(k+neighbors_index[l][m][n][2]%supercell.lattice.n_z+supercell.lattice.n_z) % supercell.lattice.n_z] \
                            [neighbors_index[l][m][n][3]]);

                            if(supercell.lattice.model_type == ModelType::Kitaev_Heisenberg) {
                                const char direction_label = kh_neighbor_direction_templates[l][m][n];
                                if(direction_label != 'x' && direction_label != 'y' && direction_label != 'z') {
                                    std::cerr << "Invalid KH direction label during neighbor materialization: "
                                              << "base site " << l
                                              << ", shell " << m
                                              << ", entry " << n
                                              << ", label \"" << direction_label << "\".\n";
                                    exit(-1);
                                }
                                supercell.site[i][j][k][l].neighbor_direction[m].push_back(direction_label);
                            } else {
                                supercell.site[i][j][k][l].neighbor_direction[m].push_back('\0');
                            }
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
