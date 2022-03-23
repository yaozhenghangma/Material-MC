#include "initialization.h"

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

void usage(char* program) {
    //TODO: Usage of the program.
    exit(1);
}

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

    // Initialize update function
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

    // Initialize the neighbors' link and energy.
    // Find the neighbors for every base sites in a cell.
    // neighbors_index[i][j][k][l]: i-th site's k-th j nearest neighbor's index in l-th direction. 4-th direction is base number.
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

        // Find link.
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

    // Set link for every sites.
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

    supercell.lattice.total_energy = supercell.energy();
    return 0;
}