#include "configure_in.h"

int ChooseHamiltonian(Supercell & supercell) {
    if(supercell.base_site.B[0] == 0.0 && supercell.base_site.B[1] == 0.0 && supercell.base_site.B[2] == 0.0) {
        if(supercell.base_site.anisotropy[0] == 0.0) {
            if(supercell.base_site.anisotropy[1] == 0.0) {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_z_anisotropy;
                }
            } else {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_y_anisotropy;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_yz_anisotropy;
                }
            }
        } else {
            if(supercell.base_site.anisotropy[1] == 0.0) {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_x_anisotropy;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_zx_anisotropy;
                }
            } else {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_xy_anisotropy;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_xyz_anisotropy;
                }
            }
        }
    } else {
        if(supercell.base_site.anisotropy[0] == 0.0) {
            if(supercell.base_site.anisotropy[1] == 0.0) {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_with_field;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_z_anisotropy_with_field;
                }
            } else {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_y_anisotropy_with_field;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_yz_anisotropy_with_field;
                }
            }
        } else {
            if(supercell.base_site.anisotropy[1] == 0.0) {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_x_anisotropy_with_field;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_zx_anisotropy_with_field;
                }
            } else {
                if(supercell.base_site.anisotropy[2] == 0.0) {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_xy_anisotropy_with_field;
                } else {
                    supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_xyz_anisotropy_with_field;
                }
            }
        }
    }
    return 0;
}

int ReadSettingFile(Supercell & supercell, MonteCarlo & monte_carlo, std::string input_file) {
    toml::table data;
    try {
        data = toml::parse_file(input_file);
        // Monte Carlo parameters
        monte_carlo.start_temperature = data["MonteCarlo"]["start_temperature"].value_or(0.0);
        monte_carlo.end_temperature = data["MonteCarlo"]["end_temperature"].value_or(0.0);
        monte_carlo.temperature_step_number = data["MonteCarlo"]["temperature_points_number"].value_or(1);
        monte_carlo.temperature_step = (monte_carlo.end_temperature-monte_carlo.start_temperature) / (monte_carlo.temperature_step_number-1);
        monte_carlo.relax_step = data["MonteCarlo"]["relaxing_steps"].value_or(1);
        monte_carlo.count_step = data["MonteCarlo"]["counting_steps"].value_or(1);
        monte_carlo.flip_number = data["MonteCarlo"]["flipping_number"].value_or(1);
        std::string tmp_string = data["MonteCarlo"]["method"].value_or("classical");
        if(tmp_string == "classical") {
            monte_carlo.methods = Methods::classical;
        } else if(tmp_string == "ptmc") {
            monte_carlo.methods = Methods::parallel_tempering;
        } else {
            monte_carlo.methods = Methods::classical;
        }
        monte_carlo.replica_exchange_step_number = data["MonteCarlo"]["exchange_step"].value_or(1);

        // Cell
        supercell.lattice.n_x = data["Lattice"]["cell_number"][0].value_or(1);
        supercell.lattice.n_y = data["Lattice"]["cell_number"][1].value_or(1);
        supercell.lattice.n_z = data["Lattice"]["cell_number"][2].value_or(1);
        supercell.lattice.tolerance_percentage = data["Lattice"]["tolerance"].value_or(0.01);
        
        // Magnetic elements
        toml::array& elements = *data.get_as<toml::array>("Elements");
        for(int i=0; i<elements.size(); i++) {
            supercell.base_site.elements.push_back(data["Elements"][i]["name"].value_or(""));
            supercell.base_site.spin_scaling.push_back(data["Elements"][i]["spin"].value_or(1.0));
            supercell.base_site.anisotropic_ratio.push_back({0.0, 0.0, 0.0});
            supercell.base_site.anisotropic_ratio[i][0] = data["Elements"][i]["anisotropic_factor"][0].value_or(1.0);
            supercell.base_site.anisotropic_ratio[i][1] = data["Elements"][i]["anisotropic_factor"][1].value_or(1.0);
            supercell.base_site.anisotropic_ratio[i][2] = data["Elements"][i]["anisotropic_factor"][2].value_or(1.0);
            auto neighbors = data["Elements"][i]["Neighbors"].as_array();
            supercell.base_site.neighbor_number.push_back(neighbors->size());
            supercell.base_site.neighbor_elements.push_back({});
            supercell.base_site.neighbor_distance_square.push_back({});
            supercell.base_site.super_exchange_parameter.push_back({});
            for(int j=0; j<neighbors->size(); j++) {
                supercell.base_site.neighbor_elements[i].push_back(data["Elements"][i]["Neighbors"][j]["name"].value_or(""));
                supercell.base_site.neighbor_distance_square[i].push_back(data["Elements"][i]["Neighbors"][j]["distance"].value_or(1.0) * data["Elements"][i]["Neighbors"][j]["distance"].value_or(1.0));
                supercell.base_site.super_exchange_parameter[i].push_back(data["Elements"][i]["Neighbors"][j]["exchange_parameter"].value_or(1.0));
            }
        }
        
        // Hamiltonian function
        supercell.base_site.B[0] = data["Hamiltonian"]["magnetic_field"][0].value_or(0.0) * MuB;
        supercell.base_site.B[1] = data["Hamiltonian"]["magnetic_field"][1].value_or(0.0) * MuB;
        supercell.base_site.B[2] = data["Hamiltonian"]["magnetic_field"][2].value_or(0.0) * MuB;
        if(supercell.base_site.B[2] != 0 || supercell.base_site.B[1] != 0 || supercell.base_site.B[0] != 0) {
            supercell.lattice.field = true;
            supercell.lattice.normalize_direction(supercell.base_site.B);
        }
        supercell.base_site.anisotropy[0] = data["Hamiltonian"]["anisotropy"][0].value_or(0.0);
        supercell.base_site.anisotropy[1] = data["Hamiltonian"]["anisotropy"][1].value_or(0.0);
        supercell.base_site.anisotropy[2] = data["Hamiltonian"]["anisotropy"][2].value_or(0.0);
        if(data["Hamiltonian"]["custom"].value_or(false)) {
            supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_custom;
        } else {
            ChooseHamiltonian(supercell);
        }
        
        // Output
        supercell.lattice.magnify_factor = data["Output"]["magnifying_factor"].value_or(1.0);

        // Initialization
        double unit_convert = PI/180.0;
        supercell.initialization.angleA[0] = data["Initialization"]["angleA"][0].value_or(0.0) * unit_convert;
        supercell.initialization.angleA[1] = data["Initialization"]["angleA"][1].value_or(0.0) * unit_convert;
        supercell.initialization.angleA[2] = data["Initialization"]["angleA"][2].value_or(0.0) * unit_convert;
        supercell.initialization.angleB[0] = data["Initialization"]["angleB"][0].value_or(0.0) * unit_convert;
        supercell.initialization.angleB[1] = data["Initialization"]["angleB"][1].value_or(0.0) * unit_convert;
        supercell.initialization.angleB[2] = data["Initialization"]["angleB"][2].value_or(0.0) * unit_convert;
        supercell.initialization.angleC[0] = data["Initialization"]["angleC"][0].value_or(0.0) * unit_convert;
        supercell.initialization.angleC[1] = data["Initialization"]["angleC"][1].value_or(0.0) * unit_convert;
        supercell.initialization.angleC[2] = data["Initialization"]["angleC"][2].value_or(0.0) * unit_convert;
        
        auto initialization_elements = data["Initialization"]["Elements"].as_array();
        
        if(initialization_elements != nullptr) {
            for(int i=0; i<initialization_elements->size(); i++) {
                auto initialization_atoms = data["Initialization"]["Elements"][i]["Atoms"].as_array();
                for(int j=0; j<initialization_atoms->size(); j++) {
                    supercell.initialization.elements.push_back(data["Initialization"]["Elements"][i]["name"].value_or(""));
                    supercell.initialization.direction.push_back({\
                    data["Initialization"]["Elements"][i]["Atoms"][j]["spin_direction"][0].value_or(1.0),\
                    data["Initialization"]["Elements"][i]["Atoms"][j]["spin_direction"][1].value_or(1.0),\
                    data["Initialization"]["Elements"][i]["Atoms"][j]["spin_direction"][2].value_or(1.0)});
                }
            }
            supercell.initialization.normalized();
        }
        
    } catch (const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        exit(-1);
    }
    return 0;
}