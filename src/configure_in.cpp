#include "configure_in.h"

/*

int ReadSettingFile(Supercell & supercell, MonteCarlo & monte_carlo, std::string input_file) {
    // Read information about enlarging and Monte Carlo from given setting file.
    std::string str;
    std::string tmp_str;
    std::ifstream in;
    in.open(input_file, std::ios::in);

    // Comment line of the file
    getline(in, str);

    // Information to control Monte Carlo simulation.
    getline(in, str); // Comment line of the simulation.
    getline(in, str); // Temperature.
    scn::scan(str, "{} {} {}", monte_carlo.start_temperature, monte_carlo.end_temperature, monte_carlo.temperature_step_number);
    monte_carlo.temperature_step = (monte_carlo.end_temperature-monte_carlo.start_temperature) / (monte_carlo.temperature_step_number-1);
    getline(in, str); // Monte Carlo steps for relaxing process.
    scn::scan(str, "{}", monte_carlo.relax_step);
    getline(in, str); // Monte Carlo steps for counting.
    scn::scan(str, "{}", monte_carlo.count_step);
    getline(in, str); // Flip number for one Monte Carlo step.
    scn::scan(str, "{}", monte_carlo.flip_number);

    // Information about lattice.
    getline(in, str); // Comment line
    getline(in, str); // Tolerance percentage
    scn::scan(str, "{}", supercell.lattice.tolerance_percentage);
    getline(in, str); // Number of cells.
    scn::scan(str, "{} {} {}", supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z);
    getline(in, str); // Hamiltonian function.
    scn::scan(str, "{}", supercell.lattice.function_choice);
    getline(in, str); // Anisotropic factor.
    if(supercell.lattice.function_choice == "Heisenberg") {
        std::cout << "Hamiltonian type: Heisenberg" << std::endl;
        scn::scan(str, "{}", supercell.base_site.anisotropic_factor_D);
    } else if(supercell.lattice.function_choice == "Heisenberg_3_axis_anisotropy") {
        std::cout << "Hamiltonian type: Heisenberg with 3 axis anisotropy" << std::endl;
        scn::scan(str, "{} {}", supercell.base_site.anisotropic_factor_D, supercell.base_site.anisotropic_factor_En);
    } else if(supercell.lattice.function_choice == "Heisenberg_external_field") {
        std::cout << "Hamiltonian type: Heisenberg with 3 axis anisotropy and external magnetic field" << std::endl;
        scn::scan(str, "{} {} {}", supercell.base_site.B[0], supercell.base_site.B[1], supercell.base_site.B[2]);
        supercell.base_site.B[0] *= MuB;
        supercell.base_site.B[1] *= MuB;
        supercell.base_site.B[2] *= MuB;
        getline(in, str);
        scn::scan(str, "{} {}", supercell.base_site.anisotropic_factor_D, supercell.base_site.anisotropic_factor_En);
    } else {
        std::cout << "Default Hamiltonian type: Heisenberg" << std::endl;
        scn::scan(str, "{}", supercell.base_site.anisotropic_factor_D);
    }
    
    getline(in, str); // Magnifying factor
    scn::scan(str, "{}", supercell.lattice.magnify_factor);

    // Information about base.
    getline(in, str); // Comment line
    getline(in, str); // Magnetic elements
    static constexpr auto pattern = ctll::fixed_string{ R"(\s+)" };
    for(auto e: ctre::split<pattern>(str)) {
        supercell.base_site.elements.push_back(std::string(e.get<0>()));
    }
    if(supercell.base_site.elements[0] == "All" || supercell.base_site.elements[0] == "all") {
        supercell.base_site.all_magnetic = true;
        supercell.base_site.elements = {};
    } else {
        supercell.base_site.all_magnetic = false;
        int i=0;
        for(; i<supercell.base_site.elements.size(); i++) {
            if(supercell.base_site.elements[i][0] == '#') {
                break;
            }
        }
        int size = supercell.base_site.elements.size();
        for(; i<size; i++) {
            supercell.base_site.elements.pop_back();
        }
    }
    std::vector<double> tmp_vector;
    std::vector<std::string> tmp_string_vector;
    int i=0;
    while(getline(in, str) && !str.empty() && str[0] != '=') { // Super-exchange parameters
        supercell.base_site.neighbor_number.push_back(0);
        scn::scan(str, "{} {}", tmp_str, supercell.base_site.neighbor_number[i]);
        supercell.base_site.super_exchange_parameter.push_back(tmp_vector);
        supercell.base_site.neighbor_distance_square.push_back(tmp_vector);
        supercell.base_site.neighbor_elements.push_back(tmp_string_vector);
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            getline(in, str);
            supercell.base_site.super_exchange_parameter[i].push_back(0);
            supercell.base_site.neighbor_distance_square[i].push_back(0);

            supercell.base_site.neighbor_elements[i].push_back(" ");
            scn::scan(str, "{} {} {}", supercell.base_site.neighbor_elements[i][j], supercell.base_site.super_exchange_parameter[i][j], supercell.base_site.neighbor_distance_square[i][j]);
            supercell.base_site.neighbor_distance_square[i][j] = supercell.base_site.neighbor_distance_square[i][j] * supercell.base_site.neighbor_distance_square[i][j];
        }

        i++;
    }
    
    getline(in, str);
    scn::scan(str, "{}", tmp_str);
    if(tmp_str[0] == 'A' || tmp_str[0] == 'a') {
        supercell.initialization.anti_ferromagnetic = true;
    } else {
        supercell.initialization.anti_ferromagnetic = false;
    }
    getline(in, str);
    scn::scan(str, "{} {} {}", \
        supercell.initialization.direction[0], \
        supercell.initialization.direction[1], \
        supercell.initialization.direction[2]);
    if(supercell.initialization.anti_ferromagnetic) {
        for(int i=0; i<supercell.base_site.elements.size(); i++) {
            supercell.initialization.anti_ferromagnetic_J.push_back({});
            getline(in, str);
            for(auto e: ctre::split<pattern>(str)) {
                tmp_str = std::string(e.get<0>());
                if (tmp_str != "" && tmp_str[0] != '#') {
                    supercell.initialization.anti_ferromagnetic_J[i].push_back(stoi(tmp_str));
                } else {
                    break;
                }
            }
        }
    }

    in.close();
    return 0;
}
*/

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
        supercell.base_site.B[1] = data["Hamiltonian"]["magnetic_field"][1].value_or(0.0) * MuB;
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