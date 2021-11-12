#include "configure_in.h"

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
    getline(in, str); // Hamiltonion function.
    scn::scan(str, "{}", supercell.lattice.function_choice);
    getline(in, str); // Anisotropic factor.
    if(supercell.lattice.function_choice == "Heisenberg") {
        std::cout << "Hamiltonion type: Heisenberg" << std::endl;
        scn::scan(str, "{}", supercell.base_site.anisotropic_factor_D);
    } else if(supercell.lattice.function_choice == "Heisenberg_3_axis_anisotropy") {
        std::cout << "Hamiltonion type: Heisenberg with 3 axis anisotropy" << std::endl;
        scn::scan(str, "{} {}", supercell.base_site.anisotropic_factor_D, supercell.base_site.anisotropic_factor_En);
    } else if(supercell.lattice.function_choice == "Heisenberg_external_field") {
        std::cout << "Hamiltonion type: Heisenberg with 3 axis anisotropy and external magnetic field" << std::endl;
        scn::scan(str, "{} {} {}", supercell.base_site.B[0], supercell.base_site.B[1], supercell.base_site.B[2]);
        supercell.base_site.B[0] *= MuB;
        supercell.base_site.B[1] *= MuB;
        supercell.base_site.B[2] *= MuB;
        getline(in, str);
        scn::scan(str, "{} {}", supercell.base_site.anisotropic_factor_D, supercell.base_site.anisotropic_factor_En);
    } else {
        std::cout << "Default Hamiltonion type: Heisenberg" << std::endl;
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