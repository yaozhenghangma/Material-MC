#include "structure_in.h"

int ReadPOSCAR(Supercell & supercell, std::string cell_structure_file) {
    // Read information about base and lattice from POSCAR(default).
    std::string str;
    static constexpr auto pattern = ctll::fixed_string{ R"(\s+)" };
    std::ifstream in;
    in.open(cell_structure_file, std::ios::in);

    getline(in, str); // Comment line.
    getline(in, str); // Lattice constant.
    scn::scan(str, "{}", supercell.lattice.scaling);

    // Vector a, b and c.
    getline(in, str); 
    scn::scan(str, "{} {} {}", supercell.lattice.a[0], supercell.lattice.a[1], supercell.lattice.a[2]);
    getline(in, str); 
    scn::scan(str, "{} {} {}", supercell.lattice.b[0], supercell.lattice.b[1], supercell.lattice.b[2]);
    getline(in, str); 
    scn::scan(str, "{} {} {}", supercell.lattice.c[0], supercell.lattice.c[1], supercell.lattice.c[2]);

    // All elements.
    std::vector<std::string> elements;
    getline(in, str);
    std::string tmp_string;
    for(auto e: ctre::split<pattern>(str)) {
        tmp_string = std::string(e.get<0>());
        if (tmp_string != "") {
            elements.push_back(tmp_string);
        }
    }

    // All elements' number.
    std::vector<int> elements_number;
    getline(in, str);
    for(auto e: ctre::split<pattern>(str)) {
        tmp_string = std::string(e.get<0>());
        if (tmp_string != "") {
            elements_number.push_back(stoi(tmp_string));
        }
    }

    // Store magnetic elements
    std::vector<double> tmp_coordinate = {0, 0, 0};
    supercell.base_site.number = 0;
    getline(in, str); // Direct of Carsitian TODO:
    {
        std::vector<std::string> magnetic_elements = supercell.base_site.elements;
        std::vector<double> spin_scaling = supercell.base_site.spin_scaling;
        std::vector<std::vector<double>> anisotropic_ratio = supercell.base_site.anisotropic_ratio;
        std::vector<int> neighbor_number = supercell.base_site.neighbor_number;
        std::vector<std::vector<std::string>> neighbor_elements = supercell.base_site.neighbor_elements;
        std::vector<std::vector<double>> neighbor_distance_square = supercell.base_site.neighbor_distance_square;
        std::vector<std::vector<double>> super_exchange_parameter = supercell.base_site.super_exchange_parameter;
        supercell.base_site.elements = {};
        supercell.base_site.coordinate = {};
        supercell.base_site.spin_scaling = {};
        supercell.base_site.anisotropic_ratio = {};
        supercell.base_site.neighbor_number = {};
        supercell.base_site.neighbor_elements = {};
        supercell.base_site.neighbor_distance_square = {};
        supercell.base_site.super_exchange_parameter = {};
        supercell.base_site.spin_initialization = {};
        int k=0;
        for(int i=0; i<elements.size(); i++) {
            if(magnetic_elements[k] == elements[i]) {
                int initialization_index = -1;
                for(int j=0; j<supercell.initialization.elements.size(); j++) {
                    if(magnetic_elements[k] == supercell.initialization.elements[j]) {
                        initialization_index = j;
                        break;
                    }
                }

                for(int j=0; j<elements_number[i]; j++) {
                    getline(in, str);
                    scn::scan(str, "{0} {1} {2}", tmp_coordinate[0], tmp_coordinate[1], tmp_coordinate[2]);
                    supercell.base_site.coordinate.push_back(tmp_coordinate);
                    supercell.base_site.elements.push_back(elements[i]);
                    supercell.base_site.spin_scaling.push_back(spin_scaling[k]);
                    supercell.base_site.anisotropic_ratio.push_back(anisotropic_ratio[k]);
                    supercell.base_site.neighbor_number.push_back(neighbor_number[k]);
                    supercell.base_site.neighbor_elements.push_back(neighbor_elements[k]);
                    supercell.base_site.neighbor_distance_square.push_back(neighbor_distance_square[k]);
                    supercell.base_site.super_exchange_parameter.push_back(super_exchange_parameter[k]);
                    if(initialization_index != -1 && j+initialization_index < supercell.initialization.elements.size() \
                    && supercell.initialization.elements[j+initialization_index] == magnetic_elements[k]) {
                        supercell.base_site.spin_initialization.push_back(supercell.initialization.direction[j+initialization_index]);
                    } else {
                        supercell.base_site.spin_initialization.push_back({1, 0, 0});
                    }
                }
                supercell.base_site.number += elements_number[i];
                k++;
            } else {
                for(int j=0; j<elements_number[i]; j++) {
                    getline(in, str);
                }
            }

            if(k>=magnetic_elements.size()) {
                break;
            }
        }
    }
    
    in.close();
    return 0;
}