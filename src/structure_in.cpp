#include "structure_in.h"

namespace {

bool ScanDouble(const std::string& line, double& value) {
    if (auto result = scn::scan<double>(line, "{}")) {
        value = result->value();
        return true;
    }
    return false;
}

bool ScanTripleSequential(const std::string& line, std::vector<double>& values) {
    if (auto result = scn::scan<double, double, double>(line, "{} {} {}")) {
        const auto parsed = result->values();
        values[0] = std::get<0>(parsed);
        values[1] = std::get<1>(parsed);
        values[2] = std::get<2>(parsed);
        return true;
    }
    return false;
}

bool ScanTriplePositional(const std::string& line, std::vector<double>& values) {
    if (auto result = scn::scan<double, double, double>(line, "{0} {1} {2}")) {
        const auto parsed = result->values();
        values[0] = std::get<0>(parsed);
        values[1] = std::get<1>(parsed);
        values[2] = std::get<2>(parsed);
        return true;
    }
    return false;
}

} // namespace

/**
 * @file structure_in.cpp
 * @brief POSCAR reader that loads lattice vectors and magnetic base sites.
 */

/**
 * @brief Parses a POSCAR file and populates lattice/base-site information.
 *
 * The parser reads lattice scaling and vectors, element names/counts, then
 * keeps only configured magnetic species in @p supercell.base_site.
 *
 * @param supercell Supercell object to populate.
 * @param cell_structure_file POSCAR file path.
 * @return int Returns 0 on completion.
 */
int ReadPOSCAR(Supercell & supercell, std::string cell_structure_file) {
    // Read lattice and site information from POSCAR.
    std::string str;
    static constexpr auto pattern = ctll::fixed_string{ R"(\s+)" };
    std::ifstream in;
    in.open(cell_structure_file, std::ios::in);

    getline(in, str); // Comment line.
    getline(in, str); // Lattice constant.
    ScanDouble(str, supercell.lattice.scaling);

    // Vector a, b and c.
    getline(in, str); 
    ScanTripleSequential(str, supercell.lattice.a);
    getline(in, str); 
    ScanTripleSequential(str, supercell.lattice.b);
    getline(in, str); 
    ScanTripleSequential(str, supercell.lattice.c);

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

    // Keep only magnetic species configured in supercell.base_site and copy species-level parameters.
    std::vector<double> tmp_coordinate = {0, 0, 0};
    supercell.base_site.number = 0;
    getline(in, str); // Coordinate mode line (Direct/Cartesian); current parser assumes fractional input.
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
                // Find optional explicit spin initialization entries for this magnetic element.
                int initialization_index = -1;
                for(int j=0; j<supercell.initialization.elements.size(); j++) {
                    if(magnetic_elements[k] == supercell.initialization.elements[j]) {
                        initialization_index = j;
                        break;
                    }
                }

                for(int j=0; j<elements_number[i]; j++) {
                    getline(in, str);
                    ScanTriplePositional(str, tmp_coordinate);
                    supercell.base_site.coordinate.push_back(tmp_coordinate);
                    supercell.base_site.elements.push_back(elements[i]);
                    supercell.base_site.spin_scaling.push_back(spin_scaling[k]);
                    supercell.base_site.anisotropic_ratio.push_back(anisotropic_ratio[k]);
                    supercell.base_site.neighbor_number.push_back(neighbor_number[k]);
                    supercell.base_site.neighbor_elements.push_back(neighbor_elements[k]);
                    supercell.base_site.neighbor_distance_square.push_back(neighbor_distance_square[k]);
                    supercell.base_site.super_exchange_parameter.push_back(super_exchange_parameter[k]);
                    // Use per-atom initialization if provided; otherwise fall back to +X.
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
                // Non-magnetic species are skipped while preserving POSCAR line consumption.
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
