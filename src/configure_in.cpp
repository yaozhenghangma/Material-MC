#include "configure_in.h"

namespace {

void TryReadHamiltonianNumber(const toml::node_view<const toml::node>& node,
                              const std::string& key,
                              const std::string& input_file,
                              double& value) {
    if(!node) {
        return;
    }

    if(const auto parsed = node.value<double>(); parsed.has_value()) {
        value = *parsed;
        return;
    }

    std::cerr << "Invalid numeric value for [Hamiltonian]." << key
              << " in " << input_file
              << ": expected integer or floating-point value, got "
              << node.type() << ".\n";
    exit(-1);
}

bool IsKhModelToken(const std::string& model_name) {
    return model_name == "Kitaev-Heisenberg" || model_name == "kitaev-heisenberg"
        || model_name == "KH" || model_name == "kh";
}

char NormalizeKhDirectionLabel(const std::string& raw_label,
                               const std::string& key_path,
                               const std::string& input_file) {
    if(raw_label.size() != 1) {
        std::cerr << "Invalid direction label for " << key_path
                  << " in " << input_file
                  << ": expected one of x/y/z, got \"" << raw_label << "\".\n";
        exit(-1);
    }

    char label = static_cast<char>(std::tolower(static_cast<unsigned char>(raw_label[0])));
    if(label != 'x' && label != 'y' && label != 'z') {
        std::cerr << "Invalid direction label for " << key_path
                  << " in " << input_file
                  << ": expected one of x/y/z, got \"" << raw_label << "\".\n";
        exit(-1);
    }
    return label;
}

void ReadKhBondTypeDirectionMapping(Supercell& supercell,
                                    const toml::table& data,
                                    const std::string& input_file,
                                    bool kh_model_selected) {
    // Keep defaults explicit for non-KH paths.
    supercell.base_site.kh_bond_type_direction = {'x', 'y', 'z'};

    if(!kh_model_selected) {
        return;
    }

    auto hamiltonian = data["Hamiltonian"];
    if(!hamiltonian) {
        std::cerr << "Missing required table [Hamiltonian] for KH model in "
                  << input_file << ".\n";
        exit(-1);
    }

    auto mapping_node = hamiltonian["BondTypeDirection"];
    if(!mapping_node) {
        std::cerr << "Missing required table [Hamiltonian.BondTypeDirection] for KH model in "
                  << input_file << ".\n";
        exit(-1);
    }

    auto mapping = mapping_node.as_table();
    if(mapping == nullptr) {
        std::cerr << "Invalid type for [Hamiltonian.BondTypeDirection] in "
                  << input_file << ": expected table.\n";
        exit(-1);
    }

    const std::array<std::string, 3> keys = {"type1", "type2", "type3"};
    supercell.base_site.kh_bond_type_direction.resize(keys.size());
    std::map<char, std::string> direction_to_key;

    for(size_t i=0; i<keys.size(); i++) {
        const std::string key_path = "[Hamiltonian.BondTypeDirection]." + keys[i];
        auto value_node = (*mapping)[keys[i]];
        if(!value_node) {
            std::cerr << "Missing required key " << key_path
                      << " in " << input_file << " for KH model.\n";
            exit(-1);
        }

        const auto value = value_node.value<std::string>();
        if(!value.has_value()) {
            std::cerr << "Invalid type for " << key_path
                      << " in " << input_file
                      << ": expected string x/y/z, got "
                      << value_node.type() << ".\n";
            exit(-1);
        }

        const char label = NormalizeKhDirectionLabel(*value, key_path, input_file);
        if(direction_to_key.find(label) != direction_to_key.end()) {
            std::cerr << "Duplicate KH direction label \"" << label << "\" in "
                      << key_path << " and " << direction_to_key[label]
                      << " in " << input_file
                      << ": type1/type2/type3 must map to unique x/y/z labels.\n";
            exit(-1);
        }

        direction_to_key[label] = key_path;
        supercell.base_site.kh_bond_type_direction[i] = label;
    }
}

void ReadKhGlobalCouplings(Supercell& supercell,
                           const toml::table& data,
                           const std::string& input_file,
                           bool kh_model_selected) {
    auto hamiltonian = data["Hamiltonian"];

    // Always reset to defaults to keep non-KH compatibility explicit.
    supercell.base_site.kh_j = 0.0;
    supercell.base_site.kh_k = 0.0;
    supercell.base_site.kh_g = 0.0;
    supercell.base_site.kh_gp = 0.0;

    if(!hamiltonian) {
        return;
    }

    // For non-KH models we intentionally ignore KH-only couplings,
    // preserving current behavior for existing inputs.
    if(!kh_model_selected) {
        return;
    }

    TryReadHamiltonianNumber(hamiltonian["J"], "J", input_file, supercell.base_site.kh_j);
    TryReadHamiltonianNumber(hamiltonian["K"], "K", input_file, supercell.base_site.kh_k);
    TryReadHamiltonianNumber(hamiltonian["G"], "G", input_file, supercell.base_site.kh_g);
    TryReadHamiltonianNumber(hamiltonian["Gp"], "Gp", input_file, supercell.base_site.kh_gp);
}

} // namespace

/**
 * @brief Select the built-in Hamiltonian variant from field/anisotropy switches.
 *
 * Selection rules:
 * - B = (0, 0, 0) uses variants without external-field terms.
 * - B != (0, 0, 0) uses variants with external-field terms.
 * - Non-zero anisotropy components (x/y/z) decide the anisotropy suffix.
 *
 * @param supercell Runtime container; this function writes
 *        supercell.lattice.hamiltonian_type in-place.
 * @return Always 0.
 */
int ChooseHamiltonian(Supercell & supercell) {
    // No external magnetic field: choose among pure Heisenberg/anisotropy variants.
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
        // External magnetic field enabled: choose corresponding "with_field" variants.
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

/**
 * @brief Parse input.toml and populate MonteCarlo/Supercell runtime parameters.
 *
 * This function is the single source of truth for user-configurable simulation
 * settings. It reads each TOML section and writes values into runtime structures.
 * Missing keys fall back to value_or(...) defaults.
 *
 * Parsed sections:
 * - [MonteCarlo]
 * - [Lattice]
 * - [[Elements]] and [[Elements.Neighbors]]
 * - [Hamiltonian]
 * - [Output]
 * - [Initialization]
 *
 * @param supercell   Lattice/model container updated in-place.
 * @param monte_carlo Monte Carlo control parameters updated in-place.
 * @param input_file  Path to TOML configuration file.
 * @return Always 0 on success. Exits process on parse failure.
 */
int ReadSettingFile(Supercell & supercell, MonteCarlo & monte_carlo, std::string input_file) {
    toml::table data;
    try {
        // Parse input TOML file.
        data = toml::parse_file(input_file);

        // [MonteCarlo] simulation schedule and sampling controls.
        monte_carlo.start_temperature = data["MonteCarlo"]["start_temperature"].value_or(0.0);
        monte_carlo.end_temperature = data["MonteCarlo"]["end_temperature"].value_or(0.0);
        monte_carlo.temperature_step_number = data["MonteCarlo"]["temperature_points_number"].value_or(1);
        // Temperature increment derived from inclusive [start, end] sampling.
        monte_carlo.temperature_step = (monte_carlo.end_temperature-monte_carlo.start_temperature) / (monte_carlo.temperature_step_number-1);
        monte_carlo.relax_step = data["MonteCarlo"]["relaxing_steps"].value_or(1);
        monte_carlo.count_step = data["MonteCarlo"]["counting_steps"].value_or(1);
        monte_carlo.flip_number = data["MonteCarlo"]["flipping_number"].value_or(1);
        std::string tmp_string = data["MonteCarlo"]["method"].value_or("classical");
        // Method is selected by first character for compact input compatibility.
        if(tmp_string[0] == 'c' || tmp_string[0] == 'C') {
            monte_carlo.methods = Methods::classical;
        } else if(tmp_string[0] == 'p' || tmp_string[0] == 'P') {
            monte_carlo.methods = Methods::parallel_tempering;
        } else {
            monte_carlo.methods = Methods::classical;
        }
        monte_carlo.replica_exchange_step_number = data["MonteCarlo"]["exchange_step"].value_or(1);

        // [Lattice] supercell replication and neighbor matching tolerance.
        supercell.lattice.n_x = data["Lattice"]["cell_number"][0].value_or(1);
        supercell.lattice.n_y = data["Lattice"]["cell_number"][1].value_or(1);
        supercell.lattice.n_z = data["Lattice"]["cell_number"][2].value_or(1);
        supercell.lattice.tolerance_percentage = data["Lattice"]["tolerance"].value_or(0.01);
        
        // [[Elements]] magnetic species definition and shell interactions.
        toml::array& elements = *data.get_as<toml::array>("Elements");
        for(int i=0; i<elements.size(); i++) {
            supercell.base_site.elements.push_back(data["Elements"][i]["name"].value_or(""));
            supercell.base_site.spin_scaling.push_back(data["Elements"][i]["spin"].value_or(1.0));
            supercell.base_site.anisotropic_ratio.push_back({0.0, 0.0, 0.0});
            supercell.base_site.anisotropic_ratio[i][0] = data["Elements"][i]["anisotropic_factor"][0].value_or(1.0);
            supercell.base_site.anisotropic_ratio[i][1] = data["Elements"][i]["anisotropic_factor"][1].value_or(1.0);
            supercell.base_site.anisotropic_ratio[i][2] = data["Elements"][i]["anisotropic_factor"][2].value_or(1.0);
            // Per-element neighbor shells: (neighbor type, distance, J).
            auto neighbors = data["Elements"][i]["Neighbors"].as_array();
            supercell.base_site.neighbor_number.push_back(neighbors->size());
            supercell.base_site.neighbor_elements.push_back({});
            supercell.base_site.neighbor_distance_square.push_back({});
            supercell.base_site.super_exchange_parameter.push_back({});
            for(int j=0; j<neighbors->size(); j++) {
                supercell.base_site.neighbor_elements[i].push_back(data["Elements"][i]["Neighbors"][j]["name"].value_or(""));
                // Store squared distance for later squared-norm neighbor checks.
                supercell.base_site.neighbor_distance_square[i].push_back(data["Elements"][i]["Neighbors"][j]["distance"].value_or(1.0) * data["Elements"][i]["Neighbors"][j]["distance"].value_or(1.0));
                supercell.base_site.super_exchange_parameter[i].push_back(data["Elements"][i]["Neighbors"][j]["exchange_parameter"].value_or(1.0));
            }
        }
        
        // [Hamiltonian] external field, anisotropy, and Hamiltonian/model dispatch.
        supercell.base_site.B[0] = data["Hamiltonian"]["magnetic_field"][0].value_or(0.0) * MuB;
        supercell.base_site.B[1] = data["Hamiltonian"]["magnetic_field"][1].value_or(0.0) * MuB;
        supercell.base_site.B[2] = data["Hamiltonian"]["magnetic_field"][2].value_or(0.0) * MuB;
        if(supercell.base_site.B[2] != 0 || supercell.base_site.B[1] != 0 || supercell.base_site.B[0] != 0) {
            supercell.lattice.field = true;
            // Cache normalized field direction for projected observables.
            supercell.lattice.normalize_direction(supercell.base_site.B);
        }
        supercell.base_site.anisotropy[0] = data["Hamiltonian"]["anisotropy"][0].value_or(0.0);
        supercell.base_site.anisotropy[1] = data["Hamiltonian"]["anisotropy"][1].value_or(0.0);
        supercell.base_site.anisotropy[2] = data["Hamiltonian"]["anisotropy"][2].value_or(0.0);
        if(data["Hamiltonian"]["custom"].value_or(false)) {
            // Use user-implemented Hamiltonian in custom/.
            supercell.lattice.hamiltonian_type = HamiltonianType::Heisenberg_custom;
        } else {
            // Auto-select built-in Hamiltonian by (field, anisotropy).
            ChooseHamiltonian(supercell);
        }
        // Select spin model family (Heisenberg, Ising, or Kitaev-Heisenberg).
        tmp_string = data["Hamiltonian"]["model"].value_or("Heisenberg");
        bool kh_model_selected = IsKhModelToken(tmp_string);
        if(!tmp_string.empty() && (tmp_string[0] == 'I' || tmp_string[0] == 'i')) {
            supercell.lattice.model_type = ModelType::Ising;
        } else if(kh_model_selected) {
            supercell.lattice.model_type = ModelType::Kitaev_Heisenberg;
        } else {
            supercell.lattice.model_type = ModelType::Heisenberg;
        }

        // KH global couplings are only consumed when KH model is selected.
        // Missing values default to 0; invalid numeric types fail fast with context.
        ReadKhGlobalCouplings(supercell, data, input_file, kh_model_selected);

        // KH bond-type to direction mapping is required for KH model only.
        ReadKhBondTypeDirectionMapping(supercell, data, input_file, kh_model_selected);

        // [Output] visualization/export controls.
        supercell.lattice.magnify_factor = data["Output"]["magnifying_factor"].value_or(1.0);
        supercell.lattice.ground_state = data["Output"]["ground_state"].value_or(false);

        // [Initialization] spin initialization angles (degrees -> radians).
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

        // Optional explicit per-atom spin directions for the primitive cell.
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