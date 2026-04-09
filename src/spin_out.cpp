#include "spin_out.h"

namespace {

std::array<double, 3> BuildCartesianFromFractional(const Supercell& supercell, const std::array<double, 3>& fractional) {
    return {
        (fractional[0]*supercell.lattice.a[0] + fractional[1]*supercell.lattice.b[0] + fractional[2]*supercell.lattice.c[0]) * supercell.lattice.magnify_factor,
        (fractional[0]*supercell.lattice.a[1] + fractional[1]*supercell.lattice.b[1] + fractional[2]*supercell.lattice.c[1]) * supercell.lattice.magnify_factor,
        (fractional[0]*supercell.lattice.a[2] + fractional[1]*supercell.lattice.b[2] + fractional[2]*supercell.lattice.c[2]) * supercell.lattice.magnify_factor
    };
}

std::array<double, 3> BuildSiteCartesian(const Supercell& supercell, int i, int j, int k, int l) {
    const std::array<double, 3> fractional = {
        supercell.base_site.coordinate[l][0] + static_cast<double>(i),
        supercell.base_site.coordinate[l][1] + static_cast<double>(j),
        supercell.base_site.coordinate[l][2] + static_cast<double>(k)
    };
    return BuildCartesianFromFractional(supercell, fractional);
}

std::array<int, 3> BondColorFromDirection(char direction) {
    if(direction == 'x') {
        return {255, 0, 0};
    }
    if(direction == 'y') {
        return {0, 255, 0};
    }
    return {0, 0, 255};
}

} // namespace

/**
 * @file spin_out.cpp
 * @brief Utilities for exporting spin configurations in XSF/VESTA format.
 */

/**
 * @brief Writes the current supercell spin configuration to an XSF file.
 *
 * Output filename: ``<spin_structure_file_prefix>.xsf``.
 *
 * File content follows standard XSF sections:
 * - ``CRYSTAL``
 * - ``PRIMVEC`` (scaled lattice vectors)
 * - ``PRIMCOORD`` (element, x, y, z, sx, sy, sz)
 *
 * @param supercell Supercell containing lattice, coordinates, and spin vectors.
 * @param spin_structure_file_prefix Prefix of the output XSF filename.
 * @return int Returns 0 on completion.
 */
int WriteSpin(Supercell & supercell, std::string spin_structure_file_prefix) {
    std::string output_file_name = spin_structure_file_prefix + ".xsf";
    auto out = fmt::output_file(output_file_name);
    out.print("CRYSTAL\n");
    out.print("PRIMVEC\n");
    out.print("{} {} {}\n", supercell.lattice.a[0]*supercell.lattice.n_x*supercell.lattice.magnify_factor, \
    supercell.lattice.a[1]*supercell.lattice.n_x*supercell.lattice.magnify_factor, \
    supercell.lattice.a[2]*supercell.lattice.n_x*supercell.lattice.magnify_factor);
    out.print("{} {} {}\n", supercell.lattice.b[0]*supercell.lattice.n_y*supercell.lattice.magnify_factor, \
    supercell.lattice.b[1]*supercell.lattice.n_y*supercell.lattice.magnify_factor, \
    supercell.lattice.b[2]*supercell.lattice.n_y*supercell.lattice.magnify_factor);
    out.print("{} {} {}\n", supercell.lattice.c[0]*supercell.lattice.n_z*supercell.lattice.magnify_factor, \
    supercell.lattice.c[1]*supercell.lattice.n_z*supercell.lattice.magnify_factor, \
    supercell.lattice.c[2]*supercell.lattice.n_z*supercell.lattice.magnify_factor);
    out.print("PRIMCOORD\n");
    out.print("{} 1\n", supercell.base_site.number*supercell.lattice.n_x*supercell.lattice.n_y*supercell.lattice.n_z);

    std::vector<double> index = {0, 0, 0};
    std::vector<double> coordinate = {0, 0, 0};
    for(int i=0; i<supercell.lattice.n_x; i++) {
        for(int j=0; j<supercell.lattice.n_y; j++) {
            for(int k=0; k<supercell.lattice.n_z; k++) {
                for(int l=0; l<supercell.base_site.number; l++) {
                    index = {supercell.base_site.coordinate[l][0] + i, supercell.base_site.coordinate[l][1] + j, supercell.base_site.coordinate[l][2] + k};
                    for(int m=0; m<3; m++) {
                        coordinate[m] = index[0]*supercell.lattice.a[m]*supercell.lattice.magnify_factor + \
                        index[1] * supercell.lattice.b[m]*supercell.lattice.magnify_factor + \
                        index[2] * supercell.lattice.c[m]*supercell.lattice.magnify_factor;
                    }
                    out.print("{} {} {} {} {} {} {}\n", supercell.base_site.elements[l], \
                    coordinate[0], \
                    coordinate[1], \
                    coordinate[2], \
                    supercell.site[i][j][k][l].spin[0], \
                    supercell.site[i][j][k][l].spin[1], \
                    supercell.site[i][j][k][l].spin[2]);
                }
            }
        }
    }
    out.close();
    return 0;
}

/**
 * @brief Writes the current supercell spin configuration to an XSF file with temperature suffix.
 *
 * Output filename: ``<spin_structure_file_prefix><T>.xsf`` where ``T`` is
 * converted by ``std::to_string``.
 *
 * The data layout is the same as the overload without temperature.
 *
 * @param supercell Supercell containing lattice, coordinates, and spin vectors.
 * @param spin_structure_file_prefix Prefix of the output XSF filename.
 * @param T Temperature value appended to the output filename.
 * @return int Returns 0 on completion.
 */
int WriteSpin(Supercell & supercell, std::string spin_structure_file_prefix, double T) {
    std::string output_file_name = spin_structure_file_prefix + std::to_string(T) + ".xsf";
    auto out = fmt::output_file(output_file_name);
    out.print("CRYSTAL\n");
    out.print("PRIMVEC\n");
    out.print("{} {} {}\n", supercell.lattice.a[0]*supercell.lattice.n_x*supercell.lattice.magnify_factor, \
    supercell.lattice.a[1]*supercell.lattice.n_x*supercell.lattice.magnify_factor, \
    supercell.lattice.a[2]*supercell.lattice.n_x*supercell.lattice.magnify_factor);
    out.print("{} {} {}\n", supercell.lattice.b[0]*supercell.lattice.n_y*supercell.lattice.magnify_factor, \
    supercell.lattice.b[1]*supercell.lattice.n_y*supercell.lattice.magnify_factor, \
    supercell.lattice.b[2]*supercell.lattice.n_y*supercell.lattice.magnify_factor);
    out.print("{} {} {}\n", supercell.lattice.c[0]*supercell.lattice.n_z*supercell.lattice.magnify_factor, \
    supercell.lattice.c[1]*supercell.lattice.n_z*supercell.lattice.magnify_factor, \
    supercell.lattice.c[2]*supercell.lattice.n_z*supercell.lattice.magnify_factor);
    out.print("PRIMCOORD\n");
    out.print("{} 1\n", supercell.base_site.number*supercell.lattice.n_x*supercell.lattice.n_y*supercell.lattice.n_z);

    std::vector<double> index = {0, 0, 0};
    std::vector<double> coordinate = {0, 0, 0};
    for(int i=0; i<supercell.lattice.n_x; i++) {
        for(int j=0; j<supercell.lattice.n_y; j++) {
            for(int k=0; k<supercell.lattice.n_z; k++) {
                for(int l=0; l<supercell.base_site.number; l++) {
                    index = {supercell.base_site.coordinate[l][0] + i, supercell.base_site.coordinate[l][1] + j, supercell.base_site.coordinate[l][2] + k};
                    for(int m=0; m<3; m++) {
                        coordinate[m] = index[0]*supercell.lattice.a[m]*supercell.lattice.magnify_factor + \
                        index[1] * supercell.lattice.b[m]*supercell.lattice.magnify_factor + \
                        index[2] * supercell.lattice.c[m]*supercell.lattice.magnify_factor;
                    }
                    out.print("{} {} {} {} {} {} {}\n", supercell.base_site.elements[l], \
                    coordinate[0], \
                    coordinate[1], \
                    coordinate[2], \
                    supercell.site[i][j][k][l].spin[0], \
                    supercell.site[i][j][k][l].spin[1], \
                    supercell.site[i][j][k][l].spin[2]);
                }
            }
        }
    }
    out.close();
    return 0;
}

int WriteVestaKhBondColor(Supercell & supercell, std::string output_file_prefix) {
    if(supercell.lattice.model_type != ModelType::Kitaev_Heisenberg) {
        return 0;
    }

    const std::string output_file_name = output_file_prefix + ".vesta";
    auto out = fmt::output_file(output_file_name);

    struct BondVector {
        std::array<double, 3> tail;
        std::array<double, 3> delta;
        char direction;
    };

    std::vector<BondVector> bond_vectors;
    std::set<std::string> visited_bonds;

    for(int i=0; i<supercell.lattice.n_x; i++) {
        for(int j=0; j<supercell.lattice.n_y; j++) {
            for(int k=0; k<supercell.lattice.n_z; k++) {
                for(int l=0; l<supercell.base_site.number; l++) {
                    const auto source_cart = BuildSiteCartesian(supercell, i, j, k, l);
                    for(int shell=0; shell<supercell.site[i][j][k][l].neighbor.size(); shell++) {
                        for(int entry=0; entry<supercell.site[i][j][k][l].neighbor[shell].size(); entry++) {
                            const char direction = supercell.site[i][j][k][l].neighbor_direction[shell][entry];
                            if(direction != 'x' && direction != 'y' && direction != 'z') {
                                continue;
                            }

                            Site* neighbor_ptr = supercell.site[i][j][k][l].neighbor[shell][entry];
                            int ni = -1;
                            int nj = -1;
                            int nk = -1;
                            int nl = -1;
                            bool found = false;
                            for(int ii=0; ii<supercell.lattice.n_x && !found; ii++) {
                                for(int jj=0; jj<supercell.lattice.n_y && !found; jj++) {
                                    for(int kk=0; kk<supercell.lattice.n_z && !found; kk++) {
                                        for(int ll=0; ll<supercell.base_site.number; ll++) {
                                            if(&supercell.site[ii][jj][kk][ll] == neighbor_ptr) {
                                                ni = ii;
                                                nj = jj;
                                                nk = kk;
                                                nl = ll;
                                                found = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }

                            if(!found) {
                                continue;
                            }

                            std::array<int, 8> left = {i, j, k, l, ni, nj, nk, nl};
                            std::array<int, 8> right = {ni, nj, nk, nl, i, j, k, l};
                            const std::array<int, 8>& key = (left < right) ? left : right;
                            std::ostringstream key_stream;
                            for(int key_index=0; key_index<8; key_index++) {
                                if(key_index != 0) {
                                    key_stream << ':';
                                }
                                key_stream << key[key_index];
                            }

                            if(visited_bonds.find(key_stream.str()) != visited_bonds.end()) {
                                continue;
                            }
                            visited_bonds.insert(key_stream.str());

                            const auto target_cart = BuildSiteCartesian(supercell, ni, nj, nk, nl);
                            bond_vectors.push_back({
                                source_cart,
                                {
                                    target_cart[0] - source_cart[0],
                                    target_cart[1] - source_cart[1],
                                    target_cart[2] - source_cart[2]
                                },
                                direction
                            });
                        }
                    }
                }
            }
        }
    }

    const std::array<double, 3> cell_a = BuildCartesianFromFractional(supercell, {static_cast<double>(supercell.lattice.n_x), 0.0, 0.0});
    const std::array<double, 3> cell_b = BuildCartesianFromFractional(supercell, {0.0, static_cast<double>(supercell.lattice.n_y), 0.0});
    const std::array<double, 3> cell_c = BuildCartesianFromFractional(supercell, {0.0, 0.0, static_cast<double>(supercell.lattice.n_z)});

    const double len_a = std::sqrt(cell_a[0]*cell_a[0] + cell_a[1]*cell_a[1] + cell_a[2]*cell_a[2]);
    const double len_b = std::sqrt(cell_b[0]*cell_b[0] + cell_b[1]*cell_b[1] + cell_b[2]*cell_b[2]);
    const double len_c = std::sqrt(cell_c[0]*cell_c[0] + cell_c[1]*cell_c[1] + cell_c[2]*cell_c[2]);

    auto clamp_unit = [](double value) {
        return std::max(-1.0, std::min(1.0, value));
    };
    const double alpha = std::acos(clamp_unit((cell_b[0]*cell_c[0] + cell_b[1]*cell_c[1] + cell_b[2]*cell_c[2]) / (len_b*len_c))) * 180.0 / M_PI;
    const double beta = std::acos(clamp_unit((cell_a[0]*cell_c[0] + cell_a[1]*cell_c[1] + cell_a[2]*cell_c[2]) / (len_a*len_c))) * 180.0 / M_PI;
    const double gamma = std::acos(clamp_unit((cell_a[0]*cell_b[0] + cell_a[1]*cell_b[1] + cell_a[2]*cell_b[2]) / (len_a*len_b))) * 180.0 / M_PI;

    out.print("#VESTA_FORMAT_VERSION 3.3.0\n");
    out.print("CRYSTAL\n");
    out.print("TITLE\n");
    out.print("{}\n", output_file_prefix);
    out.print("GROUP\n");
    out.print("1 1 P 1\n");
    out.print("SYMOP\n");
    out.print("0.000000 0.000000 0.000000 1 0 0 0 1 0 0 0 1 1\n");
    out.print("-1.0 -1.0 -1.0 0 0 0 0 0 0 0 0 0\n");
    out.print("TRANM 0\n");
    out.print("0.000000 0.000000 0.000000 1 0 0 0 1 0 0 0 1\n");
    out.print("LTRANSL\n");
    out.print("-1\n");
    out.print("0.000000 0.000000 0.000000 0.000000 0.000000 0.000000\n");
    out.print("LORIENT\n");
    out.print("-1 0 0 0 0\n");
    out.print("1.000000 0.000000 0.000000 1.000000 0.000000 0.000000\n");
    out.print("0.000000 0.000000 1.000000 0.000000 0.000000 1.000000\n");
    out.print("LMATRIX\n");
    out.print("1.000000 0.000000 0.000000 0.000000\n");
    out.print("0.000000 1.000000 0.000000 0.000000\n");
    out.print("0.000000 0.000000 1.000000 0.000000\n");
    out.print("0.000000 0.000000 0.000000 1.000000\n");
    out.print("0.000000 0.000000 0.000000\n");
    out.print("CELLP\n");
    out.print("{} {} {} {} {} {}\n", len_a, len_b, len_c, alpha, beta, gamma);
    out.print("0.000000 0.000000 0.000000 0.000000 0.000000 0.000000\n");

    out.print("STRUC\n");
    int atom_index = 1;
    std::map<std::string, int> atom_type_index;
    int next_atom_type = 1;
    for(int i=0; i<supercell.lattice.n_x; i++) {
        for(int j=0; j<supercell.lattice.n_y; j++) {
            for(int k=0; k<supercell.lattice.n_z; k++) {
                for(int l=0; l<supercell.base_site.number; l++) {
                    const std::string element = supercell.base_site.elements[l];
                    if(atom_type_index.find(element) == atom_type_index.end()) {
                        atom_type_index[element] = next_atom_type;
                        next_atom_type += 1;
                    }
                    const std::string site_name = element + std::to_string(atom_index);
                    const double fx = (supercell.base_site.coordinate[l][0] + static_cast<double>(i)) / static_cast<double>(supercell.lattice.n_x);
                    const double fy = (supercell.base_site.coordinate[l][1] + static_cast<double>(j)) / static_cast<double>(supercell.lattice.n_y);
                    const double fz = (supercell.base_site.coordinate[l][2] + static_cast<double>(k)) / static_cast<double>(supercell.lattice.n_z);

                    out.print("{} {} {} 1.0000 {} {} {} 1a 1\n", atom_index, element, site_name, fx, fy, fz);
                    out.print("0.000000 0.000000 0.000000 0.00\n");
                    atom_index += 1;
                }
            }
        }
    }
    out.print("0 0 0 0 0 0 0\n");

    out.print("THERI 0\n");
    for(int i=1; i<atom_index; i++) {
        out.print("{} A{} 1.000000\n", i, i);
    }
    out.print("0 0 0\n");

    out.print("SHAPE\n");
    out.print("0 0 0 0 0.000000 0 192 192 192 192\n");
    out.print("BOUND\n");
    out.print("-0.1 1.1 -0.1 1.1 -0.1 1.1\n");
    out.print("0 0 0 0 0\n");

    out.print("SBOND\n");
    out.print("1 XX XX 0.00000 0.00010 2 2 1 0 1 0.250 1.000 180 180 180\n");
    out.print("0 0 0 0\n");

    out.print("SITET\n");
    for(int i=1; i<atom_index; i++) {
        out.print("{} A{} 0.0001 128 128 128 128 128 128 204 0\n", i, i);
    }
    out.print("0 0 0 0 0 0\n");

    out.print("VECTR\n");
    for(int i=0; i<bond_vectors.size(); i++) {
        out.print("{} {} {} {} 0\n", i+1, bond_vectors[i].delta[0], bond_vectors[i].delta[1], bond_vectors[i].delta[2]);
        out.print("{} {} {} {} 0\n", i+1, bond_vectors[i].tail[0], bond_vectors[i].tail[1], bond_vectors[i].tail[2]);
        out.print("0 0 0 0 0\n");
    }
    out.print("0 0 0 0 0\n");

    out.print("VECTT\n");
    for(int i=0; i<bond_vectors.size(); i++) {
        const std::array<int, 3> rgb = BondColorFromDirection(bond_vectors[i].direction);
        out.print("{} 0.35 {} {} {} 0\n", i+1, rgb[0], rgb[1], rgb[2]);
    }
    out.print("0 0 0 0 0\n");

    out.print("SPLAN\n");
    out.print("0 0 0 0\n");
    out.print("LBLAT\n");
    out.print("-1\n");
    out.print("LBLSP\n");
    out.print("-1\n");
    out.print("DLATM\n");
    out.print("-1\n");
    out.print("DLBND\n");
    out.print("-1\n");
    out.print("DLPLY\n");
    out.print("-1\n");
    out.print("PLN2D\n");
    out.print("0 0 0 0\n");

    out.print("ATOMT\n");
    std::vector<std::pair<int, std::string>> atom_types;
    atom_types.reserve(atom_type_index.size());
    for(const auto& entry : atom_type_index) {
        atom_types.push_back({entry.second, entry.first});
    }
    std::sort(atom_types.begin(), atom_types.end(), [](const auto& left, const auto& right) {
        return left.first < right.first;
    });
    for(const auto& type_entry : atom_types) {
        out.print("{} {} 0.0001 128 128 128 128 128 128 204\n", type_entry.first, type_entry.second);
    }
    out.print("0 0 0 0 0 0\n");

    out.print("SCENE\n");
    out.print("1.000000 0.000000 0.000000 0.000000\n");
    out.print("0.000000 1.000000 0.000000 0.000000\n");
    out.print("0.000000 0.000000 1.000000 0.000000\n");
    out.print("0.000000 0.000000 0.000000 1.000000\n");
    out.print("0.000 0.000\n");
    out.print("0.000\n");
    out.print("1.000\n");

    out.print("HBOND 0 2\n");
    out.print("STYLE\n");
    out.print("DISPF 155743\n");
    out.print("MODEL 2 1 0\n");
    out.print("SURFS 0 1 1\n");
    out.print("SECTS 96 1\n");
    out.print("FORMS 0 1\n");
    out.print("ATOMS 0 0 1\n");
    out.print("BONDS 1\n");
    out.print("POLYS 1\n");
    out.print("VECTS 1.000000\n");

    out.print("FORMP\n");
    out.print("1 1.0 0 0 0\n");
    out.print("ATOMP\n");
    out.print("24 24 0 50 2.0 0\n");
    out.print("BONDP\n");
    out.print("1 16 0.0001 1.000 180 180 180\n");
    out.print("POLYP\n");
    out.print("204 1 1.000 180 180 180\n");

    out.print("ISURF\n");
    out.print("0 0 0 0\n");
    out.print("TEX3P\n");
    out.print("1 0.00000E+00 1.00000E+00\n");
    out.print("SECTP\n");
    out.print("1 0.00000E+00 1.00000E+00 0.00000E+00\n");
    out.print("HKLPP\n");
    out.print("192 1 1.000 255 0 255\n");
    out.print("UCOLP\n");
    out.print("0 1 1.000 0 0 0\n");
    out.print("COMPS 1\n");
    out.print("LABEL 1 12 1.000 0\n");
    out.print("PROJT 1 1.000\n");
    out.print("BKGRC\n");
    out.print("255 255 255\n");
    out.print("DPTHQ 1 -0.5000 3.5000\n");

    out.print("LIGHT0 1\n");
    out.print("1.000000 0.000000 0.000000 0.000000\n");
    out.print("0.000000 1.000000 0.000000 0.000000\n");
    out.print("0.000000 0.000000 1.000000 0.000000\n");
    out.print("0.000000 0.000000 0.000000 1.000000\n");
    out.print("0.000000 0.000000 20.000000 0.000000\n");
    out.print("0.000000 0.000000 -1.000000\n");
    out.print("196 196 196 255\n");
    out.print("178 178 178 255\n");
    out.print("255 255 255 255\n");

    for(int light_index=1; light_index<=3; light_index++) {
        out.print("LIGHT{}\n", light_index);
        out.print("1.000000 0.000000 0.000000 0.000000\n");
        out.print("0.000000 1.000000 0.000000 0.000000\n");
        out.print("0.000000 0.000000 1.000000 0.000000\n");
        out.print("0.000000 0.000000 0.000000 1.000000\n");
        out.print("0.000000 0.000000 20.000000 0.000000\n");
        out.print("0.000000 0.000000 -1.000000\n");
        out.print("0 0 0 0\n");
        out.print("0 0 0 0\n");
        out.print("0 0 0 0\n");
    }

    out.print("ATOMM\n");
    out.print("204 204 204 255\n");
    out.print("25.600\n");
    out.print("BONDM\n");
    out.print("255 255 255 255\n");
    out.print("128.000\n");
    out.print("POLYM\n");
    out.print("255 255 255 255\n");
    out.print("128.000\n");
    out.print("SURFM\n");
    out.print("0 0 0 255\n");
    out.print("128.000\n");
    out.print("FORMM\n");
    out.print("255 255 255 255\n");
    out.print("128.000\n");
    out.print("HKLPM\n");
    out.print("255 255 255 255\n");
    out.print("128.000\n");

    out.close();
    return 0;
}
