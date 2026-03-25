#include "spin_out.h"

/**
 * @file spin_out.cpp
 * @brief Utilities for exporting spin configurations in XSF format.
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