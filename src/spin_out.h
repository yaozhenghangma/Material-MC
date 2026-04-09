#ifndef SPIN_OUT
#define SPIN_OUT

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <fmt/core.h>
#include <fmt/os.h>

#include "MC_structure.h"

int WriteSpin(Supercell & supercell, std::string spin_structure_file_prefix);
int WriteSpin(Supercell & supercell, std::string spin_structure_file_prefix, double T);
int WriteVestaKhBondColor(Supercell & supercell, std::string output_file_prefix);

#endif