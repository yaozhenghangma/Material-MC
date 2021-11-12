#ifndef SPIN_OUT
#define SPIN_OUT

#include <string>

#include <fmt/core.h>
#include <fmt/os.h>

#include "MC_structure.h"

int WriteSpin(Supercell & supercell, std::string spin_structure_file_prefix, double T);

#endif