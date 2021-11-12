#ifndef STRUCTURE_IN
#define STRUCTURE_IN

#include <fstream>
#include <iostream>

#include <ctll.hpp>
#include <ctre.hpp>
#include <scn/scn.h>

#include "MC_structure.h"

int ReadPOSCAR(Supercell & supercell, std::string cell_structure_file);

#endif