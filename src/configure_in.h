#ifndef CONFIGURE_IN
#define CONFIGURE_IN

#include <fstream>
#include <iostream>
#include <string>

#include <ctll.hpp>
#include <ctre.hpp>
#include <scn/scn.h>
#include <toml++/toml.h>

#include "constants.h"
#include "MC_structure.h"

int ReadSettingFile(Supercell & supercell, MonteCarlo & monte_carlo, std::string input_file);

#endif