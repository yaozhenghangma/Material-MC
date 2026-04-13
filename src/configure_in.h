#ifndef CONFIGURE_IN
#define CONFIGURE_IN

#include <array>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <ctll.hpp>
#include <ctre.hpp>
#include <scn/scan.h>
#include <toml++/toml.h>

#include "constants.h"
#include "MC_structure.h"

int ReadSettingFile(Supercell & supercell, MonteCarlo & monte_carlo, std::string input_file);

#endif
