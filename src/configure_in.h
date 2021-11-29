#ifndef CONFIGURE_IN
#define CONFIGURE_IN

#include <fstream>
#include <numbers>
#include <iostream>
#include <string>

#include <ctll.hpp>
#include <ctre.hpp>
#include <scn/scn.h>
#include <toml++/toml.h>

#include "MC_structure.h"

const double MuB = 5.7883818012e-2;

int ReadSettingFile(Supercell & supercell, MonteCarlo & monte_carlo, std::string input_file);

#endif