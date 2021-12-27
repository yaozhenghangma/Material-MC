#ifndef OUTLOG
#define OUTLOG


#include <fmt/core.h>
#include <fmt/os.h>

#include "MC_structure.h"

int WriteLog(Supercell &, MonteCarlo &, fmt::v8::ostream &);

#endif