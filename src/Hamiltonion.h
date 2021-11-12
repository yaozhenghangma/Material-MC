#ifndef HAMILTONION
#define HAMILTONION

#include "MC_structure.h"

double Heisenberg(BaseSite & base_site, Site & site);
double Heisenberg_3_anisotropy(BaseSite & base_site, Site & site);
double Heisenberg_external_field(BaseSite & base_site, Site & site);

#endif