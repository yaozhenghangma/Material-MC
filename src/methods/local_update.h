#ifndef LOCAL_UPDATE
#define LOCAL_UPDATE

#include "../MC_structure.h"
#include "../random_function.h"

int LocalUpdateHeisenberg(Supercell & supercell, Site & one_site, double T);
int LocalUpdateIsing(Supercell & supercell, Site & one_site, double T);

#endif