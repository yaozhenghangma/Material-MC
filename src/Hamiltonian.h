#ifndef HAMILTONION
#define HAMILTONION

#include "MC_structure.h"

double Heisenberg(BaseSite & base_site, Site & site); 
double Heisenberg_xyz_anisotropy(BaseSite & base_site, Site & site); 
double Heisenberg_with_field(BaseSite & base_site, Site & site); 
double Heisenberg_xyz_anisotropy_with_field(BaseSite & base_site, Site & site);
double Heisenberg_x_anisotropy(BaseSite & base_site, Site & site); 
double Heisenberg_y_anisotropy(BaseSite & base_site, Site & site); 
double Heisenberg_z_anisotropy(BaseSite & base_site, Site & site);
double Heisenberg_xy_anisotropy(BaseSite & base_site, Site & site); 
double Heisenberg_yz_anisotropy(BaseSite & base_site, Site & site); 
double Heisenberg_zx_anisotropy(BaseSite & base_site, Site & site);
double Heisenberg_x_anisotropy_with_field(BaseSite & base_site, Site & site); 
double Heisenberg_y_anisotropy_with_field(BaseSite & base_site, Site & site); 
double Heisenberg_z_anisotropy_with_field(BaseSite & base_site, Site & site);
double Heisenberg_xy_anisotropy_with_field(BaseSite & base_site, Site & site); 
double Heisenberg_yz_anisotropy_with_field(BaseSite & base_site, Site & site); 
double Heisenberg_zx_anisotropy_with_field(BaseSite & base_site, Site & site);

#endif