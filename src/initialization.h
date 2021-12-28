#ifndef INITIALIZATION
#define INITIALIZATION

#include <string>
#include <vector>

#include "MC_structure.h"
#include "rotation.h"
#include "Hamiltonian.h"
#include "../custom/custom.h"

int ReadOptions(int argc, char** argv, std::string & cell_structure_file, std::string & input_file, std::string & output_file, std::string & spin_structure_file_prefix);
void usage(char* program);
int EnlargeCell(Supercell & supercell);
double Distance(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int>, std::vector<double>, std::vector<double>);
int AddDistance(double distance, std::vector<double> & distance_list, double tolerance_percentage);
int InitializeSupercell(Supercell & supercell);

#endif