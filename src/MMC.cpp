#include <mpi.h>
#include <iostream>
#include <random>
#include <stdlib.h>
#include <unistd.h>

#include <boost/mpi.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "MC_structure.h"
#include "structure_in.h"
#include "configure_in.h"
#include "Hamiltonion.h"
#include "spin_out.h"
#include "result_out.h"
#include "log.h"

namespace mpi = boost::mpi;

using namespace std;

const double KB = 0.08617343;

vector<int> RandomSite(int n_x, int n_y, int n_z, int base_n);
vector<double> RandomSpin(double scaling);
double RandomFloat();
int ReadOptions(int argc, char** argv, string & cell_structure_file, string & input_file, string & output_file, string & spin_structure_file_prefix);
void usage(char* program);
int EnlargeCell(Supercell & supercell);
double Distance(vector<double>, vector<double>, vector<double>, vector<int>, vector<double>, vector<double>);
int AddDistance(double distance, vector<double> & distance_list, double tolerance_percentage);
int InitializeSupercell(Supercell & supercell);
int MonteCarloRelaxing(Supercell & supercell, MonteCarlo & monte_carlo, double T);
vector<double> MonteCarloStep(Supercell & supercell, MonteCarlo & monte_carlo, double T);
int Flip(Supercell & supercell, Site & one_site, double T);


int main(int argc, char** argv) {
    mpi::environment env;
    mpi::communicator world;

    Supercell supercell;
    MonteCarlo monte_carlo;
    string cell_structure_file = "POSCAR";
    string input_file = "input.txt"; 
    string output_file = "output.txt";
    string spin_structure_file_prefix = "spin";

    shared_ptr<spdlog::logger> logger;

    // Read parameters and broadcast data using root processor
    if(world.rank() == 0) {
        // Read information from command line.
        ReadOptions(argc, argv, cell_structure_file, input_file, output_file, spin_structure_file_prefix);

        // Read information from input file.
        ReadSettingFile(supercell, monte_carlo, input_file);

        // Read information from POSCAR
        ReadPOSCAR(supercell, cell_structure_file);

        // Log file
        logger = spdlog::basic_logger_mt("basic_logger", "log.txt");
        logger->info("Successfully process input file.");
    }
    // Broadcast monte_carlo, base_site, lattice and spin_structure_file_prefix.
    broadcast(world, supercell.base_site, 0);
    broadcast(world, supercell.lattice, 0);
    broadcast(world, supercell.initialization, 0);
    broadcast(world, monte_carlo, 0);

    // Enlarge the cell with given n.
    EnlargeCell(supercell);
    InitializeSupercell(supercell);

    // Output the coordinate number
    if(world.rank() == 0) {
        logger->info("Successfully initialize the supercell.");
        WriteLog(supercell, monte_carlo, logger);
    }
    // Arrange the processors.
    int quotient = monte_carlo.temperature_step_number / world.size();
    int remainder = monte_carlo.temperature_step_number % world.size();

    // Monte Carlo
    double T;
    vector<double> result_value;
    vector<double> energy;
    vector<double> Cv;
    vector<double> moment;
    vector<double> Ki;
    vector<double> Ki_x;
    vector<double> Ki_y;
    vector<double> Ki_z;
    vector<double> moment_x;
    vector<double> moment_y;
    vector<double> moment_z;
    double energy_every_processor;
    double Cv_every_processor;
    double moment_every_processor;
    double Ki_every_processor;
    double Ki_x_every_processor;
    double Ki_y_every_processor;
    double Ki_z_every_processor;
    double moment_x_ever_processor;
    double moment_y_ever_processor;
    double moment_z_ever_processor;
    vector<double> gathered_energy;
    vector<double> gathered_Cv;
    vector<double> gathered_moment;
    vector<double> gathered_Ki;
    vector<double> gathered_Ki_x;
    vector<double> gathered_Ki_y;
    vector<double> gathered_Ki_z;
    vector<double> gathered_moment_x;
    vector<double> gathered_moment_y;
    vector<double> gathered_moment_z;
    static double one_over_number = 1.0 / (supercell.lattice.n_x * supercell.lattice.n_y * \
    supercell.lattice.n_z * supercell.base_site.number);
    for(int i=0; i<quotient+1; i++) {
        if(i<quotient || remainder !=0 ) {
            T = monte_carlo.start_temperature + (i*world.size()+world.rank())*monte_carlo.temperature_step;
            MonteCarloRelaxing(supercell, monte_carlo, T);
            result_value = MonteCarloStep(supercell, monte_carlo, T);
            WriteSpin(supercell, spin_structure_file_prefix, T);

            Cv_every_processor = (result_value[1]-result_value[0]*result_value[0])*one_over_number/(KB*T*T); //Cv
            Ki_every_processor = (result_value[3]-result_value[2]*result_value[2])/(KB*T); //Ki
            Ki_x_every_processor = (result_value[5]-result_value[4]*result_value[4])/(KB*T); //Ki
            Ki_y_every_processor = (result_value[7]-result_value[6]*result_value[6])/(KB*T); //Ki
            Ki_z_every_processor = (result_value[9]-result_value[8]*result_value[8])/(KB*T); //Ki
            energy_every_processor = result_value[0] * one_over_number; //energy
            moment_every_processor = result_value[2] * one_over_number; //moment
            moment_x_ever_processor = result_value[4] * one_over_number; //moment_x
            moment_y_ever_processor = result_value[6] * one_over_number; //moment_y
            moment_z_ever_processor = result_value[8] * one_over_number; //moment_z
            
            // Collect data form all processors.
            gather(world, energy_every_processor, gathered_energy, 0);
            gather(world, Cv_every_processor, gathered_Cv, 0);
            gather(world, moment_every_processor, gathered_moment, 0);
            gather(world, Ki_every_processor, gathered_Ki, 0);
            gather(world, moment_x_ever_processor, gathered_moment_x, 0);
            gather(world, moment_y_ever_processor, gathered_moment_y, 0);
            gather(world, moment_z_ever_processor, gathered_moment_z, 0);
            gather(world, Ki_x_every_processor, gathered_Ki_x, 0);
            gather(world, Ki_y_every_processor, gathered_Ki_y, 0);
            gather(world, Ki_z_every_processor, gathered_Ki_z, 0);

            // Store these data
            if(world.rank() == 0) {
                for(int j=0; j<world.size(); j++) {
                    energy.push_back(gathered_energy[j]);
                    Cv.push_back(gathered_Cv[j]);
                    moment.push_back(gathered_moment[j]);
                    Ki.push_back(gathered_Ki[j]);
                    moment_x.push_back(gathered_moment_x[j]);
                    moment_y.push_back(gathered_moment_y[j]);
                    moment_z.push_back(gathered_moment_z[j]);
                    Ki_x.push_back(gathered_Ki_x[j]);
                    Ki_y.push_back(gathered_Ki_y[j]);
                    Ki_z.push_back(gathered_Ki_z[j]);
                }
            }
            
        }
    }

    // Output the thermal dynamic result using root processor.
    if(world.rank() == 0) {
        logger->info("Successfully run Monte Carlo simulation.");
        WriteOutput(monte_carlo, energy, Cv, moment, Ki, moment_x, moment_y, moment_z, Ki_x, Ki_y, Ki_z, output_file);
        logger->info("Successfully output all results.");
    }
    return 0;
}

vector<int> RandomSite(int n_x, int n_y, int n_z, int base_n) {
    // Return the site index randomly.
    static random_device rd;
    static mt19937 engine(rd());
    static uniform_int_distribution<int> int_distribution_x(0, n_x-1);
    static uniform_int_distribution<int> int_distribution_y(0, n_y-1);
    static uniform_int_distribution<int> int_distribution_z(0, n_z-1);
    static uniform_int_distribution<int> int_distribution_base(0, base_n-1);

    vector<int> index = {0, 0, 0, 0};

    index[0] = int_distribution_x(engine);
    index[1] = int_distribution_y(engine);
    index[2] = int_distribution_z(engine);
    index[3] = int_distribution_base(engine);

    return index;
}

vector<double> RandomSpin(double scaling) {
    // Return a unit spin vector randomly.
    static random_device rd;
    static mt19937 engine(rd());
    static normal_distribution<double> normal{0, 1};

    double x1 = normal(engine);
    double x2 = normal(engine);
    double x3 = normal(engine);
    double factor = scaling/sqrt(x1*x1+x2*x2+x3*x3);

    return {x1*factor, x2*factor, x3*factor};
}

double RandomFloat() {
    // Return a float number between 0 and 1.
    static random_device rd;
    static mt19937 engine(rd());
    static uniform_real_distribution<double> double_distribution(0, 1);

    return double_distribution(engine);
}

int ReadOptions(int argc, char** argv, string & cell_structure_file, string & input_file, string & output_file, string & spin_structure_file_prefix) {
    // Process options from command line.
    char option;
    while ((option = getopt(argc, argv, "c:i:o:s:h")) != -1) {
        switch (option) {
            case 'c':
                cell_structure_file = optarg; break;
            case 'i':
                input_file = optarg; break;
            case 'o':
                output_file = optarg; break;
            case 's':
                spin_structure_file_prefix = optarg; break;
            case 'h':
                usage(argv[0]);
        }
    }

    return 0;
}

void usage(char* program) {
    //TODO: Usage of the program.
    exit(1);
}

int EnlargeCell(Supercell & supercell) {
    // Enlarge the system with given number and initialize the spin.
    vector<vector<vector<Site>>> site1;
    vector<vector<Site>> site2;
    vector<Site> site3;
    Site site4;

    site4.spin[0] = supercell.initialization.direction[0];
    site4.spin[1] = supercell.initialization.direction[1];
    site4.spin[2] = supercell.initialization.direction[2];

    for(int i=0; i<supercell.lattice.n_x; i++) {
        supercell.site.push_back(site1);
        for(int j=0; j<supercell.lattice.n_y; j++) {
            supercell.site[i].push_back(site2);
            for(int k=0; k<supercell.lattice.n_z; k++) {
                supercell.site[i][j].push_back(site3);
                for(int l=0; l<supercell.base_site.number; l++) {
                    supercell.site[i][j][k].push_back(site4);
                    supercell.site[i][j][k][l].spin[0] *= supercell.base_site.spin_scaling[l];
                    supercell.site[i][j][k][l].spin[1] *= supercell.base_site.spin_scaling[l];
                    supercell.site[i][j][k][l].spin[2] *= supercell.base_site.spin_scaling[l];
                    supercell.site[i][j][k][l].spin_scaling = & supercell.base_site.spin_scaling[l];
                    supercell.site[i][j][k][l].anisotropic_factor = & supercell.base_site.anisotropic_factor[l];
                    supercell.site[i][j][k][l].super_exchange_parameter = & supercell.base_site.super_exchange_parameter[l];
                    supercell.site[i][j][k][l].neighbor_number = & supercell.base_site.neighbor_number[l];
                }
            }
        }
    }

    return 0;
}

double Distance(vector<double> lattice_constant1, vector<double> lattice_constant2, vector<double> lattice_constant3, \
vector<int> index, vector<double> base_site1, vector<double> base_site2) {
    // Return squared distance of index-base_site1+base_site2.
    vector<double> total_index = {index[0]+base_site2[0]-base_site1[0], index[1]+base_site2[1]-base_site1[1], index[2]+base_site2[2]-base_site1[2]};
    double result = 0;
    for(int i=0; i<3; i++) {
        result += (total_index[0]*lattice_constant1[i] + total_index[1]*lattice_constant2[i] + total_index[2]*lattice_constant3[i]) \
        * (total_index[0]*lattice_constant1[i] + total_index[1]*lattice_constant2[i] + total_index[2]*lattice_constant3[i]);
    }
    return result;
}

int AddDistance(double distance, vector<double> & distance_list, double tolerance_percentage) {
    int s = distance_list.size()-1;
    // Check similar distance
    for(int i=0; i<s+1; i++) {
        if(distance_list[i] == 0) {
            break;
        } else if (abs(distance_list[i]-distance) < min(distance, distance_list[i]) * tolerance_percentage){
            return 0;
        }
    }

    // Add new distance and order
    if(distance_list[s] == 0 || distance_list[s] > distance) {
        distance_list[s] = distance;
        for(int i=s-1; i>=0; i--) {
            if(distance_list[i] == 0 || distance_list[i] > distance) {
                distance_list[i+1] = distance_list[i];
                distance_list[i] = distance;
            }
        }
    }
    return 0;
}

int InitializeSupercell(Supercell & supercell) {
    // Initialize Hamiltonian
    if(supercell.lattice.function_choice == "Heisenberg") {
        supercell.Hamiltonian = Heisenberg;
    } else if(supercell.lattice.function_choice == "Heisenberg_3_axis_anisotropy") {
        supercell.Hamiltonian = Heisenberg_3_anisotropy;
    } else if(supercell.lattice.function_choice == "Heisenberg_external_field") {
        supercell.Hamiltonian = Heisenberg_external_field;
    }  else {
        supercell.Hamiltonian = Heisenberg;
    }

    // Initialize the neighbors' link and energy.
    // Find the neighbors for every base sites in a cell.
    // neighbors_index[i][j][k][l]: i-th site's k-th j nearest neighbor's index in l-th direction. 4-th direction is base number.
    vector<vector<vector<vector<int>>>> neighbors_index;

    for(int i=0; i<supercell.base_site.number; i++) {
        vector<vector<vector<int>>> index1;
        neighbors_index.push_back(index1);
        for(int j=0; j<supercell.base_site.neighbor_number[i]; j++) {
            vector<vector<int>> index2;
            neighbors_index[i].push_back(index2);
        }
    }

    for(int i=0; i<supercell.base_site.number; i++) {
        double distance_square = 0;

        // Find link.
        for(int j=-supercell.base_site.neighbor_number[i]; j<supercell.base_site.neighbor_number[i]+1; j++) {
            for(int k=-supercell.base_site.neighbor_number[i]; k<supercell.base_site.neighbor_number[i]+1; k++) {
                for(int l=-supercell.base_site.neighbor_number[i]; l<supercell.base_site.neighbor_number[i]+1; l++) {
                    for(int m=0; m<supercell.base_site.number; m++) {
                        distance_square = Distance(supercell.lattice.a, supercell.lattice.b, supercell.lattice.c, {j, k, l}, \
                        supercell.base_site.coordinate[i], supercell.base_site.coordinate[m]);

                        if(distance_square == 0.0) {
                            continue;
                        } else {
                            for(int n=0; n<supercell.base_site.neighbor_number[i]; n++) {
                                if(supercell.base_site.elements[m] == supercell.base_site.neighbor_elements[i][n] \
                                && abs(distance_square - supercell.base_site.neighbor_distance_square[i][n]) \
                                < supercell.lattice.tolerance_percentage*min(distance_square, supercell.base_site.neighbor_distance_square[i][n])) {
                                    vector<int> ind = {j, k, l, m};
                                    neighbors_index[i][n].emplace_back(ind);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Set link for every sites.
    for(int i=0; i<supercell.lattice.n_x; i++) {
        for(int j=0; j<supercell.lattice.n_y; j++) {
            for(int k=0; k<supercell.lattice.n_z; k++) {
                for(int l=0; l<supercell.base_site.number; l++) {
                    vector<Site*> temp = {};
                    for(int m=0; m<supercell.base_site.neighbor_number[l]; m++) {
                        supercell.site[i][j][k][l].neighbor.push_back(temp);
                        for(int n=0; n<neighbors_index[l][m].size(); n++) {
                            supercell.site[i][j][k][l].neighbor[m].push_back( 
                            & supercell.site[(i+neighbors_index[l][m][n][0]%supercell.lattice.n_x+supercell.lattice.n_x) % supercell.lattice.n_x] \
                            [(j+neighbors_index[l][m][n][1]%supercell.lattice.n_y+supercell.lattice.n_y) % supercell.lattice.n_y] \
                            [(k+neighbors_index[l][m][n][2]%supercell.lattice.n_z+supercell.lattice.n_z) % supercell.lattice.n_z] \
                            [neighbors_index[l][m][n][3]]);
                        }
                    }
                }
            }
        }
    }

    // Initialize the spin with given configuration
    if(supercell.initialization.anti_ferromagnetic) {
        for(int i=0; i<supercell.lattice.n_x; i++) {
            for(int j=0; j<supercell.lattice.n_y; j++) {
                for(int k=0; k<supercell.lattice.n_z; k++) {
                    for(int l=0; l<supercell.base_site.number; l++) {
                        for(int m=0; m<supercell.initialization.anti_ferromagnetic_J[l].size(); m++) {
                            for(int n=0; n<supercell.site[i][j][k][l].neighbor[supercell.initialization.anti_ferromagnetic_J[l][m]-1].size(); n++) {
                                if(supercell.site[i][j][k][l].spin[0] * supercell.site[i][j][k][l].neighbor[supercell.initialization.anti_ferromagnetic_J[l][m]-1][n]->spin[0] + \
                                supercell.site[i][j][k][l].spin[1] * supercell.site[i][j][k][l].neighbor[supercell.initialization.anti_ferromagnetic_J[l][m]-1][n]->spin[1] + \
                                supercell.site[i][j][k][l].spin[2] * supercell.site[i][j][k][l].neighbor[supercell.initialization.anti_ferromagnetic_J[l][m]-1][n]->spin[2] > 0) {
                                    supercell.site[i][j][k][l].neighbor[supercell.initialization.anti_ferromagnetic_J[l][m]-1][n]->reverse_spin();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    supercell.lattice.total_energy = supercell.energy();
    return 0;
}

int MonteCarloRelaxing(Supercell & supercell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    vector<int> site_chosen;
    for(int i=0; i<monte_carlo.relax_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            Flip(supercell, supercell[site_chosen], T);
        }
    }
    
    return 0;
}

vector<double> MonteCarloStep(Supercell & supercell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    vector<int> site_chosen;
    double total_energy = 0;
    double total_energy_square = 0;
    double tmp_momentum = 0;
    double total_momentum = 0;
    double total_momentum_square = 0;
    vector<double> total_m_component = {0, 0, 0};
    vector<double> total_m_square_component = {0, 0, 0};
    vector<double> tmp_m_component;
    static double one_over_step = 1.0 / monte_carlo.count_step;
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            Flip(supercell, supercell[site_chosen], T);
        }
        total_energy += supercell.lattice.total_energy;
        total_energy_square += supercell.lattice.total_energy * supercell.lattice.total_energy;
        tmp_momentum = supercell.momentum();
        total_momentum += tmp_momentum;
        total_momentum_square += tmp_momentum * tmp_momentum;
        tmp_m_component = supercell.momentum_component();
        total_m_component[0] += tmp_m_component[0];
        total_m_component[1] += tmp_m_component[1];
        total_m_component[2] += tmp_m_component[2];
        total_m_square_component[0] += tmp_m_component[0] * tmp_m_component[0];
        total_m_square_component[1] += tmp_m_component[1] * tmp_m_component[1];
        total_m_square_component[2] += tmp_m_component[2] * tmp_m_component[2];
    }
    
    return {total_energy * one_over_step, total_energy_square * one_over_step, \
    total_momentum * one_over_step, total_momentum_square * one_over_step, \
    total_m_component[0] * one_over_step, total_m_square_component[0] * one_over_step, \
    total_m_component[1] * one_over_step, total_m_square_component[1] * one_over_step, \
    total_m_component[2] * one_over_step, total_m_square_component[2] * one_over_step};
}

int Flip(Supercell & supercell, Site & one_site, double T) {
    // Flip one spin.
    // Energy and spin before flip.
    double energy_old = supercell.Hamiltonian(supercell.base_site, one_site);
    vector<double> old_spin = one_site.spin;

    // Energy and spin after flip.
    one_site.spin = RandomSpin(*one_site.spin_scaling);
    double energy_new = supercell.Hamiltonian(supercell.base_site, one_site);
    double de = energy_new - energy_old;

    // Judge whether to flip.
    double crition = RandomFloat();
    if (crition > exp(-de/(KB*T))) {
        one_site.spin = old_spin;
    } else {
        supercell.lattice.total_energy += de;
    }

    return 0;
}