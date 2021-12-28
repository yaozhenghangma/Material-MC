#include <mpi.h>
#include <iostream>
#include <random>
#include <stdlib.h>
#include <unistd.h>

#include <boost/mpi.hpp>

#include "../custom/custom.h"
#include "MC_structure.h"
#include "structure_in.h"
#include "configure_in.h"
#include "Hamiltonian.h"
#include "spin_out.h"
#include "random_function.h"
#include "result_out.h"
#include "rotation.h"
#include "log.h"
#include "initialization.h"

namespace mpi = boost::mpi;

using namespace std;

const double KB = 0.08617343;

int MonteCarloRelaxing(Supercell & supercell, MonteCarlo & monte_carlo, double T);
vector<double> MonteCarloStep(Supercell & supercell, MonteCarlo & monte_carlo, double T);
int Flip(Supercell & supercell, Site & one_site, double T);


int main(int argc, char** argv) {
    mpi::environment env;
    mpi::communicator world;

    Supercell supercell;
    MonteCarlo monte_carlo;
    string cell_structure_file = "POSCAR";
    string input_file = "input.toml"; 
    string output_file = "output.txt";
    string spin_structure_file_prefix = "spin";

    auto logger = fmt::output_file("log.txt");

    // Read parameters and broadcast data using root processor
    if(world.rank() == 0) {
        // Read information from command line.
        ReadOptions(argc, argv, cell_structure_file, input_file, output_file, spin_structure_file_prefix);

        // Read information from input file.
        ReadSettingFile(supercell, monte_carlo, input_file);

        // Read information from POSCAR
        ReadPOSCAR(supercell, cell_structure_file);

        // Log file
        logger.print("Successfully process input file.");
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
        logger.print("Successfully initialize the supercell.");
        WriteLog(supercell, monte_carlo, logger);
        WriteSpin(supercell, "initialization_spin_configure", 0.0);
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
    vector<double> chi;
    vector<double> moment_projection;
    vector<double> chi_projection;
    double energy_every_processor;
    double Cv_every_processor;
    double moment_every_processor;
    double moment_projection_every_processor;
    double chi_projection_every_processor;
    double chi_every_processor;
    vector<double> gathered_energy;
    vector<double> gathered_Cv;
    vector<double> gathered_moment;
    vector<double> gathered_chi;
    vector<double> gathered_moment_projection;
    vector<double> gathered_chi_projection;
    static double one_over_number = 1.0 / (supercell.lattice.n_x * supercell.lattice.n_y * \
    supercell.lattice.n_z * supercell.base_site.number);
    for(int i=0; i<quotient+1; i++) {
        if(i<quotient || remainder !=0 ) {
            T = monte_carlo.start_temperature + (i*world.size()+world.rank())*monte_carlo.temperature_step;
            MonteCarloRelaxing(supercell, monte_carlo, T);
            result_value = MonteCarloStep(supercell, monte_carlo, T);
            WriteSpin(supercell, spin_structure_file_prefix, T);

            Cv_every_processor = (result_value[1]-result_value[0]*result_value[0])*one_over_number/(KB*T*T); //Cv
            chi_every_processor = (result_value[3]-result_value[2]*result_value[2])/(KB*T); //chi
            energy_every_processor = result_value[0] * one_over_number; //energy
            moment_every_processor = result_value[2] * one_over_number; //moment
            if(supercell.lattice.field) {
                chi_projection_every_processor = (result_value[5] - result_value[4]*result_value[4])/(KB*T);
                moment_projection_every_processor = result_value[4] * one_over_number;
            }
            
            // Collect data form all processors.
            gather(world, energy_every_processor, gathered_energy, 0);
            gather(world, Cv_every_processor, gathered_Cv, 0);
            gather(world, moment_every_processor, gathered_moment, 0);
            gather(world, chi_every_processor, gathered_chi, 0);
            if(supercell.lattice.field) {
                gather(world, moment_projection_every_processor, gathered_moment_projection, 0);
                gather(world, chi_projection_every_processor, gathered_chi_projection, 0);
            }

            // Store these data
            if(world.rank() == 0) {
                for(int j=0; j<world.size(); j++) {
                    energy.push_back(gathered_energy[j]);
                    Cv.push_back(gathered_Cv[j]);
                    moment.push_back(gathered_moment[j]);
                    chi.push_back(gathered_chi[j]);
                    if(supercell.lattice.field) {
                        moment_projection.push_back(gathered_moment_projection[j]);
                        chi_projection.push_back(gathered_chi_projection[j]);
                    }
                }
            }
            
        }
    }

    // Output the thermal dynamic result using root processor.
    if(world.rank() == 0) {
        logger.print("Successfully run Monte Carlo simulation.");
        WriteOutput(monte_carlo, energy, Cv, moment, chi, moment_projection, chi_projection, supercell.lattice.field, output_file);
        logger.print("Successfully output all results.");
    }
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
    double tmp_momentum_projection = 0;
    double total_momentum_projection = 0;
    double total_momentum_projection_square = 0;
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
        tmp_momentum_projection = supercell.lattice.field_direction[0] * tmp_m_component[0]
                                + supercell.lattice.field_direction[1] * tmp_m_component[1] 
                                + supercell.lattice.field_direction[2] * tmp_m_component[2];
        total_momentum_projection += tmp_momentum_projection;
        total_momentum_projection_square += tmp_momentum_projection * tmp_momentum_projection;
    }
    
    return {total_energy * one_over_step, total_energy_square * one_over_step, \
    total_momentum * one_over_step, total_momentum_square * one_over_step, \
    total_momentum_projection * one_over_step, total_momentum_projection_square * one_over_step};
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