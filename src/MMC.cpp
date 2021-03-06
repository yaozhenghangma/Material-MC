#include <ctime>
#include <mpi.h>
#include <iostream>
#include <random>
#include <stdlib.h>
#include <unistd.h>

#include <boost/mpi.hpp>

#include "MC_structure.h"
#include "structure_in.h"
#include "configure_in.h"
#include "spin_out.h"
#include "result_out.h"
#include "log.h"
#include "initialization.h"

#include "methods/classical.h"
#include "methods/parallel_tempering.h"

namespace mpi = boost::mpi;

using namespace std;

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
    time_t now = time(0);
    string dt = ctime(&now);

    // Read parameters and broadcast data using root processor
    if(world.rank() == 0) {
        logger.print("Material Monte Carlo is run on {} CPU cores.\n", world.size());
        logger.print("Start time: " + dt);
        // Read information from command line.
        ReadOptions(argc, argv, cell_structure_file, input_file, output_file, spin_structure_file_prefix);

        // Read information from input file.
        ReadSettingFile(supercell, monte_carlo, input_file);
        logger.print("Successfully process input file: " + input_file + ".\n");

        // Read information from POSCAR
        ReadPOSCAR(supercell, cell_structure_file);
        logger.print("Successfully load structure file: " + cell_structure_file + ".\n");
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
        WriteLog(supercell, monte_carlo, logger);
        WriteSpin(supercell, "structure_initialized");
    }

    vector<double> energy;
    vector<double> Cv;
    vector<double> moment;
    vector<double> chi;
    vector<double> moment_projection;
    vector<double> chi_projection;
    if(monte_carlo.methods == Methods::classical) {
        ClassicalMonteCarlo(env, world, 
            monte_carlo, supercell, spin_structure_file_prefix, 
            energy, Cv,
            moment, chi,
            moment_projection, chi_projection);
    } else {
        ParallelTemperingMonteCarlo(env, world, 
            monte_carlo, supercell, spin_structure_file_prefix, 
            energy, Cv,
            moment, chi,
            moment_projection, chi_projection);
    }

    // Output the thermal dynamic result using root processor.
    if(world.rank() == 0) {
        WriteOutput(monte_carlo, energy, Cv, moment, chi, moment_projection, chi_projection, supercell.lattice.field, output_file);
        logger.print("Successfully output all results into " + output_file + ".\n");
        now = time(0);
        dt = ctime(&now);
        logger.print("End time: " + dt);
    }
    return 0;
}