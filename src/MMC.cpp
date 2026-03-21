#include <mpi.h>
#include <iostream>
#include <random>
#include <stdlib.h>
#include <unistd.h>

#include <sstream>
#include <string>

#include <cereal/archives/binary.hpp>

#include "MC_structure.h"
#include "structure_in.h"
#include "configure_in.h"
#include "spin_out.h"
#include "result_out.h"
#include "log.h"
#include "initialization.h"

#include "methods/classical.h"
#include "methods/parallel_tempering.h"

using namespace std;

namespace mpi_cereal {
void bcast_bytes(const void* data, int size, int root, MPI_Comm comm) {
    MPI_Bcast(const_cast<void*>(data), size, MPI_BYTE, root, comm);
}

template <typename T>
void bcast_object(T& obj, int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    int size = 0;
    std::string buffer;

    if (rank == root) {
        std::ostringstream oss(std::ios::binary);
        cereal::BinaryOutputArchive archive(oss);
        archive(obj);
        buffer = oss.str();
        size = static_cast<int>(buffer.size());
    }

    MPI_Bcast(&size, 1, MPI_INT, root, comm);
    if (rank != root) {
        buffer.resize(size);
    }
    bcast_bytes(buffer.data(), size, root, comm);

    if (rank != root) {
        std::istringstream iss(buffer, std::ios::binary);
        cereal::BinaryInputArchive archive(iss);
        archive(obj);
    }
}
}

/**
 * @brief Program entry: MPI-parallel Monte Carlo simulation for materials.
 *
 * High-level flow (per README):
 * 1) Parse CLI + input.toml + POSCAR on rank 0.
 * 2) Broadcast configuration to all ranks.
 * 3) Build/initialize the supercell.
 * 4) Run selected MC method (classical or parallel tempering).
 * 5) Output logs/results on rank 0.
 */
int main(int argc, char** argv) {
    /**
     * @brief Initialize MPI environment and communicator.
     *
     * mpi::environment manages MPI init/finalize.
     * mpi::communicator world provides rank/size and collective ops.
     */
    MPI_Init(&argc, &argv);
    MPI_Comm world = MPI_COMM_WORLD;
    int world_rank = 0;
    MPI_Comm_rank(world, &world_rank);
    int world_size = 0;
    MPI_Comm_size(world, &world_size);

    /**
     * @brief Core simulation data structures.
     *
     * Supercell: lattice/structure, spin states, initialization settings.
     * MonteCarlo: MC parameters, temperature schedule, method selection, etc.
     */
    Supercell supercell;
    MonteCarlo monte_carlo;

    /**
     * @brief Default I/O file names (can be overridden by CLI).
     *
     * POSCAR: structure input (VASP-style).
     * input.toml: simulation parameters.
     * output.txt: thermodynamic results.
     * spin_*: spin configuration outputs.
     */
    string cell_structure_file = "POSCAR";
    string input_file = "input.toml";
    string output_file = "output.txt";
    string spin_structure_file_prefix = "spin";

    /**
     * @brief Log file handle (rank 0 writes status messages).
     */
    auto logger = fmt::output_file("log.txt");

    // Read parameters and broadcast data using root processor
    if(world_rank == 0) {
        /**
         * @brief Parse CLI options (may override default file names).
         */
        ReadOptions(argc, argv, cell_structure_file, input_file, output_file, spin_structure_file_prefix);

        /**
         * @brief Read simulation parameters from input.toml.
         */
        ReadSettingFile(supercell, monte_carlo, input_file);

        /**
         * @brief Read crystal structure from POSCAR into supercell.
         */
        ReadPOSCAR(supercell, cell_structure_file);

        // Log file
        logger.print("Successfully process input file.\n");
    }
    /**
     * @brief Broadcast configuration and structure data to all ranks.
     *
     * These objects are serialized by cereal and broadcast as bytes.
     */
    mpi_cereal::bcast_object(supercell.base_site, 0, world);
    mpi_cereal::bcast_object(supercell.lattice, 0, world);
    mpi_cereal::bcast_object(supercell.initialization, 0, world);
    mpi_cereal::bcast_object(monte_carlo, 0, world);

    /**
     * @brief Build the supercell and initialize spins/structure.
     */
    EnlargeCell(supercell);
    InitializeSupercell(supercell);

    // Output the coordinate number
    if(world_rank == 0) {
        /**
         * @brief Write initialization logs and initial spin structure.
         */
        logger.print("Successfully initialize the supercell.\n");
        WriteLog(supercell, monte_carlo, logger);
        WriteSpin(supercell, "structure_initialized");
    }

    /**
     * @brief Containers for thermodynamic observables.
     *
     * energy: internal energy
     * Cv: heat capacity
     * moment: magnetic moment
     * chi: susceptibility
     * *_projection: projected components (field-aligned, etc.)
     */
    vector<double> energy;
    vector<double> Cv;
    vector<double> moment;
    vector<double> chi;
    vector<double> moment_projection;
    vector<double> chi_projection;

    /**
     * @brief Dispatch to MC method chosen in input.
     *
     * Methods::classical: single-temperature Metropolis MC.
     * else: parallel tempering (replica exchange) MC.
     */
    if(monte_carlo.methods == Methods::classical) {
        ClassicalMonteCarlo(world,
            monte_carlo, supercell, spin_structure_file_prefix,
            energy, Cv,
            moment, chi,
            moment_projection, chi_projection);
    } else {
        ParallelTemperingMonteCarlo(world,
            monte_carlo, supercell, spin_structure_file_prefix,
            energy, Cv,
            moment, chi,
            moment_projection, chi_projection);
    }

    // Output the thermal dynamic result using root processor.
    if(world_rank == 0) {
        /**
         * @brief Write final results and completion log.
         */
        logger.print("Successfully run Monte Carlo simulation.\n");
        WriteOutput(monte_carlo, energy, Cv, moment, chi, moment_projection, chi_projection, supercell.lattice.field, output_file);
        logger.print("Successfully output all results.\n");
    }

    MPI_Finalize();
    return 0;
}
