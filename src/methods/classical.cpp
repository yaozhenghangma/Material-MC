#include "classical.h"

#include <string>

/**
 * @file classical.cpp
 * @brief Classical Monte Carlo workflow and MPI temperature sampling.
 */

/**
 * @brief Runs relaxation sweeps at a fixed temperature.
 *
 * @param supercell Simulation supercell.
 * @param monte_carlo Monte Carlo control parameters.
 * @param T Temperature used by local spin updates.
 * @return int Returns 0 on completion.
 */
int MonteCarloRelaxing(Supercell & supercell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    std::vector<int> site_chosen;
    for(int i=0; i<monte_carlo.relax_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            supercell.Update(supercell, supercell[site_chosen], T);
        }
    }
    
    return 0;
}

/**
 * @brief Executes measurement sweeps and accumulates thermal averages.
 *
 * Returned vector layout:
 * [0]=<E>, [1]=<E^2>, [2]=<M>, [3]=<M^2>, [4]=<M_B>, [5]=<M_B^2>.
 *
 * @param supercell Simulation supercell.
 * @param monte_carlo Monte Carlo control parameters.
 * @param T Temperature used by local spin updates.
 * @return std::vector<double> Averaged observables across count_step sweeps.
 */
std::vector<double> MonteCarloStep(Supercell & supercell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    std::vector<int> site_chosen;
    double total_energy = 0;
    double total_energy_square = 0;
    double tmp_momentum = 0;
    double total_momentum = 0;
    double total_momentum_square = 0;
    double tmp_momentum_projection = 0;
    double total_momentum_projection = 0;
    double total_momentum_projection_square = 0;
    std::vector<double> tmp_m_component;
    static double one_over_step = 1.0 / monte_carlo.count_step;
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            supercell.Update(supercell, supercell[site_chosen], T);
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

/**
 * @brief Runs measurement sweeps while tracking the minimum-energy spin state.
 *
 * In ground-state mode, the lowest-energy configuration encountered at this
 * temperature is written via WriteSpin("structure_ground_state", T).
 *
 * @param supercell Simulation supercell.
 * @param monte_carlo Monte Carlo control parameters.
 * @param T Temperature used by local spin updates.
 * @return std::vector<double> Same observable layout as MonteCarloStep().
 */
std::vector<double> MonteCarloStepGroundState(Supercell & supercell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation for outputing ground spin structure mode
    Supercell supercell_ground = supercell;
    double minimum_energy = supercell.lattice.total_energy;
    std::vector<int> site_chosen;
    double total_energy = 0;
    double total_energy_square = 0;
    double tmp_momentum = 0;
    double total_momentum = 0;
    double total_momentum_square = 0;
    double tmp_momentum_projection = 0;
    double total_momentum_projection = 0;
    double total_momentum_projection_square = 0;
    std::vector<double> tmp_m_component;
    static double one_over_step = 1.0 / monte_carlo.count_step;
    static double one_over_number = 1.0 / (supercell.lattice.n_x * supercell.lattice.n_y * \
    supercell.lattice.n_z * supercell.base_site.number);
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            supercell.Update(supercell, supercell[site_chosen], T);
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

        if(supercell.lattice.total_energy < minimum_energy) {
            supercell_ground = supercell;
            minimum_energy = supercell.lattice.total_energy;
        }
    }

    WriteSpin(supercell_ground, "structure_ground_state", T);
    WriteVestaKhBondColor(supercell_ground, "structure_ground_state_kh_bond_color_" + std::to_string(T));
    std::cout << "Minimum Energy:" << minimum_energy * one_over_number << std::endl;
    
    return {total_energy * one_over_step, total_energy_square * one_over_step, \
    total_momentum * one_over_step, total_momentum_square * one_over_step, \
    total_momentum_projection * one_over_step, total_momentum_projection_square * one_over_step};
}

/**
 * @brief Distributes temperature points across MPI ranks and gathers observables.
 *
 * Each rank samples temperatures in a strided pattern and root rank collects
 * energy, heat capacity, moment, susceptibility, and optional field-projected
 * quantities into output vectors.
 *
 * @param world MPI communicator.
 * @param monte_carlo Monte Carlo schedule and sampling controls.
 * @param supercell Simulation supercell.
 * @param spin_structure_file_prefix Unused in this implementation.
 * @param energy Output average energy per spin.
 * @param Cv Output heat capacity per spin.
 * @param moment Output average magnetic moment per spin.
 * @param chi Output magnetic susceptibility per spin.
 * @param moment_projection Output projected moment (when field enabled).
 * @param chi_projection Output projected susceptibility (when field enabled).
 * @return int Returns 0 on completion.
 */
int ClassicalMonteCarlo(MPI_Comm world,
MonteCarlo & monte_carlo, Supercell & supercell, std::string & spin_structure_file_prefix,
std::vector<double> & energy, std::vector<double> & Cv,
std::vector<double> & moment, std::vector<double> & chi,
std::vector<double> & moment_projection, std::vector<double> & chi_projection) {
    // Arrange the processors.
    int world_size = 0;
    int world_rank = 0;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &world_rank);

    // Divide temperature points across ranks as quotient + possible remainder.
    const int quotient = monte_carlo.temperature_step_number / world_size;
    const int remainder = monte_carlo.temperature_step_number % world_size;

    // Monte Carlo
    double T;
    std::vector<double> result_value;
    double energy_every_processor;
    double Cv_every_processor;
    double moment_every_processor;
    double moment_projection_every_processor;
    double chi_projection_every_processor;
    double chi_every_processor;
    std::vector<double> gathered_energy;
    std::vector<double> gathered_Cv;
    std::vector<double> gathered_moment;
    std::vector<double> gathered_chi;
    std::vector<double> gathered_moment_projection;
    std::vector<double> gathered_chi_projection;
    static double one_over_number = 1.0 / (supercell.lattice.n_x * supercell.lattice.n_y * \
    supercell.lattice.n_z * supercell.base_site.number);
    for(int i=0; i<quotient+1; i++) {
        if(i<quotient || remainder !=0 ) {
            // Rank r handles temperature index (i * world_size + r).
            T = monte_carlo.start_temperature + (i*world_size+world_rank)*monte_carlo.temperature_step;
            MonteCarloRelaxing(supercell, monte_carlo, T);
            if(supercell.lattice.ground_state && T == 0.0) {
                result_value = MonteCarloStepGroundState(supercell, monte_carlo, T);
            } else {
                result_value = MonteCarloStep(supercell, monte_carlo, T);
            }

            // Convert accumulated moments into thermodynamic observables per spin.
            Cv_every_processor = (result_value[1]-result_value[0]*result_value[0])*one_over_number/(KB*T*T); //Cv
            chi_every_processor = (result_value[3]-result_value[2]*result_value[2])*one_over_number/(KB*T); //chi
            energy_every_processor = result_value[0] * one_over_number; //energy
            moment_every_processor = result_value[2] * one_over_number; //moment
            if(supercell.lattice.field) {
                chi_projection_every_processor = (result_value[5] - result_value[4]*result_value[4])/(KB*T);
                moment_projection_every_processor = result_value[4] * one_over_number;
            }

            // Collect data form all processors.
            gathered_energy.resize(world_size);
            gathered_Cv.resize(world_size);
            gathered_moment.resize(world_size);
            gathered_chi.resize(world_size);
            if(supercell.lattice.field) {
                gathered_moment_projection.resize(world_size);
                gathered_chi_projection.resize(world_size);
            }

            MPI_Gather(&energy_every_processor, 1, MPI_DOUBLE,
                       gathered_energy.data(), 1, MPI_DOUBLE, 0, world);
            MPI_Gather(&Cv_every_processor, 1, MPI_DOUBLE,
                       gathered_Cv.data(), 1, MPI_DOUBLE, 0, world);
            MPI_Gather(&moment_every_processor, 1, MPI_DOUBLE,
                       gathered_moment.data(), 1, MPI_DOUBLE, 0, world);
            MPI_Gather(&chi_every_processor, 1, MPI_DOUBLE,
                       gathered_chi.data(), 1, MPI_DOUBLE, 0, world);
            if(supercell.lattice.field) {
                MPI_Gather(&moment_projection_every_processor, 1, MPI_DOUBLE,
                           gathered_moment_projection.data(), 1, MPI_DOUBLE, 0, world);
                MPI_Gather(&chi_projection_every_processor, 1, MPI_DOUBLE,
                           gathered_chi_projection.data(), 1, MPI_DOUBLE, 0, world);
            }

            // Store these data
            if(world_rank == 0) {
                for(int j=0; j<world_size; j++) {
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

    return 0;
}