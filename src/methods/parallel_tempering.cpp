#include "parallel_tempering.h"

int ParallelTemperingMonteCarloRelaxing(boost::mpi::environment & env, boost::mpi::communicator & world, 
Supercell & supercell, MonteCarlo & monte_carlo, double & T) {
    bool odd = true;
    static bool even_number_process = (world.size()%2 == 0);
    double energy = 0;
    double exchange_energy = 0;
    double exchange_T = 0;
    double temp_T;
    static int max_num = world.size()-1;
    bool exchange;

    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    std::vector<int> site_chosen;
    for(int i=0; i<monte_carlo.relax_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            LocalUpdate(supercell, supercell[site_chosen], T);
        }

        if(i%monte_carlo.replica_exchange_step_number == 0) {
            energy = supercell.lattice.total_energy;
            exchange_energy = supercell.lattice.total_energy;
            exchange_T = T;
            if(odd) {
                if(world.rank()%2 == 0 && world.rank() != max_num) {
                    world.recv(world.rank()+1, 0, exchange_energy);
                    world.recv(world.rank()+1, 1, exchange_T);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    world.send(world.rank()+1, 2, exchange);

                    if(exchange) {
                        world.send(world.rank()+1, 3, exchange_T);
                    }
                } else if(even_number_process || world.rank() != max_num) {
                    world.send(world.rank()-1, 0, exchange_energy);
                    world.send(world.rank()-1, 1, exchange_T);

                    world.recv(world.rank()-1, 2, exchange);

                    if(exchange) {
                        world.recv(world.rank()-1, 3, T);
                    }
                }

                if(world.rank()==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = false;
                    }
                }
                broadcast(world, odd, 0);
            } else {
                if(world.rank()%2 == 1 && world.rank() != max_num) {
                    world.recv(world.rank()+1, 0, exchange_energy);
                    world.recv(world.rank()+1, 1, exchange_T);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    world.send(world.rank()+1, 2, exchange);

                    if(exchange) {
                        world.send(world.rank()+1, 3, exchange_T);
                    }
                } else if(world.rank() != 0 && (!even_number_process || world.rank() != max_num)) {
                    world.send(world.rank()-1, 0, exchange_energy);
                    world.send(world.rank()-1, 1, exchange_T);

                    world.recv(world.rank()-1, 2, exchange);

                    if(exchange) {
                        world.recv(world.rank()-1, 3, T);
                    }
                }

                if(world.rank()==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = true;
                    }
                }
                broadcast(world, odd, 0);
            }
        }
    }
    
    return 0;
}

std::vector<double> ParallelTemperingMonteCarloStep(boost::mpi::environment & env, boost::mpi::communicator & world, 
Supercell & supercell, MonteCarlo & monte_carlo, double & T) {
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

    double exchange_total_energy;
    double exchange_total_energy_square;
    double exchange_momentum;
    double exchange_momentum_square;
    double exchange_momentum_projection;
    double exchange_momentum_projection_square;
    double energy;
    double exchange_energy;
    double exchange_T;
    double temp_T;
    static bool even_number_process = (world.size()%2 == 0);
    bool odd = true;
    bool exchange;
    static int max_num = world.size()-1;
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            LocalUpdate(supercell, supercell[site_chosen], T);
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

        if(i%monte_carlo.replica_exchange_step_number == 0) {
            energy = supercell.lattice.total_energy;
            exchange_energy = supercell.lattice.total_energy;
            exchange_T = T;
            if(odd) {
                if(world.rank()%2 == 0 && world.rank() != max_num) {
                    world.recv(world.rank()+1, 0, exchange_energy);
                    world.recv(world.rank()+1, 1, exchange_T);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    world.send(world.rank()+1, 2, exchange);

                    if(exchange) {
                        world.send(world.rank()+1, 3, exchange_T);
                        world.send(world.rank()+1, 4, total_energy);
                        world.send(world.rank()+1, 5, total_energy_square);
                        world.send(world.rank()+1, 6, total_momentum);
                        world.send(world.rank()+1, 7, total_momentum_square);
                        world.send(world.rank()+1, 8, total_momentum_projection);
                        world.send(world.rank()+1, 9, total_momentum_projection_square);

                        world.recv(world.rank()+1, 10, exchange_total_energy);
                        world.recv(world.rank()+1, 11, exchange_total_energy_square);
                        world.recv(world.rank()+1, 12, exchange_momentum);
                        world.recv(world.rank()+1, 13, exchange_momentum_square);
                        world.recv(world.rank()+1, 14, exchange_momentum_projection);
                        world.recv(world.rank()+1, 15, exchange_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(even_number_process || world.rank() != max_num) {
                    world.send(world.rank()-1, 0, exchange_energy);
                    world.send(world.rank()-1, 1, exchange_T);

                    world.recv(world.rank()-1, 2, exchange);

                    if(exchange) {
                        world.recv(world.rank()-1, 3, T);
                        world.recv(world.rank()-1, 4, exchange_total_energy);
                        world.recv(world.rank()-1, 5, exchange_total_energy_square);
                        world.recv(world.rank()-1, 6, exchange_momentum);
                        world.recv(world.rank()-1, 7, exchange_momentum_square);
                        world.recv(world.rank()-1, 8, exchange_momentum_projection);
                        world.recv(world.rank()-1, 9, exchange_momentum_projection_square);

                        world.send(world.rank()-1, 10, total_energy);
                        world.send(world.rank()-1, 11, total_energy_square);
                        world.send(world.rank()-1, 12, total_momentum);
                        world.send(world.rank()-1, 13, total_momentum_square);
                        world.send(world.rank()-1, 14, total_momentum_projection);
                        world.send(world.rank()-1, 15, total_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world.rank()==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = false;
                    }
                }
                broadcast(world, odd, 0);
            } else {
                if(world.rank()%2 == 1 && world.rank() != max_num) {
                    world.recv(world.rank()+1, 0, exchange_energy);
                    world.recv(world.rank()+1, 1, exchange_T);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    world.send(world.rank()+1, 2, exchange);

                    if(exchange) {
                        world.send(world.rank()+1, 3, exchange_T);
                        world.send(world.rank()+1, 4, total_energy);
                        world.send(world.rank()+1, 5, total_energy_square);
                        world.send(world.rank()+1, 6, total_momentum);
                        world.send(world.rank()+1, 7, total_momentum_square);
                        world.send(world.rank()+1, 8, total_momentum_projection);
                        world.send(world.rank()+1, 9, total_momentum_projection_square);

                        world.recv(world.rank()+1, 10, exchange_total_energy);
                        world.recv(world.rank()+1, 11, exchange_total_energy_square);
                        world.recv(world.rank()+1, 12, exchange_momentum);
                        world.recv(world.rank()+1, 13, exchange_momentum_square);
                        world.recv(world.rank()+1, 14, exchange_momentum_projection);
                        world.recv(world.rank()+1, 15, exchange_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(world.rank() != 0 && (!even_number_process || world.rank() != max_num)) {
                    world.send(world.rank()-1, 0, exchange_energy);
                    world.send(world.rank()-1, 1, exchange_T);

                    world.recv(world.rank()-1, 2, exchange);

                    if(exchange) {
                        world.recv(world.rank()-1, 3, T);
                        world.recv(world.rank()-1, 4, exchange_total_energy);
                        world.recv(world.rank()-1, 5, exchange_total_energy_square);
                        world.recv(world.rank()-1, 6, exchange_momentum);
                        world.recv(world.rank()-1, 7, exchange_momentum_square);
                        world.recv(world.rank()-1, 8, exchange_momentum_projection);
                        world.recv(world.rank()-1, 9, exchange_momentum_projection_square);

                        world.send(world.rank()-1, 10, total_energy);
                        world.send(world.rank()-1, 11, total_energy_square);
                        world.send(world.rank()-1, 12, total_momentum);
                        world.send(world.rank()-1, 13, total_momentum_square);
                        world.send(world.rank()-1, 14, total_momentum_projection);
                        world.send(world.rank()-1, 15, total_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world.rank()==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = true;
                    }
                }
                broadcast(world, odd, 0);
            }
        }
    }
    
    return {total_energy * one_over_step, total_energy_square * one_over_step, \
    total_momentum * one_over_step, total_momentum_square * one_over_step, \
    total_momentum_projection * one_over_step, total_momentum_projection_square * one_over_step};
}

std::vector<double> ParallelTemperingMonteCarloStepGroundState(boost::mpi::environment & env, boost::mpi::communicator & world, 
Supercell & supercell, MonteCarlo & monte_carlo, double & T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    Supercell supercell_ground = supercell;
    double minimum_energy = supercell.lattice.total_energy;
    double global_minimum_energy;
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

    double exchange_total_energy;
    double exchange_total_energy_square;
    double exchange_momentum;
    double exchange_momentum_square;
    double exchange_momentum_projection;
    double exchange_momentum_projection_square;
    double energy;
    double exchange_energy;
    double exchange_T;
    double temp_T;
    static bool even_number_process = (world.size()%2 == 0);
    bool odd = true;
    bool exchange;
    static int max_num = world.size()-1;
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            LocalUpdate(supercell, supercell[site_chosen], T);
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
            minimum_energy = supercell.lattice.total_energy;
            supercell_ground = supercell;
        }

        if(i%monte_carlo.replica_exchange_step_number == 0) {
            energy = supercell.lattice.total_energy;
            exchange_energy = supercell.lattice.total_energy;
            exchange_T = T;
            if(odd) {
                if(world.rank()%2 == 0 && world.rank() != max_num) {
                    world.recv(world.rank()+1, 0, exchange_energy);
                    world.recv(world.rank()+1, 1, exchange_T);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    world.send(world.rank()+1, 2, exchange);

                    if(exchange) {
                        world.send(world.rank()+1, 3, exchange_T);
                        world.send(world.rank()+1, 4, total_energy);
                        world.send(world.rank()+1, 5, total_energy_square);
                        world.send(world.rank()+1, 6, total_momentum);
                        world.send(world.rank()+1, 7, total_momentum_square);
                        world.send(world.rank()+1, 8, total_momentum_projection);
                        world.send(world.rank()+1, 9, total_momentum_projection_square);

                        world.recv(world.rank()+1, 10, exchange_total_energy);
                        world.recv(world.rank()+1, 11, exchange_total_energy_square);
                        world.recv(world.rank()+1, 12, exchange_momentum);
                        world.recv(world.rank()+1, 13, exchange_momentum_square);
                        world.recv(world.rank()+1, 14, exchange_momentum_projection);
                        world.recv(world.rank()+1, 15, exchange_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(even_number_process || world.rank() != max_num) {
                    world.send(world.rank()-1, 0, exchange_energy);
                    world.send(world.rank()-1, 1, exchange_T);

                    world.recv(world.rank()-1, 2, exchange);

                    if(exchange) {
                        world.recv(world.rank()-1, 3, T);
                        world.recv(world.rank()-1, 4, exchange_total_energy);
                        world.recv(world.rank()-1, 5, exchange_total_energy_square);
                        world.recv(world.rank()-1, 6, exchange_momentum);
                        world.recv(world.rank()-1, 7, exchange_momentum_square);
                        world.recv(world.rank()-1, 8, exchange_momentum_projection);
                        world.recv(world.rank()-1, 9, exchange_momentum_projection_square);

                        world.send(world.rank()-1, 10, total_energy);
                        world.send(world.rank()-1, 11, total_energy_square);
                        world.send(world.rank()-1, 12, total_momentum);
                        world.send(world.rank()-1, 13, total_momentum_square);
                        world.send(world.rank()-1, 14, total_momentum_projection);
                        world.send(world.rank()-1, 15, total_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world.rank()==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = false;
                    }
                }
                broadcast(world, odd, 0);
            } else {
                if(world.rank()%2 == 1 && world.rank() != max_num) {
                    world.recv(world.rank()+1, 0, exchange_energy);
                    world.recv(world.rank()+1, 1, exchange_T);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    world.send(world.rank()+1, 2, exchange);

                    if(exchange) {
                        world.send(world.rank()+1, 3, exchange_T);
                        world.send(world.rank()+1, 4, total_energy);
                        world.send(world.rank()+1, 5, total_energy_square);
                        world.send(world.rank()+1, 6, total_momentum);
                        world.send(world.rank()+1, 7, total_momentum_square);
                        world.send(world.rank()+1, 8, total_momentum_projection);
                        world.send(world.rank()+1, 9, total_momentum_projection_square);

                        world.recv(world.rank()+1, 10, exchange_total_energy);
                        world.recv(world.rank()+1, 11, exchange_total_energy_square);
                        world.recv(world.rank()+1, 12, exchange_momentum);
                        world.recv(world.rank()+1, 13, exchange_momentum_square);
                        world.recv(world.rank()+1, 14, exchange_momentum_projection);
                        world.recv(world.rank()+1, 15, exchange_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(world.rank() != 0 && (!even_number_process || world.rank() != max_num)) {
                    world.send(world.rank()-1, 0, exchange_energy);
                    world.send(world.rank()-1, 1, exchange_T);

                    world.recv(world.rank()-1, 2, exchange);

                    if(exchange) {
                        world.recv(world.rank()-1, 3, T);
                        world.recv(world.rank()-1, 4, exchange_total_energy);
                        world.recv(world.rank()-1, 5, exchange_total_energy_square);
                        world.recv(world.rank()-1, 6, exchange_momentum);
                        world.recv(world.rank()-1, 7, exchange_momentum_square);
                        world.recv(world.rank()-1, 8, exchange_momentum_projection);
                        world.recv(world.rank()-1, 9, exchange_momentum_projection_square);

                        world.send(world.rank()-1, 10, total_energy);
                        world.send(world.rank()-1, 11, total_energy_square);
                        world.send(world.rank()-1, 12, total_momentum);
                        world.send(world.rank()-1, 13, total_momentum_square);
                        world.send(world.rank()-1, 14, total_momentum_projection);
                        world.send(world.rank()-1, 15, total_momentum_projection_square);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world.rank()==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = true;
                    }
                }
                broadcast(world, odd, 0);
            }
        }
    }

    boost::mpi::all_reduce(world, minimum_energy, global_minimum_energy, boost::mpi::minimum<double>());
    if(minimum_energy == global_minimum_energy) {
        WriteSpin(supercell_ground, "structure_ground_state");
        std::cout << "Minimum Energy: " << global_minimum_energy << std::endl;
    }
    
    return {total_energy * one_over_step, total_energy_square * one_over_step, \
    total_momentum * one_over_step, total_momentum_square * one_over_step, \
    total_momentum_projection * one_over_step, total_momentum_projection_square * one_over_step};
}

int ParallelTemperingMonteCarlo(boost::mpi::environment & env, boost::mpi::communicator & world, 
MonteCarlo & monte_carlo, Supercell & supercell, std::string & spin_structure_file_prefix, 
std::vector<double> & energy, std::vector<double> & Cv,
std::vector<double> & moment, std::vector<double> & chi,
std::vector<double> & moment_projection, std::vector<double> & chi_projection) {
    // Arrange the processors.
    const int quotient = monte_carlo.temperature_step_number / world.size();
    const int remainder = monte_carlo.temperature_step_number % world.size();

    // Monte Carlo
    double T;
    std::vector<double> result_value;
    double energy_every_processor;
    double Cv_every_processor;
    double moment_every_processor;
    double moment_projection_every_processor;
    double chi_projection_every_processor;
    double chi_every_processor;
    std::vector<double> gathered_T;
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
            T = monte_carlo.start_temperature + (i*world.size()+world.rank())*monte_carlo.temperature_step;
            ParallelTemperingMonteCarloRelaxing(env, world, supercell, monte_carlo, T);
            if(supercell.lattice.ground_state) {
                result_value = ParallelTemperingMonteCarloStepGroundState(env, world, supercell, monte_carlo, T);
            } else {
                result_value = ParallelTemperingMonteCarloStep(env, world, supercell, monte_carlo, T);
            }

            Cv_every_processor = (result_value[1]-result_value[0]*result_value[0])*one_over_number/(KB*T*T); //Cv
            chi_every_processor = (result_value[3]-result_value[2]*result_value[2])/(KB*T); //chi
            energy_every_processor = result_value[0] * one_over_number; //energy
            moment_every_processor = result_value[2] * one_over_number; //moment
            if(supercell.lattice.field) {
                chi_projection_every_processor = (result_value[5] - result_value[4]*result_value[4])/(KB*T);
                moment_projection_every_processor = result_value[4] * one_over_number;
            }
            
            // Collect data form all processors.
            gather(world, T, gathered_T, 0);
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
                std::vector<int> index(world.size());
                for(int j=0; j<world.size(); j++) {
                    index[round((gathered_T[j]-monte_carlo.start_temperature)/monte_carlo.temperature_step)-i*world.size()] = j;
                }

                for(int j=0; j<world.size(); j++) {
                    energy.push_back(gathered_energy[index[j]]);
                    Cv.push_back(gathered_Cv[index[j]]);
                    moment.push_back(gathered_moment[index[j]]);
                    chi.push_back(gathered_chi[index[j]]);
                    if(supercell.lattice.field) {
                        moment_projection.push_back(gathered_moment_projection[index[j]]);
                        chi_projection.push_back(gathered_chi_projection[index[j]]);
                    }
                }
            }
            
        }
    }

    return 0;
}