#include "parallel_tempering.h"

int ParallelTemperingMonteCarloRelaxing(MPI_Comm world,
Supercell & supercell, MonteCarlo & monte_carlo, double & T) {
    bool odd = true;
    int world_size = 0;
    int world_rank = 0;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &world_rank);

    static bool even_number_process = (world_size%2 == 0);
    double energy = 0;
    double exchange_energy = 0;
    double exchange_T = 0;
    double temp_T;
    static int max_num = world_size-1;
    bool exchange;

    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    std::vector<int> site_chosen;
    for(int i=0; i<monte_carlo.relax_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            supercell.Update(supercell, supercell[site_chosen], T);
        }

        if(i%monte_carlo.replica_exchange_step_number == 0) {
            energy = supercell.lattice.total_energy;
            exchange_energy = supercell.lattice.total_energy;
            exchange_T = T;
            if(odd) {
                if(world_rank%2 == 0 && world_rank != max_num) {
                    MPI_Recv(&exchange_energy, 1, MPI_DOUBLE, world_rank+1, 0, world, MPI_STATUS_IGNORE);
                    MPI_Recv(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 1, world, MPI_STATUS_IGNORE);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    MPI_Send(&exchange, 1, MPI_C_BOOL, world_rank+1, 2, world);

                    if(exchange) {
                        MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 3, world);
                    }
                } else if(even_number_process || world_rank != max_num) {
                    MPI_Send(&exchange_energy, 1, MPI_DOUBLE, world_rank-1, 0, world);
                    MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank-1, 1, world);

                    MPI_Recv(&exchange, 1, MPI_C_BOOL, world_rank-1, 2, world, MPI_STATUS_IGNORE);

                    if(exchange) {
                        MPI_Recv(&T, 1, MPI_DOUBLE, world_rank-1, 3, world, MPI_STATUS_IGNORE);
                    }
                }

                if(world_rank==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = false;
                    }
                }
                MPI_Bcast(&odd, 1, MPI_C_BOOL, 0, world);
            } else {
                if(world_rank%2 == 1 && world_rank != max_num) {
                    MPI_Recv(&exchange_energy, 1, MPI_DOUBLE, world_rank+1, 0, world, MPI_STATUS_IGNORE);
                    MPI_Recv(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 1, world, MPI_STATUS_IGNORE);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    MPI_Send(&exchange, 1, MPI_C_BOOL, world_rank+1, 2, world);

                    if(exchange) {
                        MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 3, world);
                    }
                } else if(world_rank != 0 && (!even_number_process || world_rank != max_num)) {
                    MPI_Send(&exchange_energy, 1, MPI_DOUBLE, world_rank-1, 0, world);
                    MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank-1, 1, world);

                    MPI_Recv(&exchange, 1, MPI_C_BOOL, world_rank-1, 2, world, MPI_STATUS_IGNORE);

                    if(exchange) {
                        MPI_Recv(&T, 1, MPI_DOUBLE, world_rank-1, 3, world, MPI_STATUS_IGNORE);
                    }
                }

                if(world_rank==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = true;
                    }
                }
                MPI_Bcast(&odd, 1, MPI_C_BOOL, 0, world);
            }
        }
    }

    return 0;
}

std::vector<double> ParallelTemperingMonteCarloStep(MPI_Comm world,
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
    int world_size = 0;
    int world_rank = 0;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &world_rank);
    static bool even_number_process = (world_size%2 == 0);
    bool odd = true;
    bool exchange;
    static int max_num = world_size-1;
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

        if(i%monte_carlo.replica_exchange_step_number == 0) {
            energy = supercell.lattice.total_energy;
            exchange_energy = supercell.lattice.total_energy;
            exchange_T = T;
            if(odd) {
                if(world_rank%2 == 0 && world_rank != max_num) {
                    MPI_Recv(&exchange_energy, 1, MPI_DOUBLE, world_rank+1, 0, world, MPI_STATUS_IGNORE);
                    MPI_Recv(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 1, world, MPI_STATUS_IGNORE);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    MPI_Send(&exchange, 1, MPI_C_BOOL, world_rank+1, 2, world);

                    if(exchange) {
                        MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 3, world);
                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank+1, 4, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank+1, 5, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank+1, 6, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank+1, 7, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 8, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 9, world);

                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank+1, 10, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank+1, 11, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank+1, 12, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank+1, 13, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 14, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 15, world, MPI_STATUS_IGNORE);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(even_number_process || world_rank != max_num) {
                    MPI_Send(&exchange_energy, 1, MPI_DOUBLE, world_rank-1, 0, world);
                    MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank-1, 1, world);

                    MPI_Recv(&exchange, 1, MPI_C_BOOL, world_rank-1, 2, world, MPI_STATUS_IGNORE);

                    if(exchange) {
                        MPI_Recv(&T, 1, MPI_DOUBLE, world_rank-1, 3, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank-1, 4, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank-1, 5, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank-1, 6, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank-1, 7, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 8, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 9, world, MPI_STATUS_IGNORE);

                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank-1, 10, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank-1, 11, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank-1, 12, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank-1, 13, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 14, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 15, world);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world_rank==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = false;
                    }
                }
                MPI_Bcast(&odd, 1, MPI_C_BOOL, 0, world);
            } else {
                if(world_rank%2 == 1 && world_rank != max_num) {
                    MPI_Recv(&exchange_energy, 1, MPI_DOUBLE, world_rank+1, 0, world, MPI_STATUS_IGNORE);
                    MPI_Recv(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 1, world, MPI_STATUS_IGNORE);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    MPI_Send(&exchange, 1, MPI_C_BOOL, world_rank+1, 2, world);

                    if(exchange) {
                        MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 3, world);
                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank+1, 4, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank+1, 5, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank+1, 6, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank+1, 7, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 8, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 9, world);

                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank+1, 10, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank+1, 11, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank+1, 12, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank+1, 13, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 14, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 15, world, MPI_STATUS_IGNORE);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(world_rank != 0 && (!even_number_process || world_rank != max_num)) {
                    MPI_Send(&exchange_energy, 1, MPI_DOUBLE, world_rank-1, 0, world);
                    MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank-1, 1, world);

                    MPI_Recv(&exchange, 1, MPI_C_BOOL, world_rank-1, 2, world, MPI_STATUS_IGNORE);

                    if(exchange) {
                        MPI_Recv(&T, 1, MPI_DOUBLE, world_rank-1, 3, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank-1, 4, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank-1, 5, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank-1, 6, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank-1, 7, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 8, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 9, world, MPI_STATUS_IGNORE);

                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank-1, 10, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank-1, 11, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank-1, 12, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank-1, 13, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 14, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 15, world);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world_rank==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = true;
                    }
                }
                MPI_Bcast(&odd, 1, MPI_C_BOOL, 0, world);
            }
        }
    }

    return {total_energy * one_over_step, total_energy_square * one_over_step, \
    total_momentum * one_over_step, total_momentum_square * one_over_step, \
    total_momentum_projection * one_over_step, total_momentum_projection_square * one_over_step};
}

std::vector<double> ParallelTemperingMonteCarloStepGroundState(MPI_Comm world,
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
    static double one_over_number = 1.0 / (supercell.lattice.n_x * supercell.lattice.n_y * \
    supercell.lattice.n_z * supercell.base_site.number);

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
    int world_size = 0;
    int world_rank = 0;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &world_rank);
    static bool even_number_process = (world_size%2 == 0);
    bool odd = true;
    bool exchange;
    static int max_num = world_size-1;
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
            minimum_energy = supercell.lattice.total_energy;
            supercell_ground = supercell;
        }

        if(i%monte_carlo.replica_exchange_step_number == 0) {
            energy = supercell.lattice.total_energy;
            exchange_energy = supercell.lattice.total_energy;
            exchange_T = T;
            if(odd) {
                if(world_rank%2 == 0 && world_rank != max_num) {
                    MPI_Recv(&exchange_energy, 1, MPI_DOUBLE, world_rank+1, 0, world, MPI_STATUS_IGNORE);
                    MPI_Recv(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 1, world, MPI_STATUS_IGNORE);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    MPI_Send(&exchange, 1, MPI_C_BOOL, world_rank+1, 2, world);

                    if(exchange) {
                        MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 3, world);
                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank+1, 4, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank+1, 5, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank+1, 6, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank+1, 7, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 8, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 9, world);

                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank+1, 10, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank+1, 11, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank+1, 12, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank+1, 13, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 14, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 15, world, MPI_STATUS_IGNORE);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(even_number_process || world_rank != max_num) {
                    MPI_Send(&exchange_energy, 1, MPI_DOUBLE, world_rank-1, 0, world);
                    MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank-1, 1, world);

                    MPI_Recv(&exchange, 1, MPI_C_BOOL, world_rank-1, 2, world, MPI_STATUS_IGNORE);

                    if(exchange) {
                        MPI_Recv(&T, 1, MPI_DOUBLE, world_rank-1, 3, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank-1, 4, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank-1, 5, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank-1, 6, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank-1, 7, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 8, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 9, world, MPI_STATUS_IGNORE);

                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank-1, 10, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank-1, 11, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank-1, 12, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank-1, 13, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 14, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 15, world);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world_rank==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = false;
                    }
                }
                MPI_Bcast(&odd, 1, MPI_C_BOOL, 0, world);
            } else {
                if(world_rank%2 == 1 && world_rank != max_num) {
                    MPI_Recv(&exchange_energy, 1, MPI_DOUBLE, world_rank+1, 0, world, MPI_STATUS_IGNORE);
                    MPI_Recv(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 1, world, MPI_STATUS_IGNORE);

                    if (RandomFloat() > exp(-one_over_KB*(1.0/exchange_T - 1.0/T)*(energy - exchange_energy))) {
                        exchange = false;
                    } else {
                        exchange = true;
                        temp_T = T;
                        T = exchange_T;
                        exchange_T = temp_T;
                    }

                    MPI_Send(&exchange, 1, MPI_C_BOOL, world_rank+1, 2, world);

                    if(exchange) {
                        MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank+1, 3, world);
                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank+1, 4, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank+1, 5, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank+1, 6, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank+1, 7, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 8, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 9, world);

                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank+1, 10, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank+1, 11, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank+1, 12, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank+1, 13, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank+1, 14, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank+1, 15, world, MPI_STATUS_IGNORE);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                } else if(world_rank != 0 && (!even_number_process || world_rank != max_num)) {
                    MPI_Send(&exchange_energy, 1, MPI_DOUBLE, world_rank-1, 0, world);
                    MPI_Send(&exchange_T, 1, MPI_DOUBLE, world_rank-1, 1, world);

                    MPI_Recv(&exchange, 1, MPI_C_BOOL, world_rank-1, 2, world, MPI_STATUS_IGNORE);

                    if(exchange) {
                        MPI_Recv(&T, 1, MPI_DOUBLE, world_rank-1, 3, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy, 1, MPI_DOUBLE, world_rank-1, 4, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_total_energy_square, 1, MPI_DOUBLE, world_rank-1, 5, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum, 1, MPI_DOUBLE, world_rank-1, 6, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_square, 1, MPI_DOUBLE, world_rank-1, 7, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 8, world, MPI_STATUS_IGNORE);
                        MPI_Recv(&exchange_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 9, world, MPI_STATUS_IGNORE);

                        MPI_Send(&total_energy, 1, MPI_DOUBLE, world_rank-1, 10, world);
                        MPI_Send(&total_energy_square, 1, MPI_DOUBLE, world_rank-1, 11, world);
                        MPI_Send(&total_momentum, 1, MPI_DOUBLE, world_rank-1, 12, world);
                        MPI_Send(&total_momentum_square, 1, MPI_DOUBLE, world_rank-1, 13, world);
                        MPI_Send(&total_momentum_projection, 1, MPI_DOUBLE, world_rank-1, 14, world);
                        MPI_Send(&total_momentum_projection_square, 1, MPI_DOUBLE, world_rank-1, 15, world);

                        total_energy = exchange_total_energy;
                        total_energy_square = exchange_total_energy_square;
                        total_momentum = exchange_momentum;
                        total_momentum_square = exchange_momentum_square;
                        total_momentum_projection = exchange_momentum_projection;
                        total_momentum_projection_square = exchange_momentum_projection_square;
                    }
                }

                if(world_rank==0) {
                    if(RandomFloat() >= 0.5) {
                        odd = true;
                    }
                }
                MPI_Bcast(&odd, 1, MPI_C_BOOL, 0, world);
            }
        }
    }

    MPI_Allreduce(&minimum_energy, &global_minimum_energy, 1, MPI_DOUBLE, MPI_MIN, world);
    if(minimum_energy == global_minimum_energy) {
        WriteSpin(supercell_ground, "structure_ground_state");
        std::cout << "Minimum Energy: " << global_minimum_energy * one_over_number << std::endl;
    }

    return {total_energy * one_over_step, total_energy_square * one_over_step, \
    total_momentum * one_over_step, total_momentum_square * one_over_step, \
    total_momentum_projection * one_over_step, total_momentum_projection_square * one_over_step};
}

int ParallelTemperingMonteCarlo(MPI_Comm world,
MonteCarlo & monte_carlo, Supercell & supercell, std::string & spin_structure_file_prefix,
std::vector<double> & energy, std::vector<double> & Cv,
std::vector<double> & moment, std::vector<double> & chi,
std::vector<double> & moment_projection, std::vector<double> & chi_projection) {
    // Arrange the processors.
    int world_size = 0;
    int world_rank = 0;
    MPI_Comm_size(world, &world_size);
    MPI_Comm_rank(world, &world_rank);

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
            T = monte_carlo.start_temperature + (i*world_size+world_rank)*monte_carlo.temperature_step;
            ParallelTemperingMonteCarloRelaxing(world, supercell, monte_carlo, T);
            if(supercell.lattice.ground_state) {
                result_value = ParallelTemperingMonteCarloStepGroundState(world, supercell, monte_carlo, T);
            } else {
                result_value = ParallelTemperingMonteCarloStep(world, supercell, monte_carlo, T);
            }

            Cv_every_processor = (result_value[1]-result_value[0]*result_value[0])*one_over_number/(KB*T*T); //Cv
            chi_every_processor = (result_value[3]-result_value[2]*result_value[2])*one_over_number/(KB*T); //chi
            energy_every_processor = result_value[0] * one_over_number; //energy
            moment_every_processor = result_value[2] * one_over_number; //moment
            if(supercell.lattice.field) {
                chi_projection_every_processor = (result_value[5] - result_value[4]*result_value[4])/(KB*T);
                moment_projection_every_processor = result_value[4] * one_over_number;
            }

            // Collect data form all processors.
            gathered_T.resize(world_size);
            gathered_energy.resize(world_size);
            gathered_Cv.resize(world_size);
            gathered_moment.resize(world_size);
            gathered_chi.resize(world_size);
            if(supercell.lattice.field) {
                gathered_moment_projection.resize(world_size);
                gathered_chi_projection.resize(world_size);
            }

            MPI_Gather(&T, 1, MPI_DOUBLE, gathered_T.data(), 1, MPI_DOUBLE, 0, world);
            MPI_Gather(&energy_every_processor, 1, MPI_DOUBLE, gathered_energy.data(), 1, MPI_DOUBLE, 0, world);
            MPI_Gather(&Cv_every_processor, 1, MPI_DOUBLE, gathered_Cv.data(), 1, MPI_DOUBLE, 0, world);
            MPI_Gather(&moment_every_processor, 1, MPI_DOUBLE, gathered_moment.data(), 1, MPI_DOUBLE, 0, world);
            MPI_Gather(&chi_every_processor, 1, MPI_DOUBLE, gathered_chi.data(), 1, MPI_DOUBLE, 0, world);
            if(supercell.lattice.field) {
                MPI_Gather(&moment_projection_every_processor, 1, MPI_DOUBLE,
                           gathered_moment_projection.data(), 1, MPI_DOUBLE, 0, world);
                MPI_Gather(&chi_projection_every_processor, 1, MPI_DOUBLE,
                           gathered_chi_projection.data(), 1, MPI_DOUBLE, 0, world);
            }

            // Store these data
            if(world_rank == 0) {
                std::vector<int> index(world_size);
                for(int j=0; j<world_size; j++) {
                    index[round((gathered_T[j]-monte_carlo.start_temperature)/monte_carlo.temperature_step)-i*world_size] = j;
                }

                for(int j=0; j<world_size; j++) {
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
