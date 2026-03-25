#include "random_function.h"

/**
 * @file random_function.cpp
 * @brief Random sampling utilities used in Monte Carlo updates.
 */

/**
 * @brief Randomly picks one lattice site index.
 *
 * The returned vector format is {x, y, z, base_index}.
 *
 * @param n_x Lattice size along x.
 * @param n_y Lattice size along y.
 * @param n_z Lattice size along z.
 * @param base_n Number of base sites in one unit cell.
 * @return A 4-element index vector {x, y, z, base_index}.
 *
 * @note The random distributions are static and initialized on first call.
 */
std::vector<int> RandomSite(int n_x, int n_y, int n_z, int base_n) {
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::uniform_int_distribution<int> int_distribution_x(0, n_x-1);
    static std::uniform_int_distribution<int> int_distribution_y(0, n_y-1);
    static std::uniform_int_distribution<int> int_distribution_z(0, n_z-1);
    static std::uniform_int_distribution<int> int_distribution_base(0, base_n-1);

    std::vector<int> index = {0, 0, 0, 0};

    index[0] = int_distribution_x(engine);
    index[1] = int_distribution_y(engine);
    index[2] = int_distribution_z(engine);
    index[3] = int_distribution_base(engine);

    return index;
}

/**
 * @brief Generates a random 3D spin vector with a target magnitude.
 *
 * Components are sampled from a standard normal distribution,
 * then normalized and scaled by @p scaling.
 *
 * @param scaling Target vector magnitude.
 * @return A 3D vector {sx, sy, sz}.
 */
std::vector<double> RandomSpin(double scaling) {
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::normal_distribution<double> normal{0, 1};

    double x1 = normal(engine);
    double x2 = normal(engine);
    double x3 = normal(engine);
    double factor = scaling/sqrt(x1*x1+x2*x2+x3*x3);

    return {x1*factor, x2*factor, x3*factor};
}

/**
 * @brief Generates a uniform random floating-point value in [0, 1).
 *
 * @return A random double in [0, 1).
 */
double RandomFloat() {
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::uniform_real_distribution<double> double_distribution(0, 1);

    return double_distribution(engine);
}