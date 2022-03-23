#include "random_function.h"

std::vector<int> RandomSite(int n_x, int n_y, int n_z, int base_n) {
    // Return the site index randomly.
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

std::vector<double> RandomSpin(double scaling) {
    // Return a spin in 3d space randomly.
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::normal_distribution<double> normal{0, 1};

    double x1 = normal(engine);
    double x2 = normal(engine);
    double x3 = normal(engine);
    double factor = scaling/sqrt(x1*x1+x2*x2+x3*x3);

    return {x1*factor, x2*factor, x3*factor};
}

double RandomFloat() {
    // Return a float number between 0 and 1.
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::uniform_real_distribution<double> double_distribution(0, 1);

    return double_distribution(engine);
}