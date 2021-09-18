#include<iostream>
#include<vector>

using namespace std;

// Information about the lattice.
class Lattice {
public:
    // Information from POSCAR
    vector<double> a = {0, 0, 0};
    vector<double> b = {0, 0, 0};
    vector<double> c = {0, 0, 0};
    double scaling = 0;

    // Input information
    int n_x = 0;
    int n_y = 0;
    int n_z = 0;
    double anistropic_parameter = 0;
};

// Information to control Monte Carlo circling.
class MonteCarlo {
public:
    // Temperature circling.
    float start_temperature = 0;
    float end_temperature = 0;
    float temperature_step = 0;
    int temperature_step_number = 0;

    // Monte Carlo steps in relaxing process and counting process.
    int relax_step = 0;
    int count_step = 0;

    // Length of Monte Carlo steps.
    int flip_number = 0;
};

// Data of each site.
class Site {
public:
    vector<int> coordinate = {0, 0, 0};
    vector<double> spin = {0, 0, 0};

    // Information connected with chemical elements.
    char element[2];
    float spin_scaling = 0;
    double magnetic_factor = 0;
    double super_exchange_parameter[5];

    // Nearest neighbors.
    vector<Site*> first_nn;
    vector<Site*> second_nn;
    vector<Site*> third_nn;
    vector<Site*> fourth_nn;
    vector<Site*> fifth_nn;
};

int main() {
    return 0;
}