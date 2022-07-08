#include <random>
#include <tuple>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
namespace py = pybind11;

using namespace std;

int n = 5;
double kb = 1;
double J = -1.0;

double RandomFloat() {
    // Return a float number between 0 and 1.
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::uniform_real_distribution<double> double_distribution(0, 1);

    return double_distribution(engine);
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

tuple<double,double> RandomSite() {
    // Return the site index randomly.
    static std::random_device rd;
    static std::mt19937 engine(rd());
    static std::uniform_int_distribution<int> int_distribution_x(0, n-1);
    static std::uniform_int_distribution<int> int_distribution_y(0, n-1);

    int x = int_distribution_x(engine);
    int y = int_distribution_y(engine);

    return make_tuple(x, y);
}

int flip(vector<vector<vector<double> > > & lattice, double T, double & total_energy, const std::function<double(vector<vector<vector<double>>>, int, int)> &energy){
    int x, y;
    tie(x, y) = RandomSite();
    vector<double> spin = RandomSpin(1.0);
    double old_energy = energy(lattice, x, y);
    vector<double> old_spin = lattice[x][y];
    lattice[x][y] = spin;
    double new_energy = energy(lattice, x, y);
    double delta_energy = (new_energy - old_energy) * 2;
    if(exp(-delta_energy/(kb*T)) > RandomFloat()) {
        total_energy += delta_energy;
    } else {
        lattice[x][y] = old_spin;
    }
    return 0;
}

double MonteCarlo(double T, const std::function<double(vector<vector<vector<double>>>, int, int)> &energy) {
    vector<vector<vector<double>>> lattice;
    vector<vector<double>> row;
    vector<double> spin = {1, 0, 0};
    for(int i=0; i<n; i++) {
        row.push_back(spin);
    }
    for(int i=0; i<n; i++) {
        lattice.push_back(row);
    }

    int n_account = 300;
    int n_relax = 100;
    int N = 300;
    double total_energy=0;
    double average_energy=0;
    for(int i=0; i<n_relax; i++) {
        for(int j=0; j<N; j++) {
            flip(lattice, T, total_energy, energy);
        }
    }
    for(int i=0; i<n_account; i++) {
        for(int j=0; j<N; j++) {
            flip(lattice, T, total_energy, energy);
        }
    average_energy += total_energy;
    }
    average_energy /= n_account;

    return average_energy;
}

PYBIND11_MODULE(hybrid_python_function, m) {
    m.def("MonteCarlo", &MonteCarlo);
}