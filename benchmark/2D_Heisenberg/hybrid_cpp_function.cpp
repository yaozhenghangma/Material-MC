#include <random>
#include <tuple>

#include <pybind11/pybind11.h>
namespace py = pybind11;

using namespace std;

int n = 5;
double kb = 1;
double J = -1.0;

double energy(vector<vector<vector<double> > > & lattice, int x, int y){
    double energy = 0;
    if(x==0) {
        energy += J*(
            lattice[x][y][0] * lattice[x+1][y][0] +
            lattice[x][y][1] * lattice[x+1][y][1] +
            lattice[x][y][2] * lattice[x+1][y][2]);
        energy += J*(
            lattice[x][y][0] * lattice[n-1][y][0] +
            lattice[x][y][1] * lattice[n-1][y][1] +
            lattice[x][y][2] * lattice[n-1][y][2]);
    } else if(x==n-1) {
        energy += J*(
            lattice[x][y][0] * lattice[0][y][0] +
            lattice[x][y][1] * lattice[0][y][1] +
            lattice[x][y][2] * lattice[0][y][2]);
        energy += J*(
            lattice[x][y][0] * lattice[x-1][y][0] +
            lattice[x][y][1] * lattice[x-1][y][1] +
            lattice[x][y][2] * lattice[x-1][y][2]);
    } else {
        energy += J*(
            lattice[x][y][0] * lattice[x+1][y][0] +
            lattice[x][y][1] * lattice[x+1][y][1] +
            lattice[x][y][2] * lattice[x+1][y][2]);
        energy += J*(
            lattice[x][y][0] * lattice[x-1][y][0] +
            lattice[x][y][1] * lattice[x-1][y][1] +
            lattice[x][y][2] * lattice[x-1][y][2]);
    }

    if(y==0) {
        energy += J*(
            lattice[x][y][0] * lattice[x][y+1][0] +
            lattice[x][y][1] * lattice[x][y+1][1] +
            lattice[x][y][2] * lattice[x][y+1][2]);
        energy += J*(
            lattice[x][y][0] * lattice[x][n-1][0] +
            lattice[x][y][1] * lattice[x][n-1][1] +
            lattice[x][y][2] * lattice[x][n-1][2]);
    } else if(y==n-1) {
        energy += J*(
            lattice[x][y][0] * lattice[x][0][0] +
            lattice[x][y][1] * lattice[x][0][1] +
            lattice[x][y][2] * lattice[x][0][2]);
        energy += J*(
            lattice[x][y][0] * lattice[x][y-1][0] +
            lattice[x][y][1] * lattice[x][y-1][1] +
            lattice[x][y][2] * lattice[x][y-1][2]);
    } else {
        energy += J*(
            lattice[x][y][0] * lattice[x][y+1][0] +
            lattice[x][y][1] * lattice[x][y+1][1] +
            lattice[x][y][2] * lattice[x][y+1][2]);
        energy += J*(
            lattice[x][y][0] * lattice[x][y-1][0] +
            lattice[x][y][1] * lattice[x][y-1][1] +
            lattice[x][y][2] * lattice[x][y-1][2]);
    }
    return energy;
}

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

int flip(vector<vector<vector<double> > > & lattice, double T, double & total_energy){
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

double MonteCarlo(double T) {
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
            flip(lattice, T, total_energy);
        }
    }
    for(int i=0; i<n_account; i++) {
        for(int j=0; j<N; j++) {
            flip(lattice, T, total_energy);
        }
    average_energy += total_energy;
    }
    average_energy /= n_account;

    return average_energy;
}

PYBIND11_MODULE(hybrid_cpp_function, m) {
    m.def("MonteCarlo", &MonteCarlo);
}
