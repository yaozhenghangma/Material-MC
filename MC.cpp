#include<functional>
#include<math.h>
#include<iostream>
#include<vector>

using namespace std;

const double KB = 1;

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


    // Hamiltonian function to calculate energy for one site
    function<double(BaseSite &, Site &)> Hamiltonian;
};

// Information about base in the cell
class BaseSite {
public:
    // Information from POSCAR
    int number;
    vector<vector<double>> coordinate;
    vector<char[2]> elements;
};

// Data of each site.
class Site {
public:
    vector<double> spin = {0, 0, 0};

    // Information connected with chemical elements.
    char element[2];
    float spin_scaling = 0;
    double magnetic_factor = 0;
    double super_exchange_parameter[5];
    double energy = 0;

    // Nearest neighbors.
    vector<Site*> first_nn;
    vector<Site*> second_nn;
    vector<Site*> third_nn;
    vector<Site*> fourth_nn;
    vector<Site*> fifth_nn;
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

class SuperCell {
public:
//TODO: overwrite the operator []
};

int main() {
    //TODO: Read information from command line.
    //TODO: Read information from POSCAR
    //TODO: Read information from input file.
    //TODO: Enlarge the cell with given n.
    //TODO: Initialize the site information (spin and energy).
    //TODO: Monte Carlo
    //TODO: Output spin state at every temperature.
    //TODO: Output the thermal dynamic result.
    return 0;
}

int RandomInt(int n) {
    //TODO: Return a integer between 0 and n.
}

vector<double> RandomSpin() {
    //TODO: Return a spin vector randomly.
}

double RandomFloat() {
    //TODO: Return a float number between 0 and 1.
}

int ReadOptions() {
    //TODO: Process options from command line.
}

int ReadPOSCAR() {
    //TODO: Read information about base and lattice from POSCAR(default).
}

int ReadSettingFile() {
    //TODO: Read information about enlarging and Monte Carlo from given setting file.
}

int EnlargeCell() {
    //TODO: Enlarge the system with given number.
}

int MonteCarloStep(Lattice & lattice_data, BaseSite & base_data, SuperCell & site_data, MonteCarlo & monte_carlo, double T) {
    //TODO: Do Monte Carlo simulation, with given flipping number.
    int x, y, z = 0;
    int base_n = 0;
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            x = RandomInt(lattice_data.n_x);
            y = RandomInt(lattice_data.n_y);
            z = RandomInt(lattice_data.n_z);
            base_n = RandomInt(base_data.number);
        }
    }
}

int Flip(Lattice & lattice_data, BaseSite & base_data, Site & one_site, double T) {
    // Flip one spin.
    // Energy and spin before flip.
    double energy = one_site.energy;
    vector<double> old_spin = one_site.spin;

    // Energy and spin after flip.
    one_site.spin = RandomSpin();
    one_site.energy = lattice_data.Hamiltonian(base_data, one_site);
    double de = one_site.energy - energy;

    // Judge whether to flip.
    double crition = RandomFloat();
    if (crition > exp(-de/(KB*T))) {
        one_site.spin = old_spin;
        one_site.energy = energy;
    }

    return 0;
}

int WriteSpin() {
    //TODO: Output spin states of all atoms.
}

int WriteOutput() {
    //TODO: Output Monte Carlo results.
}