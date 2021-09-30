#include<functional>
#include<math.h>
#include<iostream>
#include<random>
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
    int neighbor_number = 3;

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
    vector<double> super_exchange_parameter = {};
    double energy = 0;

    // Nearest neighbors.
    vector<vector<Site*>> neighbor = {};
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
    Lattice lattice;
    BaseSite base_site;
    vector<vector<vector<vector<Site>>>> site;

    Site & operator[](vector<int> n);
    double energy();
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

Site & SuperCell::operator[](vector<int> n) {
    return this->site[n[0]][n[1]][n[2]][n[3]];
}

double SuperCell::energy() {
    double energy = 0;
    for(int i=0; i<this->lattice.n_x; i++) {
        for(int j=0; j<this->lattice.n_y; j++) {
            for(int k=0; k<this->lattice.n_z; k++) {
                for(int l=0; l<this->base_site.number; k++) {
                    energy += this->site[i][j][k][l].energy;
                }
            }
        }
    }

    return energy*0.5;
}

vector<int> RandomSite(int n_x, int n_y, int n_z, int base_n) {
    // Return the site index randomly.
    static random_device rd;
    static mt19937 engine(rd());
    static uniform_int_distribution<int> int_distribution_x(0, n_x-1);
    static uniform_int_distribution<int> int_distribution_y(0, n_y-1);
    static uniform_int_distribution<int> int_distribution_z(0, n_z-1);
    static uniform_int_distribution<int> int_distribution_base(0, base_n-1);

    vector<int> index = {0, 0, 0, 0};

    index[0] = int_distribution_x(engine);
    index[1] = int_distribution_y(engine);
    index[2] = int_distribution_z(engine);
    index[3] = int_distribution_base(engine);

    return index;
}

vector<double> RandomSpin() {
    // Return a unit spin vector randomly.
    static random_device rd;
    static mt19937 engine(rd());
    static normal_distribution<double> normal{0, 1};

    double x1 = normal(engine);
    double x2 = normal(engine);
    double x3 = normal(engine);
    double factor = 1/sqrt(x1*x1+x2*x2+x3*x3);

    return {x1*factor, x2*factor, x3*factor};
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

int MonteCarloRelaxing(SuperCell & super_cell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    vector<int> site_chosen;
    for(int i=0; i<monte_carlo.relax_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(super_cell.lattice.n_x, super_cell.lattice.n_y, super_cell.lattice.n_z, super_cell.base_site.number);
            Flip(super_cell.lattice, super_cell.base_site, super_cell[site_chosen], T);
        }
    }
    
    return 0;
}

double MonteCarloStep(SuperCell & super_cell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    vector<int> site_chosen;
    double total_energy = 0;
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(super_cell.lattice.n_x, super_cell.lattice.n_y, super_cell.lattice.n_z, super_cell.base_site.number);
            Flip(super_cell.lattice, super_cell.base_site, super_cell[site_chosen], T);
        }
        total_energy += super_cell.energy();
    }
    
    return total_energy / monte_carlo.count_step;
}

int Flip(Lattice & lattice, BaseSite & base_site, Site & one_site, double T) {
    // Flip one spin.
    // Energy and spin before flip.
    double energy = one_site.energy;
    vector<double> old_spin = one_site.spin;

    // Energy and spin after flip.
    one_site.spin = RandomSpin();
    one_site.energy = lattice.Hamiltonian(base_site, one_site);
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