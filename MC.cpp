#include<fmt/core.h>
#include<functional>
#include<math.h>
#include<iostream>
#include<random>
#include<stdlib.h>
#include<unistd.h>
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
    
    // Layered material
    bool layered = true;

    // Hamiltonian function to calculate energy for one site
    function<double(BaseSite &, Site &)> Hamiltonian;
};

// Information about base in the cell
class BaseSite {
public:
    // Information from POSCAR
    int number;
    vector<vector<double>> coordinate;
    vector<vector<double>> super_exchange_parameter;
    vector<double> spin_scaling;
    vector<double> magnetic_factor;
    vector<char[3]> elements;
};

// Data of each site.
class Site {
public:
    vector<double> spin = {0, 0, 0};

    // Information connected with chemical elements.
    //char element[3] = "  ";
    //float spin_scaling = 0;
    //double magnetic_factor = 0;
    //vector<double> super_exchange_parameter = {};
    double energy = 0;
    //TODO: store the momentum

    // Neighbors' link. For normal crystal, only variation "neighbor_ab" is used.
    vector<vector<Site*>> neighbor_ab = {};
    vector<vector<Site*>> neighbor_c = {};
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

class Supercell {
public:
    Lattice lattice;
    BaseSite base_site;
    vector<vector<vector<vector<Site>>>> site;

    Site & operator[](vector<int> n);
    double energy();
    double momentum();
};

int ReadOptions(int argc, char** argv, string & cell_structure_file, string & input_file, string & output_file, string & spin_structure_file);
int usage(char* program);

int main(int argc, char** argv) {
    // Read information from command line.
    string  cell_structure_file = "POSCAR";
    string input_file = "input.txt"; 
    string output_file = "output.txt";
    string spin_structure_file = "spin.txt";
    ReadOptions(argc, argv, cell_structure_file, input_file, output_file, spin_structure_file);

    //TODO: Read information from POSCAR
    //TODO: Read information from input file.
    //TODO: Enlarge the cell with given n.
    //TODO: Initialize the site information (spin and energy).
    //TODO: Monte Carlo
    //TODO: Output spin state at every temperature.
    //TODO: Output the thermal dynamic result.
    return 0;
}

Site & Supercell::operator[](vector<int> n) {
    return this->site[n[0]][n[1]][n[2]][n[3]];
}

double Supercell::energy() {
    double e = 0;
    for(int i=0; i<this->lattice.n_x; i++) {
        for(int j=0; j<this->lattice.n_y; j++) {
            for(int k=0; k<this->lattice.n_z; k++) {
                for(int l=0; l<this->base_site.number; k++) {
                    e += this->site[i][j][k][l].energy;
                }
            }
        }
    }

    return e*0.5;
}

double Supercell::momentum() {
    vector<double> m = {0, 0, 0};
    for(int i=0; i<this->lattice.n_x; i++) {
        for(int j=0; j<this->lattice.n_y; j++) {
            for(int k=0; k<this->lattice.n_z; k++) {
                for(int l=0; l<this->base_site.number; k++) {
                    m[0] += this->site[i][j][k][l].spin[0];
                    m[1] += this->site[i][j][k][l].spin[1];
                    m[2] += this->site[i][j][k][l].spin[2];
                }
            }
        }
    }

    return sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
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
    // Return a float number between 0 and 1.
    static random_device rd;
    static mt19937 engine(rd());
    static uniform_real_distribution<double> double_distribution(0, 1);

    return double_distribution(engine);
}

int ReadOptions(int argc, char** argv, string & cell_structure_file, string & input_file, string & output_file, string & spin_structure_file) {
    // Process options from command line.
    char option;
    while ((option = getopt(argc, argv, "c:i:o:s:h")) != -1) {
        switch (option) {
            case 'c':
                cell_structure_file = optarg; break;
            case 'i':
                input_file = optarg; break;
            case 'o':
                output_file = optarg; break;
            case 's':
                spin_structure_file = optarg; break;
            case 'h':
                usage(argv[0]);
        }
    }
}

int usage(char* program) {
    //TODO: Usage of the program.
    exit(1);
}

int ReadPOSCAR() {
    //TODO: Read information about base and lattice from POSCAR(default).
}

int ReadSettingFile() {
    //TODO: Read information about enlarging and Monte Carlo from given setting file.
}

int EnlargeCell(Supercell & supercell) {
    // Enlarge the system with given number and initialize the spin.
    vector<vector<vector<Site>>> site1;
    vector<vector<Site>> site2;
    vector<Site> site3;
    Site site4;

    site4.spin[2] = 1.0;

    for(int i=0; i<supercell.lattice.n_x; i++) {
        supercell.site.push_back(site1);
        for(int j=0; j<supercell.lattice.n_y; j++) {
            supercell.site[i].push_back(site2);
            for(int k=0; k<supercell.lattice.n_z; k++) {
                supercell.site[i][j].push_back(site3);
                for(int l=0; l<supercell.base_site.number; l++) {
                    supercell.site[i][j][k].push_back(site4);
                    supercell.site[i][j][k][l].spin[2] *= supercell.base_site.spin_scaling[l];
                }
            }
        }
    }

    return 0;
}

double Distance(vector<double> lattice_constant1, vector<double> lattice_constant2, vector<double> lattice_constant3, \
vector<int> index, vector<double> base_site1, vector<double> base_site2) {
    // Return squared distance of index-base_site1+base_site2.
    vector<double> total_index = {index[0]+base_site2[0]-base_site1[0], index[1]+base_site2[1]-base_site1[1], index[2]+base_site2[2]-base_site1[2]};
    double result = 0;
    for(int i=0; i<3; i++) {
        result += (total_index[0]*lattice_constant1[i] + total_index[1]*lattice_constant1[i] + total_index[2]*lattice_constant1[i]) \
        * (total_index[0]*lattice_constant1[i] + total_index[1]*lattice_constant1[i] + total_index[2]*lattice_constant1[i]);
    }
    return result;
}

int AddDistance(double distance, vector<double> & distance_list) {
    int s = distance_list.size()-1;
    if(distance_list[s] == 0 || distance_list[s] > distance) {
        distance_list[s] = distance;
        for(int i=s-1; i>=0; i--) {
            if(distance_list[i] == 0 || distance_list[i] > distance) {
                distance_list[i+1] = distance_list[i];
                distance_list[i] = distance;
            }
        }
    }
    return 0;
}

int InitializeSupercell(Supercell & supercell) {
    // Initialize the neighbors' link and energy.
    if(supercell.lattice.layered) {

    } else {
        // Find the neighbors for every base sites in a cell.
        // neighbors_index[i][j][k][l]: i-th site's k-th j nearest neighbor's index in l-th direction. 4-th direction is base number.
        vector<vector<vector<vector<int>>>> neighbors_index;

        for(int i=0; i<supercell.base_site.number; i++) {
            vector<vector<vector<int>>> index1;
            neighbors_index.push_back(index1);
            for(int j=0; j<supercell.lattice.neighbor_number; j++) {
                vector<vector<int>> index2;
                neighbors_index[i].push_back(index2);
            }
        }

        for(int i=0; i<supercell.base_site.number; i++) {
            // Neighbors for i-th base site.
            vector<double> distance_list(supercell.lattice.neighbor_number, 0);
            double distance_square = 0;

            // Find the distance values.
            for(int j=-supercell.lattice.neighbor_number; j<supercell.lattice.neighbor_number+1; j++) {
                for(int k=-supercell.lattice.neighbor_number+abs(j); k<supercell.lattice.neighbor_number-abs(j); k++) {
                    for(int l=-supercell.lattice.neighbor_number+abs(j)+abs(k); l<supercell.lattice.neighbor_number-abs(k)-abs(j); l++) {
                        for(int m=0; m<supercell.base_site.number; m++) {
                            distance_square = Distance(supercell.lattice.a, supercell.lattice.b, supercell.lattice.c, {i, j, k}, \
                            supercell.base_site.coordinate[i], supercell.base_site.coordinate[m]);

                            if(distance_square == 0.0) {
                                continue;
                            } else {
                                AddDistance(distance_square, distance_list);
                            }
                        }
                    }
                }
            }

            // Find link.
            for(int j=-supercell.lattice.neighbor_number; j<supercell.lattice.neighbor_number+1; j++) {
                for(int k=-supercell.lattice.neighbor_number+abs(j); k<supercell.lattice.neighbor_number-abs(j); k++) {
                    for(int l=-supercell.lattice.neighbor_number+abs(j)+abs(k); l<supercell.lattice.neighbor_number-abs(k)-abs(j); l++) {
                        for(int m=0; m<supercell.base_site.number; m++) {
                            distance_square = Distance(supercell.lattice.a, supercell.lattice.b, supercell.lattice.c, {i, j, k}, \
                            supercell.base_site.coordinate[i], supercell.base_site.coordinate[m]);

                            if(distance_square == 0.0) {
                                continue;
                            } else {
                                for(int n=0; n<supercell.lattice.neighbor_number; n++) {
                                    if(distance_square == distance_list[n]) {
                                        vector<int> ind = {j, k, l, m};
                                        neighbors_index[i][n].emplace_back(ind);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Set link for every sites.
        for(int i=0; i<supercell.lattice.n_x; i++) {
            for(int j=0; j<supercell.lattice.n_y; j++) {
                for(int k=0; k<supercell.lattice.n_z; k++) {
                    for(int l=0; l<supercell.base_site.number; l++) {
                        vector<Site*> temp = {};
                        for(int m=0; m<supercell.lattice.neighbor_number; m++) {
                            supercell.site[i][j][k][l].neighbor_ab.push_back(temp);
                            for(int n=0; n<neighbors_index[l][m].size(); n++) {
                                supercell.site[i][j][k][l].neighbor_ab[m].push_back( \ 
                                & supercell.site[(i+neighbors_index[l][m][n][0]+supercell.lattice.n_x) % supercell.lattice.n_x] \
                                [(j+neighbors_index[l][m][n][1]+supercell.lattice.n_y) % supercell.lattice.n_y] \
                                [(k+neighbors_index[l][m][n][2]+supercell.lattice.n_z) % supercell.lattice.n_z] \
                                [neighbors_index[l][m][n][3]]);
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int MonteCarloRelaxing(Supercell & supercell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    vector<int> site_chosen;
    for(int i=0; i<monte_carlo.relax_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            Flip(supercell.lattice, supercell.base_site, supercell[site_chosen], T);
        }
    }
    
    return 0;
}

vector<double> MonteCarloStep(Supercell & supercell, MonteCarlo & monte_carlo, double T) {
    // Monte Carlo simulation, with given flipping number and count number, at a specific temperature.
    vector<int> site_chosen;
    double total_energy = 0;
    double total_momentum = 0;
    static double one_over_step = 1 / monte_carlo.count_step;
    for(int i=0; i<monte_carlo.count_step; i++) {
        for(int j=0; j<monte_carlo.flip_number; j++) {
            site_chosen = RandomSite(supercell.lattice.n_x, supercell.lattice.n_y, supercell.lattice.n_z, supercell.base_site.number);
            Flip(supercell.lattice, supercell.base_site, supercell[site_chosen], T);
        }
        total_energy += supercell.energy();
        total_momentum += supercell.momentum();
    }
    
    return {total_energy * one_over_step, total_momentum * one_over_step};
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