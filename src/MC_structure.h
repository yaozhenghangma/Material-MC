#ifndef MC_STRUCTURE
#define MC_STRUCTURE

#include <functional>
#include <math.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

const double gs = 2.0;

enum class HamiltonionType {
    Heisenberg, Heisenberg_xyz_anisotropy, Heisenberg_with_field, Heisenberg_xyz_anisotropy_with_field,
    Heisenberg_x_anisotropy, Heisenberg_y_anisotropy, Heisenberg_z_anisotropy,
    Heisenberg_xy_anisotropy, Heisenberg_yz_anisotropy, Heisenberg_zx_anisotropy,
    Heisenberg_x_anisotropy_with_field, Heisenberg_y_anisotropy_with_field, Heisenberg_z_anisotropy_with_field,
    Heisenberg_xy_anisotropy_with_field, Heisenberg_yz_anisotropy_with_field, Heisenberg_zx_anisotropy_with_field
};

// Information about base in the cell
class BaseSite {
public:
    // Information from POSCAR
    int number;
    std::vector<std::vector<double>> coordinate;
    std::vector<double> spin_scaling; // Default value: 1.0.
    std::vector<std::vector<double>> anisotropic_ratio;
    std::vector<std::string> elements;

    // Input information    
    std::vector<int> neighbor_number;
    std::vector<std::vector<std::string>> neighbor_elements;
    std::vector<std::vector<double>> neighbor_distance_square;
    std::vector<std::vector<double>> super_exchange_parameter;

    //TODO: record and output coordination number

    // Anisotropy factor
    std::vector<double> anisotropy = {0, 0, 0};

    // External field
    std::vector<double> B = {0, 0, 0};
};

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, BaseSite & base_site, const unsigned int version)
{
    ar & base_site.number;
    ar & base_site.coordinate;
    ar & base_site.spin_scaling;
    ar & base_site.anisotropic_ratio;
    ar & base_site.elements;
    ar & base_site.neighbor_number;
    ar & base_site.neighbor_elements;
    ar & base_site.neighbor_distance_square;
    ar & base_site.super_exchange_parameter;
    ar & base_site.anisotropy;
    ar & base_site.B;
}

}
}

// Data of each site.
class Site {
public:
    std::vector<double> spin = {0, 0, 0};
    double * spin_scaling;
    int * neighbor_number;
    std::vector<double> * anisotropic_ratio;
    std::vector<double> * super_exchange_parameter;

    // Neighbors' link.
    std::vector<std::vector<Site*>> neighbor = {};

    int reverse_spin();
};

// Information about the lattice.
class Lattice {
public:
    // Information from POSCAR
    std::vector<double> a = {0, 0, 0};
    std::vector<double> b = {0, 0, 0};
    std::vector<double> c = {0, 0, 0};
    double scaling;

    // Input information
    int n_x;
    int n_y;
    int n_z;

    double total_energy;

    HamiltonionType hamiltonion_type = HamiltonionType::Heisenberg;

    // Maximum relative error in distance computation
    double tolerance_percentage;
    
    // Output information
    double magnify_factor = 2.0;
};

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, Lattice & lattice, const unsigned int version)
{
    ar & lattice.a;
    ar & lattice.b;
    ar & lattice.c;
    ar & lattice.scaling;
    ar & lattice.n_x;
    ar & lattice.n_y;
    ar & lattice.n_z;
    ar & lattice.total_energy;
    ar & lattice.magnify_factor;
    ar & lattice.tolerance_percentage;
    ar & lattice.hamiltonion_type;
}

}
}

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

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, MonteCarlo & monte_carlo, const unsigned int version)
{
    ar & monte_carlo.start_temperature;
    ar & monte_carlo.end_temperature;
    ar & monte_carlo.temperature_step;
    ar & monte_carlo.temperature_step_number;
    ar & monte_carlo.relax_step;
    ar & monte_carlo.count_step;
    ar & monte_carlo.flip_number;
}

}
}

class Initialization {
public:
    std::vector<double> direction = {0, 0, 1};
    bool anti_ferromagnetic = false;
    std::vector<std::vector<int>> anti_ferromagnetic_J;

    int normalized();
};

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, Initialization & initialization, const unsigned int version)
{
    ar & initialization.direction;
    ar & initialization.anti_ferromagnetic;
    ar & initialization.anti_ferromagnetic_J;
}

}
}

class Supercell {
public:
    Lattice lattice;
    BaseSite base_site;
    Initialization initialization;
    std::vector<std::vector<std::vector<std::vector<Site>>>> site;

    // Hamiltonian function to calculate energy for one site
    std::function<double(BaseSite &, Site &)> Hamiltonian;

    Site & operator[](std::vector<int> n);
    double energy();
    double momentum();
    std::vector<double> momentum_component();
};

#endif