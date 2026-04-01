#ifndef MC_STRUCTURE
#define MC_STRUCTURE

#include <functional>
#include <math.h>

#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include "constants.h"

enum class HamiltonianType {
    Heisenberg, Heisenberg_xyz_anisotropy, Heisenberg_with_field, Heisenberg_xyz_anisotropy_with_field,
    Heisenberg_x_anisotropy, Heisenberg_y_anisotropy, Heisenberg_z_anisotropy,
    Heisenberg_xy_anisotropy, Heisenberg_yz_anisotropy, Heisenberg_zx_anisotropy,
    Heisenberg_x_anisotropy_with_field, Heisenberg_y_anisotropy_with_field, Heisenberg_z_anisotropy_with_field,
    Heisenberg_xy_anisotropy_with_field, Heisenberg_yz_anisotropy_with_field, Heisenberg_zx_anisotropy_with_field,
    Heisenberg_custom
};

enum class ModelType {
    Heisenberg, Ising, Kitaev_Heisenberg
};

// Information about base in the cell
class BaseSite {
public:
    // Information from POSCAR
    int number;
    std::vector<std::vector<double>> coordinate;
    std::vector<double> spin_scaling; // Default value: 1.0.
    std::vector<std::vector<double>> spin_initialization;
    std::vector<std::vector<double>> anisotropic_ratio;
    std::vector<std::string> elements;

    // Input information
    std::vector<int> neighbor_number;
    std::vector<std::vector<std::string>> neighbor_elements;
    std::vector<std::vector<double>> neighbor_distance_square;
    std::vector<std::vector<double>> super_exchange_parameter;

    // Anisotropy factor
    std::vector<double> anisotropy = {0, 0, 0};

    // External field
    std::vector<double> B = {0, 0, 0};

    // Global Kitaev-Heisenberg couplings from [Hamiltonian].
    // These remain zero for non-KH models.
    double kh_j = 0.0;
    double kh_k = 0.0;
    double kh_g = 0.0;
    double kh_gp = 0.0;

    // KH bond-type to direction mapping from [Hamiltonian.BondTypeDirection].
    // Index: 0->type1, 1->type2, 2->type3. Values are x/y/z.
    std::vector<char> kh_bond_type_direction = {'x', 'y', 'z'};

    template <class Archive>
    void serialize(Archive & ar) {
        ar(number,
           coordinate,
           spin_scaling,
           spin_initialization,
           anisotropic_ratio,
           elements,
           neighbor_number,
           neighbor_elements,
           neighbor_distance_square,
           super_exchange_parameter,
           anisotropy,
           B,
           kh_j,
           kh_k,
           kh_g,
           kh_gp,
           kh_bond_type_direction);
    }
};

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

    // Per-edge KH bond direction labels (x/y/z), aligned with neighbor[m][n].
    // Empty/zero char means unlabeled (non-KH path).
    std::vector<std::vector<char>> neighbor_direction = {};

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

    HamiltonianType hamiltonian_type = HamiltonianType::Heisenberg;
    ModelType model_type = ModelType::Heisenberg;

    // Maximum relative error in distance computation
    double tolerance_percentage;

    // Output information
    double magnify_factor = 2.0;
    bool ground_state = false;
    bool field = false;
    std::vector<double> field_direction = {0, 0, 0};

    int normalize_direction(std::vector<double>);

    template <class Archive>
    void serialize(Archive & ar) {
        ar(a,
           b,
           c,
           scaling,
           n_x,
           n_y,
           n_z,
           total_energy,
           magnify_factor,
           tolerance_percentage,
           hamiltonian_type,
           model_type,
           ground_state,
           field,
           field_direction);
    }
};

enum class Methods {
    classical, parallel_tempering
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

    // Monte Carlo Methods
    Methods methods = Methods::classical;
    int replica_exchange_step_number = 1;

    template <class Archive>
    void serialize(Archive & ar) {
        ar(start_temperature,
           end_temperature,
           temperature_step,
           temperature_step_number,
           relax_step,
           count_step,
           flip_number,
           methods,
           replica_exchange_step_number);
    }
};

class Initialization {
public:
    std::vector<std::vector<double>> direction;
    std::vector<std::string> elements;
    std::vector<double> angleA = {0, 0, 0};
    std::vector<double> angleB = {0, 0, 0};
    std::vector<double> angleC = {0, 0, 0};

    int normalized();

    template <class Archive>
    void serialize(Archive & ar) {
        ar(direction,
           elements,
           angleA,
           angleB,
           angleC);
    }
};

class Supercell {
public:
    Lattice lattice;
    BaseSite base_site;
    Initialization initialization;
    std::vector<std::vector<std::vector<std::vector<Site>>>> site;

    // Hamiltonian function to calculate energy for one site
    std::function<double(BaseSite &, Site &)> Hamiltonian;
    std::function<double(BaseSite &, Site &)> HamiltonianBase;

    // Update function of each Monte Carlo step.
    std::function<int(Supercell &, Site &, double)> Update;

    Site & operator[](std::vector<int> n);
    double energy();
    double momentum();
    std::vector<double> momentum_component();
};

#endif