#include "result_out.h"

int WriteOutput(MonteCarlo & monte_carlo, std::vector<double> energy, std::vector<double> Cv, \
std::vector<double> moment, std::vector<double> chi, \
std::vector<double> moment_projection, std::vector<double> chi_projection, \
bool field, std::string output_file) {
    // Output Monte Carlo results.
    double T = monte_carlo.start_temperature;
    auto out = fmt::output_file(output_file);
    out.print("#T\tEnergy\tCv\tMoment\tchi");
    if(field) {
        out.print("\tMoment_B\tchi_B\n");
    } else {
        out.print("\n");
    }
    for(int i=0; i<energy.size(); i++) {
        out.print("{:.2f}\t{:.3f}\t{:.5f}\t{:.5f}\t{:.5f}\t", T, energy[i], Cv[i], moment[i], chi[i]);
        if(field) {
            out.print("{:.5f}\t{:.5f}\n", moment_projection[i], chi_projection[i]);
        } else {
            out.print("\n");
        }
        T += monte_carlo.temperature_step;
    }
    out.close();
    return 0;
}