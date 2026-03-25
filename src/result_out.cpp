#include "result_out.h"

/**
 * @file result_out.cpp
 * @brief Output utilities for Monte Carlo thermodynamic results.
 */

/**
 * @brief Writes simulation results to a tab-separated text file.
 *
 * The output always includes temperature, energy, heat capacity, moment,
 * and susceptibility. If @p field is true, projected moment/susceptibility
 * columns are appended.
 *
 * @param monte_carlo Monte Carlo configuration (temperature range and step).
 * @param energy Energy values by temperature point.
 * @param Cv Heat capacity values by temperature point.
 * @param moment Magnetic moment values by temperature point.
 * @param chi Susceptibility values by temperature point.
 * @param moment_projection Projected moment values (used when @p field is true).
 * @param chi_projection Projected susceptibility values (used when @p field is true).
 * @param field Whether to include projected field-related columns.
 * @param output_file Output file path.
 * @return int Returns 0 on completion.
 */
int WriteOutput(MonteCarlo & monte_carlo, std::vector<double> energy, std::vector<double> Cv, \
std::vector<double> moment, std::vector<double> chi, \
std::vector<double> moment_projection, std::vector<double> chi_projection, \
bool field, std::string output_file) {
    double T = monte_carlo.start_temperature;
    auto out = fmt::output_file(output_file);
    out.print("#T\tEnergy\tCv\tMoment\tchi");
    if(field) {
        out.print("\tMoment_B\tchi_B\n");
    } else {
        out.print("\n");
    }
    for(int i=0; i<energy.size(); i++) {
        out.print("{:.2f}\t{:.3f}\t{:.5f}\t{:.5f}\t{:.5f}", T, energy[i], Cv[i], moment[i], chi[i]);
        if(field) {
            out.print("\t{:.5f}\t{:.5f}\n", moment_projection[i], chi_projection[i]);
        } else {
            out.print("\n");
        }
        T += monte_carlo.temperature_step;
    }
    out.close();
    return 0;
}