#include "result_out.h"

int WriteOutput(MonteCarlo & monte_carlo, std::vector<double> energy, std::vector<double> Cv, \
std::vector<double> moment, std::vector<double> Ki, \
std::vector<double> moment_x, std::vector<double> moment_y, std::vector<double> moment_z,\
std::vector<double> Ki_x, std::vector<double> Ki_y, std::vector<double> Ki_z, \
std::string output_file) {
    // Output Monte Carlo results.
    double T = monte_carlo.start_temperature;
    auto out = fmt::output_file(output_file);
    out.print("#T\tEnergy\tCv\tMoment\tKi\tmoment_x\tKi_x\tmoment_y\tKi_y\tmoment_z\tKi_z\n");
    for(int i=0; i<energy.size(); i++) {
        out.print("{:.2f}\t{:.3f}\t{:.5f}\t{:.5f}\t{:.5f}\t", T, energy[i], Cv[i], moment[i], Ki[i]);
        out.print("{:5f}\t{:5f}\t{:5f}\t{:5f}\t{:5f}\t{:5f}\n", moment_x[i], Ki_x[i], moment_y[i], Ki_y[i], moment_z[i], Ki_z[i]);
        T += monte_carlo.temperature_step;
    }
    out.close();
    return 0;
}