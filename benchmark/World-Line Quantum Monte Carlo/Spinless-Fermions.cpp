#include <cmath>
#include <random>
#include <iostream>

using namespace std;

float t = 1.0;
float V = 0.0;
float delta_tau = 0.1;

int main(){
    vector<double> beta_list={2.0, 1.3, 1.0, 0.8, 0.7, 0.6, 0.5, 0.4};
    int number_fermions = 20;
    int number_lattice = 40;
    int relax_steps = 200;
    int Monte_Carlo_steps = 1000;
    int sweeps = 5;

    //Initialization for a time slice
    vector<int> site_slice(number_lattice, 0);
    initialize(site_slice, number_fermions, number_lattice);

    for(int i_beta=0; i_beta<8; i_beta++) {
        double energy = 0;
        vector<vector<int>> site;
        int L = int(beta_list[i_beta]/delta_tau);
        for(int i=0; i<2*L; i++) {
            site.push_back(site_slice);
        }
        //Relaxation
        for(int j=0; j<relax_steps; j++) {
            for(int k=0; k<sweeps; k++) {
                sweep(site, number_lattice, L);
            }
        }
        //Monte Carlo
        for(int j=0; j<Monte_Carlo_steps; j++) {
            for(int k=0; k<sweeps; k++) {
                sweep(site, number_lattice, L);
            }
            energy += count_energy(site);
        }
        //TODO:Output energy
    }

    return 0;
}

int initialize(vector<int> & site_slice, int number_fermions, int number_lattice){
    //Initialize site for a specific time slice randomly
    static random_device rd;
    static mt19937 engine(rd());
    static uniform_real_distribution<double> distribution(0, 1);
    for(int i=0, fermions_left=number_fermions; i<number_lattice; i++) {
        if(number_lattice-i<=fermions_left){
            site_slice[i] = 1;
            fermions_left--;
        } else if(fermions_left==0){
            continue;
        } else if(distribution(engine)>0.5){
            site_slice[i] = 1;
            fermions_left--;
        } else{
            continue;
        }
    }
    return 0;
}

int sweep(vector<vector<int>> & site, int number_lattice, int L){
    int j=0; //j for time
    int l=0; //l for lattice
    int i=0; //i=2*l or i=2*l+1
    int s=0, u=0, v=0; //parameters

    double R = 0;
    double P = 0;

    static random_device rd;
    static mt19937 engine(rd());
    static uniform_real_distribution<double> distribution(0, 1);

    static double tanh_value = tanh(delta_tau*t)*tanh(delta_tau*t);
    static double cosh_value = cosh(delta_tau*t)*cosh(delta_tau*t);
    static double cosh2_value = cosh(delta_tau*t*2.0);
    static double cosh4_value = cosh(delta_tau*t*4.0);
    static double sinh_value = sinh(delta_tau*t)*sinh(delta_tau*t);
    static double sech_value = 1.0/(cosh(delta_tau*t)*cosh(delta_tau*t));
    static double sech2_value = 1.0/cosh(delta_tau*t*2.0);
    static double exp_value = exp(delta_tau*V/2);
    static double inverse_exp_value = exp(-delta_tau*V/2.0);

    static double Psp2up1vp1 = tanh_value * exp_value;
    static double Psp2up1v0  = sech2_value * sinh_value;
    static double Psp2up1vn1 = 8.0*inverse_exp_value*sinh_value/(-1.0+8.0*cosh2_value+cosh4_value);
    static double Psp2u0vp1  = exp_value*cosh_value/(1.0+cosh_value);
    static double Psp2u0v0   = 0.5;
    static double Psp2u0vn1  = 2.0*inverse_exp_value/(3.0+cosh2_value);
    static double Psp2un1vp1 = 8.0*exp_value*cosh_value*cosh_value/(-1.0+8.0*cosh2_value+cosh4_value);
    static double Psp2un1v0  = cosh_value*sech2_value;
    static double Psp2un1vn1 = inverse_exp_value*sech_value;
    static double Psn2up1vp1 = inverse_exp_value*sech_value;
    static double Psn2up1v0  = cosh_value*sech2_value;
    static double Psn2up1vn1 = 8.0*exp_value*cosh_value*cosh_value/(-1.0+8.0*cosh2_value+cosh4_value);
    static double Psn2u0vp1  = 2*inverse_exp_value/(3.0+cosh2_value);
    static double Psn2u0v0   = 0.5;
    static double Psn2u0vn1  = exp_value*cosh_value/(1.0+cosh_value);
    static double Psn2un1vp1 = 8.0*inverse_exp_value*sinh_value/(-1.0+8.0*cosh2_value+cosh4_value);
    static double Psn2un1v0  = sech2_value*sinh_value;
    static double Psn2un1vn1 = exp_value*tanh_value;

    //loop for j=0
    for(l=0; l<number_lattice/2-1; l++){
        i = 2*l+1;
        if(i==number_lattice-1) {
            s = site[i][j] + site[i][j+1] - site[0][j] - site[0][j+1];
            u = 1 - site[0][L-1] - site[0][j+2];
            v = site[i-1][j] - site[1][j];
        } else {
            s = site[i][j] + site[i][j+1] - site[i+1][j] - site[i+1][j+1];
            u = 1 - site[i+1][L-1] - site[i+1][j+2];
            v = site[i-1][j] - site[i+2][j];
        }
        if(s==2 && u==1) {
            if(distribution(engine) < P) {

            }
        }
    }
}

double count_energy(vector<vector<int>> & site) {
    return 0;
}