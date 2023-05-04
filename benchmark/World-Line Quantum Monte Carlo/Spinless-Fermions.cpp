#include <random>
#include <iostream>

using namespace std;

int main(){
    double t = 1.0;
    double V = 0.0;
    double delta_tau = 0.1;
    vector<double> beta_list={2.0, 1.3, 1.0, 0.8, 0.7, 0.6, 0.5, 0.4};
    int number_fermions = 20;
    int number_lattice = 40;
    int relax_steps = 200;
    int Monte_Carlo_steps = 1000;
    int sweeps = 5;

    static random_device rd;
    static mt19937 engine(rd());
    static uniform_real_distribution<double> distribution(0, 1);

    //Initialization for a time slice
    vector<int> site_slice(number_lattice, 0);
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

    for(int i_beta=0; i_beta<8; i_beta++) {
        vector<vector<int>> site;
        int L = int(beta_list[i_beta]/delta_tau);
        for(int i=0; i<=2*L; i++) {
            site.push_back(site_slice);
        }
        //TODO:Relaxation
        //TODO:Monte Carlo
        //TODO:Output energy
    }

    return 0;
}