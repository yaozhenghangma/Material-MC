#include "MC_structure.h"

Site & Supercell::operator[](std::vector<int> n) {
    return this->site[n[0]][n[1]][n[2]][n[3]];
}

double Supercell::energy() {
    double e = 0;
    for(int i=0; i<this->lattice.n_x; i++) {
        for(int j=0; j<this->lattice.n_y; j++) {
            for(int k=0; k<this->lattice.n_z; k++) {
                for(int l=0; l<this->base_site.number; l++) {
                    e += this->Hamiltonian(this->base_site, this->site[i][j][k][l]);
                }
            }
        }
    }

    return e*0.5;
}

double Supercell::momentum() {
    std::vector<double> m = {0, 0, 0};
    for(int i=0; i<this->lattice.n_x; i++) {
        for(int j=0; j<this->lattice.n_y; j++) {
            for(int k=0; k<this->lattice.n_z; k++) {
                for(int l=0; l<this->base_site.number; l++) {
                    m[0] += this->site[i][j][k][l].spin[0];
                    m[1] += this->site[i][j][k][l].spin[1];
                    m[2] += this->site[i][j][k][l].spin[2];
                }
            }
        }
    }

    return sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2])*gs;
}

std::vector<double> Supercell::momentum_component() {
    double mx = 0;
    double my = 0;
    double mz = 0;
    for(int i=0; i<this->lattice.n_x; i++) {
        for(int j=0; j<this->lattice.n_y; j++) {
            for(int k=0; k<this->lattice.n_z; k++) {
                for(int l=0; l<this->base_site.number; l++) {
                    mx += this->site[i][j][k][l].spin[0];
                    my += this->site[i][j][k][l].spin[1];
                    mz += this->site[i][j][k][l].spin[2];
                }
            }
        }
    }
    return {mx*gs, my*gs, mz*gs};
}

int Initialization::normalized() {
    for(int i=0; i<this->direction.size(); i++) { 
        double norm = sqrt(1.0 / (this->direction[i][0]*this->direction[i][0] + \
        this->direction[i][1]*this->direction[i][1] + \
        this->direction[i][2]*this->direction[i][2]));
        this->direction[i][0] *= norm;
        this->direction[i][1] *= norm;
        this->direction[i][2] *= norm;
    }
    return 0;
}

int Site::reverse_spin() {
    this->spin[0] = - this->spin[0];
    this->spin[1] = - this->spin[1];
    this->spin[2] = - this->spin[2];
    return 0;
}