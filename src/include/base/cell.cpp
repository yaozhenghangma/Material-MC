//   Copyright (C) 2022  Yaozhenghang Ma
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef CELL
#define CELL

#include <string>
#include <vector>
#include <cmath>

class Atom {
    public:
    std::string symbol = "";
    std::vector<double> coordinate = {0, 0, 0};
};

class Atoms {
    public:
    std::vector<Atom> atoms;
    int number = 0;
};

class Lattice {
    public:
    std::vector<std::vector<double>> lattice = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

class Cell {
    public:
    Lattice lattice;
    Atoms atoms;
};

class Index {
    public:
    int ix = 0;
    int iy = 0;
    int iz = 0;
    int atom_index = 0;

    Index(std::vector<int> index);
};

Index::Index(std::vector<int> index) {
    this->ix = index[0];
    this->iy = index[1];
    this->iz = index[2];
    this->atom_index = index[3];
}

class SuperCell {
    private:
    std::vector<std::vector<std::vector<std::vector<Atom>>>> atoms;

    protected:
    int nx=0, ny=0, nz=0;
    int atoms_per_cell = 0;

    public:
    Lattice lattice;

    SuperCell(int nx, int ny, int nz, Cell cell);

    std::vector<int> shape();
    std::vector<double> coordinate(Index index);
    std::vector<double> displacement(Index index1, Index index2);
    double distance(Index index1, Index index2);
};

SuperCell::SuperCell(int nx, int ny, int nz, Cell cell) {
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    this->atoms_per_cell = cell.atoms.number;
    this->lattice = cell.lattice;
    for(int i=0; i<nx; i++) {
        std::vector<std::vector<std::vector<Atom>>> list_atom2;
        for(int j=0; j<ny; j++) {
            std::vector<std::vector<Atom>> list_atom1;
            for(int k=0; k<nz; k++) {
                list_atom1.push_back(cell.atoms.atoms);
            }
            list_atom2.push_back(list_atom1);
        }
        this->atoms.push_back(list_atom2);
    }
}

std::vector<int> SuperCell::shape() {
    return {this->nx, this->ny, this->nz, this->atoms_per_cell};
}

std::vector<double> SuperCell::coordinate(Index index) {
    return this->atoms[index.ix][index.iy][index.iz][index.atom_index].coordinate;
}

std::vector<double> SuperCell::displacement(Index index1, Index index2) {
    std::vector<double> coordinate1 = this->coordinate(index1);
    std::vector<double> coordinate2 = this->coordinate(index2);
    return {coordinate2[0] - coordinate1[0],
            coordinate2[1] - coordinate1[1],
            coordinate2[2] - coordinate1[2]};
}

double SuperCell::distance(Index index1, Index index2) {
    std::vector<double> displace = this->displacement(index1, index2);
    return sqrt(displace[0]*displace[0] + displace[1]*displace[1] + displace[2]*displace[2]);
}

#endif