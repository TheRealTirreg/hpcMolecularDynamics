//
// Created by Gerrit Freiwald on 19.05.2023.
//

#include "atoms.h"
#include "types.h"

Atoms::Atoms(const size_t nb_atoms, bool randomize) :
      names{}, positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms} {
    if (randomize) {
        positions.setRandom();
        velocities.setRandom();
        forces.setRandom();
    } else {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }
}

Atoms::Atoms(const Positions_t &p) :
      names{}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} {
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
}

Atoms::Atoms(const Names_t &names, const Positions_t &p) :
      names{names}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} {
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
}

Atoms::Atoms(const Positions_t &p, const Velocities_t &v) :
      names{}, positions{p}, velocities{v}, forces{3, p.cols()} {
    assert(p.cols() == v.cols());
    forces.setZero();
    masses.setOnes();
}

void Atoms::set_masses(double mass) {
    masses = mass;
}


size_t Atoms::nb_atoms() const {
    return positions.cols();
}

double Atoms::e_kin(double mass) const {
    return 0.5 * (velocities.colwise().squaredNorm() * mass).sum();
}

double Atoms::temperature(double mass, bool lj_units) const {
    // E_{kin} = 3/2 * k_B * T   with k_B being the Boltzmann constant
    // <=> T = 2 * E_{kin} / (3 * k_B)
    double k_B = 1;

    // Use for "real units" in the EAM scenario
    if (!lj_units) {
        k_B = 8.617333262 * 0.00001;  // unit: eV/K
    }

    return 2./3. * e_kin(mass) / (k_B * nb_atoms());
}