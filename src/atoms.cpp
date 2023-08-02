//
// Created by Gerrit Freiwald on 19.05.2023.
//

#include "atoms.h"
#include "types.h"

Atoms::Atoms(const size_t nb_atoms, bool randomize) :
      names{}, positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms} {
    if (randomize) {
        positions.setRandom();
        velocities.setRandom();
        forces.setRandom();
    } else {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
    }
    masses.setOnes();
}

Atoms::Atoms(const Positions_t &p) :
      names{}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
}

Atoms::Atoms(const Names_t &names, const Positions_t &p) :
      names{names}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
}

Atoms::Atoms(const Positions_t &p, const Velocities_t &v) :
      names{}, positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()} {
    assert(p.cols() == v.cols());
    forces.setZero();
    masses.setOnes();
}

void Atoms::set_masses(double mass) {
    masses *= mass;
}

void Atoms::resize(const int n) {
    positions.conservativeResize(3, n);
    velocities.conservativeResize(3, n);
    forces.conservativeResize(3, n);
    masses.conservativeResize(n);
    names.resize(n);
}

size_t Atoms::nb_atoms() const {
    return positions.cols();
}

double Atoms::e_kin() const {
    return 0.5 * (velocities.colwise().squaredNorm() * masses).sum();
}

double Atoms::local_e_kin(int nb_local) const {
    const Eigen::ArrayXd squared_norm = velocities.leftCols(nb_local).colwise().squaredNorm();
    return 0.5 * (squared_norm * masses).sum();
}

double Atoms::temperature(bool lj_units) const {
    // E_{kin} = 3/2 * k_B * T   with k_B being the Boltzmann constant
    // <=> T = 2 * E_{kin} / (3 * k_B)
    double k_B = 1;

    // Use for "real units" in the EAM scenario
    if (!lj_units) {
        k_B = 8.617333262 * 0.00001;  // unit: eV/K
    }

    return 2./3. * e_kin() / (k_B * nb_atoms());
}

double Atoms::local_temperature(int nb_local, bool lj_units) const {
    // E_{kin} = 3/2 * k_B * T   with k_B being the Boltzmann constant
    // <=> T = 2 * E_{kin} / (3 * k_B)
    double k_B = 1;

    // Use for "real units" in the EAM scenario
    if (!lj_units) {
        k_B = 8.617333262 * 0.00001;  // unit: eV/K
    }

    return 2./3. * local_e_kin(nb_local) / (k_B * nb_local);
}
