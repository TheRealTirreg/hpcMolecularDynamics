//
// Created by Gerrit Freiwald on 08.05.2023.
//

#ifndef MY_MD_CODE_ATOMS_H
#define MY_MD_CODE_ATOMS_H

#include <Eigen/Dense>
#include "types.h"

class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(const size_t nb_atoms, bool randomize = false) :
          positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms} {
        if (randomize) {
            positions.setRandom();
            velocities.setRandom();
            forces.setRandom();
        } else {
            positions.setZero();
            velocities.setZero();
            forces.setZero();
        }
    }

    explicit Atoms(const Positions_t &p) :
          positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
          positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }
};

#endif // MY_MD_CODE_ATOMS_H
