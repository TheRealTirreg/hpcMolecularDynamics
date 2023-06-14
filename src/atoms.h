//
// Created by Gerrit Freiwald on 08.05.2023.
//

#ifndef MY_MD_CODE_ATOMS_H
#define MY_MD_CODE_ATOMS_H

#include <Eigen/Dense>
#include "types.h"

class Atoms {
  public:
    Names_t names;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    explicit Atoms(size_t nb_atoms, bool randomize = false);
    explicit Atoms(const Positions_t &p);
    Atoms(const Positions_t &p, const Velocities_t &v);
    Atoms(const Names_t &names, const Positions_t &p);

    [[nodiscard]] size_t nb_atoms() const;
    [[nodiscard]] double e_kin(double mass = 1) const;
    [[nodiscard]] double temperature(double mass = 1) const;
};

#endif // MY_MD_CODE_ATOMS_H
