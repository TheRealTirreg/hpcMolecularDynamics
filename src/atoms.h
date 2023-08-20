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
    Masses_t masses;
    Stress_t stress;

    explicit Atoms(size_t nb_atoms, bool randomize = false);
    explicit Atoms(const Positions_t &p);
    Atoms(const Positions_t &p, const Velocities_t &v);
    Atoms(const Names_t &names, const Positions_t &p);

    void set_masses(double mass);
    void resize(int num);

    [[nodiscard]] size_t nb_atoms() const;
    [[nodiscard]] double e_kin() const;
    [[nodiscard]] double local_e_kin(int nb_local) const;
    [[nodiscard]] double temperature(bool lj_units = true) const;
    [[nodiscard]] double local_temperature(int nb_local, bool lj_units = true) const;
};

#endif // MY_MD_CODE_ATOMS_H
