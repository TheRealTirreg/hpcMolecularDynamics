//
// Created by Gerrit Freiwald on 4/26/23.
//

#ifndef CPP_MOLECULAR_DYNAMICS_VERLET_H
#define CPP_MOLECULAR_DYNAMICS_VERLET_H

#include <Eigen/Dense>
#include "types.h"

const int MASS_RECIPROKE = 1;

void verlet_step1(double &x, double &y, double &z,
                  double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep);

void verlet_step2(double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep);

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep);

void verlet_step2(Velocities_t &velocities, const Forces_t &forces,
                  double timestep);

#endif //CPP_MOLECULAR_DYNAMICS_VERLET_H
