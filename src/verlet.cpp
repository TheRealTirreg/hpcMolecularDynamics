//
// Created by Gerrit Freiwald on 4/26/23.
//

#include <Eigen/Dense>
#include "verlet.h"
#include "types.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep) {
    // velocities
    vx = vx + 0.5 * fx * timestep * MASS_RECIPROKE;
    vy = vy + 0.5 * fy * timestep * MASS_RECIPROKE;
    vz = vz + 0.5 * fz * timestep * MASS_RECIPROKE;

    // positions
    x = x + vx * timestep;
    y = y + vy * timestep;
    z = z + vz * timestep;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep) {
    // add force difference over time to velocities
    vx = vx + 0.5 * fx * timestep * MASS_RECIPROKE;
    vy = vy + 0.5 * fy * timestep * MASS_RECIPROKE;
    vz = vz + 0.5 * fz * timestep * MASS_RECIPROKE;
}


void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep) {
    velocities = velocities + 0.5 * forces * timestep * MASS_RECIPROKE;
    positions = positions + velocities * timestep;
}

void verlet_step2(Velocities_t &velocities, const Forces_t &forces,
                  double timestep) {
    velocities = velocities + 0.5 * forces * timestep * MASS_RECIPROKE;
}