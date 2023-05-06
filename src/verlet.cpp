//
// Created by Gerrit Freiwald on 4/26/23.
//

#include "verlet.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep) {
    vx = vx + 0.5 *
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep) {
    // ... implement Verlet step2 here ...
    //vx = vx + 0.5 * timestep * 1/m * (fx + fxplustimestep);
    //vy = vy + 0.5 * timestep * 1/m * (fy + fyplustimestep);
    //vz = vz + 0.5 * timestep * 1/m * (fz + fzplustimestep);
}