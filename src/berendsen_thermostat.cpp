//
// Created by Gerrit Freiwald on 14.06.2023.
//

#include "berendsen_thermostat.h"

void berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time, double mass, bool lj_units) {
    const double temperature_relation = temperature / atoms.temperature(mass, lj_units);
    const double lambda = std::sqrt(
        1 + (temperature_relation - 1) * timestep / relaxation_time);

    atoms.velocities = atoms.velocities * lambda;
}
