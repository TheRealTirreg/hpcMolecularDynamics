//
// Created by Gerrit Freiwald on 14.06.2023.
//

#include "berendsen_thermostat.h"

void berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time, bool lj_units) {
    const double temperature_relation = temperature / atoms.temperature(lj_units);
    const double lambda = std::sqrt(
        1 + (temperature_relation - 1) * timestep / relaxation_time);

    atoms.velocities *= lambda;
}

void berendsen_thermostat_mp(Atoms &atoms, Domain &domain, double temperature, double timestep,
                          double relaxation_time, bool lj_units) {
    const double temperature_relation = temperature / atoms.local_temperature(domain.nb_local(), lj_units);
    const double lambda = std::sqrt(
        1 + (temperature_relation - 1) * timestep / relaxation_time);

    atoms.velocities.leftCols(domain.nb_local()) *= lambda;
}
