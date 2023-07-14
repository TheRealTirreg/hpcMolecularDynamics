//
// Created by Gerrit Freiwald on 14.06.2023.
//

#ifndef MY_MD_CODE_BERENDSEN_THERMOSTAT_H
#define MY_MD_CODE_BERENDSEN_THERMOSTAT_H

#include "types.h"
#include "atoms.h"

void berendsen_thermostat(Atoms &atoms, double temperature,
                          double timestep, double relaxation_time,
                          double mass = 1, bool lj_units = true);


#endif // MY_MD_CODE_BERENDSEN_THERMOSTAT_H
