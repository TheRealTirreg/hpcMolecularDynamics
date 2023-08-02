#include "atoms.h"
#include "types.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "xyz.h"
#include "write_file.h"
#include "berendsen_thermostat.h"
#include <iostream>


/* Heat up system by rescaling velocities
 * Plot energies & temperature with energy jumps \delta Q
 * In the temperature plot, look at melting points
 * */


int main() {
    std::cout << "Starting...\n";
    std::string cluster_num = "28741";
    std::string filename = "clusters/cluster_" + cluster_num + ".xyz";
    auto [names, positions]{read_xyz(filename)};

    Atoms atoms = Atoms(Names_t(names), Positions_t(positions));
    std::cout << "num atoms: " << atoms.nb_atoms() << "\n";

    double mass_gold = 196.97;
    double mass_unit_factor = 103.6;
    double mass = mass_gold * mass_unit_factor;  // unit: g/mol
    atoms.set_masses(mass);
    double timestep = 0.5;  // unit: fs
    double total_time = 80000;  // unit: fs
    double current_time = 0;  // unit: fs

    double e_pot = 0;

    // Neighbors
    double neighbors_cutoff = 8;
    NeighborList neighbors_list{neighbors_cutoff};
    neighbors_list.update(atoms);

    // Thermostat
    double thermostat_relaxation_time = 2000 * timestep;  // 1ps for timestep 0.5fs
    double goal_temperature = 300;  // unit: K
    bool use_thermostat = true;
    double thermostat_duration = timestep * 1000;  // unit: fs

    // Temperature fitter
    double energy_increment = 190;  // unit: eV
    double wait_after_energy_injection = 200 * timestep;  // 100fs for timestep 0.5fs
    double measurement_time = wait_after_energy_injection + 1800 * timestep;  // 100fs + 900fs for timestep 0.5fs
    double measurement_time_currently = 0;  // unit: fs
    double average_temperature = 0;  // unit: K
    double steps_for_average = 0;
    double added_energy_sum = 0;  // unit: eV
    double last_avg_temperature = 0;

    std::ofstream traj("milestones/07/ovito/traj_" + cluster_num + ".xyz");
    std::ofstream energy("milestones/07/ovito/energy_" + cluster_num + ".csv");

    // Initialize forces
    e_pot = ducastelle(atoms, neighbors_list, neighbors_cutoff - 1);

    while (current_time < total_time) {
        // Verlet step 1
        Acceleration_t acceleration = atoms.forces / mass;
        verlet_step1(atoms.positions, atoms.velocities, acceleration, timestep);

        // Update forces
        neighbors_list.update(atoms);
        e_pot = ducastelle(atoms, neighbors_list, neighbors_cutoff - 1);

        // Verlet step 2
        acceleration = atoms.forces / mass;
        verlet_step2(atoms.velocities, acceleration, timestep);

        // Berendsen Thermostat (fit velocities)
        if (use_thermostat) {
            berendsen_thermostat(atoms, goal_temperature, timestep, thermostat_relaxation_time, false);
            if (current_time >= thermostat_duration) {
                use_thermostat = false;
                std::cout << "Stop using Thermostat\n";
            }
        // Alter temperature by rescaling velocities, then wait for some time, then measure over time, repeat
        } else {
            measurement_time_currently += timestep;

            if (wait_after_energy_injection < measurement_time_currently && measurement_time_currently < measurement_time) {
                average_temperature += atoms.temperature(false);
                steps_for_average++;
            }
            else if (measurement_time_currently >= measurement_time) {
                // Write to files
                write_xyz(traj, atoms);
                // TODO should be last_avg_temperature / steps_for_average or so
                write_E_T_C(energy, added_energy_sum, average_temperature / steps_for_average, energy_increment / (average_temperature / steps_for_average - last_avg_temperature));
                // std::cout << current_time << "/" << total_time << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin(mass) << "\tTotal energy: " << e_pot + atoms.e_kin(mass) << "\tTemperature: " << atoms.temperature(mass, false) << "\n";
                std::cout << current_time << "/" << total_time << "\tEnergy: " << e_pot + atoms.e_kin() << "\tAdded energy: " << added_energy_sum << "\tTemperature: " << average_temperature / steps_for_average << "\n";

                // Increment energy
                atoms.velocities *= std::sqrt(1 + energy_increment / atoms.e_kin());
                added_energy_sum += energy_increment;

                last_avg_temperature = average_temperature;
                average_temperature = 0;
                steps_for_average = 0;
                measurement_time_currently = 0;
            }
        }

        // Increment time
        current_time += timestep;
    }

    traj.close();
    std::cout << "Done simulating\n";
    return 0;
}
