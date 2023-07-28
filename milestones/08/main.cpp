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
    std::string cluster_num = "55";
    std::string filename = "external/cluster_" + cluster_num + ".xyz";
    auto [names, positions]{read_xyz(filename)};

    Atoms atoms = Atoms(Names_t(names), Positions_t(positions));
    std::cout << "num atoms: " << atoms.nb_atoms() << "\n";

    double mass_gold = 196.97;
    double mass_unit_factor = 103.6;
    double mass = mass_gold * mass_unit_factor;  // unit: g/mol
    double timestep = 0.5;  // unit: fs
    double total_time = 500000;  // unit: fs
    double current_time = 0;  // unit: fs

    double e_pot = 0;

    // Neighbors
    double neighbors_cutoff = 10;
    NeighborList neighbors_list{neighbors_cutoff};
    neighbors_list.update(atoms);

    // Thermostat
    double thermostat_relaxation_time = 500 * timestep;  // 250fs for timestep 0.5fs
    double goal_temperature = 300;  // unit: K
    bool use_thermostat = true;
    double thermostat_duration = timestep * 1000;  // unit: fs

    // Temperature fitter
    double energy_increment = 0.2;  // unit: eV
    double wait_after_energy_injection = 2000 * timestep;  // 1ps for timestep 0.5fs
    double relaxation_time = wait_after_energy_injection + 5000 * timestep;  // 1ps + 2.5ps for timestep 0.5fs
    double relaxation_time_currently = 0;
    double average_energy = 0;  // unit: eV
    double average_temperature = 0;  // unit: K
    double steps_for_average = 0;

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
            berendsen_thermostat(atoms, goal_temperature, timestep, thermostat_relaxation_time, mass, false);
            if (current_time >= thermostat_duration) {
                use_thermostat = false;
                std::cout << "Stop using Thermostat\n";
            }
        // Alter temperature by rescaling velocities
        } else {
            relaxation_time_currently += timestep;
            if (wait_after_energy_injection < relaxation_time_currently && relaxation_time_currently < relaxation_time) {
                average_energy += e_pot + atoms.e_kin(mass);
                average_temperature += atoms.temperature(mass, false);
                steps_for_average++;
            }
            else if (relaxation_time_currently >= relaxation_time) {
                // Increment energy
                atoms.velocities *= std::sqrt(1 + energy_increment / atoms.e_kin(mass));

                // Write to files
                write_xyz(traj, atoms);
                write_E_T(energy, average_energy / steps_for_average, average_temperature / steps_for_average);
                // std::cout << current_time << "/" << total_time << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin(mass) << "\tTotal energy: " << e_pot + atoms.e_kin(mass) << "\tTemperature: " << atoms.temperature(mass, false) << "\n";
                std::cout << current_time << "/" << total_time << "\tEnergy: " << average_energy / steps_for_average << "\tTemperature: " << average_temperature / steps_for_average << "\n";

                average_energy = 0;
                average_temperature = 0;
                steps_for_average = 0;
                relaxation_time_currently = 0;
            }
        }

        // Increment time
        current_time += timestep;
    }

    traj.close();
    std::cout << "Done simulating\n";
    return 0;
}
