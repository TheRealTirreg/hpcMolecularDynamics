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
    std::string cluster_num = "309";
    std::string filename = "external/cluster_" + cluster_num + ".xyz";
    auto [names, positions]{read_xyz(filename)};

    Atoms atoms = Atoms(Names_t(names), Positions_t(positions));
    std::cout << "num atoms: " << atoms.nb_atoms() << "\n";

    double mass_gold = 196.97;
    double mass_unit_factor = 103.6;
    double mass = mass_gold * mass_unit_factor;  // unit: g/mol
    double timestep = 0.5;  // unit: fs
    double total_time = 50000;  // unit: fs
    double current_time = 0;  // unit: fs

    double e_pot = 0;

    // Neighbors
    double neighbors_cutoff = 9;
    NeighborList neighbors_list{neighbors_cutoff};
    neighbors_list.update(atoms);

    // Thermostat
    double relaxation_time = 500 * timestep;  // 1ps for timestep = 0.5fs
    double goal_temperature = 300;  // unit: K
    bool use_thermostat = true;
    double thermostat_duration = timestep * 1000;  // unit: fs

    // Temperature fitter
    double energy_increment = 1;  // unit: eV
    double wait_after_energy_injection = 100 * timestep;
    double relaxation_time_currently = 0;

    int write_every_n_steps = 100;
    int write_counter = 100;

    std::ofstream traj("milestones/07/ovito/traj_" + cluster_num + ".xyz");
    std::ofstream energy("milestones/07/ovito/energy_" + cluster_num + ".csv");

    // Initialize forces
    e_pot = ducastelle(atoms, neighbors_list, neighbors_cutoff - 1);

    while (current_time < total_time) {
        // Write to file
        if (write_counter >= write_every_n_steps && relaxation_time_currently >= wait_after_energy_injection) {
            write_xyz(traj, atoms);
            write_energy(energy, e_pot, atoms.e_kin(mass), atoms.temperature(mass, false));

            // Compute total energy in the system
            std::cout << current_time << "/" << total_time << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin(mass) << "\tTotal energy: " << e_pot + atoms.e_kin(mass) << "\tTemperature: " << atoms.temperature(mass, false) << "\n";

            write_counter = 0;
        }

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
            berendsen_thermostat(atoms, goal_temperature, timestep, relaxation_time, mass, false);
            if (current_time >= thermostat_duration) {
                use_thermostat = false;
                std::cout << "Stop using Thermostat\n";
            }
        // Alter temperature by rescaling velocities
        } else {
            relaxation_time_currently += timestep;
            if (relaxation_time_currently > relaxation_time) {
                std::cout << "incrementing velocities by " << std::sqrt(1 + energy_increment / atoms.e_kin(mass)) << "\n";
                atoms.velocities *= std::sqrt(1 + energy_increment / atoms.e_kin(mass));
                relaxation_time_currently = 0;
            }
        }

        // increment time
        current_time += timestep;
        write_counter++;
    }

    traj.close();
    std::cout << "Done simulating\n";
    return 0;
}
