#include "atoms.h"
#include "types.h"
#include "neighbors.h"
#include "lj.h"
#include "verlet.h"
#include "xyz.h"
#include "write_file.h"
#include <iostream>


int main() {
    std::cout << "Starting...\n";
    auto [names, positions, velocities]
        {read_xyz_with_velocities("milestones/06/lj54.xyz")};
    Atoms atoms = Atoms(Positions_t(positions), Velocities_t(velocities));

    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    double timestep = 0.001 * std::sqrt(mass * sigma * sigma / epsilon);
    double total_time = 100 * std::sqrt(mass * sigma * sigma / epsilon);

    double e_pot = 0;

    double current_time = 0;

    double lj_cutoff = 5.0;
    double neighbors_cutoff = 5.1;

    int write_every_n_steps = 1000;
    int write_counter = 0;

    std::ofstream traj("milestones/06/ovito/traj.xyz");
    std::ofstream energy("milestones/06/ovito/energy.csv");

    // Initialize neighbors list
    NeighborList neighbors_list{neighbors_cutoff};
    neighbors_list.update(atoms);

    // Calculate energy at lj_cutoff
    const double pow_6 = std::pow(sigma / lj_cutoff, 6);
    double cutoff_energy = 4 * epsilon * (pow_6 * pow_6 - pow_6);

    // Initialize forces
    e_pot = lj_neighbors(atoms, neighbors_list, lj_cutoff, cutoff_energy, epsilon, sigma);
    std::cout << current_time << "/" << total_time << "\tForces: " << atoms.forces.col(0).transpose() << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin() << "\tTotal energy: " << e_pot + atoms.e_kin() << "\n";

    while (current_time < total_time) {
        // Write to file
        if (write_counter == write_every_n_steps) {
            write_xyz(traj, atoms);
            write_energy(energy, e_pot, atoms.e_kin());
        }

        // Verlet step 1
        Acceleration_t acceleration = atoms.forces / mass;
        verlet_step1(atoms.positions, atoms.velocities, acceleration, timestep);

        // Update forces
        neighbors_list.update(atoms);
        e_pot = lj_neighbors(atoms, neighbors_list, lj_cutoff, cutoff_energy, epsilon, sigma);

        // Verlet step 2
        acceleration = atoms.forces / mass;
        verlet_step2(atoms.velocities, acceleration, timestep);

        // Debug print
        if (write_counter == write_every_n_steps) {
            // Compute total energy in the system
            std::cout << current_time << "/" << total_time << "\tForces: " << atoms.forces.col(0).transpose() << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin() << "\tTotal energy: " << e_pot + atoms.e_kin() << "\n";
            write_counter = 0;
        }

        // increment time
        current_time += timestep;
        write_counter++;
    }

    traj.close();
    std::cout << "Done simulating\n";
    return 0;
}
