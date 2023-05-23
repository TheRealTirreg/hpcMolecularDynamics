#include "atoms.h"
#include "types.h"
#include "lj_direction_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>


int main() {
    std::cout << "Starting...\n";
    auto [names, positions, velocities]
        {read_xyz_with_velocities("milestones/04/lj54.xyz")};
    Atoms atoms = Atoms(Positions_t(positions), Velocities_t(velocities));

    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    double timestep = 0.001 * std::sqrt(mass * sigma * sigma / epsilon);
    double total_time = 100 * std::sqrt(mass * sigma * sigma / epsilon);

    double current_time = 0;

    int write_every_n_steps = 100;
    int write_counter = 0;

    std::ofstream traj("milestones/04/ovito/traj.xyz");

    // Initialize forces
    lj_direct_summation_vectorized(atoms, epsilon, sigma);

    while (current_time < total_time) {
        // Write to file
        if (write_counter == write_every_n_steps) {
            write_xyz(traj, atoms);
        }

        // Verlet step 1
        Acceleration_t acceleration = atoms.forces / mass;
        verlet_step1(atoms.positions, atoms.velocities, acceleration, timestep);

        // Update forces
        double e_pot = lj_direct_summation_vectorized(atoms, epsilon, sigma);

        // Debug print
        if (write_counter == write_every_n_steps) {
            // Compute total energy in the system
            std::cout << current_time << "/" << total_time << "\tForces: " << atoms.forces.col(0).transpose() << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin() << "\tTotal energy: " << e_pot + atoms.e_kin() << "\n";
            write_counter = 0;
        }

        // Verlet step 2
        acceleration = atoms.forces / mass;
        verlet_step2(atoms.velocities, acceleration, timestep);

        // increment time
        current_time += timestep;
        write_counter++;
    }

    traj.close();
    std::cout << "Done simulating\n";
    return 0;
}
