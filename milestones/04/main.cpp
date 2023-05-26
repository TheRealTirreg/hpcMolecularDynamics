#include "atoms.h"
#include "types.h"
#include "lj_direction_summation.h"
#include "verlet.h"
#include "xyz.h"
#include "write_file.h"
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

    double e_pot = 0;

    double current_time = 0;

    int write_every_n_steps = 1000;
    int write_counter = 0;

    std::ofstream traj("milestones/04/ovito/traj.xyz");
    std::ofstream energy("milestones/04/ovito/energy.csv");

    // Initialize forces
    e_pot = lj_direct_summation_vectorized_comparism(atoms, epsilon, sigma);
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
        e_pot = lj_direct_summation_vectorized_comparism(atoms, epsilon, sigma);

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
