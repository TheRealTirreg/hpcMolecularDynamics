#include "atoms.h"
#include "types.h"
#include "lj_direction_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>


int main() {
    /*
    size_t num_atoms = 4;
    Atoms atoms = Atoms(num_atoms, true);
     */

    auto [names, positions, velocities]
        {read_xyz_with_velocities("/home/gerri/Code/Uni/hpcMolecularDynamics/milestones/04/lj54.xyz")};
    Atoms atoms = Atoms(Positions_t(positions), Velocities_t(velocities));

    /*
    double energy{lj_direct_summation(atoms)};
    std::cout << "Total Energy: " << energy << "\n";
    */

    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    double timestep = 0.001 * std::sqrt(mass * sigma * sigma / epsilon);
    double total_time = 100 * std::sqrt(mass * sigma * sigma / epsilon);

    double current_time = 0;

    int print_every_n_steps = 1000;
    int print_counter = 0;

    while (current_time < total_time) {
        // Verlet step 1
        Acceleration_t acceleration = atoms.forces / mass;
        verlet_step1(atoms.positions, atoms.velocities, acceleration, timestep);

        // Update forces
        double total_energy = lj_direct_summation(atoms, epsilon, sigma);
        // std::cout << current_time << "%\tTotal energy: " << total_energy << "\n";
        if (print_counter == print_every_n_steps) {
            std::cout << current_time << "/" << total_time << "\tForces: " << atoms.forces.col(0).transpose() << "\n";
            print_counter = 0;
        }

        // Verlet step 2
        acceleration = atoms.forces / mass;
        verlet_step2(atoms.velocities, acceleration, timestep);

        // increment time
        current_time += timestep;
        print_counter++;
    }

    std::cout << "Done simulating\n";
    return 0;
}
