#include "atoms.h"
#include "types.h"
#include "lj_direction_summation.h"
#include "neighbors.h"
#include "lj.h"
#include "berendsen_thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include "write_file.h"
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

int main() {
    // Create grid of atoms
    std::cout << "Creating grid...\n";
    int grid_n = 1000;
    std::string grid_filename = "milestones/05/ovito/grid" + std::to_string(grid_n) + ".xyz";
    std::string file_rel_path = "../../../milestones/05/ovito/grid" + std::to_string(grid_n) + ".xyz";
    write_lattice_cube(file_rel_path, grid_n, 1);

    std::cout << "File " << file_rel_path << "\n";
    auto [names, positions, velocities]
        {read_xyz_with_velocities(file_rel_path)};

    Atoms atoms = Atoms(Positions_t(positions), Velocities_t(velocities));
    atoms.velocities = atoms.velocities.setRandom() * 0.2;
    std::cout << "num atoms: " << atoms.nb_atoms() << "\tekin: " << atoms.e_kin() << "\ttmp: " << atoms.temperature() << "\n";

#ifdef USE_MPI
    MPI_Init(&argc, &argv);

    // Retrieve process infos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    std::cout << "Starting...\n";

    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    double timestep = 0.001 * std::sqrt(mass * sigma * sigma / epsilon);
    double total_time = 100 * std::sqrt(mass * sigma * sigma / epsilon);

    double e_pot = 0;

    double current_time = 0;

    // neighbors
    double lj_cutoff = 2.5;
    double neighbors_cutoff = 2.6;
    // NeighborList neighbors_list{neighbors_cutoff};
    // neighbors_list.update(atoms);
    const double pow_6 = std::pow(sigma / lj_cutoff, 6);
    double cutoff_energy = 4 * epsilon * (pow_6 * pow_6 - pow_6);

    // Thermostat
    double relaxation_time = 2000 * timestep;
    double goal_temperature = 0.5;

    int write_every_n_steps = 100;
    int write_counter = 0;

    std::ofstream traj("milestones/05/ovito/traj.xyz");
    std::ofstream energy("milestones/05/ovito/energy.csv");

    // Initialize forces
    e_pot = lj_direct_summation(atoms, epsilon, sigma);
    // e_pot = lj_neighbors(atoms, neighbors_list, lj_cutoff, cutoff_energy, epsilon, sigma);

    while (current_time < total_time) {
        // Write to file
        if (write_counter == write_every_n_steps) {
            write_xyz(traj, atoms);
            write_energy(energy, e_pot, atoms.e_kin(), atoms.temperature());
        }

        // Verlet step 1
        Acceleration_t acceleration = atoms.forces / mass;
        verlet_step1(atoms.positions, atoms.velocities, acceleration, timestep);

        // Update forces
        e_pot = lj_direct_summation(atoms, epsilon, sigma);
        // neighbors_list.update(atoms);
        // e_pot = lj_neighbors(atoms, neighbors_list, lj_cutoff, cutoff_energy, epsilon, sigma);

        // Verlet step 2
        acceleration = atoms.forces / mass;
        verlet_step2(atoms.velocities, acceleration, timestep);

        // Debug print
        if (write_counter == write_every_n_steps) {
            // Compute total energy in the system
            std::cout << current_time << "/" << total_time << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin() << "\tTotal energy: " << e_pot + atoms.e_kin() << "\tTemperature: " << atoms.temperature() << "\n";

            write_counter = 0;
        }

        // Berendsen Thermostat (fit velocities)
        berendsen_thermostat(atoms, goal_temperature, timestep, relaxation_time);

        // increment time
        current_time += timestep;
        write_counter++;
    }

    traj.close();
    std::cout << "Done simulating\n";

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
