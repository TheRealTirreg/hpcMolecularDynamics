#include "atoms.h"
#include "types.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "xyz.h"
#include "write_file.h"
#include "berendsen_thermostat.h"
#include <iostream>

#include <mpi.h>
#include "mpi_support.h"
#include "domain.h"

void simulate(std::string cluster_num) {
    // Retrieve process infos and setup domain
    int rank; int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Using MPI with rank " << rank << " and size " << size << "\n";

    std::string filename = "external/cluster_" + cluster_num + ".xyz";
    auto [names, positions]{read_xyz(filename)};

    Atoms atoms = Atoms(Names_t(names), Positions_t(positions));
    std::cout << "num atoms: " << atoms.nb_atoms() << "\n";

    double mass_gold = 196.97;
    double mass_unit_factor = 103.6;
    double mass = mass_gold * mass_unit_factor;  // unit: g/mol
    atoms.set_masses(mass);
    double timestep = 0.5;  // unit: fs
    double total_time = 500000;  // unit: fs
    double current_time = 0;  // unit: fs

    double e_pot = 0;

    // Set up domains
    double cluster_diameter = 2 * atoms.positions.row(0).maxCoeff();
    double cluster_rim = 10;
    Domain domain{MPI_COMM_WORLD, {cluster_diameter + cluster_rim, cluster_diameter + cluster_rim, cluster_diameter + cluster_rim}, {2, 2, 1}, {0, 0, 0}};

    // Shift positions of atoms, as the center of the read cluster is 0,0,0
    atoms.positions += (cluster_diameter + cluster_rim) / 2;

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

    std::ofstream traj;
    std::ofstream energy;
    if (rank == 0) {
        traj = std::ofstream("milestones/08/ovito/traj_" + cluster_num + ".xyz");
        energy = std::ofstream("milestones/08/ovito/energy_" + cluster_num + ".csv");
        write_xyz(traj, atoms);
    }

    // Initialize forces
    e_pot = ducastelle(atoms, neighbors_list, neighbors_cutoff - 1);

    if (rank == 0) std::cout << "after ducastelle\n";
    std::cout << atoms.masses.rows() << "\n";
    std::cout << atoms.positions.cols() << "\n";
    std::cout << atoms.velocities.cols() << "\n";
    std::cout << atoms.forces.cols() << "\n";
    // while(true){}

    // Make domain ready for main loop
    domain.enable(atoms);

    domain.update_ghosts(atoms, 2 * (neighbors_cutoff - 1));
    neighbors_list.update(atoms);

    if (rank == 0) std::cout << "before while\n";

    while (current_time < total_time) {
        // Verlet step 1
        Acceleration_t acceleration = atoms.forces / mass;
        verlet_step1(atoms.positions, atoms.velocities, acceleration, timestep);

        // Exchange information between domains
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * (neighbors_cutoff - 1));

        // Update forces
        neighbors_list.update(atoms);
        e_pot = ducastelle(atoms, neighbors_list, neighbors_cutoff - 1);

        // Verlet step 2
        acceleration = atoms.forces / mass;
        verlet_step2(atoms.velocities, acceleration, timestep);

        // Berendsen Thermostat (fit velocities)

        if (int(current_time) % 1000 == 0) {
            domain.disable(atoms);
            if (rank == 0) {
                // Write to files
                write_xyz(traj, atoms);
                write_E_T(energy, average_energy / steps_for_average, average_temperature / steps_for_average);
                // std::cout << current_time << "/" << total_time << "\tPot energy: " << e_pot << "\tKin energy: " << atoms.e_kin(mass) << "\tTotal energy: " << e_pot + atoms.e_kin(mass) << "\tTemperature: " << atoms.temperature(mass, false) << "\n";
                std::cout << current_time << "/" << total_time << "\tEnergy: " << average_energy / steps_for_average << "\tTemperature: " << average_temperature / steps_for_average << "\n";
            }
            domain.enable(atoms);
            domain.update_ghosts(atoms, 2 * (neighbors_cutoff - 1));
            neighbors_list.update(atoms);
            // ducastelle?
        }

        // Increment time
        current_time += timestep;
    }

    if (rank == 0) {
        traj.close();
    }
    std::cout << "Done simulating\n";
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    simulate("55");

    MPI_Finalize();

    return 0;
}
