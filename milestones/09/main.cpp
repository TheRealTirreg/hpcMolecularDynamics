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

void simulate(std::string whisker_name) {
    // Retrieve process infos and setup domain
    int rank; int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Using MPI with rank " << rank << " and size " << size << "\n";

    std::string root = "./";
    // std::string root = "../../../";  // Needed because the location of the executable is within the cmake-build folder

    std::string filename = root + "whiskers/" + whisker_name + ".xyz";
    auto [names, positions]{read_xyz(filename)};

    Atoms atoms = Atoms(Names_t(names), Positions_t(positions));
    std::cout << "num atoms: " << atoms.nb_atoms() << "\n";

    double energy_injection = 0.0; // Unit: eV
    double mass_gold = 196.97;
    double mass_unit_factor = 103.6;
    double mass = mass_gold * mass_unit_factor;  // unit: g/mol
    atoms.set_masses(mass);
    double timestep = 0.5;  // unit: fs
    double max_total_time = 500000;  // unit: fs
    double current_time = 0;  // unit: fs

    // Set up domains
    double whisker_x = atoms.positions.row(0).maxCoeff();
    double whisker_y = atoms.positions.row(1).maxCoeff();
    double whisker_z = atoms.positions.row(2).maxCoeff();
    std::cout << "Whisker dimensions: " << whisker_x << ", " << whisker_y << ", " << whisker_z << "\n";
    double current_z = whisker_z;
    double z_increment = 5;
    double rim = 5;
    Domain domain{MPI_COMM_WORLD,
                  {whisker_x + rim, whisker_y + rim, current_z},
                  {1, 1, 1},
                  {1, 1, 1}};
    if (rank == 0) std::cout << "Setting domain size to " << whisker_x + rim << ", " << whisker_y + rim << ", " << current_z << "\n";

    // Neighbors
    double neighbors_cutoff = 10;
    NeighborList neighbors_list{neighbors_cutoff};
    neighbors_list.update(atoms);

    // Temperature fitter
    double wait_after_energy_injection = 200 * timestep;  // 100fs for timestep 0.5fs
    double measurement_time = wait_after_energy_injection + 1800 * timestep;  // 100fs + 900fs for timestep 0.5fs
    double measurement_time_currently = 0;  // unit: fs
    double avg_temperature_local = 0;  // unit: K
    double steps_for_average = 0;
    double added_energy_sum = 0;  // unit: eV
    double last_avg_temperature_total = 0;
    double avg_temperature_total;
    double e_pot_total;
    double energy_total;

    // Initialize IO
    std::ofstream traj;
    std::ofstream energy;

    // Initialize forces
    double e_pot_local = ducastelle_mp(atoms, neighbors_list, domain.nb_local(), neighbors_cutoff - 0.1);

    // Make domain ready for main loop
    domain.enable(atoms);
    // domain.scale(atoms, {whisker_x + rim, whisker_y + rim, current_z});
    // if (rank == 0) std::cout << "Scaling to " << whisker_x + rim << ", " << whisker_y + rim << ", " << current_z << "\n";
    domain.update_ghosts(atoms, 2 * (neighbors_cutoff - 0.1));
    neighbors_list.update(atoms);

    // Write initial state
    // domain.disable(atoms);
    if (rank == 0) {
        traj = std::ofstream(root + "milestones/09/ovito/traj_" + whisker_name + ".xyz");
        energy = std::ofstream(root + "milestones/09/ovito/energy_" + whisker_name + ".csv");
        write_xyz(traj, atoms);
        std::cout << "Wrote initial state\n";
    }
    // domain.update_ghosts(atoms, 2 * (neighbors_cutoff - 0.1));
    // neighbors_list.update(atoms);

    // Before the actual simulation, run the thermostat
    bool use_thermostat = true;
    double thermostat_goal_tmp = 300;
    double thermostat_relaxation_time = 500;
    double thermostat_total_time = 1000;

    // Main simulation loop
    while (current_time < max_total_time) {
        // Verlet step 1
        Acceleration_t acceleration = atoms.forces / mass;
        verlet_step1(atoms.positions, atoms.velocities, acceleration, timestep);

        // Exchange information between domains
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * (neighbors_cutoff - 0.1));

        // Update forces
        neighbors_list.update(atoms);
        e_pot_local = ducastelle_mp(atoms, neighbors_list, domain.nb_local(), neighbors_cutoff - 0.1);

        // Verlet step 2
        acceleration = atoms.forces / mass;
        verlet_step2(atoms.velocities, acceleration, timestep);

        // Berendsen Thermostat (fit velocities)
        if (use_thermostat) {
            berendsen_thermostat_mp(atoms, domain, thermostat_goal_tmp, timestep, thermostat_relaxation_time, false);
            if (current_time >= thermostat_total_time) {
                use_thermostat = false;
                if (rank == 0) std::cout << "Stop using Thermostat\n";
            }
        // Alter temperature by rescaling velocities, then wait for some time, then measure over time, repeat
        } else {
            measurement_time_currently += timestep;

            if (wait_after_energy_injection < measurement_time_currently && measurement_time_currently < measurement_time) {
                const double local_tmp = atoms.local_temperature(domain.nb_local(), false);
                // avg_temperature_total is calculated by summing the average local temperature summands
                // e.g. we have 3 atoms in 2 domains, avg_tmp_total=100K, dom(1) with 1 atom has tmp_local=60, dom(2) with 2 atoms has tmp_local=120
                //      then avg_tmp_total = (60*1/3) + (120*2/3)
                avg_temperature_local += local_tmp * domain.nb_local() / atoms.nb_atoms();
                steps_for_average++;
            }
            else if (measurement_time_currently >= measurement_time) {
                if (rank == 0) std::cout << "Disabling domains on timestep " << current_time << "/" << max_total_time << "\tAtoms: " << atoms.positions.cols() << "\n";
                domain.disable(atoms);
                if (rank == 0) std::cout << "total atoms: " << atoms.positions.cols() << "\n";

                // Calculate total system temperatures and energies
                avg_temperature_local /= steps_for_average;
                avg_temperature_total = MPI::allreduce(avg_temperature_local, MPI_SUM, MPI_COMM_WORLD);
                e_pot_total = MPI::allreduce(e_pot_local, MPI_SUM, MPI_COMM_WORLD);
                energy_total = e_pot_total + atoms.e_kin();

                // Write to files
                if (rank == 0) {
                    write_xyz(traj, atoms);
                    write_E_T(energy, energy_total, avg_temperature_total);
                    std::cout << current_time << "/" << max_total_time << "\tE_kin: " << atoms.e_kin() << "\tE_pot_total: " << e_pot_total << "\tEnergy: " << energy_total << "\tTotal added energy: " << added_energy_sum << "\tTemperature: " << avg_temperature_total << "\n";
                }

                // Increment energy
                atoms.velocities *= std::sqrt(1 + energy_injection / atoms.e_kin());
                added_energy_sum += energy_injection;

                // Rescaling domains
                // current_z += z_increment;
                // domain.scale(atoms, {whisker_x + rim, whisker_y + rim, current_z});
                // if (rank == 0) std::cout << "Scaling to " << whisker_x + rim << ", " << whisker_y + rim << ", " << current_z << "\n";

                // Reset variables and start domain separation
                last_avg_temperature_total = avg_temperature_total;
                avg_temperature_total = 0;
                steps_for_average = 0;
                measurement_time_currently = 0;

                domain.enable(atoms);
                domain.update_ghosts(atoms, 2 * (neighbors_cutoff - 0.1));
                neighbors_list.update(atoms);
            }
        }

        // Increment time
        current_time += timestep;
    }

    // Clean up
    domain.disable(atoms);
    if (rank == 0) {
        traj.close();
        energy.close();
    }
    std::cout << "Done simulating " << whisker_name << "\n";
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    simulate("whisker_small");

    std::cout << "Done.\n";

    MPI_Finalize();

    return 0;
}
