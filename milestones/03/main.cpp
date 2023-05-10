#include <iostream>
#include "verlet.h"
#include "types.h"

#ifdef USE_MPI
#include <mpi.h>
#endif


int main() {
    size_t num_atoms = 10;

    Positions_t positions(3, num_atoms);
    positions.setRandom();

    Velocities_t velocities(3, num_atoms);
    velocities.setRandom();
    velocities(0, 0) = 1;

    Forces_t forces(3, num_atoms);
    forces.setRandom();
    forces.row(0).setOnes();

    Masses_t masses(num_atoms);
    masses.setRandom();
    masses = masses.abs();  // masses cannot be negative

    Acceleration_t accelerations = forces.rowwise() / masses.transpose();

    double timestep = 1;

    size_t num_steps = 10;
    std::cout << "Starting velocity-verlet-algorithm for " << num_steps << " steps\n";

    // Below is some MPI code, try compiling with `cmake -DUSE_MPI=ON ..`
#ifdef USE_MPI
    MPI_Init(&argc, &argv);

    // Retrieve process infos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    for (size_t i = 0; i < num_steps; i++) {
        verlet_step1(positions, velocities, accelerations, timestep);
        // compute forces
        verlet_step2(velocities, accelerations, timestep);
        std::cout << "Step " << i << " done:\n"
                  << "Positions:\n" << positions << "\n\n"
                  << "Velocities:\n" << velocities << "\n\n"
                  << "Forces:\n" << forces << "\n\n"
                  << "Masses:\n" << masses.transpose() << "\n\n"
                  << "Accelerations:\n" << accelerations << "\n\n";
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
