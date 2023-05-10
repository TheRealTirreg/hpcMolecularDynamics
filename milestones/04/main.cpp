#include <iostream>
#include "verlet.h"
#include "types.h"

#ifdef USE_MPI
#include <mpi.h>
#endif


int main(int argc, char *argv[]) {
    size_t num_atoms = 10;

    Positions_t positions(3, num_atoms);
    positions.setZero();

    Velocities_t velocities(3, num_atoms);
    velocities.setZero();
    velocities(0, 0) = 1;

    Forces_t forces(3, num_atoms);
    forces.setZero();
    forces.row(0).setOnes();

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
        verlet_step1(positions, velocities, forces, timestep);
        // compute forces
        verlet_step2(velocities, forces, timestep);
        std::cout << "Step " << i << " done:\n"
                  << "Positions:\n" << positions << "\n\n"
                  << "Velocities:\n" << velocities << "\n\n"
                  << "Forces:\n" << forces << "\n\n";
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
