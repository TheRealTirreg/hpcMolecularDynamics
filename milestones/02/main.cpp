#include "verlet.h"
#include <Eigen/Dense>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif


int main(int argc, char *argv[]) {
    // initialize positions and velocities
    double x = 0;  double y = 0;  double z = 0;
    double vx = 1; double vy = 0; double vz = 0;
    double fx = 1; double fy = 1; double fz = 1;
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
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);
        // compute forces
        verlet_step2(vx, vy, vz, fx, fy, fz, timestep);
        std::cout << "Step " << i << " done:\n"
                  << x << "\t" << vx << "\t" << fx << "\t" << "\n"
                  << y << "\t" << vy << "\t" << fy << "\t" << "\n"
                  << z << "\t" << vz << "\t" << fz << "\t" << "\n";
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
