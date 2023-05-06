#include "verlet.h"
#include <Eigen/Dense>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif


int main(int argc, char *argv[]) {
    size_t num_steps = 10;
    std::cout << "Starting velocity-verlet-algorithm for " << num_steps << " steps\n";

    // Below is some MPI code, try compiling with `cmake -DUSE_MPI=ON ..`
#ifdef USE_MPI
    MPI_Init(&argc, &argv);

    // Retrieve process infos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    for (int i = 0; i < num_steps; i++) {
        std::cout << "Step " << i << "\n";
        // verlet_step1();
        // compute forces
        // verlet_step2();
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
