#include "atoms.h"
#include "types.h"
#include "lj_direction_summation.h"
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif


int main() {
    size_t num_atoms = 10;
    Atoms atoms = Atoms(num_atoms, true);

    double energy{lj_direct_summation(atoms)};

    std::cout << "Energy: " << energy << "\n";

    return 0;
}
