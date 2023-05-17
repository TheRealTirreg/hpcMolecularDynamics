#include "atoms.h"
#include "types.h"
#include "lj_direction_summation.h"
#include "xyz.h"
#include <iostream>


int main() {
    /*
    size_t num_atoms = 4;
    Atoms atoms = Atoms(num_atoms, true);
     */

    auto [names, positions, velocities]{read_xyz_with_velocities("../milestones/04/"
                                                                 "lj54.xyz")};
    Atoms atoms = Atoms(Positions_t(positions), Velocities_t(velocities));

    double energy{lj_direct_summation(atoms)};

    std::cout << "Total Energy: " << energy << "\n";

    return 0;
}
