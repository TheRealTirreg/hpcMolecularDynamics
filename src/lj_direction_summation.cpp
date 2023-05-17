//
// Created by Gerrit Freiwald on 10.05.2023.
//

#include "lj_direction_summation.h"

double lj_potential(double r, double epsilon, double sigma) {
    // Math for deriving the Lennard-Jones potential
    // (r being the positions and j describing the cartesian directions):
    // E_{pot} = 0.5 \sum_{ij} 4ε [ (σ/r_ij)^12 - (σ/r_ij)^6 ]
    // => f_j = -∇_j E_{pot} = 48ε/r^2 * (σ/r)^6 * [(σ/r)^6 - 0.5]
    double pow_6 = std::pow(sigma / r, 6);
    return 24 * epsilon * (2 * pow_6 * pow_6 - pow_6);
}

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double energy = 0;
    long idx_write = 0;
    long idx_read = 0;

    // Number of connections between atoms: n(n-1)/2
    Potential_t potentials((int)(atoms.nb_atoms() * (atoms.nb_atoms() - 1) * 0.5));

    for (long i = 0; i < (long)atoms.nb_atoms(); i++) {
        for (long j = 0; j < (long)atoms.nb_atoms(); j++) {

            // Only calculate the potential for i->j, as j->i = - i->j (Newtons third law)
            if (i < j) {
                // Get distance r_ij
                Eigen::Vector3d r1(atoms.positions.col(i));
                Eigen::Vector3d r2(atoms.positions.col(j));
                double r = (r1 - r2).norm();

                // Add lj_potential for r_ij to total potential energy
                double potential = lj_potential(r, epsilon, sigma);

                potentials(idx_write++) = potential;
                //std::cout << "In:  (" << i << ", " << j << "): " << potentials(idx_write-1) << "\n";
                energy += potential;
            }
            // Exclude i == j
            else if (i > j) {
                energy -= potentials(idx_read++);
                //std::cout << "Out: (" << i << ", " << j << "): " << -potentials(idx_read-1) << "\n";
            }
        }
    }

    return energy;
}
