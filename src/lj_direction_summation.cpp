//
// Created by Gerrit Freiwald on 10.05.2023.
//

#include "lj_direction_summation.h"


double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double total_energy = 0;
    atoms.forces.setZero();

    for (long i = 0; i < (long)atoms.nb_atoms(); i++) {
        for (long j = i + 1; j < (long)atoms.nb_atoms(); j++) {
            // Get distance r_ij
            Eigen::Vector3d r1(atoms.positions.col(i));
            Eigen::Vector3d r2(atoms.positions.col(j));
            double distance = (r2 - r1).norm();

            // Add E_{pot} for r_ij to total potential energy
            double pow_6 = std::pow(sigma / distance, 6);
            total_energy += 4 * epsilon * (pow_6 * pow_6 - pow_6);

            // Calculate forces (which is -∇_j E_{pot}).
            // i and j repel each other (Newton's 3rd Law)
            // Here, the chain rule is used for forceDirection:
            // V being pair interaction energy between i and j
            // f_k = -δE_{pot}/δr_k = \sum_{i<j} δV/δr_ij * dist_normalized_ij
            double force = 24 * epsilon / distance * (pow_6 - 2 * pow_6 * pow_6);
            Eigen::Array3d forceDirection = force * (r2 - r1).normalized();
            atoms.forces.col(i) += forceDirection;
            atoms.forces.col(j) -= forceDirection;
        }
    }

    return total_energy;
}
