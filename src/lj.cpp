//
// Created by Gerrit Freiwald on 10.05.2023.
//

#include "lj.h"


double lj_neighbors(Atoms &atoms, NeighborList &neighbors, double cutoff, double cutoff_energy, double epsilon, double sigma) {
    double total_energy = 0;
    atoms.forces.setZero();

    // Iterate over cellwise neighbored atom pairs given by the neighbors data structure
    for (auto [i, j] : neighbors) {
        // Skip mirrored pairs
        if (i >= j) {
            continue;
        }

        // Get distance r_ij
        const Eigen::Vector3d r1(atoms.positions.col(i));
        const Eigen::Vector3d r2(atoms.positions.col(j));
        const double distance = (r2 - r1).norm();

        // Skip pair if distance is too high
        if (distance >= cutoff) {
            continue;
        }

        // Add E_{pot} for r_ij to total potential energy
        const double pow_6 = std::pow(sigma / distance, 6);
        total_energy += 4 * epsilon * (pow_6 * pow_6 - pow_6);
        // We move the potential energy, so that the potential energy at the cutoff is 0.
        total_energy -= cutoff_energy;

        // Calculate forces (which is -∇_j E_{pot}).
        // i and j repel each other (Newton's 3rd Law)
        // Here, the chain rule is used for forceDirection:
        // V being pair interaction energy between i and j
        // f_k = -δE_{pot}/δr_k = \sum_{i<j} δV/δr_ij * dist_normalized_ij
        const double force = 24 * epsilon / distance * (pow_6 - 2 * pow_6 * pow_6);
        Eigen::Array3d forceDirection = force * (r2 - r1).normalized();
        atoms.forces.col(i) += forceDirection;
        atoms.forces.col(j) -= forceDirection;
    }

    return total_energy;
}
