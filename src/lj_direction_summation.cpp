//
// Created by Gerrit Freiwald on 10.05.2023.
//

#include "lj_direction_summation.h"


double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    const long n = (long)atoms.nb_atoms();
    double total_energy = 0;
    atoms.forces.setZero();

    for (long i = 0; i < n; i++) {
        for (long j = i + 1; j < n; j++) {
            // Get distance r_ij
            const Eigen::Vector3d r1(atoms.positions.col(i));
            const Eigen::Vector3d r2(atoms.positions.col(j));
            const double distance = (r2 - r1).norm();

            // Add E_{pot} for r_ij to total potential energy
            const double pow_6 = std::pow(sigma / distance, 6);
            total_energy += 4 * epsilon * (pow_6 * pow_6 - pow_6);

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
    }

    return total_energy;
}

double lj_direct_summation_vectorized(Atoms &atoms, double epsilon, double sigma) {
    double total_energy = 0;
    atoms.forces.setZero();

    const long n = (long)atoms.nb_atoms();

    for (long i = 0; i < n; i++) {
        const Eigen::Array3d r_i = atoms.positions.col(i);
        const Eigen::Array3Xd r_js = atoms.positions.rightCols(n - i - 1);

        // Compute distances between r_i and all r_js with j > i elements
        const Eigen::Array3Xd directions = r_js.colwise() - r_i;
        const Eigen::ArrayXd distances = directions.colwise().norm();

        // Compute pair interaction energies for r_ij
        const Eigen::ArrayXd pow_6 = (sigma / distances).pow(6);
        const Eigen::ArrayXd pair_energies = 4 * epsilon * (pow_6 * pow_6 - pow_6);

        // Add to total potential energy
        total_energy += pair_energies.sum();

        // Compute forces
        const Eigen::ArrayXd force_factor = 24 * epsilon / distances * (pow_6 - 2 * pow_6.pow(2));
        const Eigen::Array3Xd force_directions = directions.colwise().normalized();
        const Eigen::Array3Xd forces = force_directions.rowwise() * force_factor.transpose();

        // Update forces
        atoms.forces.rightCols(n - i - 1) -= forces;
        atoms.forces.col(i) += forces.rowwise().sum();
    }

    return total_energy;
}
