//
// Created by Gerrit Freiwald on 10.05.2023.
//

#include "lj_direction_summation.h"


double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    // Math for deriving the Lennard-Jones potential
    // (r being the positions and j describing the cartesian directions):
    // E_{pot} = 0.5 \sum_{ij} 4ε [ (σ/r_ij)^12 - (σ/r_ij)^6 ]
    // => f_j = -∇_j E_{pot} = 48ε/r^2 * (σ/r)^6 * [(σ/r)^6 - 0.5]
    // => f_j = \sum_i f_ij for all atoms i<n

    double ret = 0;
    double eps48 = 48 * epsilon;
    double sigma_pow_6 = sigma * sigma * sigma * sigma * sigma * sigma;

    for (long i = 0; i < (long)atoms.nb_atoms(); i++) {
        double x = atoms.positions(0, i);
        double y = atoms.positions(1, i);
        double z = atoms.positions(2, i);

        double x_pow_2 = x * x;
        double y_pow_2 = y * y;
        double z_pow_2 = z * z;
        double x_pow_6 = x_pow_2 * x_pow_2 * x_pow_2;
        double y_pow_6 = y_pow_2 * y_pow_2 * y_pow_2;
        double z_pow_6 = z_pow_2 * z_pow_2 * z_pow_2;
        double sigma_div_x = sigma_pow_6 / x_pow_6;
        double sigma_div_y = sigma_pow_6 / y_pow_6;
        double sigma_div_z = sigma_pow_6 / z_pow_6;

        ret += eps48 / x_pow_2 * sigma_div_x * (sigma_div_x - 0.5);
        ret += eps48 / y_pow_2 * sigma_div_y * (sigma_div_y - 0.5);
        ret += eps48 / z_pow_2 * sigma_div_z * (sigma_div_z - 0.5);
    }

    return ret;
}
