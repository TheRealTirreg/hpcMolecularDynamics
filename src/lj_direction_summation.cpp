//
// Created by Gerrit Freiwald on 10.05.2023.
//

#include "lj_direction_summation.h"


double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    // Math for deriving the Lennard-Jones potential:
    // E_{pot} = 0.5 \sum_{ij} 4ε [ (σ/r_ij)^12 - (σ/r_ij)^6 ]
    // => f_j = -∇_j E_{pot} = 48ε/r^2 * (σ/r)^6 * [(σ/r)^6 - 0.5]
    // => f_j = \sum_i f_ij


}
