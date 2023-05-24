//
// Created by Gerrit Freiwald on 10.05.2023.
//

#ifndef MY_MD_CODE_LJ_DIRECTION_SUMMATION_H
#define MY_MD_CODE_LJ_DIRECTION_SUMMATION_H

#include "types.h"
#include "atoms.h"

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

double lj_direct_summation_vectorized(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

double lj_direct_summation_vectorized_comparism(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

#endif // MY_MD_CODE_LJ_DIRECTION_SUMMATION_H
