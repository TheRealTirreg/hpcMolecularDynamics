//
// Created by Gerrit Freiwald on 10.05.2023.
//

#ifndef MY_MD_CODE_LJ_H
#define MY_MD_CODE_LJ_H

#include "types.h"
#include "atoms.h"
#include "neighbors.h"

double lj_neighbors(Atoms &atoms, NeighborList &neighbor_list, double cutoff, double cutoff_energy, double epsilon = 1.0, double sigma = 1.0);

#endif // MY_MD_CODE_LJ_H
