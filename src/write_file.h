//
// Created by Gerrit Freiwald on 22.05.2023.
//

#ifndef MY_MD_CODE_WRITE_FILE_H
#define MY_MD_CODE_WRITE_FILE_H

#include <fstream>
#include <Eigen/Dense>

void write_energy(std::ofstream &file, double epot, double ekin, double temperature);

void write_E_T(std::ofstream &file, double etotal, double tmp);

void write_E_T_Stress(std::ofstream &file, double etotal, double tmp, Eigen::Array3d &stress, double z_domain_size);

void write_lattice_cube(const std::string &filename, int n, double lattice_const_epsilon);

#endif // MY_MD_CODE_WRITE_FILE_H
