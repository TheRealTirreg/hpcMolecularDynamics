//
// Created by Gerrit Freiwald on 22.05.2023.
//

#ifndef MY_MD_CODE_WRITE_FILE_H
#define MY_MD_CODE_WRITE_FILE_H

#include <fstream>

void write_energy(std::ofstream &file, double epot, double ekin);

void write_lattice_cube(const std::string &filename, int n, double lattice_const_epsilon);

#endif // MY_MD_CODE_WRITE_FILE_H
