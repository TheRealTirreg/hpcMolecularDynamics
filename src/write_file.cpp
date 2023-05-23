//
// Created by Gerrit Freiwald on 22.05.2023.
//

#include "write_file.h"

void write_energy(std::ofstream &file, double epot, double ekin) {
    file << epot << "\t" << ekin << "\t" << epot + ekin << "\n";
}