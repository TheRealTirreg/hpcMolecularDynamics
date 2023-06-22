//
// Created by Gerrit Freiwald on 22.05.2023.
//

#include "write_file.h"
#include <cmath>

void write_energy(std::ofstream &file, double epot, double ekin) {
    file << epot << "\t" << ekin << "\t" << epot + ekin << "\n";
}

void write_lattice_cube(const std::string &filename, int n, double lattice_const_epsilon) {
    std::ofstream file(filename);
    int cube_root_n = std::floor(cbrt(n));
    int num_atoms = std::pow(cube_root_n, 3);

    file << num_atoms << "\n\n";

    for (int x = 0; x < cube_root_n; x++) {
        for (int y = 0; y < cube_root_n; y++) {
            for (int z = 0; z < cube_root_n; z++) {
                double px = x * lattice_const_epsilon;
                double py = y * lattice_const_epsilon;
                double pz = z * lattice_const_epsilon;
                double vx = 0;
                double vy = 0;
                double vz = 0;
                file << "H\t" << px << "\t" << py << "\t" << pz << "\t" << vx << "\t" << vy << "\t" << vz << "\n";
            }
        }
    }

}