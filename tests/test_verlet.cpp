//
// Created by Gerrit Freiwald on 08.05.2023.
//

#include "verlet.h"
#include <gtest/gtest.h>

TEST(Verlet, VerletStep1DoublesNoForces) {
    double x = 0;  double y = 0;  double z = 0;
    double vx = 1; double vy = 0; double vz = 0;
    double fx = 0; double fy = 0; double fz = 0;
    double timestep = 1;

    verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);

    ASSERT_DOUBLE_EQ(x, 1);
    ASSERT_DOUBLE_EQ(y, 0);
    ASSERT_DOUBLE_EQ(z, 0);

    x = 0; y = 0;
    vx = 0; vy = 1;

    verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);

    ASSERT_DOUBLE_EQ(x, 0);
    ASSERT_DOUBLE_EQ(y, 1);
    ASSERT_DOUBLE_EQ(z, 0);

    y = 0; z = 0;
    vy = 0; vz = 1;

    verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);

    ASSERT_DOUBLE_EQ(x, 0);
    ASSERT_DOUBLE_EQ(y, 0);
    ASSERT_DOUBLE_EQ(z, 1);
}

TEST(Verlet, VerletStep1DoublesWithForces) {
    double x = 0;  double y = 0;  double z = 0;
    double vx = 0; double vy = 0; double vz = 0;
    double fx = 1; double fy = 1; double fz = 1;
    double timestep = 1;

    verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);

    ASSERT_DOUBLE_EQ(x, 0.5);
    ASSERT_DOUBLE_EQ(y, 0.5);
    ASSERT_DOUBLE_EQ(z, 0.5);

    x = 0; y = 0; z = 0;
    vx = 1; vy = 0; vz = 0;

    verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);

    ASSERT_DOUBLE_EQ(x, 1.5);
    ASSERT_DOUBLE_EQ(y, 0.5);
    ASSERT_DOUBLE_EQ(z, 0.5);
}

TEST(Verlet, VerletStep2DoublesWithForces) {
    double vx = 1; double vy = 0; double vz = 0;
    double fx = 1; double fy = 1; double fz = 1;
    double timestep = 1;

    verlet_step2(vx, vy, vz, fx, fy, fz, timestep);

    ASSERT_DOUBLE_EQ(vx, 1.5);
    ASSERT_DOUBLE_EQ(vy, 0.5);
    ASSERT_DOUBLE_EQ(vz, 0.5);
}

TEST(Verlet, VerletStep1EigenNoForces) {
    size_t num_atoms = 2;

    Positions_t positions(3, num_atoms);
    positions.setZero();

    Velocities_t velocities(3, num_atoms);
    velocities.setZero();
    velocities(0, 0) = 1;

    Forces_t forces(3, num_atoms);
    forces.setZero();

    double timestep = 1;

    verlet_step1(positions, velocities, forces, timestep);

    Positions_t expectedPositions(3, num_atoms);
    expectedPositions.setZero();
    expectedPositions(0, 0) = 1;

    ASSERT_TRUE(expectedPositions.isApprox(positions));
}

TEST(Verlet, VerletStep1EigenWithForces) {
    size_t num_atoms = 2;

    Positions_t positions(3, num_atoms);
    positions.setZero();

    Velocities_t velocities(3, num_atoms);
    velocities.setZero();
    velocities(0, 0) = 1;

    Forces_t forces(3, num_atoms);
    forces.setOnes();

    double timestep = 1;

    verlet_step1(positions, velocities, forces, timestep);

    Positions_t expectedPositions(3, num_atoms);
    expectedPositions.setZero();
    expectedPositions(0, 0) = 1;
    expectedPositions += 0.5;

    ASSERT_TRUE(expectedPositions.isApprox(positions));
}

TEST(Verlet, VerletStep2EigenWithForces) {
    size_t num_atoms = 2;

    Velocities_t velocities(3, num_atoms);
    velocities.setZero();
    velocities(0, 0) = 1;

    Forces_t forces(3, num_atoms);
    forces.setOnes();

    double timestep = 1;

    verlet_step2(velocities, forces, timestep);

    Velocities_t expectedVelocities(3, num_atoms);
    expectedVelocities.setZero();
    expectedVelocities(0, 0) = 1;
    expectedVelocities += 0.5;

    ASSERT_TRUE(expectedVelocities.isApprox(velocities));
}
