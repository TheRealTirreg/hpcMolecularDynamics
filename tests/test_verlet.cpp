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

TEST(Verlet, VerletStep1EigenNoAccelerations) {
    size_t num_atoms = 2;

    Positions_t positions(3, num_atoms);
    positions.setZero();

    Velocities_t velocities(3, num_atoms);
    velocities.setZero();
    velocities(0, 0) = 1;

    Acceleration_t accelerations(3, num_atoms);
    accelerations.setZero();

    double timestep = 1;

    verlet_step1(positions, velocities, accelerations, timestep);

    Positions_t expectedPositions(3, num_atoms);
    expectedPositions.setZero();
    expectedPositions(0, 0) = 1;

    ASSERT_TRUE(expectedPositions.isApprox(positions));
}

TEST(Verlet, VerletStep1EigenWithAccelerations) {
    size_t num_atoms = 2;

    Positions_t positions(3, num_atoms);
    positions.setZero();

    Velocities_t velocities(3, num_atoms);
    velocities.setZero();
    velocities(0, 0) = 1;

    Acceleration_t accelerations(3, num_atoms);
    accelerations.setOnes();

    double timestep = 1;

    verlet_step1(positions, velocities, accelerations, timestep);

    Positions_t expectedPositions(3, num_atoms);
    expectedPositions.setZero();
    expectedPositions(0, 0) = 1;
    expectedPositions += 0.5;

    ASSERT_TRUE(expectedPositions.isApprox(positions));
}

TEST(Verlet, VerletStep2EigenWithAccelerations) {
    size_t num_atoms = 2;

    Velocities_t velocities(3, num_atoms);
    velocities.setZero();
    velocities(0, 0) = 1;

    Acceleration_t accelerations(3, num_atoms);
    accelerations.setOnes();

    double timestep = 1;

    verlet_step2(velocities, accelerations, timestep);

    Velocities_t expectedVelocities(3, num_atoms);
    expectedVelocities.setZero();
    expectedVelocities(0, 0) = 1;
    expectedVelocities += 0.5;

    ASSERT_TRUE(expectedVelocities.isApprox(velocities));
}

TEST(Verlet, AnalyticDoubles) {
    double x = 0;  double y = 0;  double z = 0;
    double vx = 1; double vy = 0; double vz = 0;
    double fx = 1; double fy = 1; double fz = 1;
    double timestep = 1; double num_timesteps = 100;

    double x0 = x; double y0 = y; double z0 = z;
    double vx0 = vx; double vy0 = vy; double vz0 = vz;
    double fx0 = fx; double fy0 = fy; double fz0 = fz;

    // Get verlet values
    for (int i = 0; i < num_timesteps; i++) {
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);
        // todo compute forces here
        verlet_step2(vx, vy, vz, fx, fy, fz, timestep);
    }

    // Get analytical values
    // v(t) = F/m*t + v0
    // r(t) = 0.5*F/m*t*t + v0*t + r0
    double t = num_timesteps * timestep;
    double analytical_vx = fx0 * MASS_RECIPROKE * t + vx0;
    double analytical_vy = fy0 * MASS_RECIPROKE * t + vy0;
    double analytical_vz = fz0 * MASS_RECIPROKE * t + vz0;

    double t_squared = t * t;
    double analytical_x = 0.5 * fx0 * MASS_RECIPROKE * t_squared + vx0 * t + x0;
    double analytical_y = 0.5 * fy0 * MASS_RECIPROKE * t_squared + vy0 * t + y0;
    double analytical_z = 0.5 * fz0 * MASS_RECIPROKE * t_squared + vz0 * t + z0;

    ASSERT_DOUBLE_EQ(x, analytical_x);
    ASSERT_DOUBLE_EQ(y, analytical_y);
    ASSERT_DOUBLE_EQ(z, analytical_z);
    ASSERT_DOUBLE_EQ(vx, analytical_vx);
    ASSERT_DOUBLE_EQ(vy, analytical_vy);
    ASSERT_DOUBLE_EQ(vz, analytical_vz);
}

TEST(Verlet, AnalyticEigen) {
    size_t num_atoms = 10;

    Positions_t positions(3, num_atoms);
    positions.setRandom();

    Velocities_t velocities(3, num_atoms);
    velocities.setRandom();

    Forces_t forces(3, num_atoms);
    forces.setRandom();

    Masses_t masses(num_atoms);
    masses.setRandom();

    Acceleration_t accelerations = forces.rowwise() / masses.transpose();

    double timestep = 1; double num_timesteps = 100;

    // copy start values
    Positions_t r0 = Positions_t(positions);
    Velocities_t v0 = Velocities_t(velocities);
    Acceleration_t a0 = Acceleration_t(accelerations);

    // Get verlet values
    for (int i = 0; i < num_timesteps; i++) {
        verlet_step1(positions, velocities, accelerations, timestep);
        // practically, compute forces here
        verlet_step2(velocities, accelerations, timestep);
    }

    // Get analytical values
    // v(t) = F/m*t + v0
    // r(t) = 0.5*F/m*t*t + v0*t + r0
    double t = num_timesteps * timestep;
    double t_squared = t * t;
    Velocities_t analytical_v = a0 * t + v0;
    Positions_t analytical_r = 0.5 * a0 * t_squared + v0 * t + r0;

    for (size_t i = 0; i < num_atoms; i++) {
        EXPECT_NEAR(positions(0, i), analytical_r(0, i), 1e-10);
        EXPECT_NEAR(positions(1, i), analytical_r(1, i), 1e-10);
        EXPECT_NEAR(positions(2, i), analytical_r(2, i), 1e-10);
    }
}
