//
// Created by Gerrit Freiwald on 6/15/23.
//

#include "atoms.h"
#include "berendsen_thermostat.h"
#include <gtest/gtest.h>

TEST(BerendsenThermostatTest, NoTemperatureChange) {
    // Create atoms with known velocities
    const size_t n = 5;

    Positions_t positions(3, n);
    positions << 0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0;

    Velocities_t velocities(3, n);
    velocities.setOnes();

    Atoms atoms{positions, velocities};

    // Set up the test environment
    double target_temperature = 5.0;
    double timestep = 0.001;
    double relaxation_time = 0.1;

    // Check the temperature before applying the thermostat
    EXPECT_NEAR(atoms.temperature(), target_temperature, 1e-6);

    // Apply the Berendsen thermostat
    berendsen_thermostat(atoms, target_temperature, timestep, relaxation_time);

    // Check that velocities are scaled correctly
    EXPECT_NEAR(atoms.velocities(0, 0), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 0), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 0), 0.997, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 1), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 1), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 1), 0.997, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 2), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 2), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 2), 0.997, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 3), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 3), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 3), 0.997, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 4), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 4), 0.997, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 4), 0.997, 1e-6);
}

TEST(BerendsenThermostatTest, TemperatureDecrease) {
    // Create atoms with known velocities
    const size_t n = 5;

    Positions_t positions(3, n);
    positions << 0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0;

    Velocities_t velocities(3, n);
    velocities.setOnes();

    Atoms atoms{positions, velocities};

    // Set up the test environment
    double target_temperature = 2.0;
    double timestep = 0.001;
    double relaxation_time = 0.1;

    // Check the temperature before applying the thermostat
    EXPECT_NEAR(atoms.temperature(), 5.0, 1e-6);

    // Apply the Berendsen thermostat
    berendsen_thermostat(atoms, target_temperature, timestep, relaxation_time);

    // Check that velocities are scaled correctly
    EXPECT_NEAR(atoms.velocities(0, 0), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 0), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 0), 0.996995486449, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 1), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 1), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 1), 0.996995486449, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 2), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 2), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 2), 0.996995486449, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 3), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 3), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 3), 0.996995486449, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 4), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 4), 0.996995486449, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 4), 0.996995486449, 1e-6);
}

TEST(BerendsenThermostatTest, TemperatureIncrease) {
    // Create atoms with known velocities
    const size_t n = 5;

    Positions_t positions(3, n);
    positions << 0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0,
        0.0, 1.0, 2.0, 3.0, 4.0;

    Velocities_t velocities(3, n);
    velocities.setOnes();

    Atoms atoms{positions, velocities};

    // Set up the test environment
    double target_temperature = 10.0;
    double timestep = 0.001;
    double relaxation_time = 0.1;

    // Check the temperature before applying the thermostat
    EXPECT_NEAR(atoms.temperature(), 5.0, 1e-6);

    // Apply the Berendsen thermostat
    berendsen_thermostat(atoms, target_temperature, timestep, relaxation_time);

    // Check that velocities are scaled correctly
    EXPECT_NEAR(atoms.velocities(0, 0), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 0), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 0), 1.00498756211, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 1), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 1), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 1), 1.00498756211, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 2), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 2), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 2), 1.00498756211, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 3), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 3), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 3), 1.00498756211, 1e-6);

    EXPECT_NEAR(atoms.velocities(0, 4), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(1, 4), 1.00498756211, 1e-6);
    EXPECT_NEAR(atoms.velocities(2, 4), 1.00498756211, 1e-6);
}
