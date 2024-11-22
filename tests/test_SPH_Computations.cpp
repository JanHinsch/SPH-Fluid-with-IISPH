#define CATCH_CONFIG_MAIN
#include "../external/catch2/catch.hpp"
#include "../include/ParticleSystem.h"
#include "../include/Grid.h"
#include "../include/Globals.h"
#include "../include/SimulationEOS.h"
#include "../include/SimulationIISPH.h"



int numParticles = 121;

void initializeParticlesForTest(ParticleSystem &particleSystem) {
    // set roh0
    density_rest = 0.01f;

    // place the particles in a box pattern (with h spacing)
    particleSystem.addBox(numParticles);

    // initialize index search arrays C and L
    particleSystem.C.clear();
    for (auto& p : particleSystem.m_particles) {
        int cellIndex = particleSystem.getCellIndex(p.position);
        particleSystem.C[cellIndex]++;

        // set mass; set density and pressure to 0
        float particleVolume = std::pow(h, 2);
        float particleMass = particleVolume * density_rest;
        p.mass = particleMass;
        p.density = 0.0f; // Initialize density to zero
        p.pressure = 0.0f; // Initialize pressure to zero
    }
    particleSystem.L = particleSystem.generateSortedList(particleSystem.m_particles);
}


TEST_CASE("SPH Computation Test (Density, Pressure Computations)", "[Simulation]") {


    ParticleSystem particleSystem(numParticles);
    initializeParticlesForTest(particleSystem);

    // Im doing Index Search on Particle in the middle of the box of particles
    Particle refParticle = particleSystem.m_particles[60];

    std::vector<Particle*> neighbours;
    // Perform the neighbor search on middle Particle
    particleSystem.iterateNeighbours(refParticle, neighbours);

    // The local density of refParticle should match the initialized restDensity since the particles are uniformly spaced
    // and the kernel sum should be approximately 1.
    SECTION("Density Computation Test") {

        float solutionDensity = SPHComputationsEOS::computeDensity(refParticle, neighbours);

        printf("neighbours: %lu \n", neighbours.size());
        printf("mass: %f \n", refParticle.mass);
        printf("solutionDensity: %f \n", solutionDensity);
        printf("density_rest: %f \n", density_rest);

        // error tolerance
        float tolerance = 0.01f;

        // REQUIRE(Approx(solutionDensity).epsilon(0.01) == density_rest);
        REQUIRE(std::abs(solutionDensity - density_rest) < tolerance);
    }

    // If the local particle density (density_i) equals the rest density , the pressure should be zero.
    SECTION("Pressure Computation Test") {

        stiffness_constant_k = 1.2f;
        // error tolerance
        float tolerance = 0.01f;


        float density_i = SPHComputationsEOS::computeDensity(refParticle, neighbours);
        // REQUIRE(Approx(density_i).epsilon(0.01) == density_rest);
        REQUIRE(std::abs(density_i - density_rest) < tolerance);

        float pressure_i = SPHComputationsEOS::computePressure(density_i);
        // REQUIRE(Approx(pressure_i).epsilon(0.1) == 0.0f);

        REQUIRE(std::abs(pressure_i - 0.0f) < tolerance);
    }
}