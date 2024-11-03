#define CATCH_CONFIG_MAIN
#include "../external/catch2/catch.hpp"
#include "../include/ParticleSystem.h"
#include "../include/Grid.h"
#include "../include/Globals.h"

int numParticles = 121;

void initializeParticlesForTest(ParticleSystem &particleSystem) {

    // place the particles in a box pattern (with h spacing)
    particleSystem.addBox(121);

    // initialize index search arrays C and L
    particleSystem.C.clear();
    for (auto& p : particleSystem.m_particles) {
        int cellIndex = particleSystem.getCellIndex(p.position);
        particleSystem.C[cellIndex]++;
    }
    particleSystem.L = particleSystem.generateSortedList(particleSystem.m_particles);
}

TEST_CASE("Neighbor Search", "[ParticleSystem]") {
    ParticleSystem particleSystem(numParticles);
    initializeParticlesForTest(particleSystem);

    // Im doing Index Search on Particle in the middle of the box of particles
    Particle refParticle = particleSystem.m_particles[60];


    SECTION("Kernel Support 2.1 * h") {
        std::vector<Particle> neighbours;

        // change the kernel support
        kernelSupport = 2.1 * h;

        // Perform the neighbor search on middle Particle
        particleSystem.iterateNeighbours(refParticle, neighbours);

        // count neighbours found
        size_t numNeighbours = neighbours.size();

        REQUIRE(numNeighbours == 13);
    }

    SECTION("Kernel Support 1.9 * h") {
        std::vector<Particle> neighbours;

        // change the kernel support
        kernelSupport = 1.9 * h;

        // Perform the neighbor search on middle Particle
        particleSystem.iterateNeighbours(refParticle, neighbours);

        // count neighbours found
        size_t numNeighbours = neighbours.size();

        REQUIRE(numNeighbours == 9);
    }

    SECTION("Kernel Support 1.1 * h") {
        std::vector<Particle> neighbours;

        // change the kernel support
        kernelSupport = 1.1 * h;

        // Perform the neighbor search on middle Particle
        particleSystem.iterateNeighbours(refParticle, neighbours);

        // count neighbours found
        size_t numNeighbours = neighbours.size();

        REQUIRE(numNeighbours == 5);
    }
}
