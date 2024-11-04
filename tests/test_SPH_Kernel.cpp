#define CATCH_CONFIG_MAIN
#include "../external/catch2/catch.hpp"
#include "../include/ParticleSystem.h"
#include "../include/Grid.h"
#include "../include/Globals.h"
#include "../include/SimulationEOS.h"



TEST_CASE("SPH Computation Test (Kernel Function)", "[Simulation]") {

    SECTION("Kernel Whole Neighborhood Ideal Sampling on all Inner Particles") {

        // place 121 particles in a box pattern (with h spacing) in a 11x11 Box
        ParticleSystem particleSystem(121);
        particleSystem.addBox(121);

        // initialize index search arrays C and L
        particleSystem.C.clear();
        for (auto& p : particleSystem.m_particles) {
            int cellIndex = particleSystem.getCellIndex(p.position);
            particleSystem.C[cellIndex]++;
        }
        particleSystem.L = particleSystem.generateSortedList(particleSystem.m_particles);

        for (auto& p : particleSystem.m_particles) {

            std::vector<Particle> neighbours;
            particleSystem.iterateNeighbours(p, neighbours);
            float kernelSum = 0.0f;

            // we only want to test inner Particles not Particle on the edge of box
            if (neighbours.size() != 9) continue;

            for (auto& particle : neighbours) {
                kernelSum += SPHComputations::kernel(p.position, particle.position);
            }

            float h_squared = std::pow(h,2);

            // error tolerance
            float tolerance = 0.01f;

            REQUIRE(std::abs(kernelSum - 1/h_squared) < tolerance);
            // REQUIRE(kernelSum == 1/h_squared);
        }
    }

    SECTION("Kernel Test particle outside kernel support") {
    // define two postions (representative for 2 particles) -> both same pos
    sf::Vector2f positionA;
    sf::Vector2f positionB;

    positionA.x = 100.0f;
    positionA.y = 100.0f;

    positionB.x = 103.0f;
    positionB.y = 103.0f;


        float W_q = SPHComputations::kernel(positionA, positionB);

        REQUIRE(W_q == 0);
    }

    SECTION("Kernel Test q == 0") {
    // define two postions (representative for 2 particles) -> both same pos
    sf::Vector2f positionA;
    sf::Vector2f positionB;

    positionA.x = 100.0f;
    positionA.y = 100.0f;

    positionB.x = 100.0f;
    positionB.y = 100.0f;


        float W_q = SPHComputations::kernel(positionA, positionB);
        W_q = (W_q*(14*M_PI*std::pow(h,2)));

        REQUIRE(W_q == 20);
    }


    SECTION("Kernel Test q == 1") {
    // define two postions (representative for 2 particles) -> with distance h on one axis
    sf::Vector2f positionA;
    sf::Vector2f positionB;

    positionA.x = 100.0f;
    positionA.y = 100.0f;

    positionB.x = 100.0f + h;
    positionB.y = 100.0f;


        float W_q = SPHComputations::kernel(positionB, positionA);
        W_q = (W_q*(14*M_PI*std::pow(h,2)));

        REQUIRE(W_q == 5);
    }

    SECTION("Kernel Test q == sqrt(2)") {
    // define two postions (representative for 2 particles) -> with distance h on both axis
    sf::Vector2f positionA;
    sf::Vector2f positionB;

    positionA.x = 100.0f;
    positionA.y = 100.0f;

    positionB.x = 100.0f + h;
    positionB.y = 100.0f + h;


        float W_q = SPHComputations::kernel(positionB, positionA);
        W_q = (W_q*(14*M_PI*std::pow(h,2)));

        REQUIRE(W_q == Approx(1.0f).epsilon(0.01)); // allow small error epsilon 1%
    }


    SECTION("Kernel Test") {

        std::vector <Particle> particles;
        int i = 0;
        // generate 9 particles with h spacing
        for (float y = 0; i < 3; y += h, i++) {
            int j = 0;
            for (float x = 0; j < 3; x += h, j++) {
                Particle p;
                p.position.x = x;
                p.position.y = y;
                particles.push_back(p);
            }
        }

        //Particle middleParticle = particles[4];
        float sum = 0.0f;
        for (auto& particle : particles) {
            float W_q = SPHComputations::kernel(particles[4].position, particle.position);
            sum += W_q;
        }

        sum = sum * std::pow(h,2);

        REQUIRE(sum == Approx(1.0).epsilon(0.01)); // allow small error epsilon 1%
    }
}   

// helper functions

float alpha() {
    return 5 / (14 * M_PI * std::pow(h, 2));
}

float beta() {
    return ((-3.0f) * std::pow((2 - sqrt(2.0f)),2));
}

TEST_CASE("SPH Computation Test (Kernel Derivative Function)", "[Simulation]") {

    SECTION("Kernel Derivative Whole Neighborhood Ideal Sampling on all Inner Particles") {

        // place 121 particles in a box pattern (with h spacing) in a 11x11 Box
        ParticleSystem particleSystem(121);
        particleSystem.addBox(121);

        // initialize index search arrays C and L
        particleSystem.C.clear();
        for (auto& p : particleSystem.m_particles) {
            int cellIndex = particleSystem.getCellIndex(p.position);
            particleSystem.C[cellIndex]++;
        }
        particleSystem.L = particleSystem.generateSortedList(particleSystem.m_particles);

        for (auto& p : particleSystem.m_particles) {

            std::vector<Particle> neighbours;
            particleSystem.iterateNeighbours(p, neighbours);
            sf::Vector2f kernelSum(0.0f, 0.0f);

            // we only want to test inner Particles not Particle on the edge of box
            if (neighbours.size() != 13) continue;

            for (auto& particle : neighbours) {
                kernelSum += SPHComputations::kernelGradient(p.position, particle.position);
            }

            REQUIRE(kernelSum == sf::Vector2f(0.0f, 0.0f));
        }
    }


    SECTION("Kernel Derivative Test") {

        // define the middle Particles with position (0,0)
        Particle middleParticle;
        middleParticle.position = sf::Vector2f(0.0f, 0.0f);

        Particle nextParticle;


        sf::Vector2f solutionVec;
        sf::Vector2f expectedVec;
        
        SECTION("Kernel Derivative grad(W((0,0) - (0,0))") {
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, middleParticle.position);

            REQUIRE(solutionVec == sf::Vector2f(0.0f, 0.0f));
        }

        SECTION("Kernel Derivative grad(W((0,0) - (h,0))") {
            nextParticle.position = sf::Vector2f(h, 0.0f);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            expectedVec = sf::Vector2f((3 * alpha() / h), 0.0f);

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }

        SECTION("Kernel Derivative grad(W((0,0) - (-h,0))") {
            nextParticle.position = sf::Vector2f(-h, 0.0f);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            // * -1 because Derivative Kernel is not symmetrical
            solutionVec.x *= -1;
            solutionVec.y *= -1;

            expectedVec = sf::Vector2f((3 * alpha() / h), 0.0f);

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }

        SECTION("Kernel Derivative grad(W((0,0) - (0,h))") {
            nextParticle.position = sf::Vector2f(0.0f, h);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            expectedVec = sf::Vector2f(0.0f, (3 * alpha() / h));

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }

        SECTION("Kernel Derivative grad(W((0,0) - (0,-h))") {
            nextParticle.position = sf::Vector2f(0.0f, -h);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            // * -1 because Derivative Kernel is not symmetrical
            solutionVec.x *= -1;
            solutionVec.y *= -1;

            expectedVec = sf::Vector2f(0.0f, (3 * alpha() / h));

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }

        SECTION("Kernel Derivative grad(W((0,0) - (h,h))") {
            nextParticle.position = sf::Vector2f(h, h);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            //expectedVec = (((-1)*(((1 / (h * sqrt(2)))) * alpha() * beta())), ((-1)*((1 / (h * sqrt(2))) * alpha() * beta()))); 
            expectedVec = alpha() * ((sf::Vector2f(0.0f, 0.0f) - sf::Vector2f(h,h))/(std::pow(h,2)*sqrt(2))) * beta();

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }

        SECTION("Kernel Derivative grad(W((0,0) - (-h,-h))") {
            nextParticle.position = sf::Vector2f(-h, -h);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            // * -1 because Derivative Kernel is not symmetrical
            solutionVec.x *= -1;
            solutionVec.y *= -1;

            //expectedVec = (((-1)*(((1 / (h * sqrt(2)))) * alpha() * beta())), ((-1)*((1 / (h * sqrt(2))) * alpha() * beta()))); 
            expectedVec = alpha() * ((sf::Vector2f(0.0f, 0.0f) - sf::Vector2f(h,h))/(std::pow(h,2)*sqrt(2))) * beta();

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }

        SECTION("Kernel Derivative grad(W((0,0) - (h,-h))") {
            nextParticle.position = sf::Vector2f(h, -h);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            //expectedVec = (((-1)*(((1 / (h * sqrt(2)))) * alpha() * beta())), ((-1)*((1 / (h * sqrt(2))) * alpha() * beta()))); 
            expectedVec = alpha() * ((sf::Vector2f(0.0f, 0.0f) - sf::Vector2f(h,-h))/(std::pow(h,2)*sqrt(2))) * beta();

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }

        SECTION("Kernel Derivative grad(W((0,0) - (-h,h))") {
            nextParticle.position = sf::Vector2f(-h, h);
            solutionVec = SPHComputations::kernelGradient(middleParticle.position, nextParticle.position);

            // * -1 because Derivative Kernel is not symmetrical
            solutionVec.x *= -1;
            solutionVec.y *= -1;

            //expectedVec = (((-1)*(((1 / (h * sqrt(2)))) * alpha() * beta())), ((-1)*((1 / (h * sqrt(2))) * alpha() * beta()))); 
            expectedVec = alpha() * ((sf::Vector2f(0.0f, 0.0f) - sf::Vector2f(h,-h))/(std::pow(h,2)*sqrt(2))) * beta();

            REQUIRE(solutionVec.x == expectedVec.x);
            REQUIRE(solutionVec.y == expectedVec.y);
        }
    }
}