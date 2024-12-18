// ParticleSystem.cpp
#include "../include/ParticleSystem.h"
#include "../include/Globals.h"
#include "../include/Grid.h"
#include "../include/SimulationEOS.h"
#include "../include/SimulationIISPH_PressureBoundaries.h"
#include "../include/UIManager.h"
#include <cstdlib>
#include <cmath>
#include <SimulationIISPH_MLSExtrapolation.h>
#include <../external/eigen-master/Eigen/Dense>

#include "ParticleSpawning.h"


ParticleSystem::ParticleSystem(unsigned int count) : m_particles(count), C(numCellsX, numCellsY), L(m_particles.size() +
                                                                                                    15 * 2 +
                                                                                                    110 * 2
                                                         // 250 * 1
                                                     ) {}

////////////////////// Simulation Part ////////////////////////////

void ParticleSystem::initiateParticles(std::vector<Particle>& particles) {

    // // write important globals to file
    writeGlobalsToFile("globals.txt");

    // addBox2(315, 524);
    // addBox(100);
    //resetSimulation();


    ///////////////////////////// Test Setting //////////////////////////////////////////

    // // bottom Boundary
    // addBoundaries(100, 420, 662); // short bottom
    // addBoundaries(200, 320, 664); // short bottom

    // left boundary
    // addBoundaries2(99, 420, 464);
    // addBoundaries2(200, 418, 264);

    // right boundary
    // addBoundaries2(99, 618, 464);
    // addBoundaries2(200, 620, 264);

    // top boundary
    // addBoundaries(100, 420, 462);

    ///////////////////////////// Half Circle with Triangle Setting /////////////////////

    // // bottom Boundary
    // addBoundaries(200, 320, 662); // short bottom
    // // addBoundaries(200, 320, 664); // short bottom
    //
    // // left boundary
    // addBoundaries2(200, 320, 264);
    // // addBoundaries2(200, 318, 264);
    //
    // // right boundary
    // addBoundaries2(200, 718, 264);
    // // addBoundaries2(200, 720, 264);
    //
    // // triangle
    // addBoundariesWithAngle(80, 520, 414, 45.0f);
    // addBoundariesWithAngle(80, 520, 414, 135.0f);
    // addBoundaries(113, 408, 527);
    //
    // // addBoundariesWithAngle(80, 520, 416, 45.0f);
    // // addBoundariesWithAngle(80, 520, 416, 135.0f);
    // // addBoundaries(113, 408, 525);
    //
    // addBox2(420, 210); // with 10000 and addbox 100


    // addBoundariesInHalfCircle({},519, 446, 199.0f, 0.0f);
    // addBoundariesInHalfCircle(400,519, 448, 199.0f, 0.0f);

    ///////////////////////////// Analyse /////////////////////

    // use 1008 particles with 12 addbox

    ParticleSpawning::addBoundaries(15, 318, 662, m_particles); // short bottom
    // addBoundaries(15, 318, 664); // short bottom
    // left boundary
    ParticleSpawning::addBoundaries2(110, 316, 444, m_particles);
    // addBoundaries2(100, 318, 464);
    // right boundary
    ParticleSpawning::addBoundaries2(110, 348, 444, m_particles);
    // addBoundaries2(100, 346, 464);
    // top boundary
    ParticleSpawning::addBoundaries(15, 318, 444, m_particles);

    ParticleSpawning::addBoxAnalyse(321, 486, 12, m_particles);



    ///////////////////////////// Rotating Circle /////////////////////

    // use 1480 Particles and 40 layer

    // addBoxAnalyse(240, 255, 40);
    //
    // // addRotatingCircle(300, 300, 100, 0.25f);
    //
    // // addRotatingCircleWithTriangles(300, 300, 100, 0.25f, 2.0f, 3, 20.0f);
    // // Center at (300, 300), radius 100, angular velocity 0.5, spacing 5, 3 triangles, height 20
    //
    // addRotatingCircleWithRectangle(
    // 300, 300, 100, 0.25f, 2.0f, 50.0f, 20.0f, 6.0f);
    // // Center: (300, 300), radius: 100, angular velocity: 0.5, spacing: 5,
    // // rectangle width: 50, height: 20, offset: 10


    //////////////////////////////////////////////////////////

    for (auto& particle : particles) {
        // Set uniform mass for all particles
        float particleRestVolume = std::pow(h, 2);
        
        float particleMass = particleRestVolume * density_rest;
        // float particleMass = 30;
        // density_rest = particleMass/particleRestVolume;


        particle.mass = particleMass;
        particle.density = 0.0f;

        particle.restVolume = particleRestVolume;
        particle.volume = 0.0f;


        //printf("Volume %f", particleVolume);
        //printf("Mass %f", particleMass);

        // count Fluid Particles
        if (!(particle.isStatic)) {
            countFluidParticles++;
        }
    }
}

// int countDir = 0; // for moveable boundary 

// Simulation Step of solver with EOS
void ParticleSystem::updateParticlesEOS(std::vector<Particle>& particles, int x_size_screen, int y_size_screen) {

    /////////// Index search
    // Calculate the total number of cells in the grid and initialize C
    C.clear();

    // Index Search
    // for neighbor search assign particle to a cell
    for (auto& p : particles) {
        int cellIndex = getCellIndex(p.position);
        C[cellIndex]++;
    }
    L = generateSortedList(particles);
    //////////////////////

    std::vector<Particle*> neighbours;

    // compute density and pressure
    for (auto& particle : particles) {
        // if (particle.isStatic) continue;
        neighbours.clear();

        iterateNeighbours(particle, neighbours);

        // quadraticNeighbourSearch(particle, neighbours);

        particle.density = SPHComputationsEOS::computeDensity(particle, neighbours);
        if (SPHComputationsEOS::isParticleCompressed(particle.density)) {
            //printf("denisty: %f", particle.density);
            particle.pressure = SPHComputationsEOS::computePressure(particle.density);
        } else {
            particle.pressure = 0;
        }
        // if (particle.isStatic == false )printf("density_i: %f", particle.density);
    }

    // compute accelerations and remove out of bounds
    int it = 0;
    for (auto& particle : particles) {
        if (particle.isStatic) continue;
        neighbours.clear();

        iterateNeighbours(particle, neighbours);

        // quadraticNeighbourSearch(particle, neighbours);

        particle.acceleration = SPHComputationsEOS::computeTotalAcceleration(particle, neighbours);

        if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
            m_particles.erase(m_particles.begin() + it);
        }
        it++;
    }

    SPHComputationsEOS::advectParticles(particles);

    SPHComputationsEOS::isCFLConditionTrue(particles);


    // Compute and write the average density to the file
    // std::string densityFile = "density_values_t_groß.txt";
    // float avgDensity = SPHComputations::getAvgDensity(particles, densityFile);

    // float cfl_number = SPHComputations::getCourantNumber(particles);
    // printf("Courant Number: %f \n", cf  l_number);

    // for moving boundary
    // if (countDir % 500 == 0) {
    //     switchDir = !switchDir;
    // }
    // countDir++;
    // moveBoundaryParticles(particles, 5.5, 5.5, switchDir);
    // rotateBoundaryParticles(particles, 320, 482, 0.9f);
}

// Simulation Step of solver with IISPH and Pressure Boundaries
void ParticleSystem::updateParticlesIISPHPressureBoundaries(std::vector<Particle>& particles, int x_size_screen, int y_size_screen) {

     /////////// Index search
     // Calculate the total number of cells in the grid and initialize C
     C.clear();

     // Index Search
     // for neighbor search assign particle to a cell
     for (auto& p : particles) {
         int cellIndex = getCellIndex(p.position);
         C[cellIndex]++;
     }
     L = generateSortedList(particles);
     ////////////////////

    for (auto& particle : particles) {
        particle.particle_neighbours.clear();
        iterateNeighbours(particle, particle.particle_neighbours);
    }

    float real_V_error_sum = 0.0f;

    // as in Pressure Boundary IISPH Paper 2018
    for (auto& particle : particles) {
        if (particle.isStatic) {
            particle.restVolume = SPHComputationsIISPH_PressureBoundaries::computeRestVolumeBoundary(particle);
            particle.volume = SPHComputationsIISPH_PressureBoundaries::computeVolumeBoundary(particle);
        }
    }

    for (auto& particle : particles) {
        if (!(particle.isStatic)) {
            particle.volume = SPHComputationsIISPH_PressureBoundaries::computeVolumeFluid(particle);
        }

        if (!(particle.isStatic)) {
            // compute real error
            // real_V_error += (particle.volume - particleRestVolume);
            // float real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
            // real_V_error_sum += real_V_error;

            float real_V_error;
            if (particle.volume > particleRestVolume) {
                real_V_error = 0.0f;
            } else {
                // real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
                real_V_error = ((particle.volume - particleRestVolume)) / particleRestVolume) * 100.0f;
                // Calculate the relative volume error
                // real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
            }
            real_V_error_sum += real_V_error;
        }
    }

    real_V_error_sum /= countFluidParticles;

    // Compute the percentage error with respect to particleRestVolume
    float percentage_V_error = std::abs(real_V_error_sum - particleRestVolume) / particleRestVolume  * 100.0f;

    for (auto& particle : particles) {
        if (particle.isStatic) continue;
        particle.predicted_velocity = SPHComputationsIISPH_PressureBoundaries::predictVelocity(particle);
    }

    for (auto& particle : particles) {
        particle.sourceTerm = SPHComputationsIISPH_PressureBoundaries::computeSourceTerm(particle);

        if (particle.isStatic) {
            particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementBoundary(particle);
        } else {
            particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(particle);
        }

        // normal pressure initialising
        // particle.pressure = 0.0f;

        particle.pressure = particle.pressure * 0.4f;

        // als test TODO: Entfernen allerding glaube ich das hier ist besser -> mach hierzu mal Iteration Analyse -> mit omega multiplizieren
        // particle.pressure = std::max(particle.sourceTerm/particle.diagonalElement, 0.0f) * omega;
    }

    int iterations_l = 0;
    float V_error = 0.0f;

    while (iterations_l < 100) {
        int it = 0;                      // For out of bounds
        V_error = 0.0f;                  // Reset volume error for this iteration

        // for each fluid sample compute pressure acc
        for (auto& particle : particles) {
            if (particle.isStatic) continue;
            particle.acceleration = SPHComputationsIISPH_PressureBoundaries::computePressureAcceleration(particle);

            // remove out of bounds particles
            if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
                m_particles.erase(m_particles.begin() + it);
            }
            it++;
        }

        for (auto& particle : particles) {
            if (particle.isStatic) {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApBoundary(particle);
            } else {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApFluid(particle);
            }

            // mit omega boundary handling
            if (particle.diagonalElement != 0.0f) {
                if (particle.isStatic) {
                    particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressureBoundaries(particle);
                } else {
                    particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressure(particle);
                }
            }

            // ohne spezielle update pressure omega handling
            // particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressure(particle);

            // volume error fluid and boundary
            float particle_volume_error = 100.0f * std::max(particle.Ap - particle.sourceTerm, 0.0f) / particleRestVolume;
            V_error += particle_volume_error;  // Accumulate total volume error

            if (std::isnan(particle.Ap)) {
                printf("AP inf");
            }

            if (std::isinf(particle.pressure)) {
                printf("Pressure inf");
            }
        }

        // Average the volume and density errors across particles
        V_error /= particles.size();

        if (V_error < 0.01f && iterations_l > 1) {
            currentIterations = iterations_l;
            // Compute and write the curren Iterattion to the file TODO file write section
            // std::string iterationFile = "currentIterations_maxPressureInit.txt";
            // SPHComputationsIISPH_PressureBoundaries::writeCurrentIterationToFile(iterationFile);

            break;
        }

        iterations_l++;
    }

    currentVolumeError = real_V_error_sum;
    // std::cout << "\rreal Error: " << real_V_error_sum << "% Time Step: " << timeStep << std::flush;
    // std::cout << "\rreal Error: " << percentage_V_error << std::flush;

    // printf("real Error: %f \n", real_V_error);
    // printf("next \n \n");

    // TODO write volume error txt
    // Compute and write the curren Iterattion to the file
    // std::string volumeErrorFile = "currentVolumeError_PressureBoundaries.txt";
    // SPHComputationsIISPH_PressureBoundaries::writeCurrentIterationToFile(volumeErrorFile);

     SPHComputationsIISPH_PressureBoundaries::advectParticles(particles);

     SPHComputationsIISPH_PressureBoundaries::isCFLConditionTrue(particles);

    // printf("Courant Number: %f \n", SPHComputationsIISPH::getCourantNumber(particles));
    // std::cout << "\rCourant Number: %f" << SPHComputationsIISPH::getCourantNumber(particles) << std::flush;

    // all for circle rotating
    // moveRotatingCircle(300, 300, 100, 0.25f, timeStep);
    // moveRotatingCircleWithTriangles(300, 300, 100, 0.25f, timeStep);
    // moveRotatingCircleWithRectangle(
    // 300, 300, 100, 0.5f, timeStep, 5.0f);

}

// Simulation Step of solver with IISPH and SPH Pressure Extrapolation (eq 3 for pressure at boundary)
void ParticleSystem::updateParticlesIISPH_Extrapolation(std::vector<Particle>& particles, int x_size_screen, int y_size_screen) {

     /////////// Index search
     // Calculate the total number of cells in the grid and initialize C
     C.clear();

     // Index Search
     // for neighbor search assign particle to a cell
     for (auto& p : particles) {
         int cellIndex = getCellIndex(p.position);
         C[cellIndex]++;
     }
     L = generateSortedList(particles);
     ////////////////////

    for (auto& particle : particles) {
        particle.particle_neighbours.clear();
        iterateNeighbours(particle, particle.particle_neighbours);
    }

    float real_V_error_sum = 0.0f;

    for (auto& particle : particles) {

        particle.volume = SPHComputationsIISPH_PressureBoundaries::computeVolumeFluid(particle);

        if (!(particle.isStatic)) {
            // compute real error
            // real_V_error += (particle.volume - particleRestVolume);
            float real_V_error;
            if (particle.volume > particleRestVolume) {
                real_V_error = 0.0f;
            } else {
                real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
            }
            real_V_error_sum += real_V_error;
        }
    }

    real_V_error_sum /= countFluidParticles;

    // Compute the percentage error with respect to particleRestVolume
    float percentage_V_error = std::abs(real_V_error_sum - particleRestVolume) / particleRestVolume  * 100.0f;

    for (auto& particle : particles) {
        if (particle.isStatic) continue;
        particle.predicted_velocity = SPHComputationsIISPH_PressureBoundaries::predictVelocity(particle);
    }

    for (auto& particle : particles) {
        particle.sourceTerm = SPHComputationsIISPH_PressureBoundaries::computeSourceTerm(particle);

        if (particle.isStatic) {
            particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementBoundary(particle);
        } else {
            particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(particle);
        }

        // normale pressure initilisierung
        // particle.pressure = 0.0f;

        // als test TODO: Entfernen allerding glaube ich das hier ist besser
        particle.pressure = std::max(particle.sourceTerm/particle.diagonalElement, 0.0f);

        // if (particle.isStatic) {
        //     particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressureBoundariesExtrapolation(particle);
        // }
    }

    int iterations_l = 0;
    float V_error = 0.0f;

    while (iterations_l < 100) {
        int it = 0;                      // For out of bounds
        V_error = 0.0f;                  // Reset volume error for this iteration

        // for each fluid sample compute pressure acc
        for (auto& particle : particles) {
            if (particle.isStatic) continue;
            particle.acceleration = SPHComputationsIISPH_PressureBoundaries::computePressureAccelerationExtrapolation(particle);

            // remove out of bounds particles
            if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
                m_particles.erase(m_particles.begin() + it);
            }
            it++;
        }

        for (auto& particle : particles) {
            if (std::abs(particle.diagonalElement) != 0.0f) {
                if (particle.isStatic) {
                    particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressureBoundariesExtrapolation(particle);
                }
            }
        }

        for (auto& particle : particles) {
            if (particle.isStatic) {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApBoundary(particle);
            } else {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApFluid(particle);
            }

            // mit extrapolation boundary handling
            if (std::abs(particle.diagonalElement) > 1e-6f) {
                if (particle.isStatic) {
                    // particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressureBoundariesExtrapolation(particle);
                } else {
                    particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressure(particle);
                }
            }

            // ohne extrapolation handling
            // particle.pressure = SPHComputationsIISPH::updatePressure(particle);

            // volume error fluid and boundary
            float particle_volume_error = 100.0f * std::max(particle.Ap - particle.sourceTerm, 0.0f) / particleRestVolume;
            V_error += particle_volume_error;  // Accumulate total volume error

            if (std::isnan(particle.Ap)) {
                printf("AP inf");
            }

            if (std::isinf(particle.pressure)) {
                printf("Pressure inf");
            }
        }

        // Average the volume and density errors across particles
        V_error /= particles.size();

        if (V_error < 0.1f && iterations_l > 1) {
            currentIterations = iterations_l;
            // printf("V_error: %f\n", V_error);
            // printf("iterations needed: %d \n", iterations_l);

            // Compute and write the curren Iterattion to the file
            // std::string iterationFile = "currentIterations_1500.txt";
            // SPHComputationsIISPH_PressureBoundaries::writeCurrentIterationToFile(iterationFile);
            break;
        }

        iterations_l++;
    }

    std::cout << "\rreal Error: " << real_V_error_sum << "% Time Step: " << timeStep << std::flush;

     SPHComputationsIISPH_PressureBoundaries::advectParticles(particles);

     SPHComputationsIISPH_PressureBoundaries::isCFLConditionTrue(particles);

    ParticleSpawning::moveRotatingCircle(300, 300, 100, 0.25f, timeStep, m_particles);
}

// Simulation Step of solver with IISPH and Mirroring for Pressure at Boundaries
void ParticleSystem::updateParticlesIISPH_Mirroring(std::vector<Particle>& particles, int x_size_screen, int y_size_screen) {

     /////////// Index search
     // Calculate the total number of cells in the grid and initialize C
     C.clear();

     // Index Search
     // for neighbor search assign particle to a cell
     for (auto& p : particles) {
         int cellIndex = getCellIndex(p.position);
         C[cellIndex]++;
     }
     L = generateSortedList(particles);
     ////////////////////

    for (auto& particle : particles) {
        particle.particle_neighbours.clear();
        iterateNeighbours(particle, particle.particle_neighbours);
    }

    float real_V_error_sum = 0.0f;

    for (auto& particle : particles) {
        particle.volume = SPHComputationsIISPH_PressureBoundaries::computeVolumeFluid(particle);

        if (!(particle.isStatic)) {
            float real_V_error;
            if (particle.volume > particleRestVolume) {
                real_V_error = 0.0f;
            } else {
                real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
            }
            real_V_error_sum += real_V_error;
        }
    }

    real_V_error_sum /= countFluidParticles;

    // Compute the percentage error with respect to particleRestVolume
    float percentage_V_error = std::abs(real_V_error_sum - particleRestVolume) / particleRestVolume  * 100.0f;

    for (auto& particle : particles) {
        if (particle.isStatic) continue;
        particle.predicted_velocity = SPHComputationsIISPH_PressureBoundaries::predictVelocity(particle);
    }

    for (auto& particle : particles) {
        particle.sourceTerm = SPHComputationsIISPH_PressureBoundaries::computeSourceTerm(particle);

        if (particle.isStatic) {
            particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementBoundary(particle);
        } else {
            particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(particle);
        }

        if (!(particle.isStatic)) {
            // normale pressure initilisierung
            // particle.pressure = 0.0f;

            particle.pressure = particle.pressure * 0.4f;

            // als test TODO: Entfernen allerding glaube ich das hier ist besser
            // particle.pressure = std::max(particle.sourceTerm/particle.diagonalElement, 0.0f);
        }
        if (particle.isStatic) particle.pressure = 0.0f;
    }

    int iterations_l = 0;
    float V_error = 0.0f;

    while (iterations_l < 100) {
        int it = 0;                      // For out of bounds
        V_error = 0.0f;                  // Reset volume error for this iteration

        // for each fluid sample compute pressure acc
        for (auto& particle : particles) {
            if (particle.isStatic) continue;
            particle.acceleration = SPHComputationsIISPH_PressureBoundaries::computePressureAccelerationMirror(particle);

            // remove out of bounds particles
            if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
                m_particles.erase(m_particles.begin() + it);
            }
            it++;
        }

        for (auto& particle : particles) {
            if (particle.isStatic) {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApBoundary(particle);
            } else {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApFluid(particle);
            }

            // mit mirroring boundary handling (no pressure calculation at boundary particles)
            if (particle.diagonalElement != 0.0f) {
                if (particle.isStatic) continue;
                particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressure(particle);
            }

            // ohne extrapolation handling
            // particle.pressure = SPHComputationsIISPH::updatePressure(particle);

            // volume error fluid and boundary
            float particle_volume_error = 100.0f * std::max(particle.Ap - particle.sourceTerm, 0.0f) / particleRestVolume;
            V_error += particle_volume_error;  // Accumulate total volume error

            if (std::isnan(particle.Ap)) {
                printf("AP inf");
            }

            if (std::isinf(particle.pressure)) {
                printf("Pressure inf");
            }
        }

        // Average the volume and density errors across particles
        V_error /= particles.size();

        if (V_error < 0.01f && iterations_l > 1) {
            currentIterations = iterations_l;
            break;
        }

        iterations_l++;
    }

    currentVolumeError = real_V_error_sum;
    // std::cout << "\rreal Error: " << real_V_error_sum << "% Time Step: " << timeStep << std::flush;
    // std::cout << "\rreal Error: " << percentage_V_error << std::flush;

    // printf("real Error: %f \n", real_V_error);
    // printf("next \n \n");

    // TODO write volume error txt
    // Compute and write the curren Iterattion to the file
    // std::string volumeErrorFile = "currentVolumeError_Mirror_Gamma_1.0.txt";
    // SPHComputationsIISPH_PressureBoundaries::writeCurrentIterationToFile(volumeErrorFile);

     SPHComputationsIISPH_PressureBoundaries::advectParticles(particles);

     SPHComputationsIISPH_PressureBoundaries::isCFLConditionTrue(particles);

    // printf("Courant Number: %f \n", SPHComputationsIISPH::getCourantNumber(particles));
    // std::cout << "\rCourant Number: %f" << SPHComputationsIISPH::getCourantNumber(particles) << std::flush;

    // moveRotatingCircle(300, 300, 100, 0.25f, timeStep);

    // all for circle rotatimg
    // moveRotatingCircleWithTriangles(300, 300, 100, 0.25f, timeStep);
    // moveRotatingCircleWithRectangle(
    // 300, 300, 100, 0.25f, timeStep, 5.0f);
}

// Simulation Step of solver with IISPH and Mirroring for Pressure at Boundaries
void ParticleSystem::updateParticlesIISPH_MLSExtrapolation(std::vector<Particle>& particles, int x_size_screen, int y_size_screen) {

     /////////// Index search
     // Calculate the total number of cells in the grid and initialize C
     C.clear();

     // Index Search
     // for neighbor search assign particle to a cell
     for (auto& p : particles) {
         int cellIndex = getCellIndex(p.position);
         C[cellIndex]++;
     }
     L = generateSortedList(particles);
     ////////////////////

    for (auto& particle : particles) {
        particle.particle_neighbours.clear();
        iterateNeighbours(particle, particle.particle_neighbours);
    }

    float real_V_error_sum = 0.0f;

    for (auto& particle : particles) {
        particle.volume = SPHComputationsIISPH_PressureBoundaries::computeVolumeFluid(particle);

        particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(particle);

        if (!(particle.isStatic)) {
            // compute real error
            // real_V_error += (particle.volume - particleRestVolume);
            // float real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
            // real_V_error_sum += real_V_error;

            float real_V_error;
            if (particle.volume > particleRestVolume) {
                real_V_error = 0.0f;
            } else {
                real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
            }
            real_V_error_sum += real_V_error;
        }
    }

    real_V_error_sum /= countFluidParticles;

    // Compute the percentage error with respect to particleRestVolume
    float percentage_V_error = std::abs(real_V_error_sum - particleRestVolume) / particleRestVolume  * 100.0f;

    for (auto& particle : particles) {
        if (particle.isStatic) continue;
        particle.predicted_velocity = SPHComputationsIISPH_PressureBoundaries::predictVelocity(particle);
    }

    for (auto& particle : particles) {
        particle.sourceTerm = SPHComputationsIISPH_PressureBoundaries::computeSourceTerm(particle);

        if (particle.isStatic) {
            // particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementBoundary(particle);
        } else {
            // ich mach das schon oben boundary braucht keins
            // particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(particle);
        }

        // normale pressure initilisierung
        particle.pressure = 0.0f;

        // als test TODO: Entfernen allerding glaube ich das hier ist besser
        // particle.pressure = std::max(particle.sourceTerm/particle.diagonalElement, 0.0f);
    }

    int iterations_l = 0;
    float V_error = 0.0f;

    while (iterations_l < 100) {
        int it = 0;                      // For out of bounds
        V_error = 0.0f;                  // Reset volume error for this iteration

        for (auto& particle : particles) {
            if (std::abs(particle.diagonalElement) > 1e-6f) {
                if (particle.isStatic) {
                    particle.pressure = SPHComputationsIISPH_MLSExtra::computePressureMLS(particle);
                }
            }
            // particle.pressure = SPHComputationsIISPH_MLSExtra::computePressureMLS(particle);
        }

        // for each fluid sample compute pressure acc
        for (auto& particle : particles) {
            if (particle.isStatic) continue;
            particle.acceleration = SPHComputationsIISPH_MLSExtra::computePressureAcceleration(particle);

            // remove out of bounds particles
            if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
                m_particles.erase(m_particles.begin() + it);
            }
            it++;
        }

        for (auto& particle : particles) {
            if (particle.isStatic) {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApBoundary(particle);
            } else {
                particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApFluid(particle);
            }

            // mit MLS boundary handling
            if (std::abs(particle.diagonalElement) > 1e-6f) {
                if (particle.isStatic) {
                    // particle.pressure = SPHComputationsIISPH_MLSExtra::computePressureMLS(particle);
                } else {
                    particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressure(particle);
                    // particle.pressure = SPHComputationsIISPH_MLSExtra::updatePressureMLSFluid(particle);
                }
            }

            // volume error fluid and boundary
            float particle_volume_error = 100.0f * std::max(particle.Ap - particle.sourceTerm, 0.0f) / particleRestVolume;
            V_error += particle_volume_error;  // Accumulate total volume error

            if (std::isnan(particle.Ap)) {
                printf("AP inf");
            }

            if (std::isinf(particle.pressure)) {
                printf("Pressure inf");
            }

            if (std::isnan(particle.pressure)) {
                printf("Pressure NAN");
            }
        }

        // Average the volume and density errors across particles
        V_error /= particles.size();

        if (V_error < 0.1f && iterations_l > 1) {
            currentIterations = iterations_l;
            // printf("V_error: %f\n", V_error);
            // printf("iterations needed: %d \n", iterations_l);

            // Compute and write the curren Iterattion to the file
            // std::string iterationFile = "currentIterations_1500.txt";
            // SPHComputationsIISPH_PressureBoundaries::writeCurrentIterationToFile(iterationFile);
            break;
        }

        iterations_l++;
    }

    std::cout << "\rreal Error: " << real_V_error_sum << "% Time Step: " << timeStep << std::flush;
    // std::cout << "\rreal Error: " << percentage_V_error << std::flush;

    // printf("real Error: %f \n", real_V_error);
    // printf("next \n \n");

     SPHComputationsIISPH_PressureBoundaries::advectParticles(particles);

     SPHComputationsIISPH_PressureBoundaries::isCFLConditionTrue(particles);

    // printf("Courant Number: %f \n", SPHComputationsIISPH::getCourantNumber(particles));
    // std::cout << "\rCourant Number: %f" << SPHComputationsIISPH::getCourantNumber(particles) << std::flush;

    ParticleSpawning::moveRotatingCircle(300, 300, 100, 0.25f, timeStep, m_particles);
}

// // Simulation Step of solver with IISPH and Mirroring for Pressure at Boundaries -> new version but does not work!
// void ParticleSystem::updateParticlesIISPH_MLSExtrapolation(std::vector<Particle>& particles, int x_size_screen, int y_size_screen) {
//
//      /////////// Index search
//      // Calculate the total number of cells in the grid and initialize C
//      C.clear();
//
//      // Index Search
//      // for neighbor search assign particle to a cell
//      for (auto& p : particles) {
//          int cellIndex = getCellIndex(p.position);
//          C[cellIndex]++;
//      }
//      L = generateSortedList(particles);
//      ////////////////////
//
//     for (auto& particle : particles) {
//         particle.particle_neighbours.clear();
//         iterateNeighbours(particle, particle.particle_neighbours);
//     }
//
//     float real_V_error_sum = 0.0f;
//
//     for (auto& particle : particles) {
//         particle.volume = SPHComputationsIISPH_PressureBoundaries::computeVolumeFluid(particle);
//
//         particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(particle);
//
//         if (!(particle.isStatic)) {
//             // compute real error
//             // real_V_error += (particle.volume - particleRestVolume);
//             float real_V_error = std::abs(particle.volume - particleRestVolume) / particleRestVolume * 100.0f;
//             real_V_error_sum += real_V_error;
//         }
//     }
//
//     real_V_error_sum /= countFluidParticles;
//
//     // Compute the percentage error with respect to particleRestVolume
//     float percentage_V_error = std::abs(real_V_error_sum - particleRestVolume) / particleRestVolume  * 100.0f;
//
//     for (auto& particle : particles) {
//         if (particle.isStatic) continue;
//         particle.predicted_velocity = SPHComputationsIISPH_PressureBoundaries::predictVelocity(particle);
//     }
//
//     for (auto& particle : particles) {
//         particle.sourceTerm = SPHComputationsIISPH_MLSExtra::computeSourceTermDensityError(particle);
//
//
//         // delete
//         // if (particle.isStatic) {
//         //     particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementBoundary(particle);
//         // } else {
//         //     particle.diagonalElement = SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(particle);
//         // }
//
//         // normale pressure initilisierung
//         particle.pressure = 0.0f;
//
//         // als test TODO: Entfernen allerding glaube ich das hier ist besser
//         // particle.pressure = std::max(particle.sourceTerm/particle.diagonalElement, 0.0f);
//
//     }
//
//     int iterations_l = 0;
//     float V_error = 0.0f;
//
//     while (iterations_l < 100) {
//         int it = 0;                      // For out of bounds
//         V_error = 0.0f;                  // Reset volume error for this iteration
//
//         for (auto& particle : particles) {
//             if (particle.diagonalElement != 0.0f) {
//                 if (particle.isStatic) {
//                     particle.pressure = SPHComputationsIISPH_MLSExtra::computePressureMLS(particle);
//                 }
//             }
//
//             // remove out of bounds particles
//             if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
//                 m_particles.erase(m_particles.begin() + it);
//             }
//             it++;
//         }
//
//         for (auto& particle : particles) {
//             if (particle.isStatic) continue;
//             particle.acceleration = SPHComputationsIISPH_MLSExtra::computePressureAcceleration(particle);
//         }
//
//         for (auto& particle : particles) {
//             if (particle.diagonalElement != 0.0f) {
//                 if (particle.isStatic) continue;
//                 particle.pressure = SPHComputationsIISPH_MLSExtra::updatePressureMLSFluid(particle);
//                 }
//             // volume error fluid and boundary
//             float particle_volume_error = 100.0f * std::max(particle.pressure - particle.sourceTerm, 0.0f) / particleRestVolume;
//             V_error += particle_volume_error;  // Accumulate total volume error
//         }
//
//
//
//         // // for each fluid sample compute pressure acc
//         // for (auto& particle : particles) {
//         //     if (particle.isStatic) continue;
//         //     particle.acceleration = SPHComputationsIISPH_PressureBoundaries::computePressureAcceleration(particle);
//         //
//         //     // remove out of bounds particles
//         //     if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
//         //         m_particles.erase(m_particles.begin() + it);
//         //     }
//         //     it++;
//         // }
//         //
//         // for (auto& particle : particles) {
//         //     if (particle.isStatic) {
//         //         particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApBoundary(particle);
//         //     } else {
//         //         particle.Ap = SPHComputationsIISPH_PressureBoundaries::computeApFluid(particle);
//         //     }
//         //
//         //     // mit mirroring boundary handling (no pressure calculation at boundary particles)
//         //     if (particle.diagonalElement != 0.0f) {
//         //         if (particle.isStatic) {
//         //             particle.pressure = SPHComputationsIISPH_MLSExtra::computePressureMLS(particle);
//         //
//         //         } else {
//         //             particle.pressure = SPHComputationsIISPH_PressureBoundaries::updatePressure(particle);
//         //         }
//         //     }
//         //
//         //     // volume error fluid and boundary
//         //     float particle_volume_error = 100.0f * std::max(particle.Ap - particle.sourceTerm, 0.0f) / particleRestVolume;
//         //     V_error += particle_volume_error;  // Accumulate total volume error
//
//         //     if (std::isnan(particle.Ap)) {
//         //         printf("AP inf");
//         //     }
//         //
//         //     if (std::isinf(particle.pressure)) {
//         //         printf("Pressure inf");
//         //     }
//         // }
//
//
//         // Average the volume and density errors across particles
//         V_error /= particles.size();
//
//         if (V_error < 0.1f && iterations_l > 1) {
//             currentIterations = iterations_l;
//             // printf("V_error: %f\n", V_error);
//             // printf("iterations needed: %d \n", iterations_l);
//
//             // Compute and write the curren Iterattion to the file
//             // std::string iterationFile = "currentIterations_1500.txt";
//             // SPHComputationsIISPH_PressureBoundaries::writeCurrentIterationToFile(iterationFile);
//             break;
//         }
//
//         iterations_l++;
//     }
//
//     std::cout << "\rreal Error: " << real_V_error_sum << "% Time Step: " << timeStep << std::flush;
//     // std::cout << "\rreal Error: " << percentage_V_error << std::flush;
//
//     // printf("real Error: %f \n", real_V_error);
//     // printf("next \n \n");
//
//      SPHComputationsIISPH_PressureBoundaries::advectParticles(particles);
//
//      SPHComputationsIISPH_PressureBoundaries::isCFLConditionTrue(particles);
//
//     // printf("Courant Number: %f \n", SPHComputationsIISPH::getCourantNumber(particles));
//     // std::cout << "\rCourant Number: %f" << SPHComputationsIISPH::getCourantNumber(particles) << std::flush;
// }







////////////////////////////////////////////////////////////////////////////////////

void ParticleSystem::resetSimulation() {
    // addBox2(439, 584);
    ParticleSpawning::addBoxAnalyse(260, 315, 12, m_particles);
}

//////////////////////// Neighbour Search ///////////////////////

// quadratic search for test
void ParticleSystem::quadraticNeighbourSearch(Particle &particle, std::vector<Particle> &m_neighbours) {
    for (auto &neighbour: m_particles) {
        if (inRadius(particle.position, neighbour.position)) {
            m_neighbours.push_back(neighbour);
        }
    }
}

// helper function to get Cell Index for Particles (get 1D Index from particles coords)
int ParticleSystem::getCellIndex(const sf::Vector2f& position) const {
    int xIndex = static_cast<int>(position.x) / cellSize;
    int yIndex = static_cast<int>(position.y) / cellSize;
    return xIndex + yIndex * numCellsX; // Assuming cell indices are row-major
}

// get Cell 1D Index from (x, y) grid coords
int ParticleSystem::getCellCoords(int x, int y) {
    return y * numCellsY + x;
}

// helper function accumulate counter for C array
void ParticleSystem::accumulateCounters(){
    int accCounter = 0;
    for (int y = 0; y < numCellsY; ++y) {
        for (int x = 0; x < numCellsX; ++x) {
            accCounter += C.get(x, y);
            C.set(x, y, accCounter);
        }
    }
}

// get sorted List of particles L (as in lecture)
std::vector<Particle*> ParticleSystem::generateSortedList(std::vector<Particle>& particles) {

    // accum C
    accumulateCounters();

    for (int i = 0; i < m_particles.size(); ++i) {
        int x = m_particles[i].position.x;
        int y =  m_particles[i].position.y;
        int cellIndex = getCellIndex(m_particles[i].position);
        C[cellIndex]--; // -1 an der Stelle für C
        L[C[cellIndex]] = &m_particles[i];
    }
    return L;
}

void ParticleSystem::iterateNeighbours(Particle& particle, std::vector<Particle*>& m_neighbours) {
    int currCell = getCellIndex(particle.position);

    std::pair<int, int> currCellCoord = C.getCoordinatesFromIndex(currCell);

    // Define the neighbor offsets 3x3 grid
    std::vector<std::pair<int, int>> neighborOffsets = {
        { -1, -1 }, { 0, -1 }, { 1, -1 },
        { -1, 0 }, { 0, 0 }, { 1, 0 },
        { -1, 1 }, { 0, 1 }, { 1, 1 }
    };
    
    for (const auto& offset : neighborOffsets) {
        std::pair<int, int> neighborCell = { currCellCoord.first + offset.first, currCellCoord.second + offset.second };

        // Ensure neighborCell is within bounds
        if (C.isValid(neighborCell.first, neighborCell.second)) {
            int neighborCellIndex = neighborCell.second * numCellsX + neighborCell.first;
            int numParticlesInCell = C[neighborCellIndex + 1] - C[neighborCellIndex];
            for (int i = 0; i < numParticlesInCell; ++i) {
                // Access the particles in this cell
                int particleIndex = C[neighborCellIndex] + i;

                //L[particleIndex]->color = sf::Color(140, 0, 255, 255); // Purple color for neighbors visual test

                // get real neighbours in radius of particle (for kernel)
                if (inRadius(particle.position, L[particleIndex]->position)) {

                    //L[particleIndex]->color = sf::Color(0, 255, 255, 255); // visual test

                    m_neighbours.push_back(L[particleIndex]); // add to neighbour list
                }
            }
        }
    }
}

// helperfunction to get particles inside a specified radius (in Abhängigkeit von global h) uses squared distance 
// (radius = kernelSupport)
bool ParticleSystem::inRadius(sf::Vector2f positionA, sf::Vector2f positionB, float radius) {
    float dx = positionA.x - positionB.x;
    float dy = positionA.y - positionB.y;
    float distance = dx * dx + dy * dy;

    return distance < (radius * radius);
}

///////////////////////////////////////////////////////////////////////