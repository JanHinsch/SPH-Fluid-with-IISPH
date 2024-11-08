// SimulationIISPH.cpp
#include "../include/SimulationEOS.h"
#include "../include/SimulationIISPH.h"
#include "../include/ParticleSystem.h"
#include "../include/Globals.h"
#include <cmath>
#include <iostream>
#include <fstream>


float SPHComputationsIISPH::kernelAlpha() {
    return 5 / (14 * M_PI * std::pow(h, 2));
}

float SPHComputationsIISPH::kernel(sf::Vector2f positionA,sf::Vector2f positionB) {
    float distance = euclideanDistance(positionA, positionB);
    float q = distance / h;

    float result = 0.0f;
    if (q < 0.0f) {
        throw std::runtime_error("Invalid distance");
    } else if (q < 1.0f) {
        result = kernelAlpha() * (std::pow((2 - q),3) - 4 * std::pow((1 - q), 3));
    } else if (q < 2.0f) {
        result = kernelAlpha() * std::pow((2 - q),3);
    } else {
        result = 0;
    }
    //printf("Kernel Value for distance %f: %f and q = %f \n", distance, result, q);
    return result;
}

sf::Vector2f SPHComputationsIISPH::kernelGradient(sf::Vector2f positionA,sf::Vector2f positionB) {
    float distance = euclideanDistance(positionA, positionB);
    float q = distance / h;
    sf::Vector2f diff = positionA - positionB;

    // if difference == 0 return (0,0)
    if (diff == sf::Vector2f(0.0f, 0.0f)) {
        return sf::Vector2f(0.0f, 0.0f);
    }

    if (q < 0.0f) {
        throw std::runtime_error("Invalid distance");
    } else if (q < 1.0f) {
        return (kernelAlpha() * (diff / (distance * h))) * (-3 * std::pow((2 - q),2) + 12 * std::pow((1- q), 2));
    }
        else if (q < 2.0f) {
        return (kernelAlpha() * (diff / (distance * h))) * (-3 * std::pow((2 - q),2));
    } else {
        return sf::Vector2f(0.0f, 0.0f);
    }
}

float SPHComputationsIISPH::computeDensity(Particle &particle_i, std::vector<Particle*> &neighbours) {
    float density_i = 0.0f;
    for (auto& particle_j : neighbours) {
        density_i += particle_j->mass * kernel(particle_i.position, particle_j->position);
    }
    return density_i;
}

sf::Vector2f SPHComputationsIISPH::predictVelocity(Particle& particle_i) {
    sf::Vector2f non_pressure_acc = computeViscosity(particle_i)
                                    + gravity;
    non_pressure_acc *= timeStep;
    return particle_i.velocity + non_pressure_acc;
}

float SPHComputationsIISPH::computeSourceTerm(Particle& particle_i) {
    float source_term = density_rest - particle_i.density;
    float fluidSum = 0.0f;
    float boundarySum = 0.0f;

    for (auto& particle_j : particle_i.particle_neighbours) {
        // if boundaries should move as well with flow change here the 0 term
        if (particle_j->isStatic) {
            boundarySum += dotProduct((particle_i.predicted_velocity - 0), kernelGradient(particle_i.position, particle_j->position))
                            * particle_j->mass;
        } else {
            fluidSum += dotProduct((particle_i.predicted_velocity - particle_j->predicted_velocity),kernelGradient(particle_i.position, particle_j->position))
                            * particle_j->mass;
        }
    }
    fluidSum *= timeStep;
    boundarySum *= timeStep;
    return source_term - fluidSum - boundarySum;
}

// test version
// float SPHComputationsIISPH::computeSourceTerm(Particle& particle_i) {
//     float source_term = density_rest - particle_i.density;
//
//     for (auto& particle_j : particle_i.particle_neighbours) {
//         sf::Vector2f velocity_difference(0.0f, 0.0f);
//
//         if (particle_j->isStatic) {
//             // Boundary particles: velocity difference with stationary velocity (0)
//             velocity_difference = particle_i.predicted_velocity;
//         } else {
//             // Fluid particles: velocity difference with neighbor's predicted velocity
//             velocity_difference = particle_i.predicted_velocity - particle_j->predicted_velocity;
//         }
//
//         // Kernel gradient between particle i and neighbor j
//         sf::Vector2f grad = kernelGradient(particle_i.position, particle_j->position);
//
//         // Dot product between velocity difference and gradient
//         float dot_product = velocity_difference.x * grad.x + velocity_difference.y * grad.y;
//
//         // Update source term
//         source_term -= timeStep * particle_j->mass * dot_product;
//     }
//
//     return source_term;
// }

// as in notesOnIISPH
float SPHComputationsIISPH::computeDiagonalElement(Particle& particle_i) {
    float firstLine = 0.0f;
    float secondLine = 0.0f;
    float thirdLine = 0.0f;
    sf::Vector2f innerSum1(0.0f, 0.0f);
    sf::Vector2f innerSum2(0.0f, 0.0f);
    sf::Vector2f innerSum3(0.0f, 0.0f);
    float timeStepSq = timeStep * timeStep;
    float densityRestSq = density_rest * density_rest;

    for (auto& particle_j : particle_i.particle_neighbours) {

        // first term
        if (particle_j->isStatic == false) {
            for (auto& particle_j_j : particle_j->particle_neighbours) {
                if (particle_j_j->isStatic == false) {
                    innerSum1 -= (particle_j_j->mass / densityRestSq) * kernelGradient(particle_i.position, particle_j_j->position);
                } else {
                    innerSum1 -= 2 * gamma_1 *(particle_j_j->mass / densityRestSq) * kernelGradient(particle_i.position, particle_j_j->position);
                }
            }
            sf::Vector2f kernelG1 = kernelGradient(particle_i.position, particle_j->position);
            float dot_product1 = innerSum1.x * kernelG1.x + innerSum1.y * kernelG1.y;
            firstLine += particle_j->mass * dot_product1 * timeStepSq;

            // second term
            innerSum2 += particle_j->mass * (particle_i.mass / densityRestSq) * kernelGradient(particle_j->position, particle_i.position);
            sf::Vector2f kernelG2 = kernelGradient(particle_i.position, particle_j->position);
            float dotProduct2 = innerSum2.x * kernelG2.x + innerSum2.x * kernelG2.y;
            secondLine += particle_j->mass * dotProduct2 * densityRestSq;

        } else {

            // third term
            for (auto& particle_j_j : particle_j->particle_neighbours) {
                if (particle_j_j->isStatic == false) {
                    innerSum3 -= (particle_j_j->mass / densityRestSq) * kernelGradient(particle_i.position, particle_j_j->position);
                } else {
                    innerSum3 -= 2 * gamma_1 *(particle_j_j->mass / densityRestSq) * kernelGradient(particle_i.position, particle_j_j->position);
                }
            }
            sf::Vector2f kernelG3 = kernelGradient(particle_i.position, particle_j->position);
            float dot_product3 = innerSum3.x * kernelG3.x + innerSum3.y * kernelG3.y;
            thirdLine += particle_j->mass * dot_product3 * timeStepSq;
        }
    }
    return firstLine + secondLine + thirdLine;
}

sf::Vector2f SPHComputationsIISPH::computePressureAcceleration(Particle &particle_i) {
    sf::Vector2f pressure_acc(0.0f, 0.0f);
    float densityRestSq = density_rest * density_rest;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic == false) {
            pressure_acc -= particle_j->mass * (particle_i.pressure/densityRestSq) + (particle_j->pressure/densityRestSq)
                            * kernelGradient(particle_i.position, particle_j->position);
        } else {
            pressure_acc -= gamma_1 * particle_j->mass * (particle_i.pressure/densityRestSq)
                            * kernelGradient(particle_i.position, particle_j->position);
        }
    }
    return pressure_acc;
}

// compute the divergence of the velocity change due to pressure acc (Laplacian)
float SPHComputationsIISPH::computeDivergence(Particle &particle_i) {
    float Ap= 0.0f;
    float timeStepSq = timeStep * timeStep;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sf::Vector2f kernelG1 = kernelGradient(particle_i.position, particle_j->position);
            float dotProduct1 =  particle_j->acceleration.x * kernelG1.x + particle_j->acceleration.y * kernelG1.y;
            Ap += particle_j->mass * timeStepSq * dotProduct1;
        } else {
            sf::Vector2f kernelG2 = kernelGradient(particle_i.position, particle_j->position);
            sf::Vector2f innerSub = particle_i.acceleration - particle_j->acceleration;
            float dotProduct2 = innerSub.x * kernelG2.x + innerSub.y * kernelG2.y;
            Ap += particle_j->mass * timeStepSq * dotProduct2;
        }
    }
    return Ap;
}

float SPHComputationsIISPH::updatePressure(Particle &particle_i) {
    float pressure = omega * ((particle_i.sourceTerm - particle_i.laplacian) / particle_i.diagonalElement) + particle_i.
                     pressure;
    if (pressure < 0) {
        return 0.0f;
    } else {
        return pressure;
    }
}


// fixed version
sf::Vector2f SPHComputationsIISPH::computeViscosity(Particle& particle_i) {
    sf::Vector2f sumFluid(0.0f, 0.0f);
    sf::Vector2f sumBoundary(0.0f, 0.0f);
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic == false) {
            sf::Vector2f v_diff = particle_i.velocity - particle_j->velocity;
            sf::Vector2f x_diff = particle_i.position - particle_j->position;

            sumFluid += (particle_j->mass / particle_j->density)
                    * ((v_diff * x_diff)/(x_diff * x_diff + 0.01f * std::pow((h),2)))
                    * kernelGradient(particle_i.position, particle_j->position);
        }
        if (particle_j->isStatic == true) {
            sf::Vector2f v_diff = particle_i.velocity - particle_j->velocity;
            sf::Vector2f x_diff = particle_i.position - particle_j->position;

            sumBoundary += (particle_j->mass / particle_j->density)
                    * ((v_diff * x_diff)/(x_diff * x_diff + 0.01f * std::pow((h),2)))
                    * kernelGradient(particle_i.position, particle_j->position);
        }
    }

    return (2 * viscosityFactor * sumFluid) + (viscosityFactor * sumBoundary);
}

sf::Vector2f SPHComputationsIISPH::computeSurfaceTension(Particle& particle_i) {
    sf::Vector2f sum(0.0f, 0.0f);

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic == false) {
            sum += kernelGradient(particle_i.position, particle_j->position);
        }
    }
    if (sum != sf::Vector2f(0.0f, 0.0f)) {
        return sum;
    } else {
        return sf::Vector2f(0.0f, 0.0f);
    }
}

sf::Vector2f SPHComputationsIISPH::computeTotalAcceleration(Particle& particle) {
    sf::Vector2f pressureAcceleration = computePressureAcceleration(particle);
    sf::Vector2f viscosityAcceleration = computeViscosity(particle);

    sf::Vector2f surfaceTension = computeSurfaceTension(particle);

    sf::Vector2f totalAcceleration = pressureAcceleration + viscosityAcceleration + gravity + surfaceTension * surfaceTensionFactor;

    return totalAcceleration;
}

void SPHComputationsIISPH::advectParticles(std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        if (particle.isStatic) {
            continue;
        }
        particle.velocity = particle.velocity + timeStep * particle.acceleration;
        // or maybe
        // particle.velocity = timeStep * particle.acceleration + particle.predicted_velocity;
        particle.position = particle.position + timeStep * particle.velocity;
    }
}

// max velo and color gradient to speed
float SPHComputationsIISPH::getMaxVelocity(std::vector<Particle>& particles) {
    float maxSpeed = 0;
    float speed = 0;
    float CFLNumber = 0;

    for (auto& particle : particles) {
        if (particle.isStatic) {
            continue;
        }

        speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    // Loop through each particle to apply the color based on the absolute speed
    for (auto& particle : particles) {
        if (particle.isStatic) continue;

        // Calculate the speed of the particle
        float speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);

        sf::Uint8 r, g, b;

        if (speed < speedThreshold1) {
            // Blue to Cyan
            float ratio = speed / speedThreshold1; // [0, speedThreshold1] -> [0, 1]
            r = 0;
            g = static_cast<sf::Uint8>(ratio * 255);
            b = 255;
        } else if (speed < speedThreshold2) {
            // Cyan to Green
            float ratio = (speed - speedThreshold1) / (speedThreshold2 - speedThreshold1); // [speedThreshold1, speedThreshold2] -> [0, 1]
            r = 0;
            g = 255;
            b = static_cast<sf::Uint8>((1 - ratio) * 255);
        } else {
            // Green to Red
            float ratio = (speed - speedThreshold2) / (maxSpeed - speedThreshold2); // [speedThreshold2, maxSpeed] -> [0, 1]
            r = static_cast<sf::Uint8>(ratio * 255);
            g = static_cast<sf::Uint8>((1 - ratio) * 255);
            b = 0;
        }

        // Set the particle's color
        particle.color = sf::Color(r, g, b, 255);
    }


    // color according to cfl number
    // for (auto& particle : particles) {
    //     if (particle.isStatic) continue;
    //     speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
    //     float normalizedSpeed = speed / maxSpeed;
    //     if (normalizedSpeed < 0.5f) {
    //     // Blue to White
    //     float ratio = normalizedSpeed * 2.0f; // [0, 0.5] -> [0, 1]
    //     sf::Uint8 r = static_cast<sf::Uint8>(ratio * 255);
    //     sf::Uint8 g = static_cast<sf::Uint8>(ratio * 255);
    //     particle.color = sf::Color(r, g, 255, 255);
    //     } else {
    //         // White to Red
    //         float ratio = (normalizedSpeed - 0.5f) * 2.0f; // [0.5, 1] -> [0, 1]
    //         sf::Uint8 r = 255;
    //         sf::Uint8 g = static_cast<sf::Uint8>((1 - ratio) * 255);
    //         particle.color = sf::Color(r, g, g, 255);
    //     }
    // }

    // for (auto& particle : particles) {
    //     if (particle.isStatic) continue;
    //     speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
    //     CFLNumber = timeStep * (h / speed);
    //     float normalizedSpeed = (speed-0)/(100-0); //min-max normalization (0-100)
    //     //printf("speed: %f", normalizedSpeed);
    //     if (normalizedSpeed < 10.0f) {
    //         float ratio = normalizedSpeed * 2.0f; // [0, 0.5] -> [0, 1]
    //         sf::Uint8 r = static_cast<sf::Uint8>(ratio * 255);
    //         sf::Uint8 g = static_cast<sf::Uint8>(ratio * 255);
    //         particle.color = sf::Color(r, g, 255, 255);
    //     }
    //     else if (normalizedSpeed >= 10.0f && speed < 60.0f) {
    //         float ratio = normalizedSpeed * 2.0f; // [0, 0.5] -> [0, 1]
    //         sf::Uint8 r = static_cast<sf::Uint8>(ratio * 255);
    //         sf::Uint8 g = static_cast<sf::Uint8>(ratio * 255);
    //         particle.color = sf::Color(r, 255, g, 255);
    //     }
    //     else {
    //         float ratio = (normalizedSpeed - 0.5f) * 2.0f; // [0.5, 1] -> [0, 1]
    //         sf::Uint8 r = 255;
    //         sf::Uint8 g = static_cast<sf::Uint8>((1 - ratio) * 255);
    //         particle.color = sf::Color(r, g, g, 255);
    //     }
    // }

    // for (auto& particle : particles) {
    //     if (particle.isStatic) continue;
    //     speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
    //     float normalizedSpeed = (speed-0)/(100-0); //min-max normalization min:0 max:100 selbstgesetzte werte
    //     if (normalizedSpeed <= 0.02f) {
    //         float ratio = normalizedSpeed / 0.05f;
    //         sf::Uint8 r = static_cast<sf::Uint8>(0 * (1 - ratio) + 0 * ratio);
    //         sf::Uint8 g = static_cast<sf::Uint8>(0 * (1 - ratio) + 255 * ratio);
    //         sf::Uint8 b = static_cast<sf::Uint8>(255 * (1 - ratio) + 255 * ratio);
    //         particle.color = sf::Color(r, g, b, 255);
    //     }
    //     else if (normalizedSpeed <= 0.2f) {
    //         float ratio = normalizedSpeed - 0.1f / 0.05f;
    //         sf::Uint8 r = static_cast<sf::Uint8>(0 * (1 - ratio) + 0 * ratio);
    //         sf::Uint8 g = static_cast<sf::Uint8>(255 * (1 - ratio) + 255 * ratio);
    //         sf::Uint8 b = static_cast<sf::Uint8>(255 * (1 - ratio) + 0 * ratio);
    //         particle.color = sf::Color(r, g, b, 255);
    //     }
    //     else if (normalizedSpeed <= 0.3f) {
    //         float ratio = normalizedSpeed - 0.3f / 0.05f;
    //         sf::Uint8 r = static_cast<sf::Uint8>(0 * (1 - ratio) + 255 * ratio);
    //         sf::Uint8 g = static_cast<sf::Uint8>(255 * (1 - ratio) + 255 * ratio);
    //         sf::Uint8 b = static_cast<sf::Uint8>(0 * (1 - ratio) + 0 * ratio);
    //         particle.color = sf::Color(r, g, b, 255);
    //     }
    //     else if (normalizedSpeed <= 0.4f) {
    //         float ratio = normalizedSpeed - 0.5f / 0.05f;
    //         sf::Uint8 r = static_cast<sf::Uint8>(255 * (1 - ratio) + 255 * ratio);
    //         sf::Uint8 g = static_cast<sf::Uint8>(255 * (1 - ratio) + 165 * ratio);
    //         sf::Uint8 b = static_cast<sf::Uint8>(0 * (1 - ratio) + 0 * ratio);
    //         particle.color = sf::Color(r, g, b, 255);
    //     }
    //     else {
    //         float ratio = normalizedSpeed - 0.7f / 0.05f;
    //         sf::Uint8 r = static_cast<sf::Uint8>(255 * (1 - ratio) + 255 * ratio);
    //         sf::Uint8 g = static_cast<sf::Uint8>(165 * (1 - ratio) + 0 * ratio);
    //         sf::Uint8 b = static_cast<sf::Uint8>(0 * (1 - ratio) + 0 * ratio);
    //         particle.color = sf::Color(r, g, b, 255);
    //     }
    // }

    return maxSpeed;
}

void SPHComputationsIISPH::isCFLConditionTrue(std::vector<Particle>& particles) {
    float lambda = 0.9f;
    float value;
    // value = lambda * (h / getMaxVelocity(particles));
    // printf("timeStep: %f  CFL: %f    maxVeloc: %f \n", timeStep, value, getMaxVelocity(particles));
    if (timeStep <= lambda * (h / getMaxVelocity(particles))) {
        CFLCondition = true;
    } else {
        CFLCondition = false;
    }
}

float SPHComputationsIISPH::getCourantNumber(std::vector<Particle>& particles) {
    float value;
    value = timeStep * (h / getMaxVelocity(particles));
    return value;
}

float SPHComputationsIISPH::getAvgDensity(std::vector<Particle>& particles, const std::string& filename) {
    float avgDensity = 0.0f;
    int count = 0;
    for (auto& particle : particles) {
        if (particle.isStatic) continue;
        if (particle.density < density_rest) continue;
        // if (particle.density > density_rest) {
        //     avgDensity += particle.density;
        //     count++;
        // }
        avgDensity += particle.density;
        count++;

    }

    avgDensity /= count;

    // Append the average density value to the text file
    std::ofstream outFile(filename, std::ios_base::app);
    if (outFile.is_open()) {
        outFile << avgDensity << "\n";
        outFile.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }

    return avgDensity;
}

// float SPHComputations::getRelativeDensityError(std::vector<Particle>& particles, const std::string& filename) {
//     float relDensity = 0.0f;
//     int count = 0;
//     for (auto& particle : particles) {
//         if (particle.isStatic) continue;
//         if (particle.density < density_rest) continue;

//         relDensity =

//     }

//     relDensity /= particles.size();

//     // Append the average density value to the text file
//     std::ofstream outFile(filename, std::ios_base::app);
//     if (outFile.is_open()) {
//         outFile << relDensity << "\n";
//         outFile.close();
//     } else {
//         std::cerr << "Unable to open file " << filename << std::endl;
//     }

//     return relDensity;
// }