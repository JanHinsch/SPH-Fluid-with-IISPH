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

// version with volume instead of density this is the OLD version of EOS
sf::Vector2f SPHComputationsIISPH::computeViscosity(Particle& particle_i) {
    sf::Vector2f sumFluid(0.0f, 0.0f);
    sf::Vector2f sumBoundary(0.0f, 0.0f);
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic == false) {
            sf::Vector2f v_diff = particle_i.velocity - particle_j->velocity;
            sf::Vector2f x_diff = particle_i.position - particle_j->position;

            sumFluid += particle_j->volume
                    * ((v_diff * x_diff)/(x_diff * x_diff + 0.01f * std::pow((h),2)))
                    * kernelGradient(particle_i.position, particle_j->position);
        }
        if (particle_j->isStatic == true) {
            sf::Vector2f v_diff = particle_i.velocity - particle_j->velocity;
            sf::Vector2f x_diff = particle_i.position - particle_j->position;

            sumBoundary += particle_j->volume
                    * ((v_diff * x_diff)/(x_diff * x_diff + 0.01f * std::pow((h),2)))
                    * kernelGradient(particle_i.position, particle_j->position);
        }
    }

    return (2 * viscosityFactor * sumFluid) + (viscosityFactor * sumBoundary);
}

// sf::Vector2f SPHComputationsIISPH::computeSurfaceTension(Particle& particle_i) {
//     sf::Vector2f sum(0.0f, 0.0f);
//
//     for (auto& particle_j : particle_i.particle_neighbours) {
//         if (particle_j->isStatic == false) {
//             sum += kernelGradient(particle_i.position, particle_j->position);
//         }
//     }
//     if (sum != sf::Vector2f(0.0f, 0.0f)) {
//         return sum;
//     } else {
//         return sf::Vector2f(0.0f, 0.0f);
//     }
// }

sf::Vector2f SPHComputationsIISPH::computeTotalAcceleration(Particle& particle) {
    sf::Vector2f pressureAcceleration = computePressureAcceleration(particle);
    sf::Vector2f viscosityAcceleration = computeViscosity(particle);

    // sf::Vector2f surfaceTension = computeSurfaceTension(particle);

    sf::Vector2f totalAcceleration = pressureAcceleration + viscosityAcceleration + gravity;

    return totalAcceleration;
}

void SPHComputationsIISPH::advectParticles(std::vector<Particle>& particles) {
    int it = 0;
    for (auto& particle : particles) {
        if (particle.isStatic) continue;

        // particle.velocity = particle.velocity + timeStep * particle.acceleration;
        // or maybe
        particle.velocity = timeStep * particle.acceleration + particle.predicted_velocity;
        particle.position = particle.position + timeStep * particle.velocity;
    }
}

//////////////////////////// Pressure Boundaries /////////////////
// eq 12
float SPHComputationsIISPH::computeRestVolumeBoundary(Particle& particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (!(particle_j->isStatic)) continue;
        sum += kernel(particle_i.position, particle_j->position);
    }
    return gamma_1/sum;
}

// // // eq 15 doesnt work for some reason see page. 4
// float SPHComputationsIISPH::computeVolumeBoundary(Particle& particle_i) {
//     float sum = 0.0f;
//     for (auto& particle_j : particle_i.particle_neighbours) {
//         if (particle_j->isStatic) continue;
//         sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position) + gamma_1 + beta;
//     }
//     return particle_i.restVolume / sum;
// }

// eq 14
float SPHComputationsIISPH::computeVolumeBoundary(Particle& particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
        } else {
            sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
        }
    }
    return particle_i.restVolume / sum + beta;
}

// eq 13
float SPHComputationsIISPH::computeVolumeFluid(Particle& particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
        } else {
            sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
        }
    }
    return particle_i.restVolume / std::max(sum, 0.0000001f);
}

// does not work currently
sf::Vector2f SPHComputationsIISPH::computeViscosityIISPH(Particle& particle_i) {
    sf::Vector2f v_diff = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f sum = sf::Vector2f(0.0f, 0.0f);
    for (auto& particle_j : particle_i.particle_neighbours) {
        v_diff = particle_j->velocity - particle_i.velocity;
        particle_j->density = particle_j->mass / particle_j->volume;
        sum += (particle_j->mass / particle_j->density) * v_diff * kernel(particle_i.position, particle_j->position);
    }
    return (viscosityFactor/timeStep) * sum;
}


sf::Vector2f SPHComputationsIISPH::predictVelocity(Particle& particle_i) {
    sf::Vector2f non_pressure_acc = computeViscosity(particle_i)
                                    + gravity;
    non_pressure_acc *= timeStep;
    return particle_i.velocity + non_pressure_acc;
}

// eq 16
float SPHComputationsIISPH::computeVelocityDivergenceFluid(Particle &particle_i) {
    float sumBoundary = 0.0f;
    float sumFluid = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sumBoundary += particle_j->volume *
                dotProduct((particle_i.predicted_velocity - 0), kernelGradient(particle_i.position, particle_j->position));
        } else {
            sumFluid += particle_j->volume *
                dotProduct((particle_i.predicted_velocity - particle_j->predicted_velocity), kernelGradient(particle_i.position, particle_j->position));
        }
    }
    return - sumFluid - sumBoundary;
}

// eq 17
float SPHComputationsIISPH::computeVelocityDivergenceBoundary(Particle &particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (!(particle_j->isStatic)) {
            sum += particle_j->volume *
                dotProduct((particle_i.predicted_velocity - particle_j->predicted_velocity), kernelGradient(particle_i.position, particle_j->position));
        }
    }
    return - sum;
}

// eq 9 and 10
float SPHComputationsIISPH::computeSourceTerm(Particle &particle_i) {
    float velocityDivergence = 0.0f;
    if (particle_i.isStatic) {
        velocityDivergence = computeVelocityDivergenceBoundary(particle_i);
    } else {
        velocityDivergence = computeVelocityDivergenceFluid(particle_i);
    }
    return 1 - (particle_i.restVolume / particle_i.volume) + timeStep * velocityDivergence;

}

// eq 19
float SPHComputationsIISPH::computeDiagonalElementFluid(Particle &particle_i) {
    sf::Vector2f sum1 = sf::Vector2f(0.0f, 0.0f);
    float sum3 = 0.0f;

    // first term
    for (auto& particle_j : particle_i.particle_neighbours) {
        sum1 += particle_j->volume * kernelGradient(particle_i.position, particle_j->position);
    }
    // calculate the square of the norm of sum1
    float normSquared1 = sum1.x * sum1.x + sum1.y * sum1.y;

    normSquared1 *= (particle_i.volume / particle_i.mass) * (timeStep * timeStep);

    // second term
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sf::Vector2f sum2 = kernelGradient(particle_i.position, particle_j->position);
        float normSquared2 = sum2.x * sum2.x + sum2.y * sum2.y;

        sum3 += particle_j->volume * (particle_j->volume / particle_j->mass) * normSquared2;
    }
    sum3 *= (timeStep * timeStep) * particle_i.volume;


    return - normSquared1 - sum3;
}

// eq 20
float SPHComputationsIISPH::computeDiagonalElementBoundary(Particle &particle_i) {
    sf::Vector2f sum2 = sf::Vector2f(0.0f, 0.0f);
    float sum3 = 0.0f;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sum2 = kernelGradient(particle_i.position, particle_j->position);
        float normSquared2 = sum2.x * sum2.x + sum2.y * sum2.y;

        sum3 += particle_j->volume * (particle_j->volume / particle_j->mass) * normSquared2;
    }
    sum3 *= (timeStep * timeStep) * particle_i.volume;

    return - sum3;
}

// eq 8
sf::Vector2f SPHComputationsIISPH::computePressureAcceleration(Particle &particle_i) {
    sf::Vector2f sum = sf::Vector2f(0.0f, 0.0f);

    for (auto& particle_j : particle_i.particle_neighbours) {
        sum += particle_j->volume * (particle_i.pressure + particle_j->pressure)
                * kernelGradient(particle_i.position, particle_j->position);
    }
    sum *= (particle_i.volume / particle_i.mass);

    return - sum;
}

// eq 9
float SPHComputationsIISPH::computeApFluid(Particle &particle_i) {
    float sumBoundary = 0.0f;
    float sumFluid = 0.0f;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sumBoundary += particle_j->volume *
                dotProduct(particle_i.acceleration, kernelGradient(particle_i.position, particle_j->position));
        } else {
            sumFluid += particle_j->volume *
                dotProduct((particle_i.acceleration - particle_j->acceleration), kernelGradient(particle_i.position, particle_j->position));
        }
    }


    return (sumBoundary + sumFluid) * (timeStep * timeStep);
}

// eq 10
float SPHComputationsIISPH::computeApBoundary(Particle &particle_i) {
    float sumFluid = 0.0f;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sumFluid += particle_j->volume *
            dotProduct(particle_j->acceleration, kernelGradient(particle_i.position, particle_j->position));
    }

    return -1 * sumFluid * (timeStep * timeStep);
}

// eq 18
float SPHComputationsIISPH::updatePressure(Particle &particle_i) {
    float pressure = particle_i.pressure + (omega * ((particle_i.sourceTerm - particle_i.Ap) / particle_i.diagonalElement));
    if (pressure < 0) {
        return 0;
    }
    return pressure;
}

// eq 18 for boundary with non constant omega
float SPHComputationsIISPH::updatePressureBoundaries(Particle &particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sum += kernel(particle_i.position, particle_j->position);
        }
    }
    float omega_boundary = 0.5f * (gamma_1/((h*h) * sum));

    float pressure = particle_i.pressure + (omega_boundary * ((particle_i.sourceTerm - particle_i.Ap) / particle_i.diagonalElement));
    if (pressure < 0) {
        return 0;
    }
    return pressure;
}


/////////////////////////////////////////////

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

    return maxSpeed;
}

void SPHComputationsIISPH::isCFLConditionTrue(std::vector<Particle>& particles) {
    float lambda = 0.9f;
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

void SPHComputationsIISPH::writeCurrentIterationToFile(const std::string& filename) {
    // Append the average density value to the text file
    std::ofstream outFile(filename, std::ios_base::app);
    if (outFile.is_open()) {
        outFile << currentIterations << "\n";
        outFile.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}