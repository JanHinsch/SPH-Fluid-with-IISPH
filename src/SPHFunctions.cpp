#include "SPHFunctions.h"
#include "../include/SimulationEOS.h"
#include "../include/ParticleSystem.h"
#include "../include/Globals.h"
#include <cmath>
#include <iostream>
#include <fstream>


float SPHFunctions::kernelAlpha() {
    return 5 / (14 * M_PI * std::pow(h, 2));
}

float SPHFunctions::kernel(sf::Vector2f positionA,sf::Vector2f positionB) {
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
    return result;
}

sf::Vector2f SPHFunctions::kernelGradient(sf::Vector2f positionA,sf::Vector2f positionB) {
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
    } else if (q < 2.0f) {
        return (kernelAlpha() * (diff / (distance * h))) * (-3 * std::pow((2 - q),2));
    } else {
        return sf::Vector2f(0.0f, 0.0f);
    }
}

void SPHFunctions::advectParticles(std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        if (particle.isStatic) {
            continue;
        }
        particle.velocity = particle.velocity + timeStep * particle.acceleration;
        particle.position = particle.position + timeStep * particle.velocity;
    }
}








//////////////////////////////////////////////////////////

// max velo and color gradient to speed
float SPHFunctions::getMaxVelocity(std::vector<Particle>& particles) {
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

// max velo and color gradient to pressure
float SPHFunctions::getMaxVelocityPressure(std::vector<Particle>& particles) {
    float maxSpeed = 0;
    float speed = 0;

    float maxPressure = 0;
    float minPressure = 0;

    for (auto& particle : particles) {
        if (particle.isStatic) {
            continue;
        }

        if (particle.pressure > maxPressure) {
            maxPressure = particle.pressure;
        }
        if (particle.pressure < minPressure) {
            minPressure = particle.pressure;
        }

        speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    // to avoid partciles disappearing when all have same pressure
    minPressure = 0.0001f;
    // set threshholds
    float speedThreshold1 = minPressure + 0.1f * (maxPressure - minPressure);
    float speedThreshold2 = minPressure + 0.3f * (maxPressure - minPressure);

    // Loop through each particle to apply the color based on the pressure
    for (auto& particle : particles) {
        sf::Uint8 r, g, b;

        if (particle.pressure < speedThreshold1) {
            // Blue to Cyan
            float ratio = particle.pressure / speedThreshold1; // [0, speedThreshold1] -> [0, 1]
            r = 0;
            g = static_cast<sf::Uint8>(ratio * 255);
            b = 255;
        } else if (particle.pressure < speedThreshold2) {
            // Cyan to Green
            float ratio = (particle.pressure - speedThreshold1) / (speedThreshold2 - speedThreshold1); // [speedThreshold1, speedThreshold2] -> [0, 1]
            r = 0;
            g = 255;
            b = static_cast<sf::Uint8>((1 - ratio) * 255);
        } else {
            // Green to Red
            float ratio = (particle.pressure - speedThreshold2) / (maxPressure - speedThreshold2); // [speedThreshold2, maxSpeed] -> [0, 1]
            r = static_cast<sf::Uint8>(ratio * 255);
            g = static_cast<sf::Uint8>((1 - ratio) * 255);
            b = 0;
        }

        // Set the particle's color
        particle.color = sf::Color(r, g, b, 255);
    }

    return maxSpeed;
}

void SPHFunctions::isCFLConditionTrue(std::vector<Particle>& particles) {
    float lambda = 0.5f;
    float check = 0.0f;
    if (pressureColors) {
        check = lambda * (h / getMaxVelocityPressure(particles));
    } else {
        check = lambda * (h / getMaxVelocity(particles));

    }
    if (timeStep <= check) {
        CFLCondition = true;
    } else {
        CFLCondition = false;
    }
    timeStep = computeAdaptiveTimeStep(check);
}

float SPHFunctions::getCourantNumber(std::vector<Particle>& particles) {
    float value;
    value = timeStep * (h / getMaxVelocity(particles));
    return value;
}

float SPHFunctions::computeAdaptiveTimeStep(float adaptiveTimeStep) {
    // Clamp the time step to a reasonable range
    adaptiveTimeStep = std::max(minTimeStep, std::min(adaptiveTimeStep, maxTimeStep));
    return adaptiveTimeStep;
}

float SPHFunctions::getAvgDensity(std::vector<Particle>& particles, const std::string& filename) {
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

void SPHFunctions::writeCurrentIterationToFile(const std::string& filename) {
    // Append the average density value to the text file
    std::ofstream outFile(filename, std::ios_base::app);
    if (outFile.is_open()) {
        outFile << currentIterations << "\n";
        outFile.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}