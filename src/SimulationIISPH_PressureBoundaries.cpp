// SimulationIISPH.cpp
#include "../include/SimulationEOS.h"
#include "../include/SimulationIISPH_PressureBoundaries.h"
#include "../include/ParticleSystem.h"
#include "../include/Globals.h"
#include <cmath>
#include <iostream>
#include <fstream>


float SPHComputationsIISPH_PressureBoundaries::kernelAlpha() {
    return 5 / (14 * M_PI * std::pow(h, 2));
}

float SPHComputationsIISPH_PressureBoundaries::kernel(sf::Vector2f positionA,sf::Vector2f positionB) {
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

sf::Vector2f SPHComputationsIISPH_PressureBoundaries::kernelGradient(sf::Vector2f positionA,sf::Vector2f positionB) {
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

// version with volume instead of density this is the OLD version of EOS TODO i changed sumBOundary factor
sf::Vector2f SPHComputationsIISPH_PressureBoundaries::computeViscosity(Particle& particle_i) {
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

    return (2 * viscosityFactor * sumFluid) + (viscosityFactor * sumBoundary * 0.0f);
}

sf::Vector2f SPHComputationsIISPH_PressureBoundaries::computeSurfaceTension(Particle& particle_i) {
    sf::Vector2f sum(0.0f, 0.0f);

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic == false) {
            sum += kernelGradient(particle_i.position, particle_j->position) * surfaceTensionFactor;
        }
    }
    if (sum != sf::Vector2f(0.0f, 0.0f)) {
        return sum;
    } else {
        return sf::Vector2f(0.0f, 0.0f);
    }
}


void SPHComputationsIISPH_PressureBoundaries::advectParticles(std::vector<Particle>& particles) {
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
float SPHComputationsIISPH_PressureBoundaries::computeRestVolumeBoundary(Particle& particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (!(particle_j->isStatic)) continue;
        sum += kernel(particle_i.position, particle_j->position);
    }
    return gamma_1/sum;
}

// eq 15 (alternative zu eq 14)
float SPHComputationsIISPH_PressureBoundaries::computeVolumeBoundary(Particle& particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
    }
    return particle_i.restVolume / (sum + gamma_2 + beta);
}

// eq 14
// float SPHComputationsIISPH_PressureBoundaries::computeVolumeBoundary(Particle& particle_i) {
//     float sum = 0.0f;
//     for (auto& particle_j : particle_i.particle_neighbours) {
//         if (particle_j->isStatic) {
//             sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
//         } else {
//             sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
//         }
//     }
//     return particle_i.restVolume / sum + beta;
// }

// eq 13
float SPHComputationsIISPH_PressureBoundaries::computeVolumeFluid(Particle& particle_i) {
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

// eq 13 mit gamma
float SPHComputationsIISPH_PressureBoundaries::computeVolumeFluidWithGamma(Particle& particle_i) {
    float sum = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position) * gamma_1;
        } else {
            sum += particle_j->restVolume * kernel(particle_i.position, particle_j->position);
        }
    }
    return particle_i.restVolume / std::max(sum, 0.0000001f);
}

// does not work currently
sf::Vector2f SPHComputationsIISPH_PressureBoundaries::computeViscosityIISPH(Particle& particle_i) {
    sf::Vector2f v_diff = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f sum = sf::Vector2f(0.0f, 0.0f);
    for (auto& particle_j : particle_i.particle_neighbours) {
        v_diff = particle_j->velocity - particle_i.velocity;
        particle_j->density = particle_j->mass / particle_j->volume;
        sum += (particle_j->mass / particle_j->density) * v_diff * kernel(particle_i.position, particle_j->position);
    }
    return (viscosityFactor/timeStep) * sum;
}


sf::Vector2f SPHComputationsIISPH_PressureBoundaries::predictVelocity(Particle& particle_i) {
    sf::Vector2f non_pressure_acc = computeViscosity(particle_i)
                                    + gravity + computeSurfaceTension(particle_i);
    non_pressure_acc *= timeStep;
    return particle_i.velocity + non_pressure_acc;
}

// eq 16
float SPHComputationsIISPH_PressureBoundaries::computeVelocityDivergenceFluid(Particle &particle_i) {
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
float SPHComputationsIISPH_PressureBoundaries::computeVelocityDivergenceBoundary(Particle &particle_i) {
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
float SPHComputationsIISPH_PressureBoundaries::computeSourceTerm(Particle &particle_i) {
    float velocityDivergence = 0.0f;
    if (particle_i.isStatic) {
        velocityDivergence = computeVelocityDivergenceBoundary(particle_i);
    } else {
        velocityDivergence = computeVelocityDivergenceFluid(particle_i);
    }
    return 1 - (particle_i.restVolume / particle_i.volume) + timeStep * velocityDivergence;

}

// eq 19
float SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementFluid(Particle &particle_i) {
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
float SPHComputationsIISPH_PressureBoundaries::computeDiagonalElementBoundary(Particle &particle_i) {
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
sf::Vector2f SPHComputationsIISPH_PressureBoundaries::computePressureAcceleration(Particle &particle_i) {
    sf::Vector2f sum = sf::Vector2f(0.0f, 0.0f);

    for (auto& particle_j : particle_i.particle_neighbours) {
        sum += particle_j->volume * (particle_i.pressure + particle_j->pressure)
                * kernelGradient(particle_i.position, particle_j->position);
    }
    sum *= (particle_i.volume / particle_i.mass);

    return - sum;
}

// eq 9
float SPHComputationsIISPH_PressureBoundaries::computeApFluid(Particle &particle_i) {
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
float SPHComputationsIISPH_PressureBoundaries::computeApBoundary(Particle &particle_i) {
    float sumFluid = 0.0f;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sumFluid += particle_j->volume *
            dotProduct(particle_j->acceleration, kernelGradient(particle_i.position, particle_j->position));
    }

    return - sumFluid * (timeStep * timeStep);
}

// eq 18
float SPHComputationsIISPH_PressureBoundaries::updatePressure(Particle &particle_i) {
    float pressure = particle_i.pressure + (omega * ((particle_i.sourceTerm - particle_i.Ap) / particle_i.diagonalElement));
    if (pressure < 0) {
        return 0;
    }
    return pressure;
}

// eq 18 for boundary with non constant omega
float SPHComputationsIISPH_PressureBoundaries::updatePressureBoundaries(Particle &particle_i) {
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

// TODO its from SPH Extrapolation (MLS short chapter in paper)
float SPHComputationsIISPH_PressureBoundaries::updatePressureBoundariesExtrapolation(Particle &particle_i) {
    float sum1 = 0.0f;
    sf::Vector2f sum2 = sf::Vector2f(0.0f, 0.0f);
    float sum3 = 0.0f;

    // for (auto& particle_j : particle_i.particle_neighbours) {
    //     if (particle_j->isStatic) continue;
    //     sum3 += kernel(particle_i.position, particle_j->position);
    //     sum1 += particle_i.pressure * sum3;
    //     sum2 += (particle_j->mass / particle_j->volume) * (particle_i.position - particle_j->position) * sum3;
    // }
    // float sum2_2 = dotProduct(gravity, sum2);
    //
    // return (sum1 + sum2_2) / sum3;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sum3 += kernel(particle_i.position, particle_j->position);
    }
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sum2 += (particle_j->mass / particle_j->volume) * (particle_i.position - particle_j->position)
                * kernel(particle_i.position, particle_j->position);
    }
    float sum2_2 = dotProduct(gravity, sum2);
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sum1 += particle_i.pressure * kernel(particle_i.position, particle_j->position);
    }
    return (sum1 + sum2_2) / (sum3);
}

// TODO its for SPH Extrapolation (MLS short chapter in paper) eq 1
sf::Vector2f SPHComputationsIISPH_PressureBoundaries::computePressureAccelerationExtrapolation(Particle &particle_i) {
    sf::Vector2f sumFluid = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f sumBoundary = sf::Vector2f(0.0f, 0.0f);
    float density_squared_i = (particle_i.mass / particle_i.volume) * (particle_i.mass / particle_i.volume);

    for (auto& particle_j : particle_i.particle_neighbours) {
        float density_squared_j = (particle_j->mass / particle_j->volume) * (particle_j->mass / particle_j->volume);
        if (particle_j->isStatic) {
            sumBoundary += particle_j->mass * ((particle_i.pressure / density_squared_i) + (particle_j->pressure / density_squared_j))
                * kernelGradient(particle_i.position, particle_j->position);
        } else {
            sumFluid += particle_j->mass * ((particle_i.pressure / density_squared_i) + (particle_j->pressure / density_squared_j))
                * kernelGradient(particle_i.position, particle_j->position);
        }
    }

    return - sumFluid - gamma_3 * sumBoundary;
}

// TODO its for Pressure Mirroring (MLS short chapter in paper) eq 2
sf::Vector2f SPHComputationsIISPH_PressureBoundaries::computePressureAccelerationMirror(Particle &particle_i) {
    sf::Vector2f sumFluid = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f sumBoundary = sf::Vector2f(0.0f, 0.0f);
    float density_squared_i = (particle_i.mass / particle_i.volume) * (particle_i.mass / particle_i.volume);

    for (auto& particle_j : particle_i.particle_neighbours) {
        float density_squared_j = (particle_j->mass / particle_j->volume) * (particle_j->mass / particle_j->volume);
        if (particle_j->isStatic) {
            sumBoundary += particle_j->mass * 2 *(particle_i.pressure / density_squared_i)
                * kernelGradient(particle_i.position, particle_j->position);
        } else {
            sumFluid += particle_j->mass * ((particle_i.pressure / density_squared_i) + (particle_j->pressure / density_squared_j))
                * kernelGradient(particle_i.position, particle_j->position);
        }
    }
    return - sumFluid - gamma_3 * sumBoundary;
}

// TODO its for Pressure Zero at boundaries
sf::Vector2f SPHComputationsIISPH_PressureBoundaries::computePressureAccelerationZero(Particle &particle_i) {
    sf::Vector2f sumFluid = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f sumBoundary = sf::Vector2f(0.0f, 0.0f);
    float density_squared_i = (particle_i.mass / particle_i.volume) * (particle_i.mass / particle_i.volume);

    for (auto& particle_j : particle_i.particle_neighbours) {
        float density_squared_j = (particle_j->mass / particle_j->volume) * (particle_j->mass / particle_j->volume);
        if (particle_j->isStatic) continue;
        sumFluid += particle_j->mass * ((particle_i.pressure / density_squared_i) + (particle_j->pressure / density_squared_j))
            * kernelGradient(particle_i.position, particle_j->position);
    }
    return - sumFluid;
}


/////////////////////////////////////////////

// max velo and color gradient to speed
float SPHComputationsIISPH_PressureBoundaries::getMaxVelocity(std::vector<Particle>& particles) {
    float maxSpeed = 0;
    float minSpeed = 0;
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

    // minSpeed = 0.0001f;
    // // set threshholds
    // float speedThreshold1 = minSpeed + 0.1f * (maxSpeed - minSpeed);
    // float speedThreshold2 = minSpeed + 0.8f * (maxSpeed - minSpeed);

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
float SPHComputationsIISPH_PressureBoundaries::getMaxVelocityPressure(std::vector<Particle>& particles) {
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

void SPHComputationsIISPH_PressureBoundaries::isCFLConditionTrue(std::vector<Particle>& particles) {
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
    if (adaptiveTimeStepping) {
        timeStep = computeAdaptiveTimeStep(check);
    }
}

float SPHComputationsIISPH_PressureBoundaries::getCourantNumber(std::vector<Particle>& particles) {
    float value;
    value = timeStep * (h / getMaxVelocity(particles));
    return value;
}

float SPHComputationsIISPH_PressureBoundaries::computeAdaptiveTimeStep(float adaptiveTimeStep) {
    // Clamp the time step to a reasonable range
    adaptiveTimeStep = std::max(minTimeStep, std::min(adaptiveTimeStep, maxTimeStep));
    return adaptiveTimeStep;
}

float SPHComputationsIISPH_PressureBoundaries::getAvgDensity(std::vector<Particle>& particles, const std::string& filename) {
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

void SPHComputationsIISPH_PressureBoundaries::writeCurrentIterationToFile(const std::string& filename) {
    // // Append the average density value to the text file
    // std::ofstream outFile(filename, std::ios_base::app);
    // if (outFile.is_open()) {
    //     outFile << currentIterations << "\n";
    //     outFile.close();
    // } else {
    //     std::cerr << "Unable to open file " << filename << std::endl;
    // }

    // changed fast to volume error
    std::ofstream outFile(filename, std::ios_base::app);
    if (outFile.is_open()) {
        outFile << currentVolumeError << "\n";
        outFile.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}