// SimulationIISPH.cpp
#include "../include/SimulationEOS.h"
#include "../include/SimulationIISPH_MLSExtrapolation.h"
#include "../include/ParticleSystem.h"
#include "../include/Globals.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <TGUI/Vector2.hpp>
#include <../external/eigen-master/Eigen/Dense>
#include <../external/eigen-master/Eigen/SVD>


float SPHComputationsIISPH_MLSExtra::kernelAlpha() {
    return 5 / (14 * M_PI * std::pow(h, 2));
}

float SPHComputationsIISPH_MLSExtra::kernel(sf::Vector2f positionA,sf::Vector2f positionB) {
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

sf::Vector2f SPHComputationsIISPH_MLSExtra::kernelGradient(sf::Vector2f positionA,sf::Vector2f positionB) {
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

// version with volume instead of density this is the OLD version of EOS
sf::Vector2f SPHComputationsIISPH_MLSExtra::computeViscosity(Particle& particle_i) {
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

// sf::Vector2f SPHComputationsIISPH_MLSExtra::computeSurfaceTension(Particle& particle_i) {
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

sf::Vector2f SPHComputationsIISPH_MLSExtra::computeTotalAcceleration(Particle& particle) {
    sf::Vector2f pressureAcceleration = computePressureAcceleration(particle);
    sf::Vector2f viscosityAcceleration = computeViscosity(particle);

    // sf::Vector2f surfaceTension = computeSurfaceTension(particle);

    sf::Vector2f totalAcceleration = pressureAcceleration + viscosityAcceleration + gravity;

    return totalAcceleration;
}

void SPHComputationsIISPH_MLSExtra::advectParticles(std::vector<Particle>& particles) {
    int it = 0;
    for (auto& particle : particles) {
        if (particle.isStatic) continue;

        // particle.velocity = particle.velocity + timeStep * particle.acceleration;
        // or maybe
        particle.velocity = timeStep * particle.acceleration + particle.predicted_velocity;
        particle.position = particle.position + timeStep * particle.velocity;
    }
}

//////////////////////////// MLS Pressure Extrapolation /////////////////

float SPHComputationsIISPH_MLSExtra::computeDensity(Particle& particle) {
    float density_i = 0.0f;
    for (auto& particle_j : particle.particle_neighbours) {
        float kernel_value = kernel(particle.position, particle_j->position);
        float contribution = particle_j->mass * kernel_value;
        density_i += contribution;
    }
    return density_i;
}

float SPHComputationsIISPH_MLSExtra::computeFactorLambda(Particle& particle_i) {
    sf::Vector2f sum1 = sf::Vector2f(0.0f, 0.0f);
    float sum2 = 0.0f;

    for (auto& particle_j : particle_i.particle_neighbours) {
        sum1 += particle_j->mass * kernelGradient(particle_i.position, particle_j->position);
    }
    float sum1NormSq = sum1.x * sum1.x + sum1.y * sum1.y;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sf:tgui::Vector2f kernelSum = kernelGradient(particle_i.position, particle_j->position);
        float kernelSumSq = kernelSum.x * kernelSum.x + kernelSum.y * kernelSum.y;

        sum2 += particle_i.mass * particle_j->mass * kernelSumSq;
    }

    return (sum1NormSq + sum2) / (- particle_i.density * particle_i.density);
}

float SPHComputationsIISPH_MLSExtra::computeSourceTermDivergenceFree(Particle& particle_i) {
    float sumFluid = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sumFluid += (particle_j->mass / particle_j->density) *
            dotProduct((particle_i.velocity - particle_j->velocity), kernelGradient(particle_i.position, particle_j->position));
    }
    sumFluid *= -1;
    return - (1/timeStep) * sumFluid;
}

float SPHComputationsIISPH_MLSExtra::computeSourceTermDensityError(Particle& particle_i) {
    float sumFluid = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sumFluid += (particle_j->mass / particle_j->density) *
            dotProduct((particle_i.predicted_velocity - particle_j->predicted_velocity), kernelGradient(particle_i.position, particle_j->position));
    }
    sumFluid *= -1;
    return (1 / timeStep*timeStep) * (1 - (density_rest / particle_i.density)) - (1 / timeStep) * sumFluid;
}

float SPHComputationsIISPH_MLSExtra::setPressureDivergenceError(Particle& particle_i) {
    float sumFluid = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sumFluid += (particle_j->mass / particle_j->density) *
            dotProduct((particle_i.acceleration - particle_j->acceleration), kernelGradient(particle_i.position, particle_j->position));
    }
    sumFluid *= -1;

    // if (particle_i.diagonalElement == 0.0f) {
    //     return 0.0f;
    // }

    float ret = particle_i.pressure + (omega / particle_i.diagonalElement) * (particle_i.sourceTerm - sumFluid);

    return ret;
}

float SPHComputationsIISPH_MLSExtra::setPressureDensityError(Particle& particle_i) {
    float sumFluid = 0.0f;
    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sumFluid += (particle_j->mass / particle_j->density) *
            dotProduct((particle_i.acceleration - particle_j->acceleration), kernelGradient(particle_i.position, particle_j->position));
    }
    sumFluid *= -1;

    // if (particle_i.diagonalElement == 0.0f) {
    //     return 0.0f;
    // }

    return std::max(0.0f,particle_i.pressure + (omega / particle_i.diagonalElement) * (particle_i.sourceTerm - sumFluid));
}

float SPHComputationsIISPH_MLSExtra::computeAlpha(Particle& particle_i) {
    float sum1 = 0.0f;
    float sum2 = 0.0f;

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sum1 += particle_j->mass / particle_j->density * kernel(particle_i.position, particle_j->position);
        sum2 += sum1 * particle_j->pressure;
    }
    return sum2 / sum1;
}

// compute db and return the basis transformation
sf::Vector2f SPHComputationsIISPH_MLSExtra::computeDb(Particle& particle_i) {
    float sum1 = 0.0f;
    sf::Vector2f sum2 = sf::Vector2f(0.0f, 0.0f);

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) continue;
        sum1 += (particle_j->mass / particle_j->density) * kernel(particle_i.position, particle_j->position);
        sum2 += sum1 * particle_j->position;
    }

    sf::Vector2f db = sum2 / sum1;

    return particle_i.position - db;
}

// Function to safely invert a 2x2 matrix using SVD
Eigen::Matrix2f safeInvertMatrix(const Eigen::Matrix2f& M) {
    // Perform SVD on the matrix
    Eigen::JacobiSVD<Eigen::Matrix2f> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // Get singular values
    Eigen::Vector2f singularValues = svd.singularValues();

    // Regularize singular values (tolerance to avoid near-zero values)
    const float tolerance = 1e-6f;
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues[i] < tolerance) {
            singularValues[i] = 0.0f; // Treat small singular values as zero
        } else {
            singularValues[i] = 1.0f / singularValues[i]; // Invert non-zero singular values
        }
    }

    // Reconstruct the pseudo-inverse
    Eigen::Matrix2f pseudoInverse = svd.matrixV() * singularValues.asDiagonal() * svd.matrixU().transpose();

    return pseudoInverse;
}

sf::Vector2f computeHyperplaneCoefficients(const Eigen::Matrix2f& M_inv, const sf::Vector2f& sum2Value) {
    // Convert sf::Vector2f to Eigen::Vector2f
    Eigen::Vector2f eigenSum2Value(sum2Value.x, sum2Value.y);

    // Perform matrix-vector multiplication
    Eigen::Vector2f eigenResult = M_inv * eigenSum2Value;

    // Convert the result back to sf::Vector2f
    sf::Vector2f result(eigenResult.x(), eigenResult.y());

    return result;
}

// returns pressure
float SPHComputationsIISPH_MLSExtra::computePressureMLS(Particle& particle_i) {
    sf::Vector2f basisVector = computeDb(particle_i);

    // Compute elements of the 2D matrix
    float sumX2 = 0.0f, sumY2 = 0.0f, sumXY = 0.0f;
    sf::Vector2f sum2Value = sf::Vector2f(0.0f, 0.0f);
    for (auto& particle_j : particle_i.particle_neighbours) {
        // Matrix M values
        if (particle_j->isStatic) continue;
        sumX2 += basisVector.x * basisVector.x;
        sumXY += basisVector.x * basisVector.y;
        sumY2 += basisVector.y * basisVector.y;

        float kernelValue = kernel(particle_i.position, particle_j->position);
        float kernelValueSum = kernelValue * (particle_j->mass / particle_j->density);
        sumX2 *= kernelValueSum;
        sumXY *= kernelValueSum;
        sumY2 *= kernelValueSum;


        // second term values
        sum2Value += basisVector * particle_j->pressure * kernelValueSum;
    }

    // inverse M
    Eigen::Matrix2f M;
    M << sumX2, sumXY,
         sumXY, sumY2;

    // Safely invert the matrix
    Eigen::Matrix2f M_inv = safeInvertMatrix(M);

    // compute hyperplane coefficients vector
    sf::Vector2f planeCoefficients_noAlpha = computeHyperplaneCoefficients(M_inv, sum2Value);

    // construct c vector (all hyperplane coefficients with alpha, beta, gamma)
    float alpha = computeAlpha(particle_i);
    sf::Vector3f planeCoefficients = sf::Vector3f(alpha, planeCoefficients_noAlpha.x, planeCoefficients_noAlpha.y);
    sf::Vector3f basisVectorwith1 = sf::Vector3f(1, basisVector.x, basisVector.y);

    return dotProduct3f(planeCoefficients, basisVectorwith1);
}

sf::Vector2f SPHComputationsIISPH_MLSExtra::computePressureAcceleration(Particle &particle_i) {
    sf::Vector2f sum = sf::Vector2f(0.0f, 0.0f);

    for (auto& particle_j : particle_i.particle_neighbours) {
        if (particle_j->isStatic) {
            sum += particle_j->mass * ((particle_i.pressure / (particle_i. density * particle_i.density))
                                        + (particle_j->pressure / (particle_j->density * particle_j->density)))
                                        * kernelGradient(particle_i.position, particle_j->position);
        } else {
            sum += particle_j->mass * ((particle_i.pressure / (particle_i. density * particle_i.density))
                                        + (particle_j->pressure / (particle_j->density * particle_j->density)))
                                        * kernelGradient(particle_i.position, particle_j->position);
        }
    }
    return - sum;
}




//////////////////////////////////////////////////////////////////////////////////////////

// max velo and color gradient to speed
float SPHComputationsIISPH_MLSExtra::getMaxVelocity(std::vector<Particle>& particles) {
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

void SPHComputationsIISPH_MLSExtra::isCFLConditionTrue(std::vector<Particle>& particles) {
    float lambda = 0.9f;
    float check = lambda * (h / getMaxVelocity(particles));
    if (timeStep <= check) {
        CFLCondition = true;
    } else {
        CFLCondition = false;
    }
    timeStep = computeAdaptiveTimeStep(check);
}

float SPHComputationsIISPH_MLSExtra::getCourantNumber(std::vector<Particle>& particles) {
    float value;
    value = timeStep * (h / getMaxVelocity(particles));
    return value;
}

float SPHComputationsIISPH_MLSExtra::computeAdaptiveTimeStep(float adaptiveTimeStep) {
    // Clamp the time step to a reasonable range
    adaptiveTimeStep = std::max(minTimeStep, std::min(adaptiveTimeStep, maxTimeStep));
    return adaptiveTimeStep;
}

float SPHComputationsIISPH_MLSExtra::getAvgDensity(std::vector<Particle>& particles, const std::string& filename) {
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

void SPHComputationsIISPH_MLSExtra::writeCurrentIterationToFile(const std::string& filename) {
    // Append the average density value to the text file
    std::ofstream outFile(filename, std::ios_base::app);
    if (outFile.is_open()) {
        outFile << currentIterations << "\n";
        outFile.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}