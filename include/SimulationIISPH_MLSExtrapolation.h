#ifndef SIMULATIONIISPH_MLSEXTRAPOLATION_H
#define SIMULATIONIISPH_MLSEXTRAPOLATION_H

#include "Particle.h"
#include <vector>
#include <SFML/Graphics.hpp>

class SPHComputationsIISPH_MLSExtra {
public:
    static float kernelAlpha();
    static float kernel(sf::Vector2f positionA, sf::Vector2f positionB);
    static sf::Vector2f kernelGradient(sf::Vector2f positionA, sf::Vector2f positionB);

    static float computeDensity(Particle& particle_i);

    static float computeFactorLambda(Particle &particle_i);

    static float computeSourceTermDivergenceFree(Particle &particle_i);

    static float computeSourceTermDensityError(Particle &particle_i);

    static float setPressureDivergenceError(Particle &particle_i);

    static float setPressureDensityError(Particle &particle_i);

    static float computeAlpha(Particle &particle_i);

    static sf::Vector2f computeDb(Particle &particle_i);

    static float computePressureMLS(Particle &particle_i);

    static float computePressureAcceleration(float density_i);

    static sf::Vector2f predictVelocity(Particle &particle_i);

    static float computeVelocityDivergenceFluid(Particle &particle_i);

    static float computeVelocityDivergenceBoundary(Particle &particle_i);

    static float computeVelocityDivergence(Particle &particle_i);

    static float computeSourceTerm(Particle &particle_i);

    static float computeDiagonalElementFluid(Particle &particle_i);

    static float computeDiagonalElementBoundary(Particle &particle_i);

    static float computeDiagonalElement(Particle &particle_i);

    static bool isParticleCompressed(float density_i);

    static float computePressure(float density_i);

    static sf::Vector2f computePressureAcceleration(Particle& particle_i);

    static float computeApFluid(Particle &particle_i);

    static float computeApBoundary(Particle &particle_i);

    static float computeDivergence(Particle &particle_i);

    static float updatePressure(Particle &particle_i);

    static float updatePressureBoundaries(Particle &particle_i);

    static float updatePressureMLSFluid(Particle &particle_i);

    static sf::Vector2f computeViscosity(Particle& particle_i);

    static sf::Vector2f computeSurfaceTension(Particle &particle_i);

    static sf::Vector2f computeTotalAcceleration(Particle& particle);

    static void advectParticles(std::vector<Particle> &particles);

    static float getMaxVelocity(std::vector<Particle>& particles);

    static float getAvgDensity(std::vector<Particle> &particles, const std::string& filename);

    static void writeCurrentIterationToFile(const std::string &filename);

    static void isCFLConditionTrue(std::vector<Particle>& particles);

    static float getCourantNumber(std::vector<Particle> &particles);

    static float computeAdaptiveTimeStep(float adaptiveTimeStep);

    static float computeRestVolumeBoundary(Particle& particle_i);

    static float computeVolumeBoundary(Particle &particle_i);

    static float computeVolumeFluid(Particle &particle_i);

    static sf::Vector2f computeViscosityIISPH(Particle &particle_i);

private:

};


#endif //SIMULATIONIISPH_MLSEXTRAPOLATION_H
