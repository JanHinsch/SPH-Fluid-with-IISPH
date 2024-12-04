// SimulationEOS.h
#ifndef SIMULATIONEOS_H
#define SIMULATIONEOS_H

#include "Particle.h"
#include <vector>
#include <SFML/Graphics.hpp>

class SPHComputationsEOS {
public:
    static float kernelAlpha();
    static float kernel(sf::Vector2f positionA, sf::Vector2f positionB);
    static sf::Vector2f kernelGradient(sf::Vector2f positionA, sf::Vector2f positionB);

    static float computeDensity(Particle& particle, std::vector<Particle*>& neighbours);

    static bool isParticleCompressed(float density_i);

    static float computePressure(float density_i);

    static sf::Vector2f computePressureAcceleration(Particle& particle, std::vector<Particle*>& neighbours, float pressure_i);

    static sf::Vector2f computeViscosity(Particle& particle_i, std::vector<Particle*>& neighbours);

    static sf::Vector2f computeSurfaceTension(Particle &particle_i, std::vector<Particle*> &neighbours);

    static sf::Vector2f computeTotalAcceleration(Particle& particle, std::vector<Particle*>& neighbours);

    static void advectParticles(std::vector<Particle> &particles);

    static float getMaxVelocity(std::vector<Particle>& particles);

    static float getMaxVelocityPressure(std::vector<Particle> &particles);

    static float getAvgDensity(std::vector<Particle> &particles, const std::string& filename);

    static void isCFLConditionTrue(std::vector<Particle>& particles);

    static float getCourantNumber(std::vector<Particle> &particles);

    static float computeAdaptiveTimeStep(float adaptiveTimeStep);

private:
    
};


#endif // SIMULATIONEOS_H