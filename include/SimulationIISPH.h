#ifndef SIMULATIONIISPH_H
#define SIMULATIONIISPH_H

#include "Particle.h"
#include <vector>
#include <SFML/Graphics.hpp>

class SPHComputationsIISPH {
public:
    static float kernelAlpha();
    static float kernel(sf::Vector2f positionA, sf::Vector2f positionB);
    static sf::Vector2f kernelGradient(sf::Vector2f positionA, sf::Vector2f positionB);

    static float computeDensity(Particle& particle_i, std::vector<Particle*>& neighbours);

    static sf::Vector2f predictVelocity(Particle &particle_i);

    static float computeSourceTerm(Particle &particle_i);

    static float computeDiagonalElement(Particle &particle_i);

    static bool isParticleCompressed(float density_i);

    static float computePressure(float density_i);

    static sf::Vector2f computePressureAcceleration(Particle& particle_i);

    static float computeDivergence(Particle &particle_i);

    static float updatePressure(Particle &particle_i);

    static sf::Vector2f computeViscosity(Particle& particle_i);

    static sf::Vector2f computeSurfaceTension(Particle &particle_i);

    static sf::Vector2f computeTotalAcceleration(Particle& particle);

    static void advectParticles(std::vector<Particle> &particles);

    static float getMaxVelocity(std::vector<Particle>& particles);

    static float getAvgDensity(std::vector<Particle> &particles, const std::string& filename);

    static void isCFLConditionTrue(std::vector<Particle>& particles);

    static float getCourantNumber(std::vector<Particle> &particles);

private:

};



#endif //SIMULATIONIISPH_H
