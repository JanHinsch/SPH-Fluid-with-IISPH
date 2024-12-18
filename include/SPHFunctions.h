#ifndef SPHFUNCTIONS_H
#define SPHFUNCTIONS_H

#include "Particle.h"
#include <vector>
#include <SFML/Graphics.hpp>


class SPHFunctions {
public:
    static float kernelAlpha();

    static float kernel(sf::Vector2f positionA, sf::Vector2f positionB);

    static sf::Vector2f kernelGradient(sf::Vector2f positionA, sf::Vector2f positionB);

    static void advectParticles(std::vector<Particle>& particles);



    //////////////
    static float getMaxVelocity(std::vector<Particle>& particles);

    static float getMaxVelocityPressure(std::vector<Particle> &particles);

    static float getAvgDensity(std::vector<Particle> &particles, const std::string& filename);

    static void isCFLConditionTrue(std::vector<Particle>& particles);

    static float getCourantNumber(std::vector<Particle> &particles);

    static float computeAdaptiveTimeStep(float adaptiveTimeStep);

    static void writeCurrentIterationToFile(const std::string& filename);


private:

};



#endif //SPHFUNCTIONS_H
