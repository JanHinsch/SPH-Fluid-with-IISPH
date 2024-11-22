#ifndef PARTICLE_H
#define PARTICLE_H

#include <SFML/System/Vector2.hpp>
#include <SFML/Graphics.hpp>

struct Particle {

    float mass;
    sf::Vector2f velocity;
    sf::Vector2f position;
    sf::Vector2f acceleration;

    // pressure boundary
    float restVolume;
    float volume;
    float Ap;

    float density;
    float pressure;

    bool isStatic = false;
    bool isMovableBoundary = false;

    sf::Color color = sf::Color(255, 255, 0, 255); // default yellow

    std::vector<Particle*> particle_neighbours;
    sf::Vector2f predicted_velocity;
    float sourceTerm;
    float diagonalElement;
    float laplacian;
};

#endif // PARTICLE_H