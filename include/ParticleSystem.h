// ParticleSystem.h
#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include "Particle.h"
#include "Grid.h"
#include "Globals.h"
#include <vector>
#include <SFML/Graphics.hpp>
#include <list>
#include <unordered_map>

class ParticleSystem {
public:
    ParticleSystem(unsigned int count);
    std::vector<Particle> m_particles;

    // Grid to manage cells
    Grid C;
    std::vector<Particle*> L; // the L array as in lecture

    // Helper functions
    int getCellIndex(const sf::Vector2f& position) const;
    int getCellCoords(int x, int y);
    void accumulateCounters();

    void iterateNeighbours(Particle& particle, std::vector<Particle*>& m_neighbours);

    static static bool inRadius(sf::Vector2f positionA, sf::Vector2f positionB, float radius=kernelSupport);

    static static void moveBoundaryParticles(std::vector<Particle>& particles, float deltaX, float deltaY, bool switchDir); // delete!!!

    static void rotateBoundaryParticles(std::vector<Particle> &particles,float centerX, float centerY, float angle);

    void addBox(int numParticles);

    void addBox2(int x, int y);

    void addBoxAnalyse(int x, int y);

    void addBoundaries(int numParticles, int x , int y);

    void addBoundaries2(int numParticles, int x, int y);

    void addMovingBoundaries(int numParticles, int x, int y);

    void addBoundariesWithAngle(int numParticles, int x, int y, float angle);

    void addBoundariesInHalfCircle(float h, int centerX, int centerY, float radius, float startAngle);

    void addBoundaries3(int numParticles, int x, int y);

    std::vector<Particle*> generateSortedList(std::vector<Particle>& particles);

    // Set start values of simulation (place particles etc.)
    void initiateParticles(std::vector<Particle>& particles);

    void updateParticlesEOS(std::vector<Particle>& particles, int x_size_screen, int y_size_screen);

    void updateParticlesIISPH(std::vector<Particle>& particles, int x_size_screen, int y_size_screen);

    void resetSimulation();

    void quadraticNeighbourSearch(Particle &particle, std::vector<Particle> &m_neighbours);

    // mouse hover
    void checkMouseHover(const sf::Vector2i& mousePosition);

};

#endif // PARTICLE_SYSTEM_H
