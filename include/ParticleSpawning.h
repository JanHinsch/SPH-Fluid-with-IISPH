#ifndef PARTICLE_SPAWNING_H
#define PARTICLE_SPAWNING_H

#include "Particle.h"
#include <vector>
#include <SFML/Graphics.hpp>
#include <../external/eigen-master/Eigen/Dense>

class ParticleSpawning {
public:
    // Move boundary particles with a given offset
    static void moveBoundaryParticles(std::vector<Particle>& particles, float deltaX, float deltaY, bool switchDir);

    // Rotate boundary particles around a center point
    static void rotateBoundaryParticles(std::vector<Particle>& particles, float centerX, float centerY, float angle);

    // Add particles in a box shape
    static void addBox(int numParticles, std::vector<Particle>& particles);

    // Add particles in a staggered box shape
    static void addBox2(int x, int y, std::vector<Particle>& particles);

    // Add particles in an analyzable staggered box
    static void addBoxAnalyse(int x, int y, int numRow, std::vector<Particle>& particles);

    // Add straight-line boundary particles
    static void addBoundaries(int numParticles, int x, int y, std::vector<Particle>& particles);

    // Add vertical boundary particles
    static void addBoundaries2(int numParticles, int x, int y, std::vector<Particle>& particles);

    // Add moving boundaries
    static void addMovingBoundaries(int numParticles, int x, int y, std::vector<Particle>& particles);

    // Add boundaries with a specific angle
    static void addBoundariesWithAngle(int numParticles, int x, int y, float angle, std::vector<Particle>& particles);

    // Add half-circle boundary particles
    static void addBoundariesInHalfCircle(float h, int centerX, int centerY, float radius, float startAngle,
                                          std::vector<Particle>& particles);

    // Add particles in a rotating circle
    static void addRotatingCircle(float centerX, float centerY, float radius, float angularVelocity,
                                  std::vector<Particle>& particles);

    // Move particles in a rotating circle
    static void moveRotatingCircle(float centerX, float centerY, float radius, float angularVelocity, float timeStep,
                                   std::vector<Particle>& particles);

    // Add a rotating circle with attached triangles
    static void addRotatingCircleWithTriangles(float centerX, float centerY, float radius, float angularVelocity, float h,
                                               int numTriangles, float triangleHeight, std::vector<Particle>& particles);

    // Move particles in a rotating circle with attached triangles
    static void moveRotatingCircleWithTriangles(float centerX, float centerY, float radius, float angularVelocity,
                                                float timeStep, std::vector<Particle>& particles);

    // Remove overlapping particles
    static void removeOverlappingParticles(float minSpacing, std::vector<Particle>& particles);

    // Add a rotating circle with attached rectangles
    static void addRotatingCircleWithRectangle(float centerX, float centerY, float radius, float angularVelocity, float h,
                                               float rectWidth, float rectHeight, float rectOffset,
                                               std::vector<Particle>& particles);

    // Move particles in a rotating circle with attached rectangles
    static void moveRotatingCircleWithRectangle(float centerX, float centerY, float radius, float angularVelocity,
                                                float timeStep, float rectOffset, std::vector<Particle>& particles);
};

#endif // PARTICLE_SPAWNING_H
