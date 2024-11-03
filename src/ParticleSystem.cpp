// ParticleSystem.cpp
#include "../include/ParticleSystem.h"
#include "../include/Globals.h"
#include "../include/Grid.h"
#include "../include/Simulation.h"
#include "../include/UIManager.h"
#include <cstdlib>
#include <cmath>
#include "ParticleSystem.h"

ParticleSystem::ParticleSystem(unsigned int count) : m_particles(count), C(numCellsX, numCellsY), L(m_particles.size() +
                                                                                                     200 * 6 +
                                                                                                     30 * 1 +
                                                                                                     80 * 4 +
                                                                                                     113 * 2 +
                                                                                                     400 * 2
                                                                                                     ) {}

////////////////////// Simulation Part ////////////////////////////

void ParticleSystem::initiateParticles(std::vector<Particle>& particles) {

    // // write important globals to file
    writeGlobalsToFile("globals.txt");

    addBox2(440, 600);
    // addBox(100);
    //resetSimulation();

    // bottom Boundary
    addBoundaries(200, 320, 662); // short bottom
    addBoundaries(200, 320, 664); // short bottom

    // left boundary
    addBoundaries2(200, 320, 264);
    addBoundaries2(200, 318, 264);


    // right boundary
    addBoundaries2(200, 718, 264);
    addBoundaries2(200, 720, 264);


    // top boundary
    // addBoundaries(128, 448, 404);

    // moving boundaries
    // addMovingBoundaries(30, 320, 452);

    // labyrinth
    // addBoundariesWithAngle(150, 320, 314, 22.0f);
    // addBoundaries(93, 320, 389);

    // addBoundariesWithAngle(150, 718, 414, 158.0f); 
    // addBoundaries(93, 534, 389);

    // triangle
    addBoundariesWithAngle(80, 520, 414, 45.0f);
    addBoundariesWithAngle(80, 520, 414, 135.0f);
    addBoundaries(113, 408, 527);

    addBoundariesWithAngle(80, 520, 416, 45.0f);
    addBoundariesWithAngle(80, 520, 416, 135.0f);
    addBoundaries(113, 408, 525);

    // addBoundariesWithAngle(20, 520, 484, 45.0f);
    // addBoundariesWithAngle(20, 520, 484, 135.0f);
    // addBoundaries(29, 492, 482);

    // addBoundariesWithAngle(20, 420, 484, 45.0f);
    // addBoundariesWithAngle(20, 420, 484, 135.0f);
    // addBoundaries(29, 392, 512);

    // addBoundariesWithAngle(20, 620, 484, 45.0f);
    // addBoundariesWithAngle(20, 620, 484, 135.0f);
    // addBoundaries(29, 592, 512);


    addBoundariesInHalfCircle(400,519, 446, 199.0f, 0.0f);
    addBoundariesInHalfCircle(400,519, 448, 199.0f, 0.0f);


    // addBoundaries(40, 480, 454);
    
    // addBoundaries2(20, 439, 414);
    // addBoundaries2(20, 480, 454);
    // addBoundaries2(20, 560, 454);
    // addBoundaries2(20, 600, 414);


    // addBoundariesWithAngle(65, 716 + h, 364, 157.0f);
    // addBoundariesWithAngle(65, 320, 364, 22.0f);

    ///////////////////////////// Analyse /////////////////////

    // addBoundaries(15, 318, 662); // short bottom
    // addBoundaries(15, 318, 664); // short bottom
    // // left boundary
    // addBoundaries2(200, 320, 264);
    // addBoundaries2(199, 318, 264);
    // // right boundary
    // addBoundaries2(200, 344, 264);
    // addBoundaries2(199, 346, 264);

    // addBoxAnalyse(320, 460);


    //////////////////////////////////////////////////////////

    for (auto& particle : particles) {
        // Set uniform mass for all particles
        float particleVolume = std::pow(h, 2);
        
        float particleMass = particleVolume * density_rest;

        particle.mass = particleMass;
        particle.density = 0.0f;

        //printf("Volume %f", particleVolume);
        //printf("Mass %f", particleMass);
    }


}

// int countDir = 0; // for moveable boundary 

void ParticleSystem::updateParticles(std::vector<Particle>& particles, int x_size_screen, int y_size_screen) {

    /////////// Index search
    // Calculate the total number of cells in the grid and initialize C
    C.clear();

    // Index Search
    // for neighbor search assign particle to a cell
    for (auto& p : particles) {        
        int cellIndex = getCellIndex(p.position);
        C[cellIndex]++;
    }
    L = generateSortedList(particles);
    //////////////////////

    std::vector<Particle> neighbours;

    // compute density and pressure
    for (auto& particle : particles) {
        // if (particle.isStatic) continue;
        neighbours.clear();

        iterateNeighbours(particle, neighbours);

        // quadraticNeighbourSearch(particle, neighbours);
        
        particle.density = SPHComputations::computeDensity(particle, neighbours);
        if (SPHComputations::isParticleCompressed(particle.density)) {
            //printf("denisty: %f", particle.density);
            particle.pressure = SPHComputations::computePressure(particle.density);
        } else {
            particle.pressure = 0;
        }
        // if (particle.isStatic == false )printf("density_i: %f", particle.density);
    }

    // compute accelerations and remove out of bounds
    int it = 0;
    for (auto& particle : particles) {
        if (particle.isStatic) continue;
        neighbours.clear();

        iterateNeighbours(particle, neighbours);

        // quadraticNeighbourSearch(particle, neighbours);

        particle.acceleration = SPHComputations::computeTotalAcceleration(particle, neighbours);

        if (particle.position.x < 20 || particle.position.x > x_size_screen - 40 || particle.position.y < 20 || particle.position.y > y_size_screen - 40) {
            m_particles.erase(m_particles.begin() + it);
        }
        it++;
    }

    SPHComputations::advectParticles(particles);

    SPHComputations::isCFLConditionTrue(particles);


    // Compute and write the average density to the file
    // std::string densityFile = "density_values_t_groß.txt";
    // float avgDensity = SPHComputations::getAvgDensity(particles, densityFile);

    // float cfl_number = SPHComputations::getCourantNumber(particles);
    // printf("Courant Number: %f \n", cf  l_number);

    // for moving boundary
    // if (countDir % 500 == 0) {
    //     switchDir = !switchDir;
    // }
    // countDir++;
    // moveBoundaryParticles(particles, 5.5, 5.5, switchDir);
    // rotateBoundaryParticles(particles, 320, 482, 0.9f);

}


void ParticleSystem::resetSimulation() {
    
    addBox2(440, 155);
    printf("Pressed Reset\n");
}

// quadratic search for test
void ParticleSystem::quadraticNeighbourSearch(Particle& particle, std::vector<Particle>& m_neighbours) {
    for (auto& neighbour : m_particles) {
        if (inRadius(particle.position, neighbour.position)) {
            m_neighbours.push_back(neighbour);
        }
    }
}

//////////////////////// Neighbour Search ///////////////////////

// helper function to get Cell Index for Particles (get 1D Index from particles coords)
int ParticleSystem::getCellIndex(const sf::Vector2f& position) const {
    int xIndex = static_cast<int>(position.x) / cellSize;
    int yIndex = static_cast<int>(position.y) / cellSize;
    return xIndex + yIndex * numCellsX; // Assuming cell indices are row-major
}

// get Cell 1D Index from (x, y) grid coords
int ParticleSystem::getCellCoords(int x, int y) {
    return y * numCellsY + x;
}

// helper function accumulate counter for C array
void ParticleSystem::accumulateCounters(){
    int accCounter = 0;
    for (int y = 0; y < numCellsY; ++y) {
        for (int x = 0; x < numCellsX; ++x) {
            accCounter += C.get(x, y);
            C.set(x, y, accCounter);
        }
    }
}

// get sorted List of particles L (as in lecture)
std::vector<Particle*> ParticleSystem::generateSortedList(std::vector<Particle>& particles) {

    // accum C
    accumulateCounters();

    for (int i = 0; i < m_particles.size(); ++i) {
        int x = m_particles[i].position.x;
        int y =  m_particles[i].position.y;
        int cellIndex = getCellIndex(m_particles[i].position);
        C[cellIndex]--; // -1 an der Stelle für C
        L[C[cellIndex]] = &m_particles[i];
    }
    return L;
}

void ParticleSystem::iterateNeighbours(Particle& particle, std::vector<Particle>& m_neighbours) {
    int currCell = getCellIndex(particle.position);

    std::pair<int, int> currCellCoord = C.getCoordinatesFromIndex(currCell);

    // Define the neighbor offsets 3x3 grid
    std::vector<std::pair<int, int>> neighborOffsets = {
        { -1, -1 }, { 0, -1 }, { 1, -1 },
        { -1, 0 }, { 0, 0 }, { 1, 0 },
        { -1, 1 }, { 0, 1 }, { 1, 1 }
    };
    
    for (const auto& offset : neighborOffsets) {
        std::pair<int, int> neighborCell = { currCellCoord.first + offset.first, currCellCoord.second + offset.second };

        // Ensure neighborCell is within bounds
        if (C.isValid(neighborCell.first, neighborCell.second)) {
            int neighborCellIndex = neighborCell.second * numCellsX + neighborCell.first;
            int numParticlesInCell = C[neighborCellIndex + 1] - C[neighborCellIndex];
            for (int i = 0; i < numParticlesInCell; ++i) {
                // Access the particles in this cell
                int particleIndex = C[neighborCellIndex] + i;

                //L[particleIndex]->color = sf::Color(140, 0, 255, 255); // Purple color for neighbors visual test

                // get real neighbours in radius of particle (for kernel)
                if (inRadius(particle.position, L[particleIndex]->position)) {

                    //L[particleIndex]->color = sf::Color(0, 255, 255, 255); // visual test

                    m_neighbours.push_back(*L[particleIndex]); // add to neighbour list
                }
            }
        }
    }
}

// helperfunction to get particles inside a specified radius (in Abhängigkeit von global h) uses squared distance 
// (radius = kernelSupport)
bool ParticleSystem::inRadius(sf::Vector2f positionA, sf::Vector2f positionB, float radius) {
    float dx = positionA.x - positionB.x;
    float dy = positionA.y - positionB.y;
    float distance = dx * dx + dy * dy;

    return distance < (radius * radius);
}

// bool ParticleSystem::inRadius(sf::Vector2f positionA, sf::Vector2f positionB, float radius) {
//     float distance = euclideanDistance(positionA, positionB);

//     return distance < (radius);
// }

///////////////////////////////////////////////////////////////////////

void ParticleSystem::moveBoundaryParticles(std::vector<Particle>& particles, float deltaX, float deltaY, bool switchDir) {
    for (auto& particle : particles) {
        if (particle.isMovableBoundary) {
            if (switchDir == true) {
                particle.position.x += deltaX * timeStep;
            } else {
                particle.position.x -= deltaX * timeStep;
            }
        }
    }
}

void ParticleSystem::rotateBoundaryParticles(std::vector<Particle>& particles, float centerX, float centerY, float angle) {
    float radian = angle * M_PI / 180.0f; // Convert angle to radians
    float cosTheta = std::cos(radian);
    float sinTheta = std::sin(radian);
    
    for (auto& particle : particles) {
        if (particle.isMovableBoundary) {
            float dx = particle.position.x - centerX;
            float dy = particle.position.y - centerY;
            
            float newX = centerX + (dx * cosTheta - dy * sinTheta);
            float newY = centerY + (dx * sinTheta + dy * cosTheta);
            
            particle.position.x = newX;
            particle.position.y = newY;
        }
    }
}

// Place Particles in a Box Shape (mostly used for testing)
void ParticleSystem::addBox(int numParticles) {
    int xc = 0;
    int x = x_size_screen/2 - 18;
    int y = 650;

    for (auto& particle : m_particles) {
        if ((xc % 11) == 0) {
            y += 1.1 * h;
            x = x_size_screen/2 - 18;
            particle.position = sf::Vector2f(x, y);
        }
        x += 1.1 * h;
        particle.position = sf::Vector2f(x, y);
        xc++;
    }
}

// void ParticleSystem::addBox(int numParticles) {

//     int start = 100;
//     int end = h * 10 + 100;
//     int particleIndex = 0;

//     for (int x = start; x < end; x += h) {
//         for (int y = start; y < end; y += h) {
//             if (particleIndex >= m_particles.size()) {
//                 return; // Stop if we've assigned all available particles
//             }
//             m_particles[particleIndex].position.x = x;
//             m_particles[particleIndex].position.y = y;

//             printf("Particle %d: Position (%d, %d)\n", particleIndex, x, y);


//             particleIndex++;
//         }
//     }
// }


void ParticleSystem::addBox2(int x, int y) {
    int xc = 0;
    int oldx = x;

    for (auto& particle : m_particles) {
        if (particle.isStatic == false) {
            particle.velocity = sf::Vector2f(0.0f, 0.0f);
            particle.acceleration = sf::Vector2f(0.0f, 0.0f);
            particle.pressure = 0.0f;
            particle.density = 0.0f;

            if ((xc % 80) == 0) {
                y += 1.45f * h;
                x = oldx + 1;
                particle.position = sf::Vector2f(x, y);
            }
            x += 1.45f * h;
            particle.position = sf::Vector2f(x, y);
            xc++;
        }
    }
}

void ParticleSystem::addBoxAnalyse(int x, int y) {
    int xc = 0;
    int oldx = x;

    for (auto& particle : m_particles) {
        if (particle.isStatic == false) {
            particle.velocity = sf::Vector2f(0.0f, 0.0f);
            particle.acceleration = sf::Vector2f(0.0f, 0.0f);
            particle.pressure = 0.0f;
            particle.density = 0.0f;

            if ((xc % 10) == 0) {
                y += h;
                x = oldx + 1;
                particle.position = sf::Vector2f(x, y);
            }
            x += h;
            particle.position = sf::Vector2f(x, y);
            xc++;
        }
    }
}

void ParticleSystem::addBoundaries(int numParticles, int x , int y) {

    for (int i = 0; i < numParticles; i++) {
        Particle particle;
        particle.isStatic = true;
        particle.color = sf::Color(255,255,0,255);

        particle.position.x = x;
        particle.position.y = y;

        m_particles.push_back(particle);

        x += h;
    }
}

void ParticleSystem::addBoundaries2(int numParticles, int x , int y) {

    for (int i = 0; i < numParticles; i++) {
        Particle particle;
        particle.isStatic = true;
        particle.color = sf::Color(255,255,0,255);

        particle.position.x = x;
        particle.position.y = y;

        m_particles.push_back(particle);

        y += h;
    }
}

void ParticleSystem::addMovingBoundaries(int numParticles, int x , int y) {
    
    for (int i = 0; i < numParticles; i++) {
        Particle particle;
        particle.isStatic = true;
        particle.isMovableBoundary = true;
        particle.color = sf::Color(255,255,0,255);

        particle.position.x = x;
        particle.position.y = y;

        m_particles.push_back(particle);

        y += h;
    }
}

void ParticleSystem::addBoundariesWithAngle(int numParticles, int x, int y, float angle) {
    // Convert angle to radians
    float radianAngle = angle * (M_PI / 180.0f);

    // Calculate the direction vector
    float dx = cos(radianAngle);
    float dy = sin(radianAngle);

    // Scale the direction vector to ensure particles are h units apart
    dx *= h;
    dy *= h;

    // Use floating-point for precise position calculations
    float currentX = static_cast<float>(x);
    float currentY = static_cast<float>(y);

    for (int i = 0; i < numParticles; i++) {
        Particle particle;
        particle.isStatic = true;
        particle.color = sf::Color(255, 255, 0, 255);

        // Set the particle position using the rounded coordinates
        particle.position.x = currentX;
        particle.position.y = currentY;

        // Add the particle to the system
        m_particles.push_back(particle);

        // Update currentX and currentY for the next particle
        currentX += dx;
        currentY += dy;
    }
}

void ParticleSystem::addBoundariesInHalfCircle(int numParticles, int centerX, int centerY, float radius, float startAngle) {
    // Calculate the angle step based on the number of particles
    float angleStep = 180.0f / (numParticles - 1);

    // Place each particle along the arc of the half-circle
    for (int i = 0; i < numParticles; i++) {
        // Calculate the angle for this particle
        float angle = startAngle + (i * angleStep);

        // Convert angle to radians
        float radianAngle = angle * (M_PI / 180.0f);

        // Calculate the position of the particle
        int x = centerX + static_cast<int>(radius * cos(radianAngle));
        int y = centerY + static_cast<int>(radius * sin(radianAngle));

        // Create and configure the particle
        Particle particle;
        particle.isStatic = true;
        particle.color = sf::Color(255, 255, 0, 255);
        particle.position.x = x;
        particle.position.y = y;

        // Add the particle to the system
        m_particles.push_back(particle);
    }
}


























//////////////////////////////////////////////////////////////////

// void ParticleSystem::checkMouseHover(const sf::Vector2i& mousePosition) {
//     sf::Vector2f mousePos(static_cast<float>(mousePosition.x), static_cast<float>(mousePosition.y));

//     for(auto& particle : m_particles) {
//         if (abs(mousePos.x - particle.position.x) < particleRadius &&
//             abs(mousePos.y - particle.position.y) < particleRadius) {
//             particle.color = sf::Color(255, 0, 0, 255);  // Change to red on hover
//             iterateNeighbours(particle);
//         } else {
//             //particle.color = sf::Color(255, 255, 0, 255);  // Reset color if not hovering
//         }
//     }
// }

void ParticleSystem::checkMouseHover(const sf::Vector2i& mousePosition) {
    sf::Vector2f mousePos(static_cast<float>(mousePosition.x), static_cast<float>(mousePosition.y));
    int mouseCellIndex = getCellIndex(mousePos);

    // TODO delete later
    std::vector<Particle> neighbours;

    // Reset color of all particles
    for (auto& particle : m_particles) {
        //particle.color = sf::Color(255, 255, 0, 255);  // Reset color
    }

    // Check particles in the cell where the mouse pointer is located
    if (C.isValid(mousePos.x / cellSize, mousePos.y / cellSize)) {
        int mouseCellCoordX = mousePos.x / cellSize;
        int mouseCellCoordY = mousePos.y / cellSize;
        int numParticlesInCell = C.get(mouseCellCoordX, mouseCellCoordY + 1) - C.get(mouseCellCoordX, mouseCellCoordY);

        for (int i = 0; i < numParticlesInCell; ++i) {
            int particleIndex = C.get(mouseCellCoordX, mouseCellCoordY) + i;
            if (inRadius(mousePos, L[particleIndex]->position, particleRadius)) {
                L[particleIndex]->color = sf::Color(255, 0, 0, 255);  // Change to red on hover
                iterateNeighbours(*L[particleIndex], neighbours);

                for (auto& p : neighbours) {
                    p.color = sf::Color(255, 255, 255, 255);
                }
                break;  // Only change color of the hovered particle
            }
        }

    }
}