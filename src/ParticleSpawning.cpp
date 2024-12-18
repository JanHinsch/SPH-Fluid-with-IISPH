#include "../include/ParticleSpawning.h"
#include "../include/ParticleSystem.h"
#include "../include/Globals.h"
#include <../external/eigen-master/Eigen/Dense>



void ParticleSpawning::moveBoundaryParticles(std::vector<Particle>& particles, float deltaX, float deltaY, bool switchDir) {
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

void ParticleSpawning::rotateBoundaryParticles(std::vector<Particle>& particles, float centerX, float centerY, float angle) {
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
void ParticleSpawning::addBox(int numParticles, std::vector<Particle>& particles) {
    int xc = 0;
    int x = x_size_screen/2 - 18;
    int y = 650;

    for (auto& particle : particles) {
        if ((xc % 11) == 0) {
            y += 1.1f * h;
            x = x_size_screen/2 - 18;
            particle.position = sf::Vector2f(x, y);
        }
        x += 1.1 * h;
        particle.position = sf::Vector2f(x, y);
        xc++;
    }
}

void ParticleSpawning::addBox2(int x, int y, std::vector<Particle>& particles) {
    int xc = 0;
    int oldx = x;

    for (auto& particle : particles) {
        if (particle.isStatic == false) {
            particle.velocity = sf::Vector2f(0.0f, 0.0f);
            particle.acceleration = sf::Vector2f(0.0f, 0.0f);
            particle.pressure = 0.0f;
            particle.density = 0.0f;

            if ((xc % 50) == 0) {
                // y += 1.45f * h;
                y += h;
                // y += h * 0.5f;
                x = oldx + 1;
                particle.position = sf::Vector2f(x, y);
            }
            // x += 1.45f * h;
            x += h;
            // x += h * 0.5f;
            particle.position = sf::Vector2f(x, y);
            xc++;
        }
    }
}

void ParticleSpawning::addBoxAnalyse(int x, int y, int numRow, std::vector<Particle>& particles) {
    int xc = 0; // Counter for particles in a row
    int yc = 0; // Counter for rows (layers)
    int oldx = x; // Store the initial x position

    for (auto& particle : particles) {
        if (particle.isStatic == false) {
            // Reset particle properties
            particle.velocity = sf::Vector2f(0.0f, 0.0f);
            particle.acceleration = sf::Vector2f(0.0f, 0.0f);
            particle.pressure = 0.0f;
            particle.density = 0.0f;

            // Start a new layer when necessary
            if ((xc % numRow) == 0) {
                y += h; // Move to the next layer
                x = oldx; // Reset x to initial value
                xc = 0; // Reset particle count in the row

                // Apply staggered offset for every second layer
                if ((yc % 2) == 1) {
                    x -= static_cast<int>(0.5f * h); // Offset by 0.5 * h
                }
                yc++; // Increment row (layer) counter
            }

            // Set particle position and increment counters
            particle.position = sf::Vector2f(x, y);
            x += h; // Increment x for the next particle
            xc++; // Increment particle counter
        }
    }
}

void ParticleSpawning::addBoundaries(int numParticles, int x , int y, std::vector<Particle>& particles) {

    for (int i = 0; i < numParticles; i++) {
        Particle particle;
        particle.isStatic = true;
        particle.color = sf::Color(255,255,0,255);

        particle.position.x = x;
        particle.position.y = y;

        particles.push_back(particle);

        x += h;
    }
}

void ParticleSpawning::addBoundaries2(int numParticles, int x , int y, std::vector<Particle>& particles) {

    for (int i = 0; i < numParticles; i++) {
        Particle particle;
        particle.isStatic = true;
        particle.color = sf::Color(255,255,0,255);

        particle.position.x = x;
        particle.position.y = y;

        particles.push_back(particle);

        y += h;
    }
}

void ParticleSpawning::addMovingBoundaries(int numParticles, int x , int y, std::vector<Particle>& particles) {

    for (int i = 0; i < numParticles; i++) {
        Particle particle;
        particle.isStatic = true;
        particle.isMovableBoundary = true;
        particle.color = sf::Color(255,255,0,255);

        particle.position.x = x;
        particle.position.y = y;

        particles.push_back(particle);

        y += h;
    }
}

void ParticleSpawning::addBoundariesWithAngle(int numParticles, int x, int y, float angle, std::vector<Particle>& particles) {
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
        particles.push_back(particle);

        // Update currentX and currentY for the next particle
        currentX += dx;
        currentY += dy;
    }
}

void ParticleSpawning::addBoundariesInHalfCircle(float h, int centerX, int centerY, float radius, float startAngle, std::vector<Particle>& particles) {
    // Convert start angle to radians
    float startRad = startAngle * (M_PI / 180.0f);
    float endRad = (startAngle + 180.0f) * (M_PI / 180.0f);

    // Calculate the total arc length of the half-circle
    float arcLength = M_PI * radius; // Half-circle has 180 degrees or π radians

    // Calculate the number of particles based on desired spacing
    int numParticles = static_cast<int>(arcLength / h) + 1;

    // Calculate the angle step based on the number of particles
    float angleStep = (endRad - startRad) / (numParticles - 1);

    // Place each particle along the arc
    for (int i = 0; i < numParticles; i++) {
        // Calculate the angle for this particle
        float radianAngle = startRad + (i * angleStep);

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
        particles.push_back(particle);
    }
}

void ParticleSpawning::addRotatingCircle(float centerX, float centerY, float radius, float angularVelocity, std::vector<Particle>& particles) {
    // Calculate the number of particles based on spacing h
    int numParticles = static_cast<int>(std::round(2 * M_PI * radius / h));
    float angleIncrement = 2 * M_PI / numParticles; // Angular step between particles

    for (int i = 0; i < numParticles; i++) {
        float angle = i * angleIncrement;

        // Calculate particle position relative to center
        Particle particle;
        particle.isStatic = true;
        particle.isMovableBoundary = true;
        particle.color = sf::Color(255, 255, 0, 255);
        particle.position.x = centerX + radius * cos(angle); // X-coordinate
        particle.position.y = centerY + radius * sin(angle); // Y-coordinate

        // Set velocity for rotation, also relative to center
        particle.velocity.x = -angularVelocity * radius * sin(angle); // Perpendicular to radius
        particle.velocity.y = angularVelocity * radius * cos(angle);

        particles.push_back(particle);
    }
}

void ParticleSpawning::moveRotatingCircle(float centerX, float centerY, float radius, float angularVelocity, float timeStep, std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        if (particle.isMovableBoundary) {
            // Calculate current angle of the particle relative to the circle center
            float dx = particle.position.x - centerX;
            float dy = particle.position.y - centerY;
            float angle = atan2(dy, dx);

            // Increment angle based on angular velocity and timestep
            angle += angularVelocity * timeStep;

            // Update particle position to new location on the circle
            particle.position.x = centerX + radius * cos(angle);
            particle.position.y = centerY + radius * sin(angle);
        }
    }
}

void ParticleSpawning::addRotatingCircleWithTriangles(
    float centerX, float centerY, float radius, float angularVelocity, float h, int numTriangles, float triangleHeight, std::vector<Particle>& particles)
{
    // Add circle particles
    int numCircleParticles = static_cast<int>(std::round(2 * M_PI * radius / h));
    float angleIncrement = 2 * M_PI / numCircleParticles;

    for (int i = 0; i < numCircleParticles; i++) {
        float angle = i * angleIncrement;

        // Circle particle
        Particle particle;
        particle.isStatic = true;
        particle.isMovableBoundary = true;
        particle.color = sf::Color(255, 255, 0, 255);
        particle.position.x = centerX + radius * cos(angle);
        particle.position.y = centerY + radius * sin(angle);

        // Velocity for rotation
        particle.velocity.x = -angularVelocity * radius * sin(angle);
        particle.velocity.y = angularVelocity * radius * cos(angle);

        particles.push_back(particle);
    }

    // Add triangles
    float triangleAngleIncrement = 2 * M_PI / numTriangles;
    for (int t = 0; t < numTriangles; t++) {
        float baseAngle = t * triangleAngleIncrement;

        // Base vertices on the circle
        float baseX1 = centerX + radius * cos(baseAngle);
        float baseY1 = centerY + radius * sin(baseAngle);
        float baseX2 = centerX + radius * cos(baseAngle + triangleAngleIncrement / 2);
        float baseY2 = centerY + radius * sin(baseAngle + triangleAngleIncrement / 2);

        // Apex of the triangle (pointing outward)
        float apexX = centerX + (radius + triangleHeight) * cos(baseAngle + triangleAngleIncrement / 4);
        float apexY = centerY + (radius + triangleHeight) * sin(baseAngle + triangleAngleIncrement / 4);

        // Add particles along the base
        int numBaseParticles = static_cast<int>(std::round(hypot(baseX2 - baseX1, baseY2 - baseY1) / h));
        for (int i = 0; i <= numBaseParticles; i++) {
            float interp = static_cast<float>(i) / numBaseParticles;
            Particle particle;
            particle.isStatic = true;
            particle.isMovableBoundary = true;
            particle.color = sf::Color(255, 0, 0, 255);
            particle.position.x = baseX1 + interp * (baseX2 - baseX1);
            particle.position.y = baseY1 + interp * (baseY2 - baseY1);

            particles.push_back(particle);
        }

        // // Add particles along the edges (base to apex)
        // int numEdgeParticles = static_cast<int>(std::round(triangleHeight / h));
        // for (int i = 0; i <= numEdgeParticles; i++) {
        //     float interp = static_cast<float>(i) / numEdgeParticles;
        //
        //     // Left edge
        //     Particle leftEdgeParticle;
        //     leftEdgeParticle.isStatic = true;
        //     leftEdgeParticle.isMovableBoundary = true;
        //     leftEdgeParticle.color = sf::Color(255, 0, 0, 255);
        //     leftEdgeParticle.position.x = baseX1 + interp * (apexX - baseX1);
        //     leftEdgeParticle.position.y = baseY1 + interp * (apexY - baseY1);
        //     particles.push_back(leftEdgeParticle);
        //
        //     // Right edge
        //     Particle rightEdgeParticle;
        //     rightEdgeParticle.isStatic = true;
        //     rightEdgeParticle.isMovableBoundary = true;
        //     rightEdgeParticle.color = sf::Color(255, 0, 0, 255);
        //     rightEdgeParticle.position.x = baseX2 + interp * (apexX - baseX2);
        //     rightEdgeParticle.position.y = baseY2 + interp * (apexY - baseY2);
        //     particles.push_back(rightEdgeParticle);
        // }
    }
    // Handle overlapping particles
    // removeOverlappingParticles(0.9f * h);
}

void ParticleSpawning::moveRotatingCircleWithTriangles(
    float centerX, float centerY, float radius, float angularVelocity, float timeStep, std::vector<Particle>& particles)
{
    for (auto& particle : particles) {
        if (particle.isMovableBoundary) {
            // Calculate current angle relative to center
            float dx = particle.position.x - centerX;
            float dy = particle.position.y - centerY;
            float angle = atan2(dy, dx);

            // Increment angle based on angular velocity
            angle += angularVelocity * timeStep;

            // Calculate distance to center (maintains radius for triangles)
            float distance = hypot(dx, dy);

            // Update position
            particle.position.x = centerX + distance * cos(angle);
            particle.position.y = centerY + distance * sin(angle);
        }
    }
}

// Function to remove overlapping particles
void ParticleSpawning::removeOverlappingParticles(float minSpacing, std::vector<Particle>& particles) {
    // for (size_t i = 0; i < particles.size(); ++i) {
    //     for (size_t j = i + 1; j < particles.size(); ++j) {
    //         if (!particles[j].isMovableBoundary) continue;
    //         // Calculate distance between particles
    //         float dx = particles[i].position.x - particles[j].position.x;
    //         float dy = particles[i].position.y - particles[j].position.y;
    //         float distance = std::sqrt(dx * dx + dy * dy);
    //
    //         // If particles are too close, remove one
    //         if (distance <= minSpacing) {
    //             particles.erase(particles.begin() + j);
    //             --j; // Adjust index after erasing
    //         }
    //     }
    // }
    particles.erase(
    std::remove_if(particles.begin(), particles.end(),
        [](const Particle& particle) {
            return particle.isMovableBoundary && std::size(particle.particle_neighbours) > 3;
        }),
    particles.end());
}

void ParticleSpawning::addRotatingCircleWithRectangle(
    float centerX, float centerY, float radius, float angularVelocity, float h,
    float rectWidth, float rectHeight, float rectOffset, std::vector<Particle>& particles)
{
    // Add circle particles
    int numCircleParticles = static_cast<int>(std::round(2 * M_PI * radius / h));
    float angleIncrement = 2 * M_PI / numCircleParticles;

    for (int i = 0; i < numCircleParticles; i++) {
        float angle = i * angleIncrement;

        // Circle particle
        Particle particle;
        particle.isStatic = true;
        particle.isMovableBoundary = true;
        particle.color = sf::Color(255, 255, 0, 255);
        particle.position.x = centerX + radius * cos(angle);
        particle.position.y = centerY + radius * sin(angle);

        // Velocity for rotation
        particle.velocity.x = -angularVelocity * radius * sin(angle);
        particle.velocity.y = angularVelocity * radius * cos(angle);

        particles.push_back(particle);
    }

    // Function to add a rectangle
    auto addRectangle = [&](float angleOffset) {
        // Rectangle's center offset from the circle's center along a radial line
        float rectCenterX = centerX + (radius + rectOffset) * cos(angleOffset);
        float rectCenterY = centerY + (radius + rectOffset) * sin(angleOffset);

        // Rectangle vertices relative to its center
        Eigen::Vector2f topLeft(-rectWidth / 2, rectHeight / 2);
        Eigen::Vector2f topRight(rectWidth / 2, rectHeight / 2);
        Eigen::Vector2f bottomLeft(-rectWidth / 2, -rectHeight / 2);
        Eigen::Vector2f bottomRight(rectWidth / 2, -rectHeight / 2);

        // Place particles along rectangle edges
        auto addEdgeParticles = [&](Eigen::Vector2f start, Eigen::Vector2f end) {
            float edgeLength = (end - start).norm();
            int numEdgeParticles = static_cast<int>(std::round(edgeLength / h));
            for (int i = 0; i <= numEdgeParticles; i++) {
                float interp = static_cast<float>(i) / numEdgeParticles;
                Eigen::Vector2f position = start + interp * (end - start);

                Particle particle;
                particle.isStatic = true;
                particle.isMovableBoundary = true;
                particle.color = sf::Color(255, 255, 0, 255);
                particle.position.x = rectCenterX + position.x();
                particle.position.y = rectCenterY + position.y();

                particles.push_back(particle);
            }
        };

        // Add edges
        addEdgeParticles(topLeft, topRight);    // Top edge
        addEdgeParticles(topRight, bottomRight); // Right edge
        addEdgeParticles(bottomRight, bottomLeft); // Bottom edge
        addEdgeParticles(bottomLeft, topLeft);    // Left edge
    };

    // Add two rectangles
    addRectangle(0);                  // First rectangle (angle = 0 radians)
    addRectangle(M_PI);               // Second rectangle (opposite side, angle = π radians)

    removeOverlappingParticles(h,particles);
}

void ParticleSpawning::moveRotatingCircleWithRectangle(
    float centerX, float centerY, float radius, float angularVelocity, float timeStep, float rectOffset, std::vector<Particle>& particles)
{
    for (auto& particle : particles) {
        if (particle.isMovableBoundary) {
            // Calculate current angle relative to the circle center
            float dx = particle.position.x - centerX;
            float dy = particle.position.y - centerY;
            float angle = atan2(dy, dx);

            // Increment angle based on angular velocity
            angle += angularVelocity * timeStep;

            // Calculate distance to the circle center
            float distance = hypot(dx, dy);

            // Update position
            particle.position.x = centerX + distance * cos(angle);
            particle.position.y = centerY + distance * sin(angle);
        }
    }
}




//////////////////////////////////////////////////////////////////

// void ParticleSpawning::checkMouseHover(const sf::Vector2i& mousePosition) {
//     sf::Vector2f mousePos(static_cast<float>(mousePosition.x), static_cast<float>(mousePosition.y));

//     for(auto& particle : particles) {
//         if (abs(mousePos.x - particle.position.x) < particleRadius &&
//             abs(mousePos.y - particle.position.y) < particleRadius) {
//             particle.color = sf::Color(255, 0, 0, 255);  // Change to red on hover
//             iterateNeighbours(particle);
//         } else {
//             //particle.color = sf::Color(255, 255, 0, 255);  // Reset color if not hovering
//         }
//     }
// }

// void ParticleSystem::checkMouseHover(const sf::Vector2i& mousePosition) {
//     sf::Vector2f mousePos(static_cast<float>(mousePosition.x), static_cast<float>(mousePosition.y));
//     int mouseCellIndex = getCellIndex(mousePos);
//
//     // TODO delete later
//     std::vector<Particle> neighbours;
//
//     // Reset color of all particles
//     for (auto& particle : particles) {
//         //particle.color = sf::Color(255, 255, 0, 255);  // Reset color
//     }
//
//     // Check particles in the cell where the mouse pointer is located
//     if (C.isValid(mousePos.x / cellSize, mousePos.y / cellSize)) {
//         int mouseCellCoordX = mousePos.x / cellSize;
//         int mouseCellCoordY = mousePos.y / cellSize;
//         int numParticlesInCell = C.get(mouseCellCoordX, mouseCellCoordY + 1) - C.get(mouseCellCoordX, mouseCellCoordY);
//
//         for (int i = 0; i < numParticlesInCell; ++i) {
//             int particleIndex = C.get(mouseCellCoordX, mouseCellCoordY) + i;
//             if (inRadius(mousePos, L[particleIndex]->position, particleRadius)) {
//                 L[particleIndex]->color = sf::Color(255, 0, 0, 255);  // Change to red on hover
//                 iterateNeighbours(*L[particleIndex], neighbours);
//
//                 for (auto& p : neighbours) {
//                     p.color = sf::Color(255, 255, 255, 255);
//                 }
//                 break;  // Only change color of the hovered particle
//             }
//         }
//
//     }
// }