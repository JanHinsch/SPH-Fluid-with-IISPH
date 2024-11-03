#include <SFML/Graphics.hpp>
#include "../include/Globals.h"
#include <cmath>
#include <fstream>
#include <iostream>

// screen size
const int x_size_screen = 1200;
const int y_size_screen = 846;

//////////////////// important globals //////////////////////////
float h = 2.0;

float kernelSupport = 2.0 * h;

float stiffness_constant_k = 200;

float surfaceTensionFactor = 20.0f;

float density_rest = 0.1f;  // Example value

float viscosityFactor = 10.4f;  // Example value, adjust based on your fluid properties also called mü

sf::Vector2f gravity = sf::Vector2f(0.0f, 9.8f);

const float timeStep = 0.01f;

// unused 
float gamma_1 = 1.0f; // apparently often between 0.0f and 1.0f
// unused
float gamma_2 = 1.0f; // apparently often > 1.0f

bool CFLCondition = true; // true when condition is met

// for neighbor search
const float cellSize = (2 * h); // of grid

const int numCellsX = x_size_screen / cellSize;
const int numCellsY = y_size_screen / cellSize;

// diameter of particles
const float particleDiameter = h;
const float particleRadius = particleDiameter / 2;

// for pausing simulation
bool simulationPaused = false;

// for moving boundries
bool switchDir = true;


// Function to calculate the Euclidean distance between two sf::Vector2f points
float euclideanDistance(const sf::Vector2f& positionA, const sf::Vector2f& positionB) {
    float dx = positionB.x - positionA.x;
    float dy = positionB.y - positionA.y;
    return sqrt(dx * dx + dy * dy);
}

void writeGlobalsToFile(const std::string& filename) {
    std::ofstream outFile(filename);
    if (outFile.is_open()) {
        outFile << "stiffness_constant_k=" << stiffness_constant_k << "\n";
        outFile << "viscosityFactor=" << viscosityFactor << "\n";
        outFile << "timeStep=" << timeStep << "\n";
        outFile << "rest_density=" << density_rest << "\n";
        outFile.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}