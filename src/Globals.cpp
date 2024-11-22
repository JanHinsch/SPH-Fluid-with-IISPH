#include <SFML/Graphics.hpp>
#include "../include/Globals.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <SimulationIISPH.h>

// screen size
const int x_size_screen = 1200;
const int y_size_screen = 846;

//////////////////// important globals //////////////////////////
float h = 2.0f;

float kernelSupport = 2.0 * h;

float stiffness_constant_k = 90000;

float surfaceTensionFactor = 0.0f;

float particleRestVolume = std::pow(h, 2);

float density_rest = 0.001f;  // Example value

float viscosityFactor = 0.0f;  // Example value, adjust based on your fluid properties also called mü

sf::Vector2f gravity = sf::Vector2f(0.0f, 9.8f);

const float timeStep = 0.007f;

bool EOS_Pressure = false;

bool IISPH_Pressure = true;

// when true print pressure color when false print velocity
bool pressureColors = true;

// as counter to how many fluid particles
int countFluidParticles = 0;

// as counter for IISPH PPE Solver Iterations
int currentIterations = 0;


float speedThreshold1 = 20.0f;  // Blue to Cyan
float speedThreshold2 = 40.0f;  // Cyan to Green

// used in IISPH
float gamma_1 = 0.7f; // keep between 0.0f and 1.0f -> best is 0.7f
// unused
float gamma_2 = 1.0f; // apparently often > 1.0f

// used in IISPH updatePressure
float omega = 0.5f;

float beta = 0.15f * h * h;

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

// for moving boundaries
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
        if (EOS_Pressure) {
            outFile << "stiffness_constant_k=" << stiffness_constant_k << "\n";
            outFile << "viscosityFactor=" << viscosityFactor << "\n";
            outFile << "timeStep=" << timeStep << "\n";
            outFile << "rest_density=" << density_rest << "\n";
            outFile.close();
        }
        if (IISPH_Pressure) {
            outFile << "countFluidParticles=" << 500 << "\n";
            outFile << "timeStep=" << timeStep << "\n";
            outFile.close();
        }
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}

// Function to generate a unique filename for each frame
std::string generateFrameFilename(int frameNumber) {
    // Use a stringstream to format the filename with zero-padding
    std::ostringstream oss;

    // The directory name is "frames/" and the filename will be "frame_XXXX.png"
    // where XXXX is the frame number, zero-padded to 4 digits.
    oss << "../frames/frame_" << std::setw(4) << std::setfill('0') << frameNumber << ".png";

    // Return the formatted string
    return oss.str();
}