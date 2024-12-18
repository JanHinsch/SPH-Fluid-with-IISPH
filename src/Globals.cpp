#include <SFML/Graphics.hpp>
#include "../include/Globals.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <SimulationIISPH_PressureBoundaries.h>

// screen size
const int x_size_screen = 1200;
const int y_size_screen = 846;

//////////////////// important globals //////////////////////////
float h = 2.0f;

float kernelSupport = 2.0 * h;

float stiffness_constant_k = 100;

float surfaceTensionFactor = 0.0f;

float particleRestVolume = std::pow(h, 2);

float density_rest = 0.01f;  // Example value

float viscosityFactor = 0.004f;  // Example value, adjust based on your fluid properties also called mü

sf::Vector2f gravity = sf::Vector2f(0.0f, 9.8f);

// adaptive timeStepping
bool adaptiveTimeStepping = false;

float timeStep = 0.007f;

const float lambda = 0.5f;

const float minTimeStep = 0.0008f;

const float maxTimeStep = 0.01;

bool EOS_Pressure = false;

bool IISPH_Pressure_Boundaries = true;

bool IISPH_Pressure_Extrapolation = false;

bool IISPH_Pressure_Mirroring = false;

bool IISPH_Pressure_Zero = false;

bool IISPH_MLS_Pressure_Extrapolation = false;

// when true print pressure color when false print velocity
bool pressureColors = true;

// as counter to how many fluid particles
int countFluidParticles = 0;

// as counter for IISPH PPE Solver Iterations
int currentIterations = 0;

// real Volume Error
float currentVolumeError;


// float speedThreshold1 = 0.0f;  // Blue to Cyan
// float speedThreshold2 = 0.0f;  // Cyan to Green

float speedThreshold1 = 20.0f;  // Blue to Cyan for velo
float speedThreshold2 = 40.0f;  // Cyan to Green for velo

// used in IISPH
float gamma_1 = 0.7f; // keep between 0.0f and 1.0f -> best is 0.7f -> used in RestVolumeBoundary calc for Pressure Boundaries
// used
float gamma_2 = 0.25f; // apparently often > 1.0f -> used in compute volumeboundary for pressureBoundaries

float gamma_3 = 1.0f; // used in pressure acc for mirroring as well as MLS Extrapolation and SPH Extrapolation

// used in IISPH updatePressure
float omega = 0.5f;

float beta = 0.15f * h * h;

bool CFLCondition = true; // true when condition is met

// for neighbor search
const float cellSize = (2 * h); // of grid

const int numCellsX = x_size_screen / cellSize;
const int numCellsY = y_size_screen / cellSize;

// diameter of particles
const float particleDiameter = 1.5f * h;
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
        if (IISPH_Pressure_Boundaries) {
            outFile << "countFluidParticles=" << 1500 << "\n";
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