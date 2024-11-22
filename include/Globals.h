// Globals.h

#ifndef GLOBALS_H
#define GLOBALS_H

extern const int x_size_screen;
extern const int y_size_screen;

extern float h;

extern float kernelSupport;

extern float stiffness_constant_k;

extern float surfaceTensionFactor;

extern float particleRestVolume;

extern float density_rest;

extern float viscosityFactor;

extern sf::Vector2f gravity;

extern const float timeStep;

extern bool EOS_Pressure;

extern bool IISPH_Pressure;

extern int countFluidParticles;

extern int currentIterations;

extern bool pressureColors;

extern float gamma_1;

extern float gamma_2;

extern float omega;

extern float beta;

extern bool CFLCondition;

extern const float cellSize;
extern const int numCellsX;
extern const int numCellsY;

extern const float particleDiameter;
extern const float particleRadius;

extern bool simulationPaused;

extern bool switchDir;

// for color gradient
extern float speedThreshold1;
extern float speedThreshold2;

////////////////////// helper functions //////////////////////

// Function to multiply sf::Vector2f by a float
inline sf::Vector2f operator*(const sf::Vector2f& vec, float scalar) {
    return sf::Vector2f(vec.x * scalar, vec.y * scalar);
}

// Function to multiply float by sf::Vector2f
inline sf::Vector2f operator*(float scalar, const sf::Vector2f& vec) {
    return sf::Vector2f(vec.x * scalar, vec.y * scalar);
}

// Function to multiply sf::Vector2f by sf::Vector2f component-wise
inline sf::Vector2f operator*(const sf::Vector2f& lhs, const sf::Vector2f& rhs) {
    return sf::Vector2f(lhs.x * rhs.x, lhs.y * rhs.y);
}

// Function to add sf::Vector2f and a double
inline sf::Vector2f operator+(const sf::Vector2f& vec, double scalar) {
    return sf::Vector2f(vec.x + static_cast<float>(scalar), vec.y + static_cast<float>(scalar));
}

// Function to add a double and sf::Vector2f
inline sf::Vector2f operator+(double scalar, const sf::Vector2f& vec) {
    return sf::Vector2f(static_cast<float>(scalar) + vec.x, static_cast<float>(scalar) + vec.y);
}

// Function to add two sf::Vector2f
inline sf::Vector2f operator+(const sf::Vector2f& lhs, const sf::Vector2f& rhs) {
    return sf::Vector2f(lhs.x + rhs.x, lhs.y + rhs.y);
}

// Function to subtract sf::Vector2f and a double
inline sf::Vector2f operator-(const sf::Vector2f& vec, double scalar) {
    return sf::Vector2f(vec.x - static_cast<float>(scalar), vec.y - static_cast<float>(scalar));
}

// Function to subtract a double and sf::Vector2f
inline sf::Vector2f operator-(double scalar, const sf::Vector2f& vec) {
    return sf::Vector2f(static_cast<float>(scalar) - vec.x, static_cast<float>(scalar) - vec.y);
}

// Function to subtract two sf::Vector2f
inline sf::Vector2f operator-(const sf::Vector2f& lhs, const sf::Vector2f& rhs) {
    return sf::Vector2f(lhs.x - rhs.x, lhs.y - rhs.y);
}

// Overload the += operator for sf::Vector2f
inline sf::Vector2f& operator+=(sf::Vector2f& lhs, const sf::Vector2f& rhs) {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    return lhs;
}

// Overload the -= operator for sf::Vector2f
inline sf::Vector2f& operator-=(sf::Vector2f& lhs, const sf::Vector2f& rhs) {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    return lhs;
}

// Function to divide sf::Vector2f by a float
inline sf::Vector2f operator/(const sf::Vector2f& vec, float scalar) {
    return sf::Vector2f(vec.x / scalar, vec.y / scalar);
}

// Function to divide two sf::Vector2f component-wise
inline sf::Vector2f operator/(const sf::Vector2f& lhs, const sf::Vector2f& rhs) {
    return sf::Vector2f(lhs.x / rhs.x, lhs.y / rhs.y);
}

// Function to calculate the dot product of two sf::Vector2f
inline float dotProduct(const sf::Vector2f& lhs, const sf::Vector2f& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

float euclideanDistance(const sf::Vector2f& positionA, const sf::Vector2f& positionB);

void writeGlobalsToFile(const std::string& filename);

std::string generateFrameFilename(int frameNumber);

#endif // GLOBALS_H