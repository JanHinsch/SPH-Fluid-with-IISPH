// Grid.cpp
#include "../include/Grid.h"

// Constructor to initialize the grid with given width and height
Grid::Grid(int width, int height) : width(width), height(height), C(width * height + 1) {}

// Function to get value at (x, y)
int Grid::get(int x, int y) const {
    if (!isValid(x, y)) throw std::out_of_range("Invalid coordinates");
    return C[y * width + x];
}

// Function to set value at (x, y)
void Grid::set(int x, int y, int value) {
    if (!isValid(x, y)) throw std::out_of_range("Invalid coordinates");
    C[y * width + x] = value;
}

// Function to check if coordinates are valid
bool Grid::isValid(int x, int y) const {
    return x >= 0 && x < width && y >= 0 && y < height;
}

// Function to clear all entries in the grid
void Grid::clear() {
    std::fill(C.begin(), C.end(), 0);
}

// Function to get a reference to the C array
void Grid::decValueAt(int x, int y) {
    C[y * width + x] -= 1;
}

// Function to get the grid cell (x, y) coordinates from the index
std::pair<int, int> Grid::getCoordinatesFromIndex(int index) const {
    if (index < 0 || index >= width * height) throw std::out_of_range("Index out of bounds");

    int x = index % width;
    int y = static_cast<int>(index / width); // Convert to integer if y is a floating-point value

    return std::make_pair(x, y);
}

// Operator to access the C array by index
int& Grid::operator[](size_t index) {
    if (index >= C.size()) throw std::out_of_range("Index out of range");
    return C[index];
}