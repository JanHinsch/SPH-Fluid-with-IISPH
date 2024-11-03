#ifndef GRID_H
#define GRID_H

#include <vector>
#include <stdexcept>

class Grid {
public:
    // Constructor to initialize the grid with given width and height
    Grid(int width, int height);

    // Function to get value at (x, y)
    int get(int x, int y) const;

    // Function to set value at (x, y)
    void set(int x, int y, int value);

    // Function to clear all entries in the grid
    void clear();

    void decValueAt(int x, int y);

    // Function to check if coordinates are valid
    bool isValid(int x, int y) const;

    std::pair<int, int> getCoordinatesFromIndex(int index) const;

    // Operator to access the C array by index
    int& operator[](size_t index);

private:
    int width, height;
    std::vector<int> C;

};

#endif // GRID_H
