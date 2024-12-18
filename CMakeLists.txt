cmake_minimum_required(VERSION 3.10)

# Project name and version
project(ParticleSystemSimulation VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find SFML package
find_package(SFML 2.6 COMPONENTS graphics window system REQUIRED)

find_package(TGUI 1.3 REQUIRED)

# Include directories
include_directories(include)

# Enable debug symbols by default for all builds
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Add the executable
add_executable(main

    src/Globals.cpp
    src/Grid.cpp
    src/ParticleSystem.cpp
        src/SimulationEOS.cpp
        src/SimulationIISPH_PressureBoundaries.cpp
    src/Visualization.cpp
    src/UIManager.cpp
    src/main.cpp
        src/SimulationIISPH_MLSExtrapolation.cpp
        src/ParticleSpawning.cpp

)


# Link SFML libraries
target_link_libraries(main sfml-graphics sfml-window sfml-system tgui)


# Add Catch2 tests
add_executable(tests_neighbour
    tests/test_neighbour_search.cpp
    src/Globals.cpp
    src/Grid.cpp
    src/ParticleSystem.cpp
        src/SimulationEOS.cpp
    src/Visualization.cpp
        src/ParticleSpawning.cpp
)

add_executable(tests_SPH_Kernel
    tests/test_SPH_Kernel.cpp
    src/Globals.cpp
    src/Grid.cpp
    src/ParticleSystem.cpp
        src/SimulationEOS.cpp
    src/Visualization.cpp
        src/ParticleSpawning.cpp
)

# add_executable(tests_SPH_Kernel_multiple_Particles
#     tests/test_SPH_Kernel_multiple_Particles.cpp
#     src/Globals.cpp
#     src/Grid.cpp
#     src/ParticleSystem.cpp
#     src/Simulation.cpp
#     src/Visualization.cpp
# )

add_executable(tests_SPH_Computations
    tests/test_SPH_Computations.cpp
    src/Globals.cpp
    src/Grid.cpp
    src/ParticleSystem.cpp
        src/SimulationEOS.cpp
        src/SimulationIISPH_PressureBoundaries.cpp
    src/Visualization.cpp
        src/ParticleSpawning.cpp
)

# Include Catch2
include_directories(external/catch2)

# Include Eigen
include_directories(external/eigen-master/Eigen)

# Link SFML libraries to tests
target_link_libraries(tests_neighbour sfml-graphics sfml-window sfml-system)
target_link_libraries(tests_SPH_Kernel sfml-graphics sfml-window sfml-system)

#target_link_libraries(tests_SPH_Kernel_multiple_Particles sfml-graphics sfml-window sfml-system)

target_link_libraries(tests_SPH_Computations sfml-graphics sfml-window sfml-system)