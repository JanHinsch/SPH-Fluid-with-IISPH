#include <SFML/Graphics.hpp>
#include <TGUI/TGUI.hpp>
#include <sstream>
#include <iomanip>
#include "../include/UIManager.h"
#include "../include/ParticleSystem.h"
#include "../include/Visualization.h"
#include "../include/SimulationEOS.h"
#include "../include/Globals.h"
sf::Clock uiClock;

const float uiUpdateInterval = 0.1f; // 100 ms update interval
int count = 0;

int main() {

    // Create the SFML window
    sf::RenderWindow window(sf::VideoMode(x_size_screen, y_size_screen), "Particle System Simulation");
    UIManager uiManager(window);

    // Create the particle system and visualization renderer
    ParticleSystem particles(1008);
    particles.initiateParticles(particles.m_particles);

    ParticleSystemRenderer renderer;
    renderer.loadTexture("../textures/white_circle.png");  // Ensure the texture path is correct
    sf::View view(sf::FloatRect(0.f, 0.f, x_size_screen, y_size_screen)); // for camera

    // Set the initial camera position (e.g., center the view on the middle of the window)
    view.setCenter(x_size_screen / 2.0f, y_size_screen / 2.0f);

    // Initialize the particles in the renderer
    sf::Vector2u textureSize = renderer.m_texture.getSize();
    renderer.setVertices(particles.m_particles, textureSize);

    // Initialize the UI
    uiManager.initializeUI(
        [&particles]() {
            simulationPaused = !simulationPaused;
        },
        [&particles]() {
            particles.resetSimulation();
        },
        [&](float value) {
            stiffness_constant_k = value;
        },
        [&](float value) {
            viscosityFactor = value;
        },
        [&](const sf::Vector2f& value) {
            gravity = value;
        },
        [&](float value) {
            surfaceTensionFactor = value;
        },
        [](int mode) {
            UIManager::changePressureMode(mode); // Properly wrap the call to changePressureMode
        },
        [&](int value) {
            gamma_3 = value;
        }
    );


    // Force an initial update to ensure everything is correctly displayed at the start
    window.setView(view);
    window.clear();
    window.draw(renderer);
    uiManager.update();
    uiManager.draw();
    window.display();

    // Main loop
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (uiManager.handleEvents(event)) {
                // If TGUI handled the event, skip further processing
                continue;
            }

            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        // Update particles with EOS
        if (!simulationPaused && EOS_Pressure) {
            particles.updateParticlesEOS(particles.m_particles, x_size_screen, y_size_screen);
        }
        // Update particles with IISPH and Pressure Boundaries
        if (!simulationPaused && IISPH_Pressure_Boundaries) {
            particles.updateParticlesIISPHPressureBoundaries(particles.m_particles, x_size_screen, y_size_screen);
        }
        // Update particles with IISPH and SPH Extrapolation at Boundaries
        if (!simulationPaused && IISPH_Pressure_Extrapolation) {
            particles.updateParticlesIISPH_Extrapolation(particles.m_particles, x_size_screen, y_size_screen);
        }
        // Update particles with IISPH and Mirror at Boundaries (no pressure calculation at boundary)
        if (!simulationPaused && IISPH_Pressure_Mirroring) {
            particles.updateParticlesIISPH_Mirroring(particles.m_particles, x_size_screen, y_size_screen);
        }
        // Update particles with DFSPH and MLS Extrapolation for Boundaries
        if (!simulationPaused && IISPH_MLS_Pressure_Extrapolation) {
            particles.updateParticlesIISPH_MLSExtrapolation(particles.m_particles, x_size_screen, y_size_screen);
        }

        // Update vertices based on new particle positions
        renderer.setVertices(particles.m_particles, textureSize);

        // Handle camera updates
        renderer.updateCamera(window, view, event, particles);


        // Update UI periodically
        if (uiClock.getElapsedTime().asSeconds() > uiUpdateInterval) {
            uiManager.update();
            uiClock.restart();
        }

        // Rendering
        window.clear();
        window.setView(view);
        window.draw(renderer);
        uiManager.update();
        uiManager.draw();
        window.display();


        // Save each frame
        // sf::Image screenshot = window.capture();
        // screenshot.saveToFile(generateFrameFilename(count));


        // Exit the loop
        count++;
        // std::cout << "\rcount: " << count << "/1000 " << std::flush;
        // if (count == 1000) {
        // std::cout << "Simulation complete. Exiting..." << std::endl;
        //     break;
        // }

        // quick bugfix
        if (count==1) {
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);
            view.zoom(0.9f);

            view.move(-330, 130);
        }

        // if (count==5000) {
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //     view.zoom(0.9f);
        //
        //     view.move(0, 80);
        // }
    }

    return 0;
}
