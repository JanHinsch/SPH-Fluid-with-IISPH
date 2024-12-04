#include <SFML/Graphics.hpp>
#include <TGUI/TGUI.hpp>
#include <sstream>
#include <iomanip>
#include "../include/UIManager.h"
#include "../include/ParticleSystem.h"
#include "../include/Visualization.h"
#include "../include/SimulationEOS.h"
#include "../include/Globals.h"

int count = 0;

int main() {

    // Create the SFML window
    sf::RenderWindow window(sf::VideoMode(x_size_screen, y_size_screen), "Particle System Simulation");
    UIManager uiManager(window);

    // Create the particle system and visualization renderer
    ParticleSystem particles(1000);
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
        [&](int value) {
            surfaceTensionFactor = value;
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
            particles.updateParticlesIISPH(particles.m_particles, x_size_screen, y_size_screen);
        }
        // Update particles with IISPH and Extrapolation at Boundaries
        if (!simulationPaused && IISPH_Pressure_Extrapolation) {
            particles.updateParticlesIISPH_Extrapolation(particles.m_particles, x_size_screen, y_size_screen);
        }
        // Update particles with DFSPH and MLS Extrapolation for Boundaries
        if (!simulationPaused && DFSPH_Pressure) {
            particles.updateParticlesDFSPH(particles.m_particles, x_size_screen, y_size_screen);
        }

        // Update vertices based on new particle positions
        renderer.setVertices(particles.m_particles, textureSize);

        // Handle camera updates
        renderer.updateCamera(window, view, event, particles);

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
        // std::cout << "\rcount: " << count << "/500 " << std::flush;
        // if (count == 500) {
        //     std::cout << "Simulation complete. Exiting..." << std::endl;
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

            view.move(-300, 130);
        }
    }

    return 0;
}
