// Visualization.cpp
#include "../include/Visualization.h"
#include "../include/Globals.h"


void ParticleSystemRenderer::setVertices(std::vector<Particle>& particles, sf::Vector2u textureSize) {
    m_vertices.setPrimitiveType(sf::Triangles);
    m_vertices.resize(particles.size() * 3);

    float height = sqrt(3.0f) / 2.0f * 2.0f * particleRadius; // Height of equilateral triangle


    for (std::size_t i = 0; i < particles.size(); ++i) {
        sf::Vertex* tri = &m_vertices[i * 3];

        tri[0].position = particles[i].position + sf::Vector2f(-particleRadius, height / 3.0f);
        tri[1].position = particles[i].position + sf::Vector2f(particleRadius, height / 3.0f);
        tri[2].position = particles[i].position + sf::Vector2f(0, -2 * height / 3.0f);


        tri[0].texCoords = sf::Vector2f(0, 0);
        tri[1].texCoords = sf::Vector2f(textureSize.x, 0);
        tri[2].texCoords = sf::Vector2f(textureSize.x / 2, textureSize.y);
        
        tri[0].color = particles[i].color;
        tri[1].color = particles[i].color;
        tri[2].color = particles[i].color;
    }
}

void ParticleSystemRenderer::loadTexture(const std::string& texturePath) {
    if (!m_texture.loadFromFile(texturePath)) {
        // handle error...
    }
    m_texture.setSmooth(true);
}

void ParticleSystemRenderer::draw(sf::RenderTarget& target, sf::RenderStates states) const {
    states.transform *= getTransform();
    states.texture = &m_texture;
    target.draw(m_vertices, states);
}

void ParticleSystemRenderer::updateCamera(sf::RenderWindow& window, sf::View& view, sf::Event& event, ParticleSystem& particleSystem) {
    while (window.pollEvent(event)) {
        if (event.type == sf::Event::Closed)
            window.close();

        // Camera controls
        if (event.type == sf::Event::KeyPressed) {
            switch (event.key.code) {
                case sf::Keyboard::Left:
                    view.move(-10, 0);
                    break;
                case sf::Keyboard::Right:
                    view.move(10, 0);
                    break;
                case sf::Keyboard::Up:
                    view.move(0, -10);
                    break;
                case sf::Keyboard::Down:
                    view.move(0, 10);
                    break;
                case sf::Keyboard::A:  // Zoom out
                    view.zoom(1.1f);
                    break;
                case sf::Keyboard::B:  // Zoom in
                    view.zoom(0.9f);
                    break;
                case sf::Keyboard::Space:  // Toggle simulation pause
                    simulationPaused = !simulationPaused;
                    break;
            }
        }

        // Convert mouse position to world coordinates (to correct for transformations)
        sf::Vector2i pixelPos = sf::Mouse::getPosition(window);
        sf::Vector2f worldPos = window.mapPixelToCoords(pixelPos);

        // Check mouse hover
        // particleSystem.checkMouseHover(static_cast<sf::Vector2i>(worldPos)); // visual test
    }
}

void ParticleSystemRenderer::drawGrid(sf::RenderWindow& window, int x_size_screen, int y_size_screen, int cellSize) {
    sf::Color gridColor = sf::Color(100, 100, 100, 100);  // Light gray color

    // Draw vertical grid lines
    for (int x = 0; x <= x_size_screen; x += cellSize) {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(static_cast<float>(x), 0), gridColor),
            sf::Vertex(sf::Vector2f(static_cast<float>(x), static_cast<float>(y_size_screen)), gridColor)
        };
        window.draw(line, 2, sf::Lines);
    }

    // Draw horizontal grid lines
    for (int y = 0; y <= y_size_screen; y += cellSize) {
        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(0, static_cast<float>(y)), gridColor),
            sf::Vertex(sf::Vector2f(static_cast<float>(x_size_screen), static_cast<float>(y)), gridColor)
        };
        window.draw(line, 2, sf::Lines);
    }
}