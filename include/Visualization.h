// Visualization.h
#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <SFML/Graphics.hpp>
#include "Particle.h"
#include "ParticleSystem.h"
#include <vector>

class ParticleSystemRenderer : public sf::Drawable, public sf::Transformable {
public:
    sf::VertexArray m_vertices;
    sf::Texture m_texture;

    void setVertices(std::vector<Particle>& particles, sf::Vector2u textureSize);
    void loadTexture(const std::string& texturePath);
    void updateCamera(sf::RenderWindow& window, sf::View& view, sf::Event& event, ParticleSystem& particleSystem);

    void drawGrid(sf::RenderWindow& window, int x_size_screen, int y_size_screen, int cellSize);


private:
    void draw(sf::RenderTarget& target, sf::RenderStates states) const override;
};

#endif // VISUALIZATION_H
