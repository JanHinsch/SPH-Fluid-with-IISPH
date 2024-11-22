#pragma once

#include <SFML/Graphics.hpp>
#include <TGUI/TGUI.hpp>
#include <TGUI/Backend/SFML-Graphics.hpp>
#include <TGUI/Core.hpp>
#include <functional>

class UIManager {
public:
    UIManager(sf::RenderWindow& window);
    bool handleEvents(const sf::Event& event);
    void update();
    void draw();

    void initializeUI(const std::function<void()>& pauseResumeCallback,
                      const std::function<void()>& resetCallback,
                      const std::function<void(float)>& stiffnessChanged,
                      const std::function<void(float)>& viscosityChanged,
                      const std::function<void(const sf::Vector2f&)>& gravityChanged,
                      const std::function<void(float)>& gamma2Changed);

    void updateStiffnessLabel(float value);

    void updateViscosityLabel(float value);

    void updateGravityYLabel(float value);

    void updateGamma2Label(float value);

    void updateCFLConditionLabel();

    void updateCurrentIterationsLabel(int value);

private:
    tgui::Gui m_gui;
    std::function<void()> m_pauseResumeCallback;
    std::function<void()> m_resetCallback;
    std::function<void(float)> m_stiffnessChanged;
    std::function<void(float)> m_viscosityChanged;
    std::function<void(const sf::Vector2f&)> m_gravityChanged;
    std::function<void(float)> m_gamma2Changed;

    tgui::Slider::Ptr m_stiffnessSlider;
    tgui::Slider::Ptr m_viscositySlider;
    tgui::Slider::Ptr m_gravityYSlider;
    tgui::Slider::Ptr m_gamma2Slider;

    std::shared_ptr<tgui::Label> m_stiffnessLabel;
    std::shared_ptr<tgui::Label> m_viscosityLabel;
    std::shared_ptr<tgui::Label> m_gravityYLabel;
    std::shared_ptr<tgui::Label> m_gamma2Label;
    std::shared_ptr<tgui::Label> m_cflConditionLabel;
    std::shared_ptr<tgui::Label> m_currentIterationsLabel;
};
