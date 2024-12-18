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

    void initializeUI(const std::function<void()> &pauseResumeCallback,
                      const std::function<void()> &resetCallback,
                      const std::function<void(float)> &stiffnessChanged,
                      const std::function<void(float)> &viscosityChanged,
                      const std::function<void(const sf::Vector2f &)> &gravityChanged,
                      const std::function<void(float)> &surfaceChanged,
                      const std::function<void(int)> &pressureModeChanged,
                      const std::function<void(float)> &gamma1Changed);

    void updateStiffnessLabel(float value);

    void updateViscosityLabel(float value);

    void updateGravityYLabel(float value);

    void updateSurfaceLabel(float value) const;

    void updateCFLConditionLabel();

    void updateVolumeErrorLabel(float value);

    void updateCurrentIterationsLabel(int value);

    static void changePressureMode(int value);

    void updateGamma1Label(float value) const;

private:
    tgui::Gui m_gui;
    std::function<void()> m_pauseResumeCallback;
    std::function<void()> m_resetCallback;
    std::function<void(float)> m_stiffnessChanged;
    std::function<void(float)> m_viscosityChanged;
    std::function<void(const sf::Vector2f&)> m_gravityChanged;
    std::function<void(float)> m_surfaceChanged;
    std::function<void(float)> m_gamma1Changed;

    tgui::Slider::Ptr m_stiffnessSlider;
    tgui::Slider::Ptr m_viscositySlider;
    tgui::Slider::Ptr m_gravityYSlider;
    tgui::Slider::Ptr m_SurfaceTensionSlider;

    tgui::EditBox::Ptr m_gamma1EditBox;

    std::shared_ptr<tgui::Label> m_stiffnessLabel;
    std::shared_ptr<tgui::Label> m_viscosityLabel;
    std::shared_ptr<tgui::Label> m_gravityYLabel;
    std::shared_ptr<tgui::Label> m_SurfaceTensionLabel;
    std::shared_ptr<tgui::Label> m_cflConditionLabel;
    std::shared_ptr<tgui::Label> m_volumeErrorLabel;
    std::shared_ptr<tgui::Label> m_currentIterationsLabel;
    std::shared_ptr<tgui::Label> m_gamma1Label;
};
