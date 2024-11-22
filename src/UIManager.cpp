#include "../include/UIManager.h"
#include "../include/Globals.h"


UIManager::UIManager(sf::RenderWindow& window)
    : m_gui(window) {
}

bool UIManager::handleEvents(const sf::Event& event) {
    if (m_gui.handleEvent(event)) {
        return true;
    }
    return false;
}

void UIManager::update() {
    if (m_stiffnessSlider) {
        m_stiffnessChanged(m_stiffnessSlider->getValue());
    }
    if (m_viscositySlider) {
        m_viscosityChanged(m_viscositySlider->getValue());
    }
    if (m_gravityYSlider) {
        m_gravityChanged(sf::Vector2f(0.0f, m_gravityYSlider->getValue()));
    }
    if (m_gamma2Slider) {
        m_gamma2Changed(m_gamma2Slider->getValue());
    }
    // Update CFL condition label
    updateCFLConditionLabel();

    // update current Iterations
    updateCurrentIterationsLabel(currentIterations);
}

void UIManager::draw() {
    m_gui.draw();
}

void UIManager::initializeUI(const std::function<void()>& pauseResumeCallback,
                             const std::function<void()>& resetCallback,
                             const std::function<void(float)>& stiffnessChanged,
                             const std::function<void(float)>& viscosityChanged,
                             const std::function<void(const sf::Vector2f&)>& gravityChanged,
                             const std::function<void(float)>& gamma2Changed) {
    m_pauseResumeCallback = pauseResumeCallback;
    m_resetCallback = resetCallback;
    m_stiffnessChanged = stiffnessChanged;
    m_viscosityChanged = viscosityChanged;
    m_gravityChanged = gravityChanged;
    m_gamma2Changed = gamma2Changed;

    auto button = tgui::Button::create("Pause/Resume");
    button->setPosition(10, 10);
    button->setSize(100, 30);
    button->onPress([this]() {
        m_pauseResumeCallback();
    });
    m_gui.add(button);

    auto resetButton = tgui::Button::create("Reset");
    resetButton->setPosition(110, 10);
    resetButton->setSize(100, 30);
    resetButton->onPress([this]() {
        m_resetCallback();
    });
    m_gui.add(resetButton);

    // Stiffness Slider
    if (EOS_Pressure) {
        m_stiffnessSlider = tgui::Slider::create();
        m_stiffnessSlider->setPosition(10, 50);
        m_stiffnessSlider->setSize(200, 20);
        m_stiffnessSlider->setMinimum(0.0f);
        m_stiffnessSlider->setMaximum(2000.0f);
        m_stiffnessSlider->setValue(stiffness_constant_k);
        m_stiffnessSlider->onValueChange([this](float value) {
            m_stiffnessChanged(value);
            updateStiffnessLabel(value);
        });
        m_gui.add(m_stiffnessSlider);
    }

    // Stiffness Label
    if (EOS_Pressure) {
        m_stiffnessLabel = tgui::Label::create();
        m_stiffnessLabel->setPosition(220, 50);
        m_stiffnessLabel->setTextSize(16);
        updateStiffnessLabel(stiffness_constant_k);
        m_gui.add(m_stiffnessLabel);
    }

    // Viscosity Slider
    m_viscositySlider = tgui::Slider::create();
    m_viscositySlider->setPosition(10, 80);
    m_viscositySlider->setSize(200, 20);
    m_viscositySlider->setMinimum(20.0f);
    m_viscositySlider->setMaximum(50.0f);
    m_viscositySlider->setValue(viscosityFactor);
    m_viscositySlider->onValueChange([this](float value) {
        m_viscosityChanged(value);
        updateViscosityLabel(value);
    });
    m_gui.add(m_viscositySlider);

    // Viscosity Label
    m_viscosityLabel = tgui::Label::create();
    m_viscosityLabel->setPosition(220, 80);
    m_viscosityLabel->setTextSize(16);
    updateViscosityLabel(viscosityFactor);
    m_gui.add(m_viscosityLabel);

    // Gravity Y Slider
    m_gravityYSlider = tgui::Slider::create();
    m_gravityYSlider->setPosition(10, 110);
    m_gravityYSlider->setSize(200, 20);
    m_gravityYSlider->setMinimum(-20.0f);
    m_gravityYSlider->setMaximum(20.0f);
    m_gravityYSlider->setValue(gravity.y);
    m_gravityYSlider->onValueChange([this](float value) {
        m_gravityChanged(sf::Vector2f(gravity.x, value));
        updateGravityYLabel(value);
    });
    m_gui.add(m_gravityYSlider);

    // Gravity Y Label
    m_gravityYLabel = tgui::Label::create();
    m_gravityYLabel->setPosition(220, 110);
    m_gravityYLabel->setTextSize(16);
    updateGravityYLabel(gravity.y);
    m_gui.add(m_gravityYLabel);

    // Gamma 2 Slider -> now Surface Tension Factor
    m_gamma2Slider = tgui::Slider::create();
    m_gamma2Slider->setPosition(10, 140);
    m_gamma2Slider->setSize(200, 20);
    m_gamma2Slider->setMinimum(0);
    m_gamma2Slider->setMaximum(1000);
    m_gamma2Slider->setValue(surfaceTensionFactor);
    m_gamma2Slider->onValueChange([this](int value) {
        m_gamma2Changed(value);
        updateGamma2Label(value);
    });
    m_gui.add(m_gamma2Slider);

    // Gamma 2 Label  -> i changed this fast to suraceTensionFactor
    m_gamma2Label = tgui::Label::create();
    m_gamma2Label->setPosition(220, 140);
    m_gamma2Label->setTextSize(16);
    updateGamma2Label(surfaceTensionFactor);
    m_gui.add(m_gamma2Label);

    // CFL Label
    m_cflConditionLabel = tgui::Label::create();
    m_cflConditionLabel->setPosition(10, 170);
    m_cflConditionLabel->setTextSize(16);
    updateCFLConditionLabel();
    m_gui.add(m_cflConditionLabel);

    // Conditional UI for IISPH
    if (IISPH_Pressure) {
        // Add Current Iterations Label under the CFL Condition
        m_currentIterationsLabel = tgui::Label::create();
        m_currentIterationsLabel->setPosition(10, 200); // Adjust position if needed
        m_currentIterationsLabel->setTextSize(16);
        updateCurrentIterationsLabel(currentIterations);
        m_gui.add(m_currentIterationsLabel);
    }
}

void UIManager::updateStiffnessLabel(float value) {
    m_stiffnessLabel->setText("Stiffness: " + std::to_string(static_cast<int>(value)));
    m_stiffnessLabel->getRenderer()->setTextColor(sf::Color(255,255,255));
}

void UIManager::updateViscosityLabel(float value) {
    m_viscosityLabel->setText("Viscosity: " + std::to_string(static_cast<int>(value)));
    m_viscosityLabel->getRenderer()->setTextColor(sf::Color(255,255,255));
}

void UIManager::updateGravityYLabel(float value) {
    m_gravityYLabel->setText("Gravity: " + std::to_string(static_cast<int>(value)));
    m_gravityYLabel->getRenderer()->setTextColor(sf::Color(255,255,255));
}

void UIManager::updateGamma2Label(float value) {
    m_gamma2Label->setText("Surface Tension Factor: " + std::to_string(static_cast<int>(value)));
    m_gamma2Label->getRenderer()->setTextColor(sf::Color(255,255,255));
}

void UIManager::updateCFLConditionLabel() {
    if (CFLCondition == true) {
        m_cflConditionLabel->setText("CFL Condition: Yes");
        m_cflConditionLabel->getRenderer()->setTextColor(sf::Color::Green);
    } else {
        m_cflConditionLabel->setText("CFL Condition: No");
        m_cflConditionLabel->getRenderer()->setTextColor(sf::Color::Red);
    }
}

void UIManager::updateCurrentIterationsLabel(int value) {
    if (m_currentIterationsLabel) {
        m_currentIterationsLabel->setText("Current Iterations: " + std::to_string(value));
        m_currentIterationsLabel->getRenderer()->setTextColor(sf::Color(255, 255, 255));
    }
}