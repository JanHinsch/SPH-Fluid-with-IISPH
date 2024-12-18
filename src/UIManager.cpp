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
    if (m_SurfaceTensionSlider) {
        m_surfaceChanged(m_SurfaceTensionSlider->getValue());
    }
    // Update CFL condition label
    updateCFLConditionLabel();

    // update current Iterations
    updateCurrentIterationsLabel(currentIterations);

    // Update Volume Error Label
    updateVolumeErrorLabel(currentVolumeError);
}

void UIManager::draw() {
    m_gui.draw();
}

void UIManager::initializeUI(const std::function<void()>& pauseResumeCallback,
                             const std::function<void()>& resetCallback,
                             const std::function<void(float)>& stiffnessChanged,
                             const std::function<void(float)>& viscosityChanged,
                             const std::function<void(const sf::Vector2f&)>& gravityChanged,
                             const std::function<void(float)>& SurfaceChanged,
                             const std::function<void(int)>& pressureModeChanged,
                             const std::function<void(float)>& gamma1Changed) { // Added pressureModeChanged callback
    m_pauseResumeCallback = pauseResumeCallback;
    m_resetCallback = resetCallback;
    m_stiffnessChanged = stiffnessChanged;
    m_viscosityChanged = viscosityChanged;
    m_gravityChanged = gravityChanged;
    m_surfaceChanged = SurfaceChanged;
    m_gamma1Changed = gamma1Changed;

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

    // Pressure Mode Toggle
    auto pressureModeComboBox = tgui::ComboBox::create();
    pressureModeComboBox->setPosition(10, 230); // Adjust the position as needed
    pressureModeComboBox->setSize(200, 30);
    pressureModeComboBox->addItem("EOS Pressure");
    pressureModeComboBox->addItem("Pressure Boundaries");
    pressureModeComboBox->addItem("MLS Pressure Extrapolation");
    pressureModeComboBox->addItem("SPH Extrapolation");
    pressureModeComboBox->addItem("Mirroring");
    // Default selected mode
    if (EOS_Pressure) pressureModeComboBox->setSelectedItem("EOS Pressure");
    if (IISPH_Pressure_Boundaries) pressureModeComboBox->setSelectedItem("Pressure Boundaries");
    if (IISPH_MLS_Pressure_Extrapolation) pressureModeComboBox->setSelectedItem("MLS Pressure Extrapolation");
    if (IISPH_Pressure_Extrapolation) pressureModeComboBox->setSelectedItem("SPH Extrapolation");
    if (IISPH_Pressure_Mirroring) pressureModeComboBox->setSelectedItem("Mirroring");
    pressureModeComboBox->onItemSelect([this, pressureModeChanged](const tgui::String& selectedItem) {
        int mode = 0; // Map selected item to mode
        if (selectedItem == "EOS Pressure") mode = 0;
        else if (selectedItem == "Pressure Boundaries") mode = 1;
        else if (selectedItem == "MLS Pressure Extrapolation") mode = 2;
        else if (selectedItem == "SPH Extrapolation") mode = 3;
        else if (selectedItem == "Mirroring") mode = 4;
        pressureModeChanged(mode); // Call the callback with the selected mode
    });
    m_gui.add(pressureModeComboBox);

    // Gamma1 EditBox instead of Slider
    m_gamma1EditBox = tgui::EditBox::create();
    m_gamma1EditBox->setPosition(10, 270); // Position below pressure mode combo box
    m_gamma1EditBox->setSize(200, 30);
    m_gamma1EditBox->setDefaultText("Enter Gamma1 value");
    m_gamma1EditBox->setText(std::to_string(gamma_3));
    m_gamma1EditBox->setMaximumCharacters(3);
    m_gamma1EditBox->onTextChange([this](const tgui::String& text) {
    // Convert TGUI string to std::string
    std::string str = text.toStdString();
    try {
        float value = std::stof(str); // Convert string to float
        m_gamma1Changed(value); // Call callback with the entered value
        updateGamma1Label(value); // Update the label to reflect the new value
    } catch (const std::invalid_argument&) {
        // Handle invalid input (e.g., non-numeric input)
        std::cout << "Invalid input for Gamma." << std::endl;
    }
});
    m_gui.add(m_gamma1EditBox);

    // Gamma1 Label
    m_gamma1Label = tgui::Label::create();
    m_gamma1Label->setPosition(220, 270);
    m_gamma1Label->setTextSize(16);
    updateGamma1Label(gamma_3);
    m_gui.add(m_gamma1Label);

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

    // Surface Tension Factor Slider
    m_SurfaceTensionSlider = tgui::Slider::create();
    m_SurfaceTensionSlider->setPosition(10, 140);
    m_SurfaceTensionSlider->setSize(200, 20);
    m_SurfaceTensionSlider->setMinimum(0);
    m_SurfaceTensionSlider->setMaximum(1000);
    m_SurfaceTensionSlider->setValue(surfaceTensionFactor);
    m_SurfaceTensionSlider->onValueChange([this](int value) {
        updateSurfaceLabel(value);
    });
    m_gui.add(m_SurfaceTensionSlider);

    // suraceTensionFactor Label
    m_SurfaceTensionLabel = tgui::Label::create();
    m_SurfaceTensionLabel->setPosition(220, 140);
    m_SurfaceTensionLabel->setTextSize(16);
    updateSurfaceLabel(surfaceTensionFactor);
    m_gui.add(m_SurfaceTensionLabel);

    // CFL Label
    m_cflConditionLabel = tgui::Label::create();
    m_cflConditionLabel->setPosition(10, 170);
    m_cflConditionLabel->setTextSize(16);
    updateCFLConditionLabel();
    m_gui.add(m_cflConditionLabel);

    // Volume Error Label
    m_volumeErrorLabel = tgui::Label::create();
    m_volumeErrorLabel->setPosition(10, 310);
    m_volumeErrorLabel->setTextSize(16);
    updateVolumeErrorLabel(currentVolumeError);
    m_gui.add(m_volumeErrorLabel);


    // Conditional UI for IISPH
    if (IISPH_Pressure_Boundaries || IISPH_Pressure_Extrapolation || IISPH_Pressure_Mirroring || IISPH_Pressure_Zero
        || IISPH_MLS_Pressure_Extrapolation) {
        // Add Current Iterations Label under the CFL Condition
        m_currentIterationsLabel = tgui::Label::create();
        m_currentIterationsLabel->setPosition(10, 200);
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

void UIManager::updateSurfaceLabel(float value) const {
    m_SurfaceTensionLabel->setText("Surface Tension Factor: " + std::to_string(static_cast<int>(value)));
    m_SurfaceTensionLabel->getRenderer()->setTextColor(sf::Color(255,255,255));
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

void UIManager::updateVolumeErrorLabel(float value) {
    if (CFLCondition == true) {
        m_volumeErrorLabel->setText("Error: " + std::to_string(value) + "%");
        m_volumeErrorLabel->getRenderer()->setTextColor(sf::Color::White);
    }
}

void UIManager::updateCurrentIterationsLabel(int value) {
    if (m_currentIterationsLabel) {
        m_currentIterationsLabel->setText("Current Iterations: " + std::to_string(value));
        m_currentIterationsLabel->getRenderer()->setTextColor(sf::Color(255, 255, 255));
    }
}

// Define a function to handle pressure mode changes
void UIManager::changePressureMode(int mode) {
    EOS_Pressure = (mode == 0);
    IISPH_Pressure_Boundaries = (mode == 1);
    IISPH_MLS_Pressure_Extrapolation = (mode == 2);
    IISPH_Pressure_Extrapolation = (mode == 3);
    IISPH_Pressure_Mirroring = (mode == 4);

    // Optional: Print mode for debugging
    std::cout << "Pressure mode changed to: " << mode << std::endl;
};

void UIManager::updateGamma1Label(float value) const {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << value;
    m_gamma1Label->setText("Gamma");
    m_gamma1Label->getRenderer()->setTextColor(sf::Color(255, 255, 255));
}
