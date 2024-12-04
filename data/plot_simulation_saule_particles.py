# File: plot_simulation_saule_particles.py

import matplotlib.pyplot as plt
import numpy as np

def read_file(file_path):
    """Reads a file and returns its contents as a list of lines."""
    try:
        with open(file_path, 'r') as file:
            return [line.strip() for line in file.readlines()]
    except FileNotFoundError:
        print(f"Error: {file_path} not found.")
        return []

def calculate_average(data):
    """Calculates the average of a list of numerical values."""
    return np.mean(data) if data else 0

def main():
    # File names and labels for the data files
    files = [
        ("currentIterations_normalPressureInit.txt", "Normal Pressure Init"),
        ("currentIterations_maxPressureInit.txt", "Particles Max Pressure Init")
    ]
    
    # Path to the directory containing the data files
    base_path = "/Users/jan/Desktop/SPH_Fluid_with_IISPH/cmake-build-debug/"
    
    # Colors for the different lines
    colors = ['blue', 'green']
    
    # Store averages for each file
    averages = []

    plt.figure(figsize=(12, 8))

    # Read and plot data from each file
    for idx, (file_name, label) in enumerate(files):
        file_path = base_path + file_name
        iteration_data = read_file(file_path)
        try:
            # Convert each line into an integer (assumes data in the files are numbers)
            iteration_values = [int(float(value)) for value in iteration_data]
        except ValueError:
            print(f"Error: Non-numeric data found in {file_name}.")
            continue

        # Define the X-axis as simulation steps
        x_values = list(range(1, len(iteration_values) + 1))

        # Plot the data with a unique color
        plt.plot(x_values, iteration_values, marker='o', linestyle='-', color=colors[idx], label=label)

        # Calculate and store the average value
        avg_value = calculate_average(iteration_values)
        averages.append((label, avg_value))

    # Plot settings
    plt.xlabel('Simulation Steps')
    plt.ylabel('Needed Iterations')
    plt.title('Pillar Test for Different Pressure Initialization')
    plt.grid(True)
    plt.legend(loc='upper right')

    # Create average information for each file
    average_text = "\n".join([f"Avg Iterations for {label}: {avg:.2f}" for label, avg in averages])
    
    # Add a text box with averages in the lower left corner of the plot
    plt.text(0.02, 0.98, average_text, fontsize=10, verticalalignment='top', horizontalalignment='left',
             transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.8))

    # Show the plot
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
