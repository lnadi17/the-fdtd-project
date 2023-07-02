# the-fdtd-project

## Running the Simulations

To run the simulation, execute the MATLAB code provided. Both simulations, 1D and 2D show the electric field amplitudes, but it can easily be changed to the magnetic fields. 
Anyone can customize the parameters, experiment with different configurations and observe the effects on wave propagation.

## fdtd_1d.m
This MATLAB code performs a 1D Finite-Difference Time-Domain (FDTD) simulation to model the propagation of electromagnetic waves in a lossy dielectric medium. The code is based on Example 1.7 from the book "Sullivan, Dennis M. Electromagnetic Simulation Using the FDTD Method." It provides customization options to modify various parameters and control the simulation process from the code.

### 1D Simulation Features and Customizations

1. **Grid Configuration:** The grid parameters, such as physical size, number of cells, and cell size, can be customized to suit different simulation scenarios.

2. **Physical Parameters:** The code allows modification of physical parameters, such as permittivity and permeability, to simulate different materials and media.

3. **Source Selection:** Users can choose between a Gaussian pulse source or a harmonic source. Parameters such as pulse width, offset, and frequency can be adjusted.

4. **Source Types:** The simulation offers different source types, including soft, hard, or directional sources. Users can uncomment the desired source type based on their simulation requirements.

5. **Metal Placement:** By uncommenting a specific line in the code (Line 104), users can place a thin metal object within the simulation space to observe its impact on wave propagation.

6. **User Interaction:** The code provides an optional feature to prompt the user to continue or stop the simulation at a specific step. A dialog box will be displayed, allowing the user to control the simulation flow.

7. **Visualization Control:** Users can adjust the visualization step size to display the electric (or magnetic) field visualization at specific intervals.

## Acknowledgments

This code is based on the concepts and examples presented in the book "Sullivan, Dennis M. Electromagnetic Simulation Using the FDTD Method."

## License

This project is licensed under the [GPL-3.0 License](LICENSE).
