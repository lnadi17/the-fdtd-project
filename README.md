# the-fdtd-project

This repository contains MATLAB code for performing 1D and 2D Finite-Difference Time-Domain (FDTD) simulations. The simulations model the propagation of electromagnetic waves in various scenarios using the FDTD method. Different customization options are available to modify the simulation parameters and observe the effects on wave propagation.

## 1D Simulation (fdtd_1d.m)
This MATLAB code performs a 1D Finite-Difference Time-Domain (FDTD) simulation to model the propagation of electromagnetic waves in a lossy dielectric medium. The code is based on Example 1.7 from the book "Sullivan, Dennis M. Electromagnetic Simulation Using the FDTD Method." It provides customization options to modify various parameters and control the simulation process from the code.

## 2D Simulation (fdtd_2d.m)
The second FDTD simulation is a 2D simulation that models the propagation of electromagnetic waves in a rectangular domain from two coupled circular E-field sources using the same FDTD approach. The simulation parameters and customization options are outlined below.

## Running the Simulations

To run the simulation, execute the MATLAB code provided. Both simulations, 1D and 2D show the electric field amplitudes, but it can easily be changed to the magnetic fields. 
Anyone can customize the parameters, experiment with different configurations and observe the effects on wave propagation.

## 1D Simulation Features and Customizations

- **Grid Configuration:** The grid parameters, such as physical size, number of cells, and cell size, can be customized to suit different simulation scenarios.

- **Physical Parameters:** The code allows modification of physical parameters, such as permittivity and permeability, to simulate different materials and media.

- **Source Selection:** Users can choose between a Gaussian pulse source or a harmonic source. Parameters such as pulse width, offset, and frequency can be adjusted.

- **Source Types:** The simulation offers different source types, including soft, hard, or directional sources. Users can uncomment the desired source type based on their simulation requirements.

- **Metal Placement:** By uncommenting a specific line in the code (Line 104), users can place a thin metal object within the simulation space to observe its impact on wave propagation.

- **User Interaction:** The code provides an optional feature to prompt the user to continue or stop the simulation at a specific step. A dialog box will be displayed, allowing the user to control the simulation flow.

- **Visualization Control:** Users can adjust the visualization step size to display the electric (or magnetic) field visualization at specific intervals.

## 2D Simulation Features and Customizations
- **Grid Configuration:** Grid parameters can be configured similar to the 1D simulation.
  
- **Physical Parameters:** Current simulation assumes free space permittivity and permeability, but this can be easily changed if desired.
  
- **Source Selection:** Source is somewhat hardcoded to have a circular shape, and the previous point source implementation was removed from the code. However, parameters as source frequency, spacial phase variation, phase difference between sources, radius, width and other parameters can be customized.
  
- **User Interaction:** This feature works just the same as the 1D simulation counterpart.
  
- **Perfectly Matched Layer:** A Perfectly Matched Layer (PML) is implemented to absorb outgoing waves effectively.
  
- **Detector:** The 2D FDTD simulation includes a detector feature that captures the average power direction during the simulation. It accumulates the absolute values of the Dz field at each grid point. After the simulation, the detector values are normalized and shown as a graph.

## Acknowledgments

This project was created during an Applied Electromagnetics course and implements concepts inspired by the [EMPossible](https://empossible.net/academics/emp5304/). This code is also based on the concepts and examples presented in the book "Sullivan, Dennis M. Electromagnetic Simulation Using the FDTD Method."

## License

This project is licensed under the [GPL-3.0 License](LICENSE).
