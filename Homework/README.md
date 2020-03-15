# LMECA2660: Homework

Simulation of the convection-diffusion equation in a 1-D domain with an initial
gaussian condition.


### How to use ?
All the simulation parameters are stored in the file src/solver.h:
- You can choose from 5 differend FD schemes for the convective term du/dx;
    this is done by changing the ```FD_SCHEME``` variable.
- You can change the spatial discretization;
    this is done by changing the ```SPACE_DISCR``` variable.
- You can do a purely convective or convective-diffusive simulation;
    this is done by changing the ```DIFF_COEFF``` variable.
- ... and many more possibilities, which are all in the solver.h header file.

### Display the results
The C code does not display the results, instead it writes them in various files.
To view the results, a Matlab script is provided in matlab/plot_data.m. You can use browse through it to view the different viewing options.
