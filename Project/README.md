To compile this project your need:

1. GCC compiler
2. PETSc library >= 3.12.4
3. LaTeX env.
4. Python >= 3.6 with numpy, matplotlib and numexpr packages installed
5. Run the command `cd build && cmake ..`

run `python3 main.py -h` in the project directory to see all the possibilies.

# Examples

## Build and run the project.c
`python3 main.py -r`

## Building project.c and generating airfoils dimensionless plots
### 7-c zone
`python3 main.py --plot mesh_w -no_save -save_frames --limits 0.33 --adim_unit_value 0.046635 --adim_unit_symbol c -adim`

### 2-c zone
`python3 main.py --plot mesh_w -no_save -save_frames --limits 0.094 --adim_unit_value 0.046635 --adim_unit_symbol c -adim`

### Making the video from the frames
`python3 main.py --make_video mesh_w`
