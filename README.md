This is the repository for the first challenge of the course: Advanced Programming for Scientific computing (Prof. Formaggia).

This program consists in a first degree solver for finding the roots of a function. The linear search problem can be chosen to be tackled with 3 methods (for each method use the code on the right to represent it in the data input):
- Exponential - 1
- Inverse - 2
- Armijo - 0 or any other positive integer

## Parameters

The file `solver/main/parameters_input.dat` is the parameter file that is given as input to the program. To use it, change the value of the parameters and save the file. To have functions of more or less than 2 variables, add or take out lines for the initial coordinates scheme already present in the file. 

So for example, implementi a function taking 3 variables, I would need to add the following lines to the file:
- `x0_3 = ... # initial value of 3rd coordinate`
- `dfun_3 = ... # partial derivative w.r.t. 3rd variable` 

## Custom functions
At the moment of writing, the function and its gradient are NOT yet read from the `solver/main/parameters_input.dat` file using `muparserx`. This is still to be implemented. 
**To change the function to be minimized**, head over in the `solver/main/` directory and change the `main.cpp`. At the start of the file there are two functions representing the analytical function to be mimized and its gradient vector function. 

## Running the solver
To run the program head over head over in the `solver/main/` directory and type `make` to generate the executable. At that point run it with: `./main`.
To delete object files and executable then run `make clean`.