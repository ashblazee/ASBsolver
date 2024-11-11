ASBsolver is a numerical solver developed for solving ordinary and partial differential equations.

A prerequisite to using the solver would be to understand the first principle of calculus and discretization methods 
Since ASBsolver employs the Finite Difference Method, it is necessary to understand how differential equations are 
discretized using this method. https://en.wikipedia.org/wiki/Finite_difference_method

To use the solver, a finite difference equation must be first obtained, and the value of the next timestep must be made 
as the subject of the equation.
This will be the function that is passed into the solver function

An object is instantiated from the FiniteDifferenceMethod class and its boundary and intial conditions are set, using
the setter methods.

The step size, skip size and range are also specified.

The skip size specifies how often the data is recorded and stored in memory. This is done so as to conserve memory,
because it is not necessary to obtain data every timestep especially if the timestep is very small.

More than 1 initial, boundary and range can be specified, but the number of specifications has to be appropriate to the 
nature of the equation, or the solution might be inaccurate or the program would result in a runtime error.

IMPORTANT: For Ordinary Differential Equations, the boundary conditions are IGNORED by the solver and only the 
initial conditions are computed. An invalid initial condition given will result in a runtime error
Range is added in the form of a 2-element array

The FiniteDifferenceMethod class provides ODE and PDE solvers of different nature and dimensions, and the appropriate
function is to be chosen for each equation by the user.

