# Compressible_Flow_Relations

This repository contains three functions that can be used to solve compressible flow relations in MATLAB.  They are based directly on the code from the [*Compressible Aerodynamics Calculator*](http://www.dept.aoe.vt.edu/~devenpor/aoe3114/calc.html), created by William J. Devenport from the Department of Aerospace and Ocean Engineering at Virginia Tech.

The code has been formatted for use in MATLAB.  Usage information can be found in the description at the beginning of the individual functions.  This can also be accessed by typing ```help NORMAL_SHOCK```, for example.  The default output is a structure of all the solution variables.  If you only need one solution variable, this can also be specified as the optional argument.

## Example

Let's use the normal shock function in this example.  The other two functions operate in exactly the same way, just with different input and output arguments.  In this example, we have an upstream Mach number of 2.4, and a specific heat ratio (gamma) of 1.4.  First, let's say we want all the output variables at the same time.  We will type:

```
sol = NORMAL_SHOCK(2.4,'M1',1.4);
```

The output will be a structure with all the variables as individual names labeled according to the OUTPUTS section in the description.  So if we want to access the static pressure ratio in the solution structure, we can type the following:

```
sol.P2P1
```

The output will be the number 6.5533.

Now let's say we only wanted the output of the static pressure ratio for the same inputs (i.e. we didn't want the entire solution structure).  We can simply specify the output we want in the following way:

```
sol = NORMAL_SHOCK(2.4,'M1',1.4,'P2P1');
```

This formulation of the function call gives us the value ```sol = 6.5533```.

