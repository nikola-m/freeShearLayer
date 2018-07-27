# freeShearLayer

A CFD code based on the Finite Volume Method discretisation of Navier-Stokes equations for simulation of compressible shear layer.

The code was originally written by Andrei Chernousov back in early 2000's and published on his website (still available [here](http://www.geocities.ws/MotorCity/Pit/9939/freecfd.htm) ). The basic code had several versions, because each simulation case (mesh, initial and boundary conditions, etc.) was individually hardcoded. There were no efforts as I am aware to further develop these so this repository has intention of changing that. We have already made some basic changes to the code, on the side of the algorithm itself and regarding the organisation of the project.

Further we will include other cases and make code more general.

For now, the algorithm is based on Finite Volume Method discretisation of Navier-Stokes equations on Cartesian grids. Piecewise parabolic method (Woodward-Collela JCP 1984) is used for reconstruction of face values of variables. Explicit TVD Runge-Kutta methods are used for time-steping. Approximate-Riemann solvers are used for inviscid fluxes. Ghost cell approach is used to impose boundary conditions (good choice for the explicit method). Postprocesing files are written in Tecplot format. Turbulence is modelled using LES models, namely Smagorinsky with wall damping by Van Driest function is used, where applicable. This will be the focus of upgrade in the future as we plan to implement more LES models.

### Programing language

Code is written in ANSI C language.

### Parallelisation

Available soon, using MPI.

### Basic usage:

```
cd src
mkdir obj
make
cd ..
./run
```

Input file is given with .ini extension and can be easily modified.

Simulation snapshot of the Vorticity magnitude isosurface:

![Free Shear Layer Vorticity Snapshot](https://github.com/nikola-m/freeShearLayer/blob/master/vorticity-nsteps3000.png)


