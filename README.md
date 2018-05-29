Finite volume adjoint solver, 1D Euler
Author: SH Bryngelson

Solves the strong-conservation form of the 1d euler equations
and their continuous adjoint using a finite-volume scheme with Lax-Friedrichs fluxes
and either WENO3, WENO5, or no reconstruction. Forward euler and SSP RK3
time integration schemes are supplied.
