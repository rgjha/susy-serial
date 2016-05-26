# susy-serial

## Installation of the program

It is very easy to perform the installation and execution of `SUSY_LATTICE`.  Below we provide the necessary steps on \*nix systems.
1. Download and unpack the code.
2. Change the directory to `SUSY_LATTICE`.
3. Edit utilities.h to set the number of supercharges by (un)defining Q16, and lattice lengths LX, LY, LZ and T.  Note: setting e.g. LX=1 allows study of corresponding dimensionally reduced model.  Number of colors is given by NCOLOR.  The compilation parameter FULLMATRIX allows one to compute the Pfaffian and fermion eigenvalues using full matrix code (for small lattices).
4. Compile the code:
```
g++ -O *.cpp -o SUSY_LATTICE -llapack -lblas
```
Setting the parameter SCALAREIG=0 in utiliities.h allows compilation without the LAPACK and BLAS libraries at the expense of having no measurements of the (scalar) eigenvalues.
5. Modify the input parameters located in file `parameters`
6. Type `./SUSY_LATTICE &> log &` to run the code.

The output of the code is written to the following files in the working directory:
1. `cgs`: Average number of conjugate gradient (CG) iterations (see `MCG_solver.cpp`).
2. `config`: File to read in containing the site, link and plaquette field configurations from a previous run (see `read_in.cpp`).
3. `corrlines`: Correlation function between Polyakov loops as function of spatial separation (see `corrlines.cpp`).
4. `data`: Boson (1st column) and fermion (2nd column)  contributions to the total action (see `measure.cpp`).
5. `config`: Site, link and plaquette field configurations stored as ASCII (see `write_out.cpp`).
6. `eigenvalues`:  Eigenvalues of `\cUbar_a(x)\cU_a(x)`. `N` real numbers for each lattice point `x` and direction `a` (see `measure.cpp`).
7. `hmc_test`: `e^{-\Delta H}` from HMC test (see `update.cpp`).
8. `lines_s`: Spatial Wilson line (see `measure.cpp`).
9. `lines_t`: Temporal Wilson line, a.k.a. the Polyakov loop (see `measure.cpp`).
10. `log`: Log file.
11. `loops`: Wilson loops (see `loop.cpp`).
12. `scalars`: `Tr[\cUbar_a(x)\cU_a(x)]` (see `measure.cpp`).
13. `ulines_s`: The spatial Wilson line computed using the unitary part of the link (see `measure.cpp`).
14. `ulines_t`: The Polyakov loop computed using the unitary part of the link (see `measure.cpp`).

## The list of files in `SUSY_LATTICE` library

This is the list of files included in `SUSY_LATTICE` with a brief description of their purpose.
1. `action.cpp`: Compute the total action, fermionic and bosonic.
2. `corrlines.cpp`: Finds the traced product of the link matrices at various lattice sites.
3. `evolve_fields.cpp`: Leapfrog evolution algorithm. Also stores the fermion and boson forces for the next iteration. Ratio of fermion to boson time steps controlled by parameter `STEPS`.
4. `fermion_forces.cpp`: Computes the fermion kick to gauge link force.
5. `force.cpp`: Bosonic and pseudo-fermionic contribution to the force.
6. `kinetic_energy.cpp`: Computes the kinetic energy term in the Hamiltonian.
7. `line.cpp`: Computes the Polyakov lines.
8. `loop.cpp`: Computes the Wilson loops.
9. `matrix.cpp`: Builds the fermion matrix (sparse and full forms) and also computes the Pfaffian of the fermion operator.
10. `MCG_solver.cpp`: multi-mass CG solver needed for RHMC algorithm.
11. `measure.cpp`: Performs measurements on field configurations.  Writes out scalar eigenvalues, Polyakov/Wilson loops and the action.
12. `my_gen.cpp`: Computes the $SU(N)$ generator matrices.
13. `obs.cpp`: Computes fermion and gauge actions. Also returns the unitary piece of the complex link field.
14. `read_in.cpp`: Reads in the previously generated field configurations -- file `config`.
15. `read_param.cpp`: Reads in the simulation parameters from a data file called `parameters`.
16. `setup.cpp`: Contains the partial fraction coefficients necessary to represent fractional power of fermion operator -- used by Remez algorithm.
17. `sym.cpp`: The main program - performs warm up on field configurations and commences measurement sweeps once the configurations are warmed up.
18. `unit.cpp`: Extracts the unitary piece of the complex gauge links.
19. `update.cpp`: Updates the field configurations based on RHMC algorithm.
20. `update_gauge_momenta.cpp`: Called by `evolve_fields()` to change gauge momenta.
21. `update_gauge_field.cp`: Called by `evolve_fields()`.
22. `utilities.cpp`: Utility functions. Contains constructors for site, link, plaquette fields, gauge fields, twist fermions etc.  Edit to change number of supercharges and size of lattice dimensions.
23. `write_out.cpp`: Writes out the values of gauge and twist fermion fields.

## A sample input parameter file for `SUSY_LATTICE`

There is a sample input parameter file called `parameters` located in the `SUSY_LATTICE` folder.

These are the definitions of the parameters:
1. `SWEEPS`: Total number of Monte Carlo time steps intended for taking measurement steps.
2. `THERM`: Total number of Monte Carlo time steps intended for thermalizing the field configurations.
3. `GAP`: The gap between measurement steps.
4. `LAMBDA`: The `t Hooft coupling.
5. `BETA`: Inverse temperature.
6. `BMASS`: A bosonic mass parameter
7. `FMASS`: A fermionic mass paramater. Need to set non-zero for PBC.
8. `DT`: The time step put in the integrator for leapfrog evolution.
9. `READIN`: Determines whether to read in the previously generated field configurations or not.  The program will read in the previous configurations if `READIN` is set to `1`.
