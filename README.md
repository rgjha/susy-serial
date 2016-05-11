# susy_serial

Installation of the program

It is very easy to perform the installation and execution of SUSY\_LATTICE. 
Below we provide the necessary steps on Unix or Linux systems.

[1.] Download the code from Github and unpack it. 
[2.] Change the directory to SUSY\_LATTICE.
[3.] Edit utilities.h to set number of supercharges Q16, 
and lattice lengths LX,LY,LZ and T. 
Note: setting e.g. LX=1 allows
study of corresponding dimensionally reduced model. 
Number of colors is given by NCOLOR.
The compilation parameter FULLMATRIX
allows one to compute the Pfaffian and fermion
eigenvalues using full matrix code -- only for small lattices ...).
[4.] Compile the code (g++ -O *.cpp -o SUSY\_LATTICE -llapack -lblas). Setting the parameter SCALAREIG=0 in utiliities.h allows compilation
without the LAPACK and BLAS libraries at the 
expense of having no measurements of the ( scalar) eigenvalues 
[4.] Modify the input parameters located in file {\tt parameters} 
[5.] Type ./SUSY\_LATTICE $>$\& log \& to run the code.


The output of the code produces the following files in the running directory: 

[1.] {\tt cgs:} Average number of conjugate gradient (CG) iterations. (See {\tt MCG\_solver.cpp}).
[2.] {\tt config:} File to read in containing the site, link and plaquette field configurations from a previous run. (See {\tt read\_in.cpp}.)
[3.] {\tt corrlines:} Correlation function between temporal Polyakov lines as function of spatial separation (See {\tt corrlines.cpp}.)
[4.] {\tt data:} Boson (1st column) and fermion (2nd column)  contributions to the total action. (See {\tt measure.cpp}.)
[5.] {\tt config:} Site, link and plaquette field configurations stored as ASCII (See {\tt write\_out.cpp}.)
[6.] {\tt eigenvalues:}  Eigenvalues of $\cU^\dagger_a(\bx)\cU_a(\bx)$. $N$ real numbers for each lattice point $\bx$ and direction $\ba$ (See {\tt measure.cpp}.)
[7.] {\tt hmc\_test:} $e^{-\Delta H}$ from HMC test (See {\tt update.cpp}.)
[8.] {\tt lines\_s:} Spatial Polyakov line. (See {\tt measure.cpp}.)
[9.] {\tt lines\_t:} Temporal Polyakov line. (See {\tt measure.cpp}.)
[10.] {\tt log:} Log file.
[11.] {\tt loops:} Wilson loops. (See {\tt loop.cpp}.)
[12.] {\tt scalars:} ${\rm Tr} \Big(\cU^\dagger_a(\bx)\cU_a(\bx)\Big)$. (See {\tt measure.cpp}.)
[13.] {\tt ulines\_s:} The spatial Polyakov line computed using the unitary part of the link (See {\tt measure.cpp}.)
[14.] {\tt ulines\_t:} The temporal Polyakov line computed using the unitary part of the link (See {\tt measure.cpp}.)


%--------------------------------------------
\section{The list of files in SUSY\_LATTICE library}
%--------------------------------------------
We list the files included in SUSY\_LATTICE with a brief description of their purpose.

[1.]{\tt action.cpp:} Compute the total action - fermionic and bosonic.
%
[2.]{\tt corrlines.cpp:} Finds the traced product of the link matrices at various lattice sites.
%
[3.]{\tt evolve\_fields.cpp:} Leapfrog evolution algorithm. Also stores the fermion and boson forces for the next iteration. Ratio of fermion to boson
time steps controlled by parameter {\tt STEPS}.
%
[4.]{\tt fermion\_forces.cpp:} Computes the fermion kick to gauge link force.
%
[5.]{\tt force.cpp:} Bosonic and pseudo-fermionic contribution to the force.
%
[6.]{\tt kinetic\_energy.cpp:} Computes the kinetic energy term in the Hamiltonian.
%
[7.]{\tt line.cpp:} Computes the Polyakov lines.
%
[8.]{\tt loop.cpp:} Computes the Wilson loops.
%
[9.]{\tt matrix.cpp:} Builds the fermion matrix (sparse and full forms) and also computes the Pfaffian of the fermion operator.
%
[10.]{\tt MCG\_solver.cpp:} multi-mass CG solver needed for RHMC algorithm.
%
[11.]{\tt measure.cpp:} Performs measurements on field configurations. Writes out scalar eigenvalues, Polyakov/Wilson loops and the action.
%
[12.]{\tt my\_gen.cpp:} Computes the $SU(N)$ generator matrices.
%
[13.]{\tt obs.cpp:} Computes fermion and gauge actions. Also returns the unitary piece of the complex link field.
%
[14.]{\tt read\_in.cpp:} Reads in the previously generated field configurations - file {\tt config}.
%
[15.]{\tt read\_param.cpp:} Reads in the simulation parameters from a data file called {\tt parameters}.
%
[16.]{\tt setup.cpp:} Contains the partial fraction coefficients necessary to represent fractional power of fermion operator - used  by Remez algorithm.
%
[17.]{\tt sym.cpp:} The main program - performs warm up on field configurations and commences measurement sweeps once the configurations are warmed up.
%
[18.]{\tt unit.cpp:} Extracts the unitary piece of the complex gauge links.
%
[19.]{\tt update.cpp:} Updates the field configurations based on RHMC algorithm.
%
[20.]{\tt update\_gauge\_momenta.cpp:} Called by {\tt evolve\_fields()} to change gauge momenta.
%
[21.]{\tt update\_gauge\_field.cp:} Called by {\tt evolve\_fields()}.
%
[22.]{\tt utilities.cpp:} Utility functions. Contains constructors for site, link, plaquette fields, gauge fields, twist fermions etc. Edit to change
number of supercharges and size of lattice dimensions.
%
[23.]{\tt write\_out.cpp:} Writes out the values of gauge and twist fermion fields.
%

%--------------------------------------------------------
\section{A sample input parameter file for SUSY\_LATTICE}
%--------------------------------------------------------

This is a sample input parameter file called {\tt parameters} located in the SUSY\_LATTICE folder.

These are the definitions of the parameters:

[1] {\tt SWEEPS:} Total number of Monte Carlo time steps intended for taking measurement steps.
[2.] {\tt THERM:} Total number of Monte Carlo time steps intended for thermalizing the field configurations.
[3.] {\tt GAP:} The gap between measurement steps.
[4.] {\tt LAMBDA:} The `t Hooft coupling.
[5.] {\tt BETA:} Inverse temperature.
[6.] {\tt BMASS:} A bosonic mass parameter
[7.] {\tt FMASS:} A fermionic mass paramater. Need to set non-zero for PBC.
[8.] {\tt DT:} The time step put in the integrator for leapfrog evolution.
[9.] {\tt READIN:} Determines whether to read in the previously generated field configurations or not. 
The program will read in the previous configurations if {\tt READIN} is set to ${\tt 1}$. 

