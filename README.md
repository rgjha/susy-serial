Installation of the program

It is very easy to perform the installation and execution of SUSY_LATTICE. Below we provide the necessary steps on Unix or 
Linux systems.

[1.] Download the code from Github and unpack it. 

[2.] Change the directory to SUSY_LATTICE. 

[3.] Edit utilities.h to set number of supercharges Q16, and lattice lengths LX,LY,LZ and T. 
Note: setting e.g. LX=1 allows study of corresponding dimensionally reduced model. Number of colors is given by NCOLOR. 
The compilation parameter FULLMATRIX allows one to compute the Pfaffian and fermion eigenvalues using full matrix code -- 
only for small lattices ...). 

[4.] Compile the code (g++ -O *.cpp -o SUSY_LATTICE -llapack -lblas). Setting the parameter 
SCALAREIG=0 in utiliities.h allows compilation without the LAPACK and BLAS libraries at the expense of having no measurements
of the ( scalar) eigenvalues [4.] Modify the input parameters located in file {\tt parameters} 

[5.] Type ./SUSY_LATTICE $>$& log & to run the code.
