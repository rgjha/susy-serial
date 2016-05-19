#include "read_param.h"

#include <stdio.h>
#include <sys/time.h>

unsigned int random_seed()
{

 unsigned int seed;
 struct timeval tv;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
   gettimeofday(&tv,0);
   seed = tv.tv_sec + tv.tv_usec;
   printf("Got seed %u from gettimeofday()\n",seed);
 } else {
   fread(&seed,sizeof(seed),1,devrandom);
   printf("Got seed %u from /dev/random\n",seed);
   fclose(devrandom);
 }

 return(seed);

}



void read_param(void){
double LAMBDA,BETA;
ifstream f_in("parameters");
if(!f_in.good()){
	cout << "\ncan't open file parameters to read data!\n";
	exit(1);}
<<<<<<< HEAD
f_in>>SWEEPS>>THERM>>GAP>>BETA>>G>>BMASS>>C1>>C2>>DT>>READIN >>SWEEPNO ;

// assume spatial size very large for this scaling
LAMBDA=1.0;
=======
f_in>>SWEEPS>>THERM>>GAP>>LAMBDA>>G>>BMASS>>C2>>DT>>READIN >>SWEEPNO>>BETA ;

// assume spatial size very large for this scaling
//BETA=1.0;
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0


if(D==4){
KAPPA=(NCOLOR*0.5)/LAMBDA;

<<<<<<< HEAD
if((LX==1)&&(LY==1)&&(LZ==1)&&(T!=1)){
=======
if((LX==1)&&(LY==1)&&(LZ!=1)&&(T!=1)){
KAPPA = KAPPA*(T*T)/(BETA*BETA);}

// Added by Raghav // 
/*if((LX==1)&&(LY==1)&&(LZ==1)&&(T!=1)){
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
KAPPA=KAPPA*(T*T*T)/(BETA*BETA*BETA);}

if((LX==1)&&(LY==1)&&(LZ!=1)&&(T!=1)){
KAPPA=KAPPA*(T*T)/(BETA*BETA);}

if((LX==1)&&(LY!=1)&&(LZ!=1)&&(T!=1)){
<<<<<<< HEAD
KAPPA=KAPPA*T/BETA;}
=======
KAPPA=KAPPA*T/BETA;}*/
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

}

if(D==2){
KAPPA=(NCOLOR*0.5*T*T)/(LAMBDA*BETA*BETA);

if((LX==1)&&(T!=1)){
KAPPA=KAPPA*T/BETA;}

}

<<<<<<< HEAD
BMASS/=T;
=======

>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

TRAJECTORY_LENGTH=(int)(0.5/DT);


if((FERMIONS==1)&&(NUMLINK==5)){cout << "16 supercharge theory \n";}
if((FERMIONS==1)&&(NUMLINK==2)){cout << "4 supercharge theory\n";}
if((FERMIONS==0)&&(NUMLINK==5)){cout << "Quenched 16 supercharge theory \n";}
if((FERMIONS==0)&&(NUMLINK==2)){cout << "Quenched 4 supercharge theory\n";}


cout << "Number of colors " << NCOLOR <<  "\n";
cout << "Temporal extent " << T << "\n";
if(D==2){
cout << "Spatial extent " << LX << "\n";}
if(D==4){
cout << "Spatial extent " << LX << "\t" << LY << "\t" << LZ << "\n";}

cout << "Dimensionless t'Hooft coupling " << LAMBDA << "\n";
<<<<<<< HEAD
cout << "Lattice Coupling " << KAPPA << "\n";
cout << "Boson Mass " << BMASS << "\n";
cout << "C1 coeff " << C1 << "\n";
cout << "C2 coeff " << C2 << "\n";
cout << "Coupling to det " << G << "\n";
=======
cout << "C2 coeff " << C2 << "\n";
cout << "Coupling to det " << G << "\n";
cout << "Konishi mass " << BMASS << "\n";
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
cout << "Thermalization sweeps " << THERM << "\n";
cout << "Number of sweeps " << SWEEPS << "\n";
cout << "Gap between measurements " << GAP << "\n";
cout << "Time step in leapfrog eqs " << DT << "\n";
cout << "Trajectory length " << TRAJECTORY_LENGTH << "\n";
<<<<<<< HEAD
=======
cout << "Inverse temperature (BETA)" << BETA << "\n"; 
cout << "Kappa (dimensionless lattice governing constant) " << KAPPA << "\n"; 
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
cout << "Minimax approx degree " << DEGREE << "\n";
cout << "Reading initial config: (1 for yes, 0 for no) " << READIN << "\n";
//cout << "Old sweep number: " << OLDSWEEPNO << "\n";

if (PBC==1.0) {cout << "periodic temporal bc for fermions" << "\n";}
else{cout << "antiperiodic temporal bc for fermions" << "\n";}

<<<<<<< HEAD
=======
//G=G/sqrt(KAPPA);
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

srand(random_seed());
//srand(0);
setup();

return;
}
