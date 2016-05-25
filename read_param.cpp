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

<<<<<<< HEAD
void read_param(void){
double LAMBDA,BETA,GAMMA;
=======


void read_param(void){
double LAMBDA,BETA;
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88
ifstream f_in("parameters");
if(!f_in.good()){
	cout << "\ncan't open file parameters to read data!\n";
	exit(1);}
<<<<<<< HEAD
f_in>>SWEEPS>>THERM>>GAP>>BETA>>GAMMA>>DT>>READIN >>SWEEPNO ;

// BETA is Rt (BETA in units of root lambda)
// GAMMA is MASS in units of root lambda
// assume spatial size very large for this scaling
LAMBDA=1.0;
G=0.0;
=======
f_in>>SWEEPS>>THERM>>GAP>>BETA>>G>>BMASS>>C1>>C2>>DT>>READIN >>SWEEPNO ;

// assume spatial size very large for this scaling
LAMBDA=1.0;

>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88

if(D==4){
KAPPA=(NCOLOR*0.5)/LAMBDA;

if((LX==1)&&(LY==1)&&(LZ==1)&&(T!=1)){
KAPPA=KAPPA*(T*T*T)/(BETA*BETA*BETA);}

if((LX==1)&&(LY==1)&&(LZ!=1)&&(T!=1)){
KAPPA=KAPPA*(T*T)/(BETA*BETA);}

if((LX==1)&&(LY!=1)&&(LZ!=1)&&(T!=1)){
KAPPA=KAPPA*T/BETA;}

}

if(D==2){
KAPPA=(NCOLOR*0.5*T*T)/(LAMBDA*BETA*BETA);

if((LX==1)&&(T!=1)){
KAPPA=KAPPA*T/BETA;}

}

<<<<<<< HEAD
BMASS=GAMMA*BETA/T;    // Added : May 24 by Simon // 
=======
BMASS/=T;
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88

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

<<<<<<< HEAD
cout << "Inverse temperature in units of root lambda " << BETA << "\n";
cout << "Lattice Coupling " << KAPPA << "\n";
cout << "Mass in units of root lambda " << GAMMA << "\n";
cout << "Lattice scalar mass squared " << KAPPA*BMASS*BMASS << "\n";
=======
cout << "Dimensionless t'Hooft coupling " << LAMBDA << "\n";
cout << "Lattice Coupling " << KAPPA << "\n";
cout << "Boson Mass " << BMASS << "\n";
cout << "C1 coeff " << C1 << "\n";
cout << "C2 coeff " << C2 << "\n";
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88
cout << "Coupling to det " << G << "\n";
cout << "Thermalization sweeps " << THERM << "\n";
cout << "Number of sweeps " << SWEEPS << "\n";
cout << "Gap between measurements " << GAP << "\n";
cout << "Time step in leapfrog eqs " << DT << "\n";
cout << "Trajectory length " << TRAJECTORY_LENGTH << "\n";
cout << "Minimax approx degree " << DEGREE << "\n";
cout << "Reading initial config: (1 for yes, 0 for no) " << READIN << "\n";
//cout << "Old sweep number: " << OLDSWEEPNO << "\n";

if (PBC==1.0) {cout << "periodic temporal bc for fermions" << "\n";}
else{cout << "antiperiodic temporal bc for fermions" << "\n";}


srand(random_seed());
//srand(0);
setup();

return;
}
