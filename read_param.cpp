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
f_in>>SWEEPS>>THERM>>GAP>>LAMBDA>>BMASS>>DT>>READIN ;

// assume spatial size very large for this scaling
BETA=1.0;

if(D==4){
KAPPA=(NCOLOR*0.5)/LAMBDA;
}

if(D==2){
KAPPA=(NCOLOR*0.5*T*T)/(LAMBDA*BETA*BETA);
}


TRAJECTORY_LENGTH=(int)(1.0/DT);


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
cout << "Konishi mass " << BMASS << "\n";
cout << "Thermalization sweeps " << THERM << "\n";
cout << "Number of sweeps " << SWEEPS << "\n";
cout << "Gap between measurements " << GAP << "\n";
cout << "Time step in leapfrog eqs " << DT << "\n";
cout << "Trajectory length " << TRAJECTORY_LENGTH << "\n";
cout << "Remez approx degree " << DEGREE << "\n";
cout << "Reading initial config: (1 for yes, 0 for no) " << READIN << "\n";

if (PBC==1.0) {cout << "periodic temporal bc for fermions" << "\n";}
else{cout << "antiperiodic temporal bc for fermions" << "\n";}
if (TWIST==1) {cout << "twisted bcs" << endl;}

srand(random_seed());
//srand(0);
setup();

return;
}
