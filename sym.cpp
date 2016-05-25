#include "sym.h"

int SWEEPS,GAP,THERM,SEED,READIN,OLDSWEEPNO;
double KAPPA,DT,TIME,G;

double amp[DEGREE],shift[DEGREE],ampdeg;
Umatrix Lambda[NUMGEN];
<<<<<<< HEAD
double LARGECUT,SMALLCUT,BMASS,MASS;
=======
double LARGECUT,SMALLCUT,BMASS,C1,C2,MASS;
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88
int TRAJECTORY_LENGTH;
double perm[NUMLINK][NUMLINK][NUMLINK][NUMLINK][NUMLINK];
int side[D],SIMULATING,SWEEPNO;
int Lattice_Map[D];
int num_in_row[LEN],BLOCK_MEASURE,BLOCKING;

int main(int argc, char *argv[]){
int sweep;
Gauge_Field U,Up;
Twist_Fermion F;

#ifdef GPU
int gpuid=atoi(argv[1]);
cudaSetDevice(gpuid);
cout << "using GPU " << gpuid << endl;
#endif

read_param();

if(READIN){
read_in(U,F);
}
else{
U=Gauge_Field(0);
F=Twist_Fermion(1);
}

cout << "Warming up" << "\n" << flush;
DT=DT/2;
cout << "DT is " << DT << "\n" << flush;
for(sweep=1;sweep<=THERM/4;sweep++){
clock_t time= clock();
update(U,F);
cout << "sweep time is " << float(clock()-time)/CLOCKS_PER_SEC << endl;
write_out(U,F,0);
}
DT=DT*2;
cout << "DT is " << DT << "\n" << flush;
for(sweep=1;sweep<=(3*THERM)/4;sweep++){
clock_t time=clock();
update(U,F);
cout << "sweep time is " << float(clock()-time)/CLOCKS_PER_SEC << endl;
write_out(U,F,0);
}

cout << "Commencing measurement sweeps" << "\n" << flush;


for(sweep=SWEEPNO+1;sweep<=SWEEPS;sweep++){
clock_t time=clock();
update(U,F);
cout << "sweep time is " << float(clock()-time)/CLOCKS_PER_SEC << endl;
	
	//  measure config
        cout << "sweep no. " << sweep << "\n" << flush;

	if(sweep%GAP==0){
        measure(U,F,sweep);
//        write_out(U,F,sweep);
        write_out(U,F,0);}
        

}
return(0);


}
