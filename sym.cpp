#include "sym.h"

<<<<<<< HEAD
=======

// This calls : update_o.cpp, measure.cpp, write_out.cpp //  

>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
int SWEEPS,GAP,THERM,SEED,READIN,OLDSWEEPNO;
double KAPPA,DT,TIME,G;

double amp[DEGREE],shift[DEGREE],ampdeg;
Umatrix Lambda[NUMGEN];
<<<<<<< HEAD
double LARGECUT,SMALLCUT,BMASS,C1,C2,MASS;
int TRAJECTORY_LENGTH;
=======
double LARGECUT,SMALLCUT,BMASS,C2,MASS;
int TRAJECTORY_LENGTH,TOTALNONZEROES;
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
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

<<<<<<< HEAD
read_param();

if(READIN){
read_in(U,F);
}
else{
U=Gauge_Field(0);
=======
read_param();            // READ PARAMETERS // 

if(READIN){              
read_in(U,F);
}
else{
U=Gauge_Field(0);        // IF READIN = 0  i.e not reading configs // 
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
F=Twist_Fermion(1);
}

cout << "Warming up" << "\n" << flush;
<<<<<<< HEAD
DT=DT/2;
=======
DT=DT/10;
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
cout << "DT is " << DT << "\n" << flush;
for(sweep=1;sweep<=THERM/4;sweep++){
clock_t time= clock();
update(U,F);
cout << "sweep time is " << float(clock()-time)/CLOCKS_PER_SEC << endl;
write_out(U,F,0);
}
<<<<<<< HEAD
DT=DT*2;
=======
DT=DT*10;
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
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
