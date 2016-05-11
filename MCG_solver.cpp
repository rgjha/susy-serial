#include "MCG_solver.h"
#include "gpusolver.h"

#include <ctime>
// multimass CG solver 
// can return solution to (M^daggerM +beta_i) x=b for all shifts beta_i

Complex rn[LEN],bn[LEN],plist[DEGREE][LEN],pn[LEN];
Complex tn[LEN],sn[LEN],soln[DEGREE][LEN], solnGPU[DEGREE][LEN];
Complex m[LEN*NONZEROES];
int col[LEN*NONZEROES],row[LEN+1];

void MCG_solver(const Gauge_Field &U,
const Twist_Fermion &b, double shift[], Twist_Fermion sol[DEGREE], 
Twist_Fermion psol[DEGREE]){

static ofstream f_cgs;

Adjoint_Links V;
Twist_Fermion p,t,s;
Lattice_Vector x;
int sites,A;
double alpha1,alpha2,beta0,beta1,rrtmp,rrtmp2,resid,psdot;

double alphalist[DEGREE],betalist[DEGREE],xi2[DEGREE];
double xi0[DEGREE],xi1[DEGREE];

int count2;
static int av_count2=0,no_calls=0,first_time=1;
int i,n;
double dummy;

	if(first_time){
	f_cgs.open("cgs");
	if(f_cgs.bad()){ 
	cout << "failed to open cgs file\n" << flush ;}

        first_time=0;}


compute_Adjoint_Links(U,V);
MASS=0.0;

build_vector(b,bn);
build_sparse_matrix(V,U,m,col,row);


#ifdef GPU
clock_t begin_time = clock();
gpusolver(m,col,row,bn,shift,solnGPU);
    for(n=0;n<DEGREE;n++){
        for(i=0;i<LEN;i++){
            soln[n][i]=solnGPU[n][i];}}
cout << "GPU_TIME_MCG: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
#endif

#ifndef GPU
no_calls++;
count2=0;
clock_t begin_time = clock();

// initialize solver

beta0=1;
alpha1=0.0;

for(n=0;n<DEGREE;n++){

for(i=0;i<LEN;i++){
soln[n][i]=Complex();
plist[n][i]=bn[i];}

alphalist[n]=alpha1;
betalist[n]=beta0;
xi0[n]=beta0;
xi1[n]=beta0;
xi2[n]=beta0;
}

for(i=0;i<LEN;i++){
rn[i]=bn[i];
pn[i]=bn[i];}

do{

sparse_mult(m,col,row,1,pn,tn);
sparse_mult(m,col,row,-1,tn,sn);
for(i=0;i<LEN;i++){
sn[i]=sn[i]+MASS*MASS*pn[i];}

rrtmp=dot(rn,rn).real();
psdot=dot(sn,pn).real();

beta1=-rrtmp/psdot;

for(n=0;n<DEGREE;n++){
xi2[n]=(xi1[n]*xi0[n]*beta0)/
 (
beta1*alpha1*(xi0[n]-xi1[n])+
xi0[n]*beta0*(1-shift[n]*beta1)
 );
if(fabs(xi2[n])<1.0e-50){xi2[n]=1.0e-50;}

betalist[n]=beta1*xi2[n]/xi1[n];
}


for(i=0;i<LEN;i++){
rn[i]=rn[i]+beta1*sn[i];}

for(n=0;n<DEGREE;n++){

for(i=0;i<LEN;i++){
soln[n][i]=soln[n][i]-betalist[n]*plist[n][i];}
}


rrtmp2=dot(rn,rn).real();

alpha2=rrtmp2/rrtmp;

for(n=0;n<DEGREE;n++){
alphalist[n]=alpha2*(xi2[n]/xi1[n])*(betalist[n]/beta1);
}

for(i=0;i<LEN;i++){
pn[i]=rn[i]+alpha2*pn[i];}

for(n=0;n<DEGREE;n++){
for(i=0;i<LEN;i++){
plist[n][i]=xi2[n]*rn[i]+alphalist[n]*plist[n][i];}
}

resid=sqrt(rrtmp2/(LEN));

for(n=0;n<DEGREE;n++){
xi0[n]=xi1[n];
xi1[n]=xi2[n];
}

alpha1=alpha2;
beta0=beta1;
count2++;

//cout << "residual is " << resid << "\n" << flush;
}
while((resid>CG_RESIDUAL)&&(count2<(2*LEN)));
cout << "CPU_TIME_MCG " <<  double(clock()-begin_time)/CLOCKS_PER_SEC << endl; 

#endif
    
av_count2+=count2;


if(0){
cout<<"exiting residual is " << resid << " at " 
<< count2 << " iterations\n" <<
flush;}


if(no_calls%10==0){
#ifndef GPU
cout << "average number of CG iterations " <<
(double)av_count2/(no_calls) << "\n" << flush;

f_cgs  << (double)av_count2/(no_calls) << "\n" << flush;

no_calls=0;

av_count2=0;
#endif
}



// extract vector

for(n=0;n<DEGREE;n++){
for(i=0;i<LEN;i++){
pn[i]=soln[n][i];
}
sparse_mult(m,col,row,1,pn,tn);
sol[n]=extract_vector(pn);
psol[n]=extract_vector(tn);
}

// check solutions

for(n=0;n<DEGREE;n++){
for(i=0;i<LEN;i++){
pn[i]=soln[n][i];}

sparse_mult(m,col,row,1,pn,tn);
sparse_mult(m,col,row,-1,tn,sn);
for(i=0;i<LEN;i++){
sn[i]=sn[i]+MASS*MASS*pn[i];}

for(i=0;i<LEN;i++){
sn[i]=sn[i]+shift[n]*pn[i];}

double t=0.0;
for(i=0;i<LEN;i++){
t+=((bn[i]-sn[i])*conjug(bn[i]-sn[i])).real();}

if(sqrt(t)/LEN>0.0000001){cout << "poor soln\n" << sqrt(t)/LEN << endl;}

}



return;
}
