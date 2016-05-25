
 #include "CG_solver.h"
#include "gpusolver2.h"



// CG solver 
void CG_solver(const Gauge_Field &U, 
	       const Twist_Fermion &rhs, Twist_Fermion &sol){

int count2,i;
static int av_count2=0,no_calls=0;
double alpha1,alpha2,beta0,beta1,rrtmp,rrtmp2,resid,psdot;
Adjoint_Links V;
Complex rn[LEN],bn[LEN],pn[LEN];
Complex tn[LEN],sn[LEN],soln[LEN], solnGPU[LEN];
Complex m[LEN*NONZEROES];
int col[LEN*NONZEROES],row[LEN+1];

no_calls++;
compute_Adjoint_Links(U,V);

count2=0;

// initial guess


build_vector(rhs,bn);
build_sparse_matrix(V,m,col,row);

sparse_mult(m,col,row,-1,bn,sn);
for(i=0;i<LEN;i++){
bn[i]=sn[i];}

#ifdef GPU
clock_t begin_time = clock();
gpusolver2(m,col,row,bn,solnGPU);
for(i=0;i<LEN;i++){
soln[i]=solnGPU[i];}
<<<<<<< HEAD
//cout << "GPU_TIME_CG: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
=======
cout << "GPU_TIME_CG: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << endl;
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88
#endif

#ifndef GPU
no_calls++;
count2=0;
clock_t begin_time = clock();

// initialize solver

beta0=1;
alpha1=0.0;


for(i=0;i<LEN;i++){
soln[i]=Complex();
}


for(i=0;i<LEN;i++){
rn[i]=bn[i];
pn[i]=bn[i];}

do{

sparse_mult(m,col,row,1,pn,tn);
sparse_mult(m,col,row,-1,tn,sn);


rrtmp=dot(rn,rn).real();
psdot=dot(sn,pn).real();

beta1=-rrtmp/psdot;

for(i=0;i<LEN;i++){
rn[i]=rn[i]+beta1*sn[i];}


for(i=0;i<LEN;i++){
soln[i]=soln[i]-beta1*pn[i];}

rrtmp2=dot(rn,rn).real();

alpha2=rrtmp2/rrtmp;

for(i=0;i<LEN;i++){
pn[i]=rn[i]+alpha2*pn[i];}

resid=sqrt(rrtmp2/(LEN));

alpha1=alpha2;
beta0=beta1;
count2++;

//cout << "residual is " << resid << "\n" << flush;
}
while((resid>CG_RESIDUAL)&&(count2<(2*LEN)));
<<<<<<< HEAD
//cout << "CPU_TIME_CG " <<  double(clock()-begin_time)/CLOCKS_PER_SEC << endl; 
=======
cout << "CPU_TIME_CG " <<  double(clock()-begin_time)/CLOCKS_PER_SEC << endl; 
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88

#endif
    
av_count2+=count2;


if(0){
cout<<"exiting residual in propagator " << resid << " at " 
<< count2 << " iterations\n" <<
flush;}


if(no_calls%10==0){
#ifndef GPU
<<<<<<< HEAD
//cout << "average number of CG iterations " <<
//(double)av_count2/(no_calls) << "\n" << flush;
=======
cout << "average number of CG iterations " <<
(double)av_count2/(no_calls) << "\n" << flush;
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88

no_calls=0;

av_count2=0;
#endif
}



// extract vector
sol=extract_vector(soln);


// check solutions


for(i=0;i<LEN;i++){
pn[i]=soln[i];}

sparse_mult(m,col,row,1,pn,tn);
sparse_mult(m,col,row,-1,tn,sn);

for(i=0;i<LEN;i++){
sn[i]=sn[i]+MASS*MASS*pn[i];}

double t=0.0;
for(i=0;i<LEN;i++){
t+=((bn[i]-sn[i])*conjug(bn[i]-sn[i])).real();}

if(sqrt(t)/LEN>0.0000001){cout << "poor soln\n" << sqrt(t)/LEN << endl;}




return;
}
