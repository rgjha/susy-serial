#include "MCG_solver.h"

#include <ctime>
// multimass CG solver 
// can return solution to (M^daggerM +beta_i) x=b for all shifts beta_i

Twist_Fermion plist[DEGREE];
Twist_Fermion r,s,t,p;

void MCG_solver(const Gauge_Field &U,
const Twist_Fermion &b, double shift[], Twist_Fermion sol[DEGREE], 
Twist_Fermion psol[DEGREE]){

static ofstream f_cgs;

Lattice_Vector x;
int n,sites,count;
double alpha1,alpha2,beta0,beta1,rrdot,rrdot2,resid,psdot;

double alphalist[DEGREE],betalist[DEGREE],xi2[DEGREE];
double xi0[DEGREE],xi1[DEGREE];

static int av_count=0,no_calls=0,first_time=1;

	if(first_time){
	f_cgs.open("cgs");
	if(f_cgs.bad()){ 
	cout << "failed to open cgs file\n" << flush ;}

    first_time=0;}

count=0;
no_calls++;
    
beta0=1.0;
alpha1=0.0;

for(n=0;n<DEGREE;n++){
sol[n]=Twist_Fermion();
plist[n]=b;

alphalist[n]=alpha1;
betalist[n]=beta0;
xi0[n]=beta0;
xi1[n]=beta0;
xi2[n]=beta0;
}

r=b;
p=b;

check_trace(p);

double MASS=0.0;

do{

t=Fermion_op(U,p);
s=Adj_Fermion_op(U,t);

//cout << "count is " << count << endl;
//check_trace(s);

s=s+MASS*MASS*p;

rrdot=Tr(Adj(r)*r).real();
psdot=Tr(Adj(s)*p).real();

beta1=-rrdot/psdot;

for(int n=0;n<DEGREE;n++){
xi2[n]=(xi1[n]*xi0[n]*beta0)/
(
beta1*alpha1*(xi0[n]-xi1[n])+
xi0[n]*beta0*(1-shift[n]*beta1)
);
if(fabs(xi2[n])<1.0e-50) {xi2[n]=1.0e-50;}

betalist[n]=beta1*xi2[n]/xi1[n];
}

r=r+beta1*s;

for(n=0;n<DEGREE;n++){
sol[n]=sol[n]-betalist[n]*plist[n];
}

rrdot2=Tr(Adj(r)*r).real();

alpha2=rrdot2/rrdot;

for(n=0;n<DEGREE;n++){
alphalist[n]=alpha2*(xi2[n]/xi1[n])*(betalist[n]/beta1);
}

p=r+alpha2*p;

for(n=0;n<DEGREE;n++){
plist[n]=xi2[n]*r+alphalist[n]*plist[n];
}

resid=sqrt(rrdot2/LEN);

for(n=0;n<DEGREE;n++){
xi0[n]=xi1[n];
xi1[n]=xi2[n];
}
    

alpha1=alpha2;
beta0=beta1;

count++;
//cout << "residual is " << resid << "\n" << flush;
}
while((resid>CG_RESIDUAL)&&(count<(10*LEN)));

av_count+=count;


if(0){
cout<<"exiting residual is " << resid << " at " 
<< count << " iterations\n" <<
flush;}



if(no_calls%100==0){
cout << "average number of CG iterations " <<
(double)av_count/(no_calls) << "\n" << flush;
no_calls=0;
av_count=0;
}

// check solutions


for(n=0;n<DEGREE;n++){

t=Fermion_op(U,sol[n]);
psol[n]=t;
s=Adj_Fermion_op(U,t);

s=s+(MASS*MASS+shift[n])*sol[n];

double tt;
tt=Tr(Adj(s-b)*(s-b)).real();

if(sqrt(tt)/LEN>0.0000001){cout << "poor soln\n" << sqrt(tt)/LEN << endl;}

}


return;
}
