#include "action.h"

double action(const Gauge_Field &U, const Twist_Fermion F){
Lattice_Vector x,e_mu,e_nu;
Gauge_Field Udag;
Complex dum,d,trace;
Umatrix dummy;
double act_s=0.0, act_F=0.0;
int mu,nu,site; 
Twist_Fermion sol[DEGREE],psol[DEGREE];
Umatrix Fmunu,DmuUmu;

for(int n=0;n<DEGREE;n++){
sol[n]=Twist_Fermion();
psol[n]=Twist_Fermion();
}

Udag=Adj(U);

// usual pieces
site=0;
while(loop_over_lattice(x,site)){
DmuUmu=Umatrix();
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
DmuUmu=DmuUmu+U.get(x,mu)*Udag.get(x,mu)-
              Udag.get(x,-e_mu,mu)*U.get(x,-e_mu,mu);}
    
act_s=act_s+0.5*Tr(DmuUmu*DmuUmu).real();
}

// Konishi mass term
site=0;
dum=Complex();
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
dummy=U.get(x,mu)*Udag.get(x,mu)-Umatrix(1);
dum=dum+Tr(dummy*dummy);
}}

act_s=act_s+(BMASS*BMASS)*dum.real();

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
for(nu=mu+1;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
Fmunu=
U.get(x,mu)*U.get(x,e_mu,nu)-
U.get(x,nu)*U.get(x,e_nu,mu);
act_s=act_s+2.0*Tr(Fmunu*Adj(Fmunu)).real();

}
}
}

act_s=KAPPA*act_s;

// pseudofermion contribution

if(FERMIONS){

act_F=FF*ampdeg*Tr(Adj(F)*F).real();

MCG_solver(U,F,shift,sol,psol);

for(int n=0;n<DEGREE;n++){
act_F=act_F+FF*amp[n]*Tr(Adj(F)*sol[n]).real();}
}


//cout << "act_F is " << act_F << "\n" << flush;

return(act_s+act_F);
}
