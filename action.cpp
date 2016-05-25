#include "action.h"

double action(const Gauge_Field &U, const Twist_Fermion F){
Lattice_Vector x,e_mu,e_nu;
Gauge_Field Udag;
double dum, act, act1;
Complex d,trace;
Umatrix p,udum;
Umatrix dummy;

double act_s=0.0, act_F=0.0;
int mu,nu,site; 
Twist_Fermion sol[DEGREE],psol[DEGREE];
Umatrix DmuUmu,Fmunu;

Udag=Adj(U);
  
// boson action 

//SIMON: beta term modification
site=0;
while(loop_over_lattice(x,site)){
DmuUmu=Umatrix();
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
DmuUmu=DmuUmu+U.get(x,mu)*Udag.get(x,mu)-Udag.get(x-e_mu,mu)*U.get(x-e_mu,mu);
}


act_s=act_s+0.5*C2*Tr(DmuUmu*DmuUmu).real();


}

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
for(nu=mu+1;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
Fmunu=U.get(x,mu)*U.get(x+e_mu,nu)-U.get(x,nu)*U.get(x+e_nu,mu);
act_s=act_s+2.0*Tr(Fmunu*Adj(Fmunu)).real();


}
}
}

cout << "  GAUGE  " << KAPPA*act_s << "\n" << flush;
act=KAPPA*act_s; 


//act_s=KAPPA*act_s;

/* mass term for U(1) mode*/
site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){

udum=U.get(x,mu)*Udag.get(x,mu)-Umatrix(1);
act_s+= BMASS*BMASS*Tr(udum*udum).real();



//dum=(1.0/NCOLOR)*Tr(Udag.get(x,mu)*U.get(x,mu)).real()-1.0;
//act_s=act_s+(BMASS*BMASS)*dum*dum;
}}

cout << "  BMASS  " << KAPPA*act_s - act<< "\n" << flush;
act1 = KAPPA*act_s; 

act_s=KAPPA*act_s;


site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=mu+1;nu<NUMLINK;nu++){
e_mu=Lattice_Vector(mu);
e_nu=Lattice_Vector(nu);
p=U.get(x,mu)*U.get(x+e_mu,nu)*Udag.get(x+e_nu,mu)*Udag.get(x,nu);
d=det(p)-Complex(1.0,0.0);
act_s=act_s+G*(d*conjug(d)).real();
}}}


if(act_s == act1){
cout << "  DET  " << 0.00 << "\n" << flush;	
}
else {
cout << "   DET  " << act_s - act - act1 << "\n" << flush;	
}


//act_s=KAPPA*act_s;

// pseudofermion contribution


if(FERMIONS){

act_F=ampdeg*(Cjg(F)*F).real();

MCG_solver(U,F,shift,sol,psol);

for(int n=0;n<DEGREE;n++){
act_F=act_F+amp[n]*(Cjg(F)*sol[n]).real();}
}


cout << "  FERMIONS (EXCEPT MOM) " << act_F << "\n" << flush;
cout << "####################################" << endl ; 

return(act_s+act_F);
}
