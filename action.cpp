#include "action.h"

double action(const Gauge_Field &U, const Twist_Fermion F){
Lattice_Vector x,e_mu,e_nu;
Gauge_Field Udag;
<<<<<<< HEAD
double dum;
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
DmuUmu=DmuUmu+(1.0+C1)*U.get(x,mu)*Udag.get(x,mu)-Udag.get(x-e_mu,mu)*U.get(x-e_mu,mu);
trace=Tr(U.get(x,mu)*Udag.get(x,mu));
DmuUmu=DmuUmu-(C1/NCOLOR)*trace*Umatrix(1);
}


act_s=act_s+0.5*C2*Tr(DmuUmu*DmuUmu).real();

}
=======
Complex dum,d,trace;
Umatrix p,dummy;
USite_Field DmuUmu;
UPlaq_Field P;
double act_s=0.0, act_F=0.0;
int mu,nu,site; 
Twist_Fermion sol[DEGREE],psol[DEGREE];
Umatrix Fmunu;


Udag=Adj(U);


P=Plaq(U);
site=0;
while(loop_over_lattice(x,site)){
dum=Complex();
for(mu=0;mu<NUMLINK;mu++){
for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
dum=dum+(det(P.get(x,mu,nu))-Complex(1.0,0.0));}}
DmuUmu.set(x,G*dum*Umatrix(1));
}

// usual piece
site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
DmuUmu.set(x,DmuUmu.get(x)+U.get(x,mu)*Udag.get(x,mu)-
                           Udag.get(x-e_mu,mu)*U.get(x-e_mu,mu));}
    
act_s=act_s+0.5*C2*Tr(DmuUmu.get(x)*DmuUmu.get(x)).real();
}


// Konishi mass term - single trace operator - 'eig'  // 


site=0;
dum=Complex();
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){        
dummy=U.get(x,mu)*Udag.get(x,mu)-Umatrix(1);
dum=dum+Tr(dummy*dummy);
}}

act_s=act_s+(BMASS*BMASS)*dum.real();


// Below : Double trace mass term - "trace" // 

/* site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
dum=(1.0/NCOLOR)*Tr(Udag.get(x,mu)*U.get(x,mu)).real()-1.0;
act_s=act_s+(BMASS*BMASS)*dum*dum;
}}  */ 

>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0


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

<<<<<<< HEAD
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
=======
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

act_s=KAPPA*act_s;


<<<<<<< HEAD
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

//act_s=KAPPA*act_s;

// pseudofermion contribution
=======
// Now, the fermions // 

// Pseudofermion contribution
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0


if(FERMIONS){

act_F=ampdeg*(Cjg(F)*F).real();

MCG_solver(U,F,shift,sol,psol);

for(int n=0;n<DEGREE;n++){
act_F=act_F+amp[n]*(Cjg(F)*sol[n]).real();}
}


<<<<<<< HEAD
//cout << "act_F is " << act_F << "\n" << flush;
=======
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

return(act_s+act_F);
}
