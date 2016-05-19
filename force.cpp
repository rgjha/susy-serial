#include "force.h"

void force(const Gauge_Field &U, Gauge_Field &f_U, 
const Twist_Fermion &F, Twist_Fermion &f_F, int fermion){
Lattice_Vector x,e_mu,e_nu;
int sites,mu,nu,i,j,a;
Twist_Fermion sol[DEGREE],psol[DEGREE];
Gauge_Field Udag,utmp,U2,UdU,ft;
USite_Field DmuUmu;
<<<<<<< HEAD
UPlaq_Field Fmunu;
Twist_Fermion ptmp,stmp;
Umatrix tmp,tmp2,ftmp;
Complex trace;

// gauge force - boson contribution
=======
UPlaq_Field P,Fmunu;
Twist_Fermion ptmp,stmp;
Umatrix tmp2,ftmp;
Complex d,dum,tmp;
Scalar_Plaquette ZZ,WW;

>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

Udag=Adj(U);
DmuUmu=USite_Field();
Fmunu=UPlaq_Field();

f_U=Gauge_Field();

<<<<<<< HEAD
//SIMON: beta term modification
=======
Udag=Adj(U);
   

P=Plaq(U);
sites=0;
while(loop_over_lattice(x,sites)){
dum=Complex();
for(mu=0;mu<NUMLINK;mu++){
for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
dum=dum+(det(P.get(x,mu,nu))-Complex(1.0,0.0));}}

DmuUmu.set(x,G*dum*Umatrix(1));
}


>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
<<<<<<< HEAD
trace=Tr(U.get(x,mu)*Udag.get(x,mu));
DmuUmu.set(x,DmuUmu.get(x)+(1.0+C1)*U.get(x,mu)*Udag.get(x,mu)-
Udag.get(x-e_mu,mu)*U.get(x-e_mu,mu)-(C1/NCOLOR)*trace*Umatrix(1));}
=======
DmuUmu.set(x,DmuUmu.get(x)+U.get(x,mu)*Udag.get(x,mu)-
                           Udag.get(x-e_mu,mu)*U.get(x-e_mu,mu));}
}

// Derivative of det term // 


Scalar_Plaquette Z;
Complex W;
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
W=det(P.get(x,mu,nu))-Complex(1.0,0.0);
WW.set(x,mu,nu,W);
Z.set(x,mu,nu,Tr(DmuUmu.get(x)));
}}}

ZZ=mydiff(Z,WW);
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){

tmp=Complex();
for(nu=0;nu<NUMLINK;nu++){
if(nu==mu) continue;
tmp=tmp+ZZ.get(x,mu,nu);
}

f_U.set(x,mu,f_U.get(x,mu)+G*tmp*inverse(U.get(x,mu)));
}
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

<<<<<<< HEAD
f_U.set(x,mu,f_U.get(x,mu)+(1.0+C1)*U.get(x,mu)*Udag.get(x,mu)*DmuUmu.get(x));
f_U.set(x,mu,f_U.get(x,mu)-U.get(x,mu)*DmuUmu.get(x+e_mu)*Udag.get(x,mu));
f_U.set(x,mu,f_U.get(x,mu)-(C1/NCOLOR)*Tr(DmuUmu.get(x))*Udag.get(x,mu));
=======
f_U.set(x,mu,f_U.get(x,mu)+Udag.get(x,mu)*DmuUmu.get(x));
f_U.set(x,mu,f_U.get(x,mu)-DmuUmu.get(x+e_mu)*Udag.get(x,mu));
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
f_U.set(x,mu,C2*f_U.get(x,mu));
}

}

			   
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
for(nu=mu+1;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
Fmunu.set(x,mu,nu,U.get(x,mu)*U.get(x+e_mu,nu)-U.get(x,nu)*U.get(x+e_nu,mu));
Fmunu.set(x,nu,mu,-1.0*Fmunu.get(x,mu,nu));
}
}
}


sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

for(nu=0;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
if(mu==nu) continue;

<<<<<<< HEAD
f_U.set(x,mu,f_U.get(x,mu)+2.0*U.get(x,mu)*U.get(x+e_mu,nu)*Adj(Fmunu.get(x,mu,nu)));
f_U.set(x,mu,f_U.get(x,mu)-2.0*U.get(x,mu)*Adj(Fmunu.get(x-e_nu,mu,nu))*U.get(x-e_nu,nu));
=======
f_U.set(x,mu,f_U.get(x,mu)+2.0*U.get(x+e_mu,nu)*Adj(Fmunu.get(x,mu,nu)));
f_U.set(x,mu,f_U.get(x,mu)-2.0*Adj(Fmunu.get(x-e_nu,mu,nu))*U.get(x-e_nu,nu));
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}}}

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
f_U.set(x,mu,KAPPA*f_U.get(x,mu));
}}

<<<<<<< HEAD
// U(1) mass term
=======


// Konishi mass term
// Force contribution for the single trace mass term // 
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
UdU.set(x,mu,U.get(x,mu)*Udag.get(x,mu)-Umatrix(1));
}}

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){

f_U.set(x,mu,f_U.get(x,mu)+
<<<<<<< HEAD
(2*KAPPA*BMASS*BMASS)*U.get(x,mu)*Udag.get(x,mu)*UdU.get(x,mu));
}}

force_det(U,ft);

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
f_U.set(x,mu,f_U.get(x,mu)+ft.get(x,mu));
=======
(2*KAPPA*BMASS*BMASS)*Udag.get(x,mu)*UdU.get(x,mu));
}}


// Below is the force term corresponding to double trace mass term // 

/* sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
UdU.set(x,mu,U.get(x,mu)*Udag.get(x,mu)-Umatrix(1));
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}}

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
<<<<<<< HEAD
=======

f_U.set(x,mu,f_U.get(x,mu)+
(2*KAPPA*BMASS*BMASS/(NCOLOR*NCOLOR))*Tr(UdU.get(x,mu)).real()*
Udag.get(x,mu));}} 

*/ 



sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
f_U.set(x,mu,-1.0*Adj(f_U.get(x,mu)));
}}

// SU(N) case

if(NUMGEN==(NCOLOR*NCOLOR-1)){
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
f_U.set(x,mu,f_U.get(x,mu)-(1.0/NCOLOR)*Tr(f_U.get(x,mu))*Umatrix(1));
}
}
}


<<<<<<< HEAD
// add in contributions from pseudofermions
// use partial fraction approx to inverse square root of operator
=======
// Add in contributions from pseudofermions
// Use partial fraction approx to inverse square root of operator
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

f_F=Twist_Fermion();

if(FERMIONS && (fermion==1)){
f_U=Gauge_Field();
<<<<<<< HEAD
//cout << "fermion forces " << endl;
=======
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0


f_F=-ampdeg*Cjg(F);

MCG_solver(U,F,shift,sol,psol);

for(int n=0;n<DEGREE;n++){

stmp=sol[n];
ptmp=psol[n];

fermion_forces(U,utmp,stmp,ptmp);

f_F=f_F-amp[n]*Cjg(sol[n]);

<<<<<<< HEAD
// add in kick from fermion effective action
=======
// Add in kick from fermion effective action // 
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

sites=0;
while(loop_over_lattice(x,sites)){
for(int mu=0;mu<NUMLINK;mu++){
f_U.set(x,mu,f_U.get(x,mu)+amp[n]*utmp.get(x,mu));
}}



}

f_F=Cjg(f_F);

}


return;
}
