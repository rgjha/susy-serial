#include "fermion_forces.h"

void fermion_forces(const Gauge_Field &U, Gauge_Field &f_U, 
const Twist_Fermion &s, const Twist_Fermion &p){
Lattice_Vector x,e_mu,e_nu,e_a,e_b,e_c;
int sites,mu,nu,a,b,c,d,e,l,m;
Umatrix tmp;
Complex tmp2,tt,tmp3;
UPlaq_Field P;
Gauge_Field Udag;

Udag=Adj(U);
f_U=Gauge_Field();

// compute fermion kick to gauge link force
// f_U=Adj(Ms).D_U M(U,Ub) s - Adj[Adj(Ms). D_Ub M(U,Ub) s]
 
// DU on chi_{munu}D_mu(U)psi_nu 

sites=0;
while(loop_over_lattice(x,sites)){

for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

tmp=Umatrix();

for(nu=0;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
if (mu==nu) continue;

for(a=0;a<NUMGEN;a++){
for(b=0;b<NUMGEN;b++){
tmp=tmp-
conjug(p.getC().get(x,mu,nu).get(a))*s.getL().get(x+e_mu,nu).get(b)*BC(x,e_mu)*
Lambda[b]*Lambda[a]+
conjug(p.getC().get(x-e_nu,mu,nu).get(a))*s.getL().get(x-e_nu,nu).get(b)*
Lambda[a]*Lambda[b];
}}}

f_U.set(x,mu,f_U.get(x,mu)+1.0*tmp);

}

}

sites=0;
while(loop_over_lattice(x,sites)){

for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

tmp=Umatrix();

for(nu=0;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
if(mu==nu) continue;

for(a=0;a<NUMGEN;a++){
for(b=0;b<NUMGEN;b++){
tmp=tmp-
conjug(p.getL().get(x-e_nu,nu).get(a))*s.getC().get(x-e_nu,mu,nu).get(b)*
Lambda[b]*Lambda[a]+
conjug(p.getL().get(x+e_mu,nu).get(a))*BC(x,e_mu)*s.getC().get(x,mu,nu).get(b)*
Lambda[a]*Lambda[b];
}}
}

f_U.set(x,mu,f_U.get(x,mu)+1.0*tmp);
}

}

// DUb on psi_muDb_mu(U)eta
  

sites=0;
while(loop_over_lattice(x,sites)){


for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

tmp=Umatrix();
for(a=0;a<NUMGEN;a++){
for(b=0;b<NUMGEN;b++){
tmp=tmp-
conjug(p.getS().get(x).get(a))*s.getL().get(x,mu).get(b)*
Lambda[a]*Lambda[b]+
conjug(p.getS().get(x+e_mu).get(a))*BC(x,e_mu)*s.getL().get(x,mu).get(b)*
Lambda[b]*Lambda[a];
}}


f_U.set(x,mu,f_U.get(x,mu)+0.5*Adj(tmp));
}

}

sites=0;
while(loop_over_lattice(x,sites)){

for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

tmp=Umatrix();
for(a=0;a<NUMGEN;a++){
for(b=0;b<NUMGEN;b++){
tmp=tmp-
conjug(p.getL().get(x,mu).get(a))*s.getS().get(x+e_mu).get(b)*BC(x,e_mu)*
Lambda[a]*Lambda[b]+
conjug(p.getL().get(x,mu).get(a))*s.getS().get(x).get(b)*
Lambda[b]*Lambda[a];
}}


f_U.set(x,mu,f_U.get(x,mu)+0.5*Adj(tmp));
}

}

// chi_ab D_c chi_de epsilon_abcde -- Q-closed piece

if(NUMLINK==5){

sites=0;
while(loop_over_lattice(x,sites)){

for(c=0;c<NUMLINK;c++){
e_c=Lattice_Vector(c);
tmp=Umatrix();

for(d=0;d<NUMLINK;d++){
if(d==c) continue;
for(e=d+1;e<NUMLINK;e++){
if(e==c) continue;

for(a=0;a<NUMLINK;a++){
if((a==d)||(a==e)||(a==c)) continue;
e_a=Lattice_Vector(a);
for(b=a+1;b<NUMLINK;b++){
if((b==d)||(b==e)||(b==c)) continue;
e_b=Lattice_Vector(b);

for(l=0;l<NUMGEN;l++){
for(m=0;m<NUMGEN;m++){
tmp=tmp+(
conjug(p.getC().get(x+e_a+e_b+e_c,d,e).get(l))*s.getC().get(x+e_c,a,b).get(m)*
BC(x,e_a,e_b,e_c)*BC(x,e_c)*Lambda[l]*Lambda[m]-
conjug(p.getC().get(x+e_c,d,e).get(l))*s.getC().get(x-e_a-e_b,a,b).get(m)*
BC(x,-e_a,-e_b)*BC(x,e_c)*Lambda[m]*Lambda[l])*
perm[d][e][c][a][b];
}}}}
}}

f_U.set(x,c,f_U.get(x,c)-0.5*Adj(tmp));
}

}

sites=0;
while(loop_over_lattice(x,sites)){


for(c=0;c<NUMLINK;c++){
e_c=Lattice_Vector(c);
tmp=Umatrix();

for(d=0;d<NUMLINK;d++){
if(d==c) continue;
for(e=d+1;e<NUMLINK;e++){
if(e==c) continue;

for(a=0;a<NUMLINK;a++){
if((a==d)||(a==e)||(a==c)) continue;
e_a=Lattice_Vector(a);
for(b=a+1;b<NUMLINK;b++){
if((b==d)||(b==e)||(b==c)) continue;
e_b=Lattice_Vector(b);

for(l=0;l<NUMGEN;l++){
for(m=0;m<NUMGEN;m++){
tmp=tmp+(
conjug(p.getC().get(x-e_a-e_b,a,b).get(l))*s.getC().get(x+e_c,d,e).get(m)*
BC(x,-e_a,-e_b)*BC(x,e_c)*Lambda[l]*Lambda[m]-
conjug(p.getC().get(x+e_c,a,b).get(l))*s.getC().get(x+e_a+e_b+e_c,d,e).get(m)*
BC(x,e_a,e_b,e_c)*BC(x,e_c)*Lambda[m]*Lambda[l])*
perm[a][b][c][d][e];
}}}}
}}

f_U.set(x,c,f_U.get(x,c)-0.5*Adj(tmp));
}

}

}

// SIMON: susy det term modification


if(GO){

P=Plaq(U);
Scalar_Plaquette W,Z,F,Y;
Scalar_Plaquette G1[NUMGEN];

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
tmp2=det(P.get(x,mu,nu))-Complex(1.0,0.0);
W.set(x,mu,nu,tmp2);
Z.set(x,mu,nu,(Complex(1.0,0.0)+tmp2));
}}}
// first part of fermion force

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
e_nu=Lattice_Vector(nu);

for(b=0;b<NUMGEN;b++){
tmp2=conjug(p.getS().get(x).get(NUMGEN-1))*s.getL().get(x,mu).get(b)*
Z.get(x,mu,nu)+
conjug(p.getS().get(x-e_nu).get(NUMGEN-1))*s.getL().get(x,mu).get(b)*
Z.get(x-e_nu,nu,mu)*BC(x,-e_nu)-
conjug(p.getL().get(x,mu).get(b))*s.getS().get(x).get(NUMGEN-1)*
Z.get(x,mu,nu)-
conjug(p.getL().get(x,mu).get(b))*s.getS().get(x-e_nu).get(NUMGEN-1)*
Z.get(x-e_nu,nu,mu)*BC(x,-e_nu);

G1[b].set(x,mu,nu,G*0.5*Complex(0.0,1.0)*sqrt((double)NCOLOR)*tmp2);
}
}}}

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
tmp=Umatrix();

for(b=0;b<NUMGEN;b++){
for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
tmp=tmp-G1[b].get(x,mu,nu)*inverse(U.get(x,mu))*Lambda[b]*inverse(U.get(x,mu));
}}


f_U.set(x,mu,f_U.get(x,mu)-tmp);
}}

// second part

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
tmp2=Complex();

for(b=0;b<NUMGEN;b++){
tmp2=tmp2+
conjug(p.getS().get(x).get(NUMGEN-1))*Tr(inverse(U.get(x,mu))*
Lambda[b])*s.getL().get(x,mu).get(b)+
conjug(p.getS().get(x).get(NUMGEN-1))*Tr(inverse(U.get(x+e_mu,nu))*
Lambda[b])*s.getL().get(x+e_mu,nu).get(b)*BC(x,e_mu)-
conjug(p.getL().get(x,mu).get(b))*Tr(inverse(U.get(x,mu))*
Lambda[b])*s.getS().get(x).get(NUMGEN-1)-
conjug(p.getL().get(x+e_mu,nu).get(b))*Tr(inverse(U.get(x+e_mu,nu))*
Lambda[b])*s.getS().get(x).get(NUMGEN-1)*BC(x,e_mu);
}

F.set(x,mu,nu,0.5*G*Complex(0.0,1.0)*sqrt((double)NCOLOR)*tmp2);

}}}


Z=mydiff(F,W);
Y=mydiff2(F,W);

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){

tmp2=Complex();
tmp3=Complex();
for(nu=0;nu<NUMLINK;nu++){
if(mu==nu) continue;
tmp2=tmp2+Z.get(x,mu,nu);
tmp3=tmp3+Y.get(x,mu,nu);}

f_U.set(x,mu,f_U.get(x,mu)-tmp2*inverse(U.get(x,mu)));
f_U.set(x,mu,f_U.get(x,mu)-Adj(tmp3*inverse(Udag.get(x,mu))));
}}


}

// SU(N) case

if(NUMGEN==(NCOLOR*NCOLOR-1)){
sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
f_U.set(x,mu,f_U.get(x,mu)-(1.0/NCOLOR)*Tr(f_U.get(x,mu))*Umatrix(1));
}
}
}


sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
f_U.set(x,mu,-1.0*Adj(f_U.get(x,mu)));
}
}

return;
}
