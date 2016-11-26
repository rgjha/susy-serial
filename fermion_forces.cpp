#include "fermion_forces.h"

void fermion_forces(const Gauge_Field &U, Gauge_Field &f_U, 
const Twist_Fermion &s, const Twist_Fermion &p){
Lattice_Vector x,e_mu,e_nu,e_a,e_b,e_c;
int sites,mu,nu,a,b,c,d,e;
Umatrix tmp;
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

tmp=tmp+
s.getL().get(x,e_mu,nu)*Adj(p.getC().get(x,mu,nu))-
Adj(p.getC().get(x,-e_nu,mu,nu))*s.getL().get(x,-e_nu,nu);
}

f_U.set(x,mu,f_U.get(x,mu)+U.get(x,mu)*tmp);
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

tmp=tmp-
Adj(p.getL().get(x,e_mu,nu))*s.getC().get(x,mu,nu)+
s.getC().get(x,-e_nu,mu,nu)*Adj(p.getL().get(x,-e_nu,nu));
}

f_U.set(x,mu,f_U.get(x,mu)+U.get(x,mu)*tmp);
}
}

// DUb on psi_muDb_mu(U)eta

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

tmp=Umatrix();

tmp=tmp+
Adj(p.getS().get(x))*s.getL().get(x,mu)-
s.getL().get(x,mu)*Adj(p.getS().get(x,e_mu));

f_U.set(x,mu,f_U.get(x,mu)+0.5*Adj(tmp*Udag.get(x,mu)));
}}


sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

tmp=Umatrix();

tmp=tmp-
s.getS().get(x)*Adj(p.getL().get(x,mu))+
Adj(p.getL().get(x,mu))*s.getS().get(x,e_mu);

f_U.set(x,mu,f_U.get(x,mu)+0.5*Adj(tmp*Udag.get(x,mu)));

}}


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

tmp=tmp+(
Adj(p.getC().get(x,e_a,e_b,e_c,d,e))*s.getC().get(x,e_c,a,b)-
s.getC().get(x,-e_a,-e_b,a,b)*Adj(p.getC().get(x,e_c,d,e))
)*perm[d][e][c][a][b];

}}}}


f_U.set(x,c,f_U.get(x,c)+0.5*Adj(tmp*Udag.get(x,c)));
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

tmp=tmp-(
s.getC().get(x,e_a,e_b,e_c,d,e)*Adj(p.getC().get(x,e_c,a,b))-
Adj(p.getC().get(x,-e_a,-e_b,a,b))*s.getC().get(x,e_c,d,e)
)*perm[a][b][c][d][e];

}}}}

f_U.set(x,c,f_U.get(x,c)+0.5*Adj(tmp*Udag.get(x,c)));
}
}

// end NUMLINK==5 block
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
