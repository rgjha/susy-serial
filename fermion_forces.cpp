#include "fermion_forces.h"

void fermion_forces(const Gauge_Field &U, Gauge_Field &f_U, 
const Twist_Fermion &s, const Twist_Fermion &p){
Lattice_Vector x,e_mu,e_nu,e_a,e_b,e_c;
int sites,mu,nu,a,b,c,d,e,l,m;
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

for(a=0;a<NUMGEN;a++){
for(b=0;b<NUMGEN;b++){
tmp=tmp-
conjug(p.getC().get(x,mu,nu).get(a))*s.getL().get(x+e_mu,nu).get(b)*BC(x,e_mu)*
U.get(x,mu)*Lambda[b]*Lambda[a]+
conjug(p.getC().get(x-e_nu,mu,nu).get(a))*s.getL().get(x-e_nu,nu).get(b)*
U.get(x,mu)*Lambda[a]*Lambda[b];
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
U.get(x,mu)*Lambda[b]*Lambda[a]+
conjug(p.getL().get(x+e_mu,nu).get(a))*BC(x,e_mu)*s.getC().get(x,mu,nu).get(b)*
U.get(x,mu)*Lambda[a]*Lambda[b];
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
Lambda[a]*Lambda[b]*Udag.get(x,mu)+
conjug(p.getS().get(x+e_mu).get(a))*BC(x,e_mu)*s.getL().get(x,mu).get(b)*
Lambda[b]*Lambda[a]*Udag.get(x,mu);
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
Lambda[a]*Lambda[b]*Udag.get(x,mu)+
conjug(p.getL().get(x,mu).get(a))*s.getS().get(x).get(b)*
Lambda[b]*Lambda[a]*Udag.get(x,mu);
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
perm[d][e][c][a][b]*Udag.get(x,c);
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
perm[a][b][c][d][e]*Udag.get(x,c);
}}}}
}}

f_U.set(x,c,f_U.get(x,c)-0.5*Adj(tmp));
}

}

}

// SIMON: beta term modification

if(NUMGEN==(NCOLOR*NCOLOR)){

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){

tmp=Umatrix();
for(a=0;a<NUMGEN-1;a++){
for(b=0;b<NUMGEN;b++){
tmp=tmp-conjug(p.getS().get(x).get(a))*s.getL().get(x,mu).get(b)*
Lambda[a]*Lambda[b];
}}

f_U.set(x,mu,f_U.get(x,mu));
}}


sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){

tmp=Umatrix();
for(a=0;a<NUMGEN-1;a++){
for(b=0;b<NUMGEN;b++){
tmp=tmp-conjug(p.getL().get(x,mu).get(b))*s.getS().get(x).get(a)*
Lambda[a]*Lambda[b];
}}

f_U.set(x,mu,f_U.get(x,mu));
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
