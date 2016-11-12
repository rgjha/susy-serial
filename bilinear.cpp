#include "bilinear.h"

void bilinear(const Gauge_Field &U, const Twist_Fermion &F){
Lattice_Vector x,y;
int mu,nu,n,a,b,sources,sites,t;
Twist_Fermion sol,rhs; 
Site_Field Stmp;
Link_Field Ltmp;
Complex g[T],geta[T];
Afield Atmp;
Adjoint_Links V;
Gauge_Field Ubar;
static int first_time=1;
static ofstream f_bilin,f_sk,f_eta;

if(first_time==1){
f_bilin.open("bilinear",ios::app);
if(f_bilin.bad()){
cerr << "failed to open bilinear file" << "\n";exit(1);}
f_eta.open("etapsi",ios::app);
if(f_eta.bad()){
cerr << "failed to open etapsi file" << "\n";exit(1);}
f_sk.open("superK",ios::app);
if(f_sk.bad()){
cerr << "failed to open superK file" << "\n";exit(1);}
first_time=0;}

for(a=0;a<T;a++){
g[a]=Complex();
geta[a]=Complex();
}

Ubar=Adj(U);
compute_Adjoint_Links(U,V);
    
// example: Tr(chi F) source and sink

int SOURCES=10;
Complex cond=Complex();

sites=0;
double dum=0.0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
dum=dum+Tr(Adj(U.get(x,a))*U.get(x,a)).real();}
}
dum/=(SITES*NUMLINK*NCOLOR);
    
//random source pt
for(int sources=0;sources<SOURCES;sources++){
for(int n=0;n<(D-1);n++){
x.set(n,(int)(rand()/(double)RAND_MAX*LX));}
x.set(D-1,(int)(rand()/(double)RAND_MAX*T));

// do not source the trace piece of eta
for(a=0;a<NUMGEN-1;a++){

//set up eta source
Atmp=Afield();
Atmp.set(a,Complex(1.0,0.0));
Stmp=Site_Field();
Stmp.set(x,Atmp);
rhs=Twist_Fermion();
rhs.setS(Stmp);

//compute propagator
CG_solver(U,rhs,sol);

//loop over psi sink and fold in gauge link
for(mu=0;mu<NUMLINK;mu++){
for(b=0;b<NUMGEN;b++){
cond=cond+sol.getL().get(x,mu).get(b)*conjug(V.get(x,mu).get(a,b));}}

sites=0;
while(loop_over_lattice(y,sites)){
for(b=0;b<NUMGEN-1;b++){
t=(x.get(D-1)-y.get(D-1)+T)%T;
geta[t]=geta[t]+sol.getS().get(y).get(b)*Tr(Lambda[a]*U.get(x,mu)*Ubar.get(x,mu))*
                                         Tr(Lambda[b]*U.get(y,nu)*Ubar.get(y,nu));
}}
}

// now superpartner of Konishi
// <\psi_a^A(x)\psi_b^B(y)><Tr(T_A Ubar_a(x))Tr(T_B Ubar_b(y))>
    
for(a=0;a<NUMGEN;a++){
for(mu=0;mu<NUMLINK;mu++){
Atmp=Afield();
Atmp.set(a,Complex(1.0,0.0));
Ltmp=Link_Field();
Ltmp.set(x,mu,Atmp);
rhs=Twist_Fermion();
rhs.setL(Ltmp);
        
CG_solver(U,rhs,sol);
for(b=0;b<NUMGEN;b++){
for(nu=0;nu<NUMLINK;nu++){
    
sites=0;
while(loop_over_lattice(y,sites)){
t=(y.get(D-1)-x.get(D-1)+T)%T;
g[t]=g[t]+sol.getL().get(y,nu).get(b)*
     Tr(Lambda[a]*Ubar.get(x,mu)*(Ubar.get(x,mu)*U.get(x,mu)-dum*Umatrix(1)))*
     Tr(Lambda[b]*Ubar.get(y,nu)*(Ubar.get(y,nu)*U.get(y,nu)-dum*Umatrix(1)));
}
    
}
}
}}
    
}

f_bilin << (1.0/SOURCES)*cond << endl;
for(t=0;t<T;t++){
f_sk << t << "\t" << g[t] << endl;
f_eta << t << "\t" << geta[t] << endl;
}
return;
}
