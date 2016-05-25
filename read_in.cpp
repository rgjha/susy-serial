#include "read_in.h"

void read_in(Gauge_Field &u, Twist_Fermion &F){
ifstream f_read;
int LDUM, TDUM;
double BETADUM, NCOLORDUM, DTDUM;
Afield dummy2;
Umatrix dummy;
Lattice_Vector x;
int site,mu,nu;
Site_Field s;
Link_Field l;
Plaq_Field p;

f_read.open("config");
if(f_read.bad()){
cout << "error opening config file to read\n" << flush;}

f_read >> LDUM >> TDUM >> BETADUM >> DTDUM >> NCOLORDUM; 
if ((LDUM!=LX) || (TDUM!=T)) {cout << "wrong size lattice read in - abort\n";}
cout << "config coupling is " << BETADUM << "\n";
cout << "time step used for input config " << DTDUM << "\n" << flush;
cout << "number of colors " << NCOLORDUM << "\n" << flush;

site=0;
while(loop_over_lattice(x,site)){
for(int mu=0;mu<NUMLINK;mu++){
f_read >> dummy;
u.set(x,mu,dummy);}
}


site=0;
while(loop_over_lattice(x,site)){
f_read >> dummy2;
s.set(x,dummy2);
}
F.setS(s);


site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
f_read >> dummy2;
l.set(x,mu,dummy2);
}
}
F.setL(l);

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
p.set(x,mu,mu,Afield());
for(nu=mu+1;nu<NUMLINK;nu++){
f_read >> dummy2;
p.set(x,mu,nu,dummy2);
p.set(x,nu,mu,-1.0*dummy2);
}
}
}
F.setC(p);

f_read.close();
return;
}

