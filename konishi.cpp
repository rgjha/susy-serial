#include "konishi.h"

double konishi(const Gauge_Field &U){
int sites,a;
Lattice_Vector x;
double dum,K;
Gauge_Field B;
static int count=0;

if(count==0){cout << "bare lattice" << endl;}
if(count==1){cout << "block lattice" << endl;}
sites=0;
dum=0.0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
dum=dum+Tr(U.get(x,a)*Adj(U.get(x,a))).real();}
}
dum/=(SITES*NUMLINK*NCOLOR);


sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
B.set(x,a,(U.get(x,a)*Adj(U.get(x,a)))-dum*Umatrix(1));
}}


sites=0;
K=0.0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
K+=Tr(B.get(x,a)*B.get(x,a)).real();
}
}

count++;
count=count%2;

cout << "K is " << K/SITES << endl;

return(K/SITES);
}

