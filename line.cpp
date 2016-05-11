#include "line.h"

// computes Polyakov line
Complex line(const Gauge_Field  &U, const int mu){

Lattice_Vector x,y,e_mu;
Umatrix prod;
int site,t,M;
Complex poly=Complex();


if(mu==(D-2)) M=LZ;       // 0 is for x, 1 for y, 2 for z, 3 for T // 
if(mu==(D-1)) M=T;

site=0;
while(loop_over_lattice(x,site)){
prod=Umatrix(1);

e_mu=Lattice_Vector(mu);

y=x;
for(t=1;t<=M;t++){
prod=prod*(U.get(y,mu));
y=y+e_mu;
}

poly=poly+(1.0/NCOLOR)*Tr(prod);
}

return((1.0/SITES)*poly);
}
