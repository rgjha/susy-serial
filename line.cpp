#include "line.h"

// computes Polyakov line
Complex line(const Gauge_Field  &U, const int mu){

Lattice_Vector x,y,e_mu;
Umatrix prod;
int site,t,M;
Complex poly=Complex();


if(mu==0) M=LX;
if(mu==(D-1)) M=T;

site=0;
while(loop_over_lattice(x,site)){
prod=Umatrix(1);

e_mu=Lattice_Vector(mu);

y=x;
prod=U.get(x,mu);
for(t=1;t<M;t++){
prod=prod*(U.get(y,e_mu,mu));
y=y+e_mu;
}

poly=poly+(1.0/NCOLOR)*Tr(prod);
}

return((1.0/SITES)*poly);
}
