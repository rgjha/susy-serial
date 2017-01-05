#include "gauge_transform.h"

// Simon's email around ~ Dec 20. This forgot one of the projections back to traceless matrices for the gauge fields  //

Gauge_Field gauge_transform(const Gauge_Field &U){
Lattice_Vector x,e_mu;
int site,mu;
Gauge_Field G,V;
Umatrix tmp;

site=0;
while(loop_over_lattice(x,site)){
tmp=Umatrix();
for(mu=0;mu<(NCOLOR*NCOLOR-1);mu++){
tmp=tmp+Lambda[mu]*gasdev();}
G.set(x,0,exp(tmp));
}

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
tmp=G.get(x,0)*U.get(x,mu)*Adj(G.get(x+e_mu,0));
V.set(x,mu,tmp);
}}

return(V);
}
