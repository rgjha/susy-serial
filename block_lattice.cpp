#include "block_lattice.h"

void block_lattice(const Gauge_Field &U, Gauge_Field &Up){
int sites=0,mu;
Lattice_Vector x,e_mu;

Up=Gauge_Field();

while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
Up.set(x,mu,(U.get(x,mu)*U.get(x+e_mu,mu)));
}}

return;
}
