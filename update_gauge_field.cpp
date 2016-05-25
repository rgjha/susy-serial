#include "update_gauge_field.h"

void update_gauge_field(Gauge_Field &U, const Gauge_Field &p_U, double dt){
Lattice_Vector x;
int site=0;

while(loop_over_lattice(x,site)){

for(int mu=0;mu<NUMLINK;mu++){
U.set(x,mu,exp(dt*p_U.get(x,mu))*U.get(x,mu));}

}

return;
}
