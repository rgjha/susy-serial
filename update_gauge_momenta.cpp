#include "update_gauge_momenta.h"

void update_gauge_momenta(Gauge_Field &p_U, const Gauge_Field f_U, double dt){
Lattice_Vector x;
int site=0;

while(loop_over_lattice(x,site)){

for(int mu=0;mu<NUMLINK;mu++){
p_U.set(x,mu,p_U.get(x,mu)+dt*f_U.get(x,mu));}

}
return;
}
