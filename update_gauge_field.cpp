#include "update_gauge_field.h"

void update_gauge_field(Gauge_Field &U, const Gauge_Field &p_U, double dt){
Lattice_Vector x;
int site=0;

while(loop_over_lattice(x,site)){

for(int mu=0;mu<NUMLINK;mu++){
<<<<<<< HEAD
U.set(x,mu,exp(dt*p_U.get(x,mu))*U.get(x,mu));}
=======
U.set(x,mu,U.get(x,mu)+dt*p_U.get(x,mu));}
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

}

return;
}
