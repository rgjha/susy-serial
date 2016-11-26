#include "kinetic_energy.h"

double kinetic_energy(const Gauge_Field &p_U,
const Twist_Fermion &p_F){
Complex dum=Complex();
int sites,mu;
Lattice_Vector x;

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
dum=dum+Tr(Adj(p_U.get(x,mu))*p_U.get(x,mu));
}
}  
 
 
dum=dum+Tr(Adj(p_F)*p_F);
   
return(dum.real());
	
} 
