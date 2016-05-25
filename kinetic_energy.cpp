#include "kinetic_energy.h"

double kinetic_energy(const Gauge_Field &p_U,
const Twist_Fermion &p_F){
Complex dum=Complex();
int sites,mu,mid1;
Lattice_Vector x;

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
mid1=0; 
dum=dum+Tr(Adj(p_U.get(x,mu))*p_U.get(x,mu));
}
}  

mid1=dum.real(); 
cout << "  BOSON_MOM " << mid1 << "\n";
 
dum=dum+Cjg(p_F)*p_F;
cout << "  FERMIONS_MOM " << dum.real() - mid1 << "\n";
   
return(dum.real());
	
} 
