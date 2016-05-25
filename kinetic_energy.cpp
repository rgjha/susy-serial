#include "kinetic_energy.h"

double kinetic_energy(const Gauge_Field &p_U,
const Twist_Fermion &p_F){
Complex dum=Complex();
<<<<<<< HEAD
int sites,mu,mid1;
=======
int sites,mu;
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88
Lattice_Vector x;

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
<<<<<<< HEAD
mid1=0; 
dum=dum+Tr(Adj(p_U.get(x,mu))*p_U.get(x,mu));
}
}  

mid1=dum.real(); 
cout << "  BOSON_MOM " << mid1 << "\n";
 
dum=dum+Cjg(p_F)*p_F;
cout << "  FERMIONS_MOM " << dum.real() - mid1 << "\n";
=======
dum=dum+Tr(Adj(p_U.get(x,mu))*p_U.get(x,mu));
}
}  
 
 
dum=dum+Cjg(p_F)*p_F;
>>>>>>> 233423c79c47c3999f05183e0ce9d46165517c88
   
return(dum.real());
	
} 
