#include "scalars.h"

void scalars(const Gauge_Field &U, double &t1, double &t2){
int sites,mu;
Lattice_Vector x;

sites=0;
	t1=0.0;t2=0.0;
	while(loop_over_lattice(x,sites)){
    
	for(int mu=2;mu<(D-1);mu++){
	t1+=Tr(U.get(x,mu)*Adj(U.get(x,mu))).real();
	}
    
        t2+=Tr(U.get(x,D-1)*Adj(U.get(x,D-1))).real();
    
	}

       t1=t1/(SITES*NCOLOR);
       t2=t2/(SITES*NCOLOR);

return;
}

