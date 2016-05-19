#include "loop.h"
// computes wilson loops

void loop(const Gauge_Field  &U, double wilson[LX][T]){
int site=0,r,t,R,M,mu,nu;
Lattice_Vector x,y,e_mu,e_nu;
Umatrix prod,prod2;

site=0;

for(R=1;R<=(LX/2);R++)
for(M=1;M<=(T/2);M++){

wilson[R][M]=0.0;
//cout << "working on R,M " << R << "\t" << M << "\n" << flush;
site=0;
while(loop_over_lattice(x,site)){

for(mu=0;mu<(D-1);mu++)
for(nu=D-1;nu<D;nu++){

e_mu=Lattice_Vector(mu);
e_nu=Lattice_Vector(nu);

prod=Umatrix(1);
y=x;

for(r=1;r<=R;r++){
prod=prod*(U.get(y,mu));
y=y+e_mu;
}
for(t=1;t<=M;t++){
prod=prod*(U.get(y,nu));
y=y+e_nu;
}
y=y-e_mu;
for(r=1;r<=R;r++){
prod=prod*Adj(U.get(y,mu));
y=y-e_mu;
}
y=y+e_mu-e_nu;
for(t=1;t<=M;t++){
prod=prod*Adj(U.get(y,nu));
y=y-e_nu;
}

<<<<<<< HEAD
wilson[R][M]=wilson[R][M]+(1.0/NCOLOR)*Tr(prod).real();
=======
wilson[R][M]=wilson[R][M]+(1.0/NCOLOR)*Tr(prod).norm();
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}
}

//cout << "wilson loop is " << wilson[R][M] << "\n" << flush;
}

return;
}
