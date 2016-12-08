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
for(nu=mu+1;nu<(D-1);nu++){

e_mu=Lattice_Vector(mu);
e_nu=Lattice_Vector(nu);

y=x;
prod=U.get(x,mu);
for(r=1;r<R;r++){
prod=prod*U.get(y,e_mu,mu);
y=y+e_mu;
}

y=y+e_mu;
prod=prod*U.get(y,e_nu,nu);
for(t=1;t<M;t++){
prod=prod*U.get(y,e_nu,nu);
y=y+e_nu;
}

y=x;
prod2=U.get(x,nu);
for(r=1;r<R;r++){
prod2=prod2*U.get(y,e_nu,nu);
y=y+e_nu;
}

y=y+e_nu;
prod2=U.get(y,e_mu,mu);
for(t=1;t<M;t++){
prod2=prod2*U.get(y,e_mu,mu);
y=y+e_mu;
}

prod=prod*Adj(prod2);

wilson[R][M]=wilson[R][M]+(1.0/NCOLOR)*Tr(prod).norm();
}
}

//cout << "wilson loop is " << wilson[R][M] << "\n" << flush;
}

return;
}
