#include "unit.h"

void unit(const Gauge_Field &U, Gauge_Field &U2){
Lattice_Vector x,e_mu;
int mu,site,i,j;
double
b[2*NCOLOR],dummy[2],dummy2[2*NCOLOR*NCOLOR],work[4*NCOLOR],work2[4*NCOLOR];
double at[2*NCOLOR*NCOLOR];
double re,im;
Umatrix diag, sim, prod, unitary,av,tmp;


if(SMALLEIG){

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){

av=U.get(x,mu)*Adj(U.get(x,mu));

int ok,c1,c2,c3,c6;
char c4,c5;


for(i=0;i<NCOLOR;i++){
for(j=0;j<NCOLOR;j++){
at[2*(j+NCOLOR*i)]=av.get(j,i).real();
at[2*(j+NCOLOR*i)+1]=av.get(j,i).imag();
}
}

c1=NCOLOR;
c2=2*NCOLOR;
c3=1;
c4='N';
c5='V';
c6=NCOLOR;

zgeev_(&c4,&c5,&c1,at,&c1,b,dummy,&c3,dummy2,&c6,work,&c2,work2,&ok);

diag=Umatrix();
for(i=0;i<2*NCOLOR;i=i+2){
//cout << "b is " << b[i] << "\t";
diag.set(i/2,i/2,Complex(1.0/sqrt(b[i]),0.0));}

//cout << "\n";

for(i=0;i<NCOLOR;i++){
for(j=0;j<NCOLOR;j++){
re=dummy2[2*(j+NCOLOR*i)];
im=dummy2[2*(j+NCOLOR*i)+1];
sim.set(j,i,Complex(re,im));
}}


prod=sim*diag*Adj(sim);
unitary=prod*U.get(x,mu);
/*prod=unitary*Adj(unitary);

for(i=0;i<NCOLOR;i++){
for(j=0;j<NCOLOR;j++){
cout << prod.get(i,j) << "\t";}
cout << "\n";
}
cout << "\n";*/

//cout << "final matrix unitary ?" << "\n";
prod=unitary*Adj(unitary);
if(fabs(1.0/NCOLOR*Tr(prod).real()-1.0)>0.0000001){
cout << "error in getting unitary piece" << "\n";}

U2.set(x,mu,unitary);
}
}
}

return;
}
