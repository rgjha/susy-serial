#include "divdet.h"

void divdet(const Gauge_Field &U, Gauge_Field &U2){
Lattice_Vector x;
int mu,site;
Complex tmp;
double r,theta,scale;

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
tmp=det(U.get(x,mu));
r=tmp.norm();
theta=atan(tmp.imag()/tmp.real());
scale=pow(r,(-1.0/NCOLOR));
theta=-1.0*theta/NCOLOR;
tmp=scale*Complex(cos(theta),sin(theta));

U2.set(x,mu,tmp*U.get(x,mu));
}}


return;
}
