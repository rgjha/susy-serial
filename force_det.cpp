#include "force_det.h"

void force_det(const Gauge_Field &U, Gauge_Field &f_U){
Lattice_Vector x,e_b,e_c;
int site,a,b,c;
Gauge_Field Udag;
UPlaq_Field adjplaq=UPlaq_Field();
UPlaq_Field detplaq=UPlaq_Field();
Umatrix tmp,plaq;

Udag=Adj(U);
f_U=Gauge_Field();

site=0;
while(loop_over_lattice(x,site)){
for(c=0;c<NUMLINK;c++){
e_c=Lattice_Vector(c);
for(b=c+1;b<NUMLINK;b++){
e_b=Lattice_Vector(b);

plaq=U.get(x,c)*U.get(x+e_c,b)*Udag.get(x+e_b,c)*Udag.get(x,b);
adjplaq.set(x,c,b,adjugate(plaq));
adjplaq.set(x,b,c,Adj(adjplaq.get(x,c,b)));

// store determinants for each plaq as diagonal elements of a BPlaq_Field ...
tmp=Umatrix();
for(a=0;a<NCOLOR;a++){
tmp.set(a,a,det(plaq));}

detplaq.set(x,c,b,tmp);
detplaq.set(x,b,c,Adj(detplaq.get(x,c,b)));

}}}

site=0;
while(loop_over_lattice(x,site)){
for(c=0;c<NUMLINK;c++){
e_c=Lattice_Vector(c);
tmp=Umatrix();

for(b=0;b<NUMLINK;b++){
if(b==c) continue;
e_b=Lattice_Vector(b);


// derivs of det P
tmp=tmp+U.get(x+e_c,b)*Udag.get(x+e_b,c)*Udag.get(x,b)*adjplaq.get(x,c,b)*
(detplaq.get(x,b,c).get(0,0)-Complex(1.0,0.0));
tmp=tmp+Udag.get(x+e_c-e_b,b)*Udag.get(x-e_b,c)*adjplaq.get(x-e_b,b,c)*U.get(x-e_b,b)*
(detplaq.get(x-e_b,c,b).get(0,0)-Complex(1.0,0.0));

}

f_U.set(x,c,G*tmp);
}}


return;
}
