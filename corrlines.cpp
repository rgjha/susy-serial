#include "corrlines.h"

void corrlines(const Gauge_Field &U){

static int first_time=1;
int site=0,site2=0,t,num[1000],dist2,nzero,i;
Complex corr[1000],ncorr[1000];
double distance[1000];
static ofstream f_line2;
Lattice_Vector x,y,e_mu;
USite_Field prod;

 
  	if(first_time==1){
        f_line2.open("corrlines",ios::app);
	if(f_line2.bad()){
	cerr << "failed to open corrlines file" << "\n";exit(1);}
	first_time=0;}



for(i=0;i<1000;i++){
corr[i]=Complex();
num[i]=0;
}

prod=USite_Field(0);
site=0;
while(loop_over_lattice(x,site)){

e_mu=Lattice_Vector(D-1);


y=x;
for(t=1;t<=T;t++){
prod.set(x,prod.get(x)*U.get(y,D-1));
y=y+e_mu;
}

}
	
site=0;
while(loop_over_lattice(x,site)){
site2=0;
while(loop_over_lattice(y,site2)){
if(x.get(D-1)==y.get(D-1)){
dist2=length(x-y);
num[dist2]++;
corr[dist2]=corr[dist2]+(1.0/(NCOLOR*NCOLOR))*
Tr(prod.get(x))*Tr(Adj(prod.get(y)));
}}}

nzero=0;
for(i=0;i<1000;i++){
if(num[i]>0){
ncorr[nzero]=(1.0/num[i])*corr[i];
distance[nzero]=sqrt(1.0*i);nzero++;}
}

for(i=0;i<nzero;i++){
f_line2 << distance[i] << "\t" <<  ncorr[i] << "\n";}
f_line2 << flush;

return;
}
