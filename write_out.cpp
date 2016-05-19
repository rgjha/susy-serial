#include "write_out.h"

void write_out(const Gauge_Field &u, const
Twist_Fermion &F, const int num){
static ofstream f_write;
static int first_time=1;
Lattice_Vector x;
int site,mu,nu;
char s[100];

sprintf(s,"config%d",num);

//cout << "config name " << s << endl;

if(num>0){
f_write.open(s);}
else{
f_write.open("config");}

if(f_write.bad()){
cout << "error opening config file\n" << flush;}
f_write << LX << "\t" << T << "\t" << KAPPA << "\t" << DT << "\t" << NCOLOR << "\n";

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
f_write << u.get(x,mu) << "\n";
}}
f_write << "\n" << flush;
f_write.close();

return;


if(num==0){
site=0;
while(loop_over_lattice(x,site)){
f_write << F.getS().get(x)<< "\n";
}


site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
f_write <<  F.getL().get(x,mu) << "\n";
}
}

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=mu+1;nu<NUMLINK;nu++){
f_write << F.getC().get(x,mu,nu) << "\n";
}
}
}
}

f_write.close();
return;
}

