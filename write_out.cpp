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
<<<<<<< HEAD
f_write << LX << "\t" << T << "\t" << KAPPA << "\t" << DT << "\t" << NCOLOR << "\n";
=======
f_write << setprecision(10) << LX << "\t" << T << "\t" << KAPPA << "\t" << DT << "\t" << NCOLOR << "\n";
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
<<<<<<< HEAD
f_write << u.get(x,mu) << "\n";
}}
f_write << "\n" << flush;
f_write.close();

return;
=======
f_write << setprecision(PREC) << u.get(x,mu) << "\n";
}}
f_write << "\n" << flush;
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0


if(num==0){
site=0;
while(loop_over_lattice(x,site)){
<<<<<<< HEAD
f_write << F.getS().get(x)<< "\n";
=======
f_write << setprecision(PREC) << F.getS().get(x)<< "\n";  // THis is for the site //  Line 1339-1349 of utilities // 
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}


site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
<<<<<<< HEAD
f_write <<  F.getL().get(x,mu) << "\n";
=======
f_write <<  setprecision(PREC) << F.getL().get(x,mu) << "\n";   // This is for link fields // 
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}
}

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=mu+1;nu<NUMLINK;nu++){
<<<<<<< HEAD
f_write << F.getC().get(x,mu,nu) << "\n";
=======
f_write << setprecision(PREC) << F.getC().get(x,mu,nu) << "\n";   // This is for the plaquette or contour C // 
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}
}
}
}

f_write.close();
return;
}

