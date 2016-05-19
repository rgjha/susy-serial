#include "konishi.h"

double konishi(const Gauge_Field &U){
<<<<<<< HEAD
int sites,a;
Lattice_Vector x;
double dum,K;
Gauge_Field B;
static int count=0;

if(count==0){cout << "bare lattice" << endl;}
if(count==1){cout << "block lattice" << endl;}
sites=0;
=======
int sites,a,b,mu,nu;
double dum,K,KK[4][4],KK2[4][4];
Gauge_Field B;
static int count=0,first_time=1;
static ofstream f_k1,f_k2;
Lattice_Vector x,y,e_mu;
USite_Field prod;
double P[NUMLINK][NUMLINK];
 
  	if(first_time==1){
    f_k1.open("sug1",ios::app);
	if(f_k1.bad()){
	cerr << "failed to open sug1 file" << "\n";exit(1);}

    f_k2.open("sug2",ios::app);
    if(f_k2.bad()){
    cerr << "failed to open sug2 file" << "\n";exit(1);}

	first_time=0;}

P[0][0]=1/sqrt(2.0);
P[1][0]=1/sqrt(6.0);
P[2][0]=1/sqrt(12.0);
P[3][0]=1/sqrt(20.0);
P[4][0]=1/sqrt(5.0);

P[0][1]=-1/sqrt(2.0);
P[1][1]=1/sqrt(6.0);
P[2][1]=1/sqrt(12.0);
P[3][1]=1/sqrt(20.0);
P[4][1]=1/sqrt(5.0);

P[0][2]=0.0;
P[1][2]=-2/sqrt(6.0);
P[2][2]=1/sqrt(12.0);
P[3][2]=1/sqrt(20.0);
P[4][2]=1/sqrt(5.0);

P[0][3]=0.0;
P[1][3]=0.0;
P[2][3]=-3/sqrt(12.0);
P[3][3]=1/sqrt(20.0);
P[4][3]=1/sqrt(5.0);

P[0][4]=0.0;
P[1][4]=0.0;
P[2][4]=0.0;
P[3][4]=-4/sqrt(20.0);
P[4][4]=1/sqrt(5.0);


if(count==0){cout << "bare lattice" << endl;}
if(count==1){cout << "block lattice" << endl;}
/*sites=0;
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
dum=0.0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
dum=dum+Tr(U.get(x,a)*Adj(U.get(x,a))).real();}
}
<<<<<<< HEAD
dum/=(SITES*NUMLINK*NCOLOR);
=======
dum/=(SITES*NUMLINK*NCOLOR);*/
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0


sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
<<<<<<< HEAD
B.set(x,a,(U.get(x,a)*Adj(U.get(x,a)))-dum*Umatrix(1));
=======
B.set(x,a,U.get(x,a)*Adj(U.get(x,a))-Umatrix(1));
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}}


sites=0;
K=0.0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
<<<<<<< HEAD
K+=Tr(B.get(x,a)*B.get(x,a)).real();
}
=======
for(b=0;b<NUMLINK;b++){
K+=Tr(B.get(x,a)*B.get(x,b)).real();
}}
}

for(mu=0;mu<4;mu++){
for(nu=0;nu<4;nu++){
KK2[mu][nu]=0.0;}}

sites=0;
while(loop_over_lattice(x,sites)){

for(mu=0;mu<4;mu++){
for(nu=0;nu<4;nu++){
KK[mu][nu]=0.0;

for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
KK[mu][nu]=KK[mu][nu]+0.5*(P[mu][a]*P[nu][b]+P[nu][a]*P[mu][b])*Tr(B.get(x,a)*B.get(x,b)).real();
}}

}}

dum=0.0;
for(mu=0;mu<4;mu++){
dum=dum+KK[mu][mu];}

//cout << "trace is " << dum << endl;

for(mu=0;mu<4;mu++){
for(nu=0;nu<4;nu++){
if(mu==nu){
KK[mu][nu]=KK[mu][nu]-0.25*dum;}
}}

for(mu=0;mu<4;mu++){
for(nu=0;nu<4;nu++){
KK2[mu][nu]=KK2[mu][nu]+KK[mu][nu];}}

}


if(count==0){
for(mu=0;mu<4;mu++){
for(nu=0;nu<4;nu++){
f_k1 << KK2[mu][nu]/SITES << endl;}}
}
if(count==1){
for(mu=0;mu<4;mu++){
for(nu=0;nu<4;nu++){
f_k2 << KK2[mu][nu]/(16*SITES) << endl;}}
>>>>>>> f33135b5861f274b44c622ee0ce6ebc81e898eb0
}

count++;
count=count%2;

cout << "K is " << K/SITES << endl;

return(K/SITES);
}

