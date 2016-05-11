#include "correlators.h"

void correlators(const Gauge_Field &U){
Lattice_Vector x,y,z;
int a,b,sites,pt,j,i,t; 
double dum;
Gauge_Field B;
double vx[D],vy[D],vyb[D];
double e[D][D];

Complex X[SITES][NUMLINK][NUMLINK], K[SITES];
Complex XS[T][NUMLINK][NUMLINK],KS[T];
Complex GK[T], GS[T],GX[T][NUMLINK][NUMLINK];
static int first_time=1;
double d[1000],r,Kcorr[1000],Scorr[1000];
int total=0,id,num[1000]; 
static ofstream f_pic,f_pic2,f_SUGRA,f_K,f_K2,f_S2,f_num;

double KF[2*SITES+1],TF[2*SITES+1],XF[2*SITES+1][NUMLINK][NUMLINK];
int isign,ndim=D,nn[D+1];

nn[1]=T;
nn[2]=LX;
nn[3]=LX;
nn[4]=LX;

if(first_time==1){

f_SUGRA.open("corrSUGRA",ios::app);
if(f_SUGRA.bad()){
cerr << "failed to open corrSUGRA file" << "\n";exit(1);}

f_K.open("corrK",ios::app);
if(f_K.bad()){
cerr << "failed to open corrK file" << endl;}

f_K2.open("cK",ios::app);
if(f_K2.bad()){
cerr << "failed to open cK file" << endl;}

f_S2.open("cS",ios::app);
if(f_S2.bad()){
cerr << "failed to open cS file" << endl;}

f_num.open("num");
if(f_num.bad()){
cerr << "failed to open num file" << endl;}

f_pic.open("lattice");
if(f_pic.bad()){
cerr << "failed to open lattice file"<< endl;}

f_pic2.open("latticeb");
if(f_pic2.bad()){
cerr << "failed to open latticeb file"<<endl;}

first_time=0;
}

for(i=0;i<1000;i++){
num[i]=0;
d[i]=0.0;
Scorr[i]=0.0;
Kcorr[i]=0.0;
}

// assign basis vectors
e[0][0]=1.0/sqrt(2.0);
e[1][0]=1.0/sqrt(6.0);
e[2][0]=1.0/sqrt(12.0);
e[3][0]=1.0/sqrt(20.0);

e[0][1]=-1.0/sqrt(2.0);
e[1][1]=1.0/sqrt(6.0);
e[2][1]=1.0/sqrt(12.0);
e[3][1]=1.0/sqrt(20.0);

e[0][2]=0.0;
e[1][2]=-2.0/sqrt(6.0);
e[2][2]=1.0/sqrt(12.0);
e[3][2]=1.0/sqrt(20.0);

e[0][3]=0.0;
e[1][3]=0.0;
e[2][3]=-3.0/sqrt(12.0);
e[3][3]=1.0/sqrt(20.0);


//computing SUGRA and Konishi operators///
//////////////////////////////////////////


sites=0;
dum=0.0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
dum=dum+Tr(Adj(U.get(x,a))*U.get(x,a)).real();}
}
dum/=(SITES*NUMLINK*NCOLOR);


sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
B.set(x,a,(U.get(x,a)*Adj(U.get(x,a)))-dum*Umatrix(1));
}}


sites=0;
while(loop_over_lattice(x,sites)){ 
pt=lat_pack(x);

for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
X[pt][a][b]=Tr(B.get(x,a)*B.get(x,b));
}}
}


for(pt=0;pt<SITES;pt++){
K[pt]=Complex();
for(a=0;a<NUMLINK;a++){
K[pt]=K[pt]+X[pt][a][a];
}}

for(pt=0;pt<SITES;pt++){
for(a=0;a<NUMLINK;a++){
X[pt][a][a]=X[pt][a][a]-(1.0/NUMLINK)*K[pt];}
}

// first timesliced averaged correlators 

for(t=0;t<T;t++){
for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
XS[t][a][b]=Complex();
}}
KS[t]=Complex();
}

for(sites=0;sites<SITES;sites++){
lat_unpack(sites,x);
t=x.get(D-1);

for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
XS[t][a][b]=XS[t][a][b]+(1.0/(LX*LX*LX))*X[sites][a][b];}}

KS[t]=KS[t]+(1.0/(LX*LX*LX))*K[sites];
}

// compute timeslice correlator
for(t=0;t<T;t++){
GK[t]=Complex();
for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
GX[t][a][b]=Complex();
}}
}

for(sites=0;sites<T;sites++){
for(pt=0;pt<T;pt++){
t=((sites-pt)+T)%T;

for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
GX[t][a][b]=GX[t][a][b]+XS[sites][a][b]*XS[pt][a][b];
}}

GK[t]=GK[t]+KS[sites]*KS[pt];
}
}

// average the SUGRA correlators
for(t=0;t<T;t++){
GS[t]=Complex();
for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
GS[t]=GS[t]+GX[t][a][b];}}
GS[t]=(1.0/double(NUMLINK*NUMLINK*T))*GS[t];
GK[t]=GK[t]*(1.0/double(T));
}

for(t=0;t<T;t++){
f_K << t << "\t" << GK[t] << endl;

f_SUGRA << t << "\t" << GS[t] << endl;
}

// Use FFT to handle convolution implicit in fully averaged
// bosonic correlator

// map to offset arrays needed for FFT
t=1;
for(pt=0;pt<SITES;pt++){
KF[t]=K[pt].real();
KF[t+1]=0.0;
for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
XF[t][a][b]=X[pt][a][b].real();
XF[t+1][a][b]=0.0;}}
t=t+2;
}

// Konishi first
isign=1;
fourn(KF,nn,ndim,isign);

for(t=1;t<=(2*SITES);t+=2){
KF[t]=KF[t]*KF[t]+KF[t+1]*KF[t+1];
KF[t+1]=0.0;}

isign=-1;
fourn(KF,nn,ndim,isign);

pt=0;
for(t=1;t<=(2*SITES);t+=2){
K[pt]=(1.0/(double)SITES)*Complex(KF[t],0.0);pt++;
}

// now SUGRA FFTs

for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){

// temp array for fixed (a,b)
for(pt=1;pt<=(2*SITES);pt++){
TF[pt]=XF[pt][a][b];}

isign=1;
fourn(TF,nn,ndim,isign);

for(t=1;t<=(2*SITES);t+=2){
TF[t]=TF[t]*TF[t]+TF[t+1]*TF[t+1];
TF[t+1]=0.0;}

isign=-1;
fourn(TF,nn,ndim,isign);

pt=0;
for(t=1;t<=(2*SITES);t+=2){
X[pt][a][b]=(1.0/(double)SITES)*Complex(TF[t],0.0);pt++;
}

}}

// now look at rotationally averaged correlators
static int first=1;

x.set(0,0);
x.set(1,0);
x.set(2,0);
x.set(3,0);
pt=lat_pack(x);
for(t=0;t<SITES;t++){
lat_unpack(t,y);

r=distR(x,y);
id=check(d,r,total);
Kcorr[id]=Kcorr[id]+K[t].real();
num[id]=num[id]+1;

for(a=0;a<NUMLINK;a++){
for(b=0;b<NUMLINK;b++){
Scorr[id]=Scorr[id]+X[t][a][b].real();
}}
Scorr[id]=Scorr[id]*(1.0/(NUMLINK*NUMLINK));

}

if(first){
for(a=0;a<total;a++){
f_num << d[a] << "\t" << num[a]<< endl;}
first=0;
}

for(a=0;a<total;a++){
f_K2 << d[a] << "\t" << Kcorr[a]/(double)SITES << endl;
f_S2 << d[a] << "\t" << Scorr[a]/(double)SITES << endl;
}


return;
}


//function for |r|
double distR(Lattice_Vector &x, Lattice_Vector &y){  
double r=0.0;
int i,j;
Lattice_Vector zdum;
double delta[D][D];

int dum;
for(i=0;i<(D-1);i++){
dum=x.get(i)-y.get(i);
if(dum>(LX/2)) dum=dum-LX;
if(dum<(-LX/2)) dum=dum+LX;
zdum.set(i,dum);}


dum=x.get(D-1)-y.get(D-1);
if(dum>(T/2)) dum=dum-T;
if(dum<(-T/2)) dum=dum+T;
zdum.set(D-1,dum);

for(i=0;i<D;i++){
for(j=0;j<D;j++){
delta[i][j]=0.0;
if(i==j) delta[i][i]=1.0;
}}

zdum=x-y;
for(i=0;i<D;i++){
for(j=0;j<D;j++){
r+=(delta[i][j]-(1.0/(double)NUMLINK))*zdum.get(i)*zdum.get(j);
}}

r=sqrt(r);
return(r);
}

/*
double distR(Lattice_Vector &x, Lattice_Vector &y){  

int mu,nu,alpha0,alpha1,alpha2,alpha3,alpha[D],j;
double r,dsq,dsqmin;
Lattice_Vector z;
double kron[D][D];
    //cout << "in dist" << endl;
    for(mu=0;mu<D;mu++){
        for(nu=0;nu<D;nu++){
            kron[mu][nu]=0.0;
            if(mu==nu) {kron[mu][nu]=1.0;}
        }
    }
    
// compute distance between x and y while varying y over all neighboring copies of
// fundamental domain of torus
// take min distance as computed using a4* (real) coordinates
    
dsqmin=100000.0;
    
        for(alpha0=-1;alpha0<=1;alpha0++){
        for(alpha1=-1;alpha1<=1;alpha1++){
        for(alpha2=-1;alpha2<=1;alpha2++){
        for(alpha3=-1;alpha3<=1;alpha3++){
        alpha[0]=alpha0;
        alpha[1]=alpha1;
        alpha[2]=alpha2;
        alpha[3]=alpha3;    
            dsq=0.0;
            for(mu=0;mu<D;mu++){
                for(nu=0;nu<D;nu++){
                    dsq=dsq+(x.get(mu)-y.get(mu)-alpha[mu]*LX)*(x.get(nu)-y.get(nu)-alpha[nu]*LX)*
                    (kron[mu][nu]-1.0/NUMLINK);
            }
        }
            if(dsq<dsqmin){dsqmin=dsq;}
            //cout << "dsq is " << dsq << endl;
}}}}

 
r=sqrt(dsqmin);
return(r);
}
*/

//function to compute the check for uniqueness for r
int check(double d[],double r, int &total){

int i,id;
int found=0;

// now check if in fact distance r has been seen before
for(i=0;i<total;i++){
if(fabs(d[i]-r)<0.000001){
// if so return old index of this distance in array d[]
id=i;
found=1;
}
if(found) break;
}

// otherwise new distance - add to d[] and increment total
if(!found){
id=total;
d[total]=r;
total++;}

return(id);

}

