
#include "matrix.h"

int lat_pack(const Lattice_Vector &x){
int i,site=0;

for(i=0;i<D;i++){
site=site+Lattice_Map[i]*x.get(i);
}

return(site);
}

void lat_unpack(const int n, Lattice_Vector &x){
int current,i;

current=n;
for(i=D-1;i>=0;i--){
x.set(i,current/Lattice_Map[i]);
current=current-Lattice_Map[i]*x.get(i);
}

return;
}

int pack(const int flavor, const Lattice_Vector &x, const int color){
int i;
i=(SITES*NUMGEN)*flavor+NUMGEN*lat_pack(x)+color;

return(i);
}

void unpack(const int i, int &flavor, Lattice_Vector &x, int &color){
int tmp;
flavor=i/(SITES*NUMGEN);
tmp=i-flavor*(SITES*NUMGEN);
lat_unpack(tmp/NUMGEN,x);
color=tmp-(tmp/NUMGEN)*NUMGEN;
return;
}


int flavor(const int mu, const int nu){
int dum;

if(mu==0){
dum=nu;}
if(mu==1){
dum=5+nu-2;}
if(mu==2){
dum=8+nu-3;}
if(mu==3){
dum=10;}
return(dum);
}


#ifdef FULLMATRIX
void full_fermion_op(const Gauge_Field &U, Complex M[LEN][LEN]){
int sites,a,b,j1,j2,i,mu,nu,l,m,d,e,c,index,num_chis;
Lattice_Vector x,e_mu,e_nu,e_a,e_b,e_c;
Gauge_Field Udag,Udtr;
Umatrix Gam, Gam2;

Udag=Adj(U);


for(j1=0;j1<LEN;j1++){
for(j2=0;j2<LEN;j2++){
M[j1][j2]=Complex();
}}


// eta dmubar psi_mu

num_chis=(NUMLINK*(NUMLINK-1))/2;
//cout<< "number of chis " << num_chis << "\n" << flush;

sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMGEN;a++){

i=pack(0,x,a);

for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
for(b=0;b<NUMGEN;b++){

j1=pack(num_chis+1+mu,x,b);
j2=pack(num_chis+1+mu,x-e_mu,b);
Gam=twist(x,-e_mu);
M[i][j1]=M[i][j1]+Tr(Lambda[b]*Udag.get(x,mu)*Lambda[a]);
M[i][j2]=M[i][j2]-Tr(Lambda[a]*Udag.get(x,-e_mu,mu)*Gam*Lambda[b]*Adj(Gam));
}
}

}
}
// psi_mu dbarT eta

sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMGEN;a++){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
i=pack(mu+num_chis+1,x,a);

for(b=0;b<NUMGEN;b++){
j1=pack(0,x+e_mu,b);
j2=pack(0,x,b);
Gam=twist(x,e_mu);
M[i][j1]=M[i][j1]+Tr(Gam*Lambda[b]*Adj(Gam)*Udag.get(x,mu)*Lambda[a]);
M[i][j2]=M[i][j2]-Tr(Lambda[a]*Udag.get(x,mu)*Lambda[b]);
}}}}


// chi_{munu} dmu psi_nu


sites=0;
while(loop_over_lattice(x,sites)){

for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

for(nu=mu+1;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);

index=flavor(mu,nu);
//cout << "index is " << index << "\n" << flush;
for(a=0;a<NUMGEN;a++){
i=pack(index,x,a);

for(b=0;b<NUMGEN;b++){

j1=pack(num_chis+1+nu,x+e_mu,b);
j2=pack(num_chis+1+nu,x,b);
Gam=twist(x,e_mu);
M[i][j1]=M[i][j1]+Tr(Lambda[a]*U.get(x,mu)*Gam*Lambda[b]*Adj(Gam));
M[i][j2]=M[i][j2]-Tr(Lambda[b]*U.get(x,e_nu,mu)*Lambda[a]);
j1=pack(num_chis+1+mu,x+e_nu,b);
j2=pack(num_chis+1+mu,x,b);
Gam=twist(x,e_nu);
M[i][j1]=M[i][j1]-Tr(Lambda[a]*U.get(x,nu)*Gam*Lambda[b]*Adj(Gam));
M[i][j2]=M[i][j2]+Tr(Lambda[b]*U.get(x,e_mu,nu)*Lambda[a]);

}}}}}
// psi_nu dmu chi_{munu}


sites=0;
while(loop_over_lattice(x,sites)){

for(nu=0;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);

for(a=0;a<NUMGEN;a++){

for(mu=nu+1;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

index=flavor(nu,mu);
//cout << "index is " << index << "\n" << flush;
for(b=0;b<NUMGEN;b++){

i=pack(num_chis+1+nu,x,a);
j1=pack(index,x,b);
j2=pack(index,x-e_mu,b);
Gam=twist(x,-e_mu);
M[i][j1]=M[i][j1]-Tr(Lambda[a]*U.get(x,e_nu,mu)*Lambda[b]);
M[i][j2]=M[i][j2]+Tr(Gam*Lambda[b]*Adj(Gam)*U.get(x,-e_mu,mu)*Lambda[a]);
i=pack(num_chis+1+mu,x,a);
j1=pack(index,x,b);
j2=pack(index,x-e_nu,b);
Gam=twist(x,-e_nu);
M[i][j1]=M[i][j1]+Tr(Lambda[a]*U.get(x,e_mu,nu)*Lambda[b]);
M[i][j2]=M[i][j2]-Tr(Gam*Lambda[b]*Adj(Gam)*U.get(x,-e_nu,nu)*Lambda[a]);
}
}
}
}}


// Q-closed term

if(NUMLINK==5){
sites=0;
while(loop_over_lattice(x,sites)){
for(d=0;d<NUMLINK;d++){
for(e=d+1;e<NUMLINK;e++){

for(l=0;l<NUMGEN;l++){
i=pack(flavor(d,e),x,l);

for(c=0;c<NUMLINK;c++){
if((c==d)||(c==e)) continue;

e_c=Lattice_Vector(c);

for(a=0;a<NUMLINK;a++){
e_a=Lattice_Vector(a);
if((a==c)||(a==d)||(a==e)) continue;
for(b=a+1;b<NUMLINK;b++){
e_b=Lattice_Vector(b);
if((b==c)||(b==d)||(b==e)) continue;

for(m=0;m<NUMGEN;m++){

j1=pack(flavor(a,b),x-e_a-e_b,m);
j2=pack(flavor(a,b),x-e_a-e_b-e_c,m);
Gam=twist(x,-e_a,-e_b);
Gam2=twist(x,-e_a,-e_b,-e_c);

M[i][j1]=M[i][j1]+0.5*perm[d][e][c][a][b]*
Tr(Gam*Lambda[m]*Adj(Gam)*Udag.get(x,-e_a,-e_b,-e_c,c)*Lambda[l]);
M[i][j2]=M[i][j2]-0.5*perm[d][e][c][a][b]*
Tr(Lambda[l]*Udag.get(x,-e_c,c)*Gam2*Lambda[m]*Adj(Gam2));

}}
}

}}}}}


sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
e_a=Lattice_Vector(a);
for(b=a+1;b<NUMLINK;b++){
e_b=Lattice_Vector(b);

for(l=0;l<NUMGEN;l++){
i=pack(flavor(a,b),x,l);

for(c=0;c<NUMLINK;c++){
if((c==a)||(c==b)) continue;

e_c=Lattice_Vector(c);

for(d=0;d<NUMLINK;d++){

if((d==c)||(d==a)||(d==b)) continue;
for(e=d+1;e<NUMLINK;e++){

if((e==c)||(e==a)||(e==b)) continue;
for(m=0;m<NUMGEN;m++){

j1=pack(flavor(d,e),x+e_a+e_b+e_c,m);
j2=pack(flavor(d,e),x+e_a+e_b,m);
Gam=twist(x,e_a,e_b,e_c);
Gam2=twist(x,e_a,e_b);

M[i][j1]=M[i][j1]+0.5*perm[a][b][c][d][e]*
Tr(Gam*Lambda[m]*Adj(Gam)*Udag.get(x,e_a,e_b,c)*Lambda[l]);
M[i][j2]=M[i][j2]-0.5*perm[a][b][c][d][e]*
Tr(Lambda[l]*Udag.get(x,-e_c,c)*Gam2*Lambda[m]*Adj(Gam2));

}}
}

}}}}}
}

//test antisymmetric nature
for(j1=0;j1<LEN;j1++){
for(j2=j1+1;j2<LEN;j2++){
//cout << "M[j1][j2] is " << M[j1][j2] << "\n" <<flush;
Complex dum=M[j1][j2]+M[j2][j1];
if(dum.norm()>0.00000001) {
cout << "problem " << j1 << "\t" << j2 << "\t" << 
M[j1][j2] << "\t" << M[j2][j1] << "\n" << flush;}
}
}


return;
}


double cnorm(const Complex &c){
double dum;
dum=c.real()*c.real()+c.imag()*c.imag();
return(sqrt(dum));
}

void eigenvalues(Complex M[LEN][LEN]){
static int firsttime=1;
static ofstream f_feigen,f_det;
int i,j,ok,c1,c2,c3;
char c4;
Complex temp[2*LEN+1];
double b[2*LEN],dummy[2],work[4*LEN];
double at[2*LEN*LEN];

if(firsttime){
        f_feigen.open("feigen",ios::app);
        if(f_feigen.bad()){
        cout << "failed to open feigen file\n" << flush ;}
        firsttime=0;
}

if(1){
for(i=0;i<LEN;i++){
for(j=0;j<LEN;j++){
at[2*(j+LEN*i)]=M[j][i].real();
at[2*(j+LEN*i)+1]=M[j][i].imag();
}
}

c1=LEN;
c2=2*LEN;
c3=1;
c4='N';

zgeev_(&c4,&c4,&c1,at,&c1,b,dummy,&c3,dummy,&c3,work,&c2,work,&ok);

for(j=0;j<2*LEN;j=j+2){
temp[j/2+1]=Complex(b[j],b[j+1]);}

ceigsrt(temp,LEN);
int SUB=0;

double arg[LEN+1];
Complex phase[LEN+1];

for(i=1;i<=(LEN-SUB);i++){
arg[i]=0.0;
for(j=1;j<=i;j++){
arg[i]=arg[i]+atan2(temp[j].imag(),temp[j].real());
}
phase[i]=Complex(cos(arg[i]),sin(arg[i]));}


Complex tmp;
int count=0;
for(i=1;i<=(LEN);i++){
if(temp[i].norm()<0.000001) {count++;}
//f_feigen << i << "\t" << phase[i] << endl;
f_feigen << temp[i] << endl;}

f_feigen << "\n";

cout << "number of zeromodes  is " << count << endl;


}


return;
}

Complex Pfaffian(Complex M[LEN][LEN])
{
	int i,j,k,jpiv,interchange=1;
        static int firsttime=1;
	double pivot,totalangle,angle,cosine,sine,mag;
	Complex dum[LEN],scale,f;
//        cout << "in Pfaffian\n" << flush;
        static ofstream f_pfaff;
       
        if(firsttime){ 
        f_pfaff.open("pfaffian",ios::app);
        if(f_pfaff.bad()){
        cout << "failed to open pfaffian file\n" << flush ;}
        firsttime=0;
}
// loop over all rows in steps of 2
	for(i=0;i<LEN-2;i+=2){

// first row i:	
// find col whose ith component is biggest to use as pivot
	pivot=cnorm(M[i][i+1]);
	jpiv=i+1;
	for(j=i+2;j<LEN;j++){
	if(cnorm(M[i][j])>pivot){pivot=cnorm(M[i][j]);jpiv=j;}}

// interchange col(i+1) with col(jpiv)
        for(j=0;j<LEN;j++){
	dum[j]=M[j][i+1];}
	for(j=0;j<LEN;j++){
	M[j][i+1]=M[j][jpiv];
	M[j][jpiv]=dum[j];}
// interchange row(i+1) with row(jpiv)
	for(j=0;j<LEN;j++){
	dum[j]=M[i+1][j];}
	for(j=0;j<LEN;j++){
	M[i+1][j]=M[jpiv][j];
	M[jpiv][j]=dum[j];}
	
	if(jpiv!=i+1) interchange*=(-1);
	
// using this zero progressively elements of row M[i][j], j=i+2...LEN-1
        for(j=i+2;j<LEN;j++){
	scale=M[i][j]/M[i][i+1];
	
	for(k=0;k<LEN;k++){
	M[k][j]=M[k][j]-scale*M[k][i+1];}
// zero out elements along corresponding column M[j][i] too
	for(k=0;k<LEN;k++){
	M[j][k]=M[j][k]-scale*M[i+1][k];}
	}

	
// next row i+1;
        
	// using this zero progressively elements M[i][j], j=i+2...LEN-1
        for(j=i+2;j<LEN;j++){
	scale=M[i+1][j]/M[i+1][i];
	
	for(k=0;k<LEN;k++){
	M[k][j]=M[k][j]-scale*M[k][i];}
// zero out elements along corresponding column too
	for(k=0;k<LEN;k++){
	M[j][k]=M[j][k]-scale*M[i][k];}
	}
	
	}
	
	
	f=Complex(1.0,0.0);
        mag=0.0;
	totalangle=0.0;
	for(i=0;i<LEN;i+=2){
        cosine=M[i][i+1].real()/M[i][i+1].norm();
        sine=M[i][i+1].imag()/M[i][i+1].norm();
        if ((cosine>0.0) && (sine>0.0)) angle=acos(cosine);
        if ((cosine<0.0) && (sine>0.0)) angle=acos(cosine);
        if ((cosine>0.0) && (sine<0.0)) angle=4*acos(0.0)-1.0*acos(cosine);
        if ((cosine<0.0) && (sine<0.0)) angle=4*acos(0.0)-1.0*acos(cosine);
        mag+=log(M[i][i+1].norm());
        totalangle+=angle;
	f=f*M[i][i+1];
	}
        
        cout << "Pfaffian is " << Complex(cos(totalangle)*interchange,sin(totalangle)*interchange) << endl;

        f_pfaff << mag << "\t" << cos(totalangle)*interchange << "\t" 
                << sin(totalangle)*interchange << "\n" << flush;
return (f*interchange);
}
#endif
