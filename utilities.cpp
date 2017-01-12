#include "utilities.h"

// Complex type
Complex::Complex(void){re=0.0;im=0.0;}
Complex::Complex(double x, double y){re=x;im=y;}
double Complex::real(void) const{return(re);}
double Complex::imag(void) const{return(im);}
double Complex::norm(void){return(sqrt(re*re+im*im));}
void Complex::print(void) const {cout << "("<< re << ", " << im << ")";}

ostream& operator<<(ostream& out,Complex c){
  out<<c.real()<<"\t"<<c.imag();
  return out;}
istream& operator>>(istream& in,Complex & c){
  double x,y;
  in>>x>>y;
  c=Complex(x,y);
  return in;}


Complex operator /(const Complex &o1, const Complex &o2){
  Complex dum;
  double norm;
  norm=o2.real()*o2.real()+o2.imag()*o2.imag();
  dum=Complex((o1.real()*o2.real()+o1.imag()*o2.imag())/norm,
    (o1.imag()*o2.real()-o1.real()*o2.imag())/norm);
  return(dum);}
Complex pow(const Complex &o1, const int o2){
  Complex c(1,0);
  for (int i=0;i<o2;i++) c=c*o1;
  return c;}



// unitary matrix type Umatrix
Umatrix::Umatrix(void){
  for(int i=0;i<NCOLOR;i++){
    for(int j=0;j<NCOLOR;j++) {mat[i][j]=Complex();}}
  }

Umatrix::Umatrix(int k){
         if(k==1){
         for(int i=0;i<NCOLOR;i++){
          for(int j=0;j<NCOLOR;j++) {mat[i][j]=Complex();}
    }
   for(int n=0;n<NCOLOR;n++){
   mat[n][n]=Complex(1.0,0.0);}
   }
   else {cout << "wrong Umatrix constructor\n" << flush; }

}

Umatrix::Umatrix(Complex m[NCOLOR][NCOLOR]){
  for(int i=0;i<NCOLOR;i++)
    for(int j=0;j<NCOLOR;j++) mat[i][j]=m[i][j];}


Complex Umatrix::get(int i, int j) const {return(mat[i][j]);}
void Umatrix::set(int i, int j, const Complex o){mat[i][j]=o;}
void Umatrix::print(void){
  for(int i=0;i<NCOLOR;i++){
    for(int j=0;j<NCOLOR;j++){mat[i][j].print();cout << "\t";}
    cout << "\n";}
    }
Umatrix Adj(const Umatrix &u){
        Umatrix res;
        for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++){
  res.set(i,j,conjug(u.get(j,i)));}
  return(res);}


Umatrix Trans(const Umatrix &u){
        Umatrix res;
        for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++){
  res.set(i,j,u.get(j,i));}
  return(res);}


ostream& operator<<(ostream& out,Umatrix s){
for(int i=0;i<NCOLOR;i++)
for(int j=0;j<NCOLOR;j++){
  out<<s.get(i,j)<<'\t';}
  return out;}
istream& operator>>(istream& in, Umatrix & s){
  Complex v[NCOLOR][NCOLOR];
  for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++){
  in>>v[i][j];}
  s=Umatrix(v);
  return in;
   }


Umatrix operator *(const Umatrix &o1, const Umatrix &o2){
         Umatrix r;
   Complex dum;
   for(int i=0;i<NCOLOR;i++)
   for(int j=0;j<NCOLOR;j++){
   dum=Complex();
   for(int k=0;k<NCOLOR;k++){
   dum=dum+o1.get(i,k)*o2.get(k,j);}
   r.set(i,j,dum);}
   return(r);}

        Umatrix operator *(const Umatrix &o1, const Complex &o2){
        Umatrix dum;
  for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++){
  dum.set(i,j,o1.get(i,j)*o2);}
  return(dum);}
        Umatrix operator *(const Complex &o2, const Umatrix &o1){
  Umatrix dum;
  for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++){
  dum.set(i,j,o1.get(i,j)*o2);}
  return(dum);}
  Umatrix operator *(const Umatrix &o1, const double o2){
        Umatrix dum;
  for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++){
  dum.set(i,j,o1.get(i,j)*o2);}
  return(dum);}
  Umatrix operator *(const double o2, const Umatrix &o1){
        Umatrix dum;
  for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++){
  dum.set(i,j,o1.get(i,j)*o2);}
  return(dum);}
  Umatrix operator +(const Umatrix &x, const Umatrix &y){
  Umatrix dum;
  for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++)
  dum.set(i,j,x.get(i,j)+y.get(i,j));
  return(dum);
  }
        Umatrix operator -(const Umatrix &x, const Umatrix &y){
  Umatrix dum;
  for(int i=0;i<NCOLOR;i++)
  for(int j=0;j<NCOLOR;j++)
  dum.set(i,j,x.get(i,j)-y.get(i,j));
  return(dum);
  }


    Umatrix exp(const Umatrix &u){
  Umatrix c,del,prod;
  double fac=1.0;
        int i=1;
        prod=Umatrix(1);
        c=Umatrix(1);
        static int sum=0,counter=0;

  do{
        fac=fac*(double)i;
        prod=prod*u;
        del=prod*(1.0/fac);
        c=c+del;
        i++;}
        while(sqrt(Tr(del*Adj(del)).real())>GAUGETOL);

        sum+=i;
        counter++;
        if(counter==100000){
        cout << "mean no. of terms in exp() "
        << (double)sum/counter << "\n" << flush;
        counter=0;sum=0;}

  return(c);}



Complex Tr(const Umatrix &o){
         Complex dum=Complex();
   for(int i=0;i<NCOLOR;i++)
   dum=dum+o.get(i,i);
   return(dum);}

Umatrix gaussU(void){
  int i;
  Umatrix tmp=Umatrix();
  for(i=0;i<NUMGEN;i++){
  tmp=tmp+1.0/sqrt(2.0)*Complex(gasdev(),gasdev())*Lambda[i];}
  return(tmp);
  }


// adjoint field matrix type Afield
Afield::Afield(void){
  for(int i=0;i<NUMGEN;i++)
    {afield[i]=Complex();}
  }

Afield::Afield(int c){
if(c==1){
for(int i=0;i<NUMGEN;i++){
afield[i]=Complex(1.0/sqrt(2.0)*gasdev(),1.0/sqrt(2.0)*gasdev());}
}
else{cout << "error initializing afield\n"; exit(1);}
return;
}

Afield::Afield(Complex m[NUMGEN]){
  for(int i=0;i<NUMGEN;i++)
    afield[i]=m[i];}

Afield gaussA(void){
Afield dum=Afield();

for(int i=0;i<NUMGEN;i++){
dum.set(i,Complex(1.0/sqrt(2.0)*gasdev(),1.0/sqrt(2.0)*gasdev()));
}

return(dum);
}

Complex Afield::get(int i) const {return(afield[i]);}
void Afield::set(int i, const Complex o){afield[i]=o;}
void Afield::print(void){
  for(int i=0;i<NUMGEN;i++){
    afield[i].print();cout << "\t";}
    cout << "\n";
    }

// nb using antihermitian basis for Lambdas
Afield Cjg(const Afield &u){
        Afield res;
        for(int i=0;i<NUMGEN;i++){
  res.set(i,conjug(u.get(i)));}
  return(res);}

ostream& operator<<(ostream& out,Afield s){
for(int i=0;i<NUMGEN;i++)
{
  out<<s.get(i)<<'\t';}
  return out;}

istream& operator>>(istream& in, Afield & s){
  Complex v[NUMGEN];
  for(int j=0;j<NUMGEN;j++){
  in>>v[j];}
  s=Afield(v);
  return in;
  }


        Afield operator *(const Afield &o1, const Complex &o2){
        Afield dum;
  for(int i=0;i<NUMGEN;i++)
  {
  dum.set(i,o1.get(i)*o2);}
  return(dum);}
        Afield operator *(const Complex &o2, const Afield &o1){
  Afield dum;
  for(int i=0;i<NUMGEN;i++)
  {
  dum.set(i,o1.get(i)*o2);}
  return(dum);}
  Afield operator *(const Afield &o1, const double o2){
        Afield dum;
  for(int i=0;i<NUMGEN;i++)
  {
  dum.set(i,o1.get(i)*o2);}
  return(dum);}
  Afield operator *(const double o2, const Afield &o1){
        Afield dum;
  for(int i=0;i<NUMGEN;i++)
  {
  dum.set(i,o1.get(i)*o2);}
  return(dum);}
  Afield operator +(const Afield &x, const Afield &y){
  Afield dum;
  for(int i=0;i<NUMGEN;i++){
  dum.set(i,x.get(i)+y.get(i));}
  return(dum);
  }
        Afield operator -(const Afield &x, const Afield &y){
  Afield dum;
  for(int i=0;i<NUMGEN;i++){
  dum.set(i,x.get(i)-y.get(i));}
  return(dum);
  }



 Complex operator *(const Afield &A, const Afield &B){
 Complex dum=Complex();
 for(int i=0;i<NUMGEN;i++){
 dum=dum+A.get(i)*B.get(i);}
 return(dum);
 }


Lattice_Vector::Lattice_Vector(void){for(int i=0;i<D;i++)coords[i]=0;}
Lattice_Vector::Lattice_Vector(int mu){
if(mu<D){
for(int i=0;i<D;i++){
coords[i]=0;}
coords[mu]=1;
return;}

for(int i=0;i<D;i++){
coords[i]=-1;
}
//cout << "here ...\n";
return;
}

void Lattice_Vector::set(int i, int a){
coords[i]=a;
return;}
int Lattice_Vector::get(int i) const{
return(coords[i]);}
void Lattice_Vector::print(void) const {
for(int i=0;i<D;i++) {cout << coords[i] << "\t";} cout << "\n";}


Lattice_Vector operator +(const Lattice_Vector &x, const Lattice_Vector &y){
Lattice_Vector dum;
for(int i=0;i<D;i++){
if((x.get(i)+y.get(i))>=0){
dum.set(i,(x.get(i)+y.get(i))%side[i]);}
if((x.get(i)+y.get(i))<0){
dum.set(i,(x.get(i)+y.get(i)+side[i])%side[i]);}
}


return(dum);}

Lattice_Vector operator -(const Lattice_Vector &x, const Lattice_Vector &y){
Lattice_Vector dum;
for(int i=0;i<D;i++){
if((x.get(i)-y.get(i))>=0){
dum.set(i,(x.get(i)-y.get(i))%side[i]);}
if((x.get(i)-y.get(i))<0){
dum.set(i,(x.get(i)-y.get(i)+side[i])%side[i]);}
}


return(dum);}


Lattice_Vector operator -(const Lattice_Vector &x){
Lattice_Vector dum;
for(int i=0;i<D;i++)
dum.set(i,-1*x.get(i));
return(dum);
}

int length(const Lattice_Vector &x){
int dum=0,tmp;
for(int i=0;i<(D-1);i++){
tmp=x.get(i);
if(tmp>(LX/2)){tmp=tmp-LX/2;}
dum+=tmp*tmp;}

tmp=x.get(D-1);
if(tmp>(T/2)){tmp=tmp-T/2;}
dum+=tmp*tmp;

return(dum);
}

double BC(const Lattice_Vector &x, const Lattice_Vector &y){

if(x.get(D-1)+y.get(D-1)<0) return(PBC);
if(x.get(D-1)+y.get(D-1)>(T-1))return(PBC);

return(1.0);
}

double BC(const Lattice_Vector &x, const Lattice_Vector &y,
const Lattice_Vector &z){

if(x.get(D-1)+y.get(D-1)+z.get(D-1)<0) return(PBC);
if(x.get(D-1)+y.get(D-1)+z.get(D-1)>(T-1))return(PBC);

return(1.0);
}

double BC(const Lattice_Vector &x, const Lattice_Vector &y,
const Lattice_Vector &z, const Lattice_Vector &w){

if(x.get(D-1)+y.get(D-1)+z.get(D-1)+w.get(D-1)<0) return(PBC);
if(x.get(D-1)+y.get(D-1)+z.get(D-1)+w.get(D-1)>(T-1))return(PBC);

return(1.0);
}

int blocklatsite(Lattice_Vector &x ){
int i,c;

c=0;
for(i=0;i<D;i++){
if(x.get(i)%2==0){c++;}
}

if(c!=4) {return(0);}
else
{return(1);}
}

int loop_over_lattice(Lattice_Vector &x, int &site){
int test,i,current;


current=site;
for(i=D-1;i>=0;i--){
x.set(i,current/Lattice_Map[i]);
current=current-Lattice_Map[i]*x.get(i);
}

if(current!=0){cout << "error in loop_over_lattice" << "\n";}

if(site==SITES)
test=1;
else
test=0;

site++;
return(!test);
}



Gauge_Field::Gauge_Field(void){
for(int i=0;i<SITES;i++)
for(int j=0;j<NUMLINK;j++){
link[i][j]=Umatrix();}
}

Gauge_Field::Gauge_Field(int hot){
if(hot==2){
for(int i=0;i<SITES;i++){
for(int j=0;j<NUMLINK;j++){
link[i][j]=gaussU();}}
return;
}

if(hot==0){
for(int i=0;i<SITES;i++){
for(int j=0;j<NUMLINK;j++){
link[i][j]=exp(0.01*gaussU());}
}
return;
}


cout << "error in gauge field constructor " << "\n" << flush;

return;
}


Umatrix Gauge_Field::get(const Lattice_Vector &x, const int mu) const{
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(link[site][mu]);
}

void Gauge_Field::set(const Lattice_Vector &x, const int mu, const Umatrix &u){
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

link[site][mu]=u;
return;
}


Gauge_Field Adj(const Gauge_Field & U){
int site;
Gauge_Field W;
Lattice_Vector x;
site=0;
while(loop_over_lattice(x,site)){
for(int mu=0;mu<NUMLINK;mu++){
W.set(x,mu,Adj(U.get(x,mu)));
}}
return W;
}

Gauge_Field Transpose(const Gauge_Field & U){
int site;
Gauge_Field W;
Lattice_Vector x;
site=0;
while(loop_over_lattice(x,site)){
for(int mu=0;mu<NUMLINK;mu++){
W.set(x,mu,Trans(U.get(x,mu)));
}}
return W;
}


USite_Field::USite_Field(void){
for(int i=0;i<SITES;i++){
points[i]=Umatrix();}
return;
}

USite_Field::USite_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
points[i]=gaussU();
}
return;}

if(c==0){
for(int i=0;i<SITES;i++){
points[i]=Umatrix(1);}

return;}

cout << "error in site constructor\n" << flush;

}

Umatrix USite_Field::get(const Lattice_Vector &x) const{
int site=0,i;

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(points[site]);
}

void USite_Field::set(const Lattice_Vector &x, const Umatrix &u){
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

points[site]=u;
return;
}

void USite_Field::print(void){
cout << "site field values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;
cout << points[i] << "\n" << flush;}
return;
}


USite_Field operator +(const USite_Field &s1, const USite_Field &s2){
int sites=0;
Lattice_Vector x;
USite_Field dum=USite_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,s1.get(x)+s2.get(x));
}

return(dum);
}


USite_Field operator -(const USite_Field &s1, const USite_Field &s2){
int sites=0;
Lattice_Vector x;
USite_Field dum=USite_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,s1.get(x)-s2.get(x));
}

return(dum);
}

USite_Field operator *(const double o, const USite_Field &s){
int sites=0;
Lattice_Vector x;
USite_Field dum=USite_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,o*s.get(x));
}

return(dum);
}


USite_Field operator *(const Complex &o, const USite_Field &s){
int sites=0;
Lattice_Vector x;
USite_Field dum=USite_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,o*s.get(x));
}

return(dum);
}

Umatrix operator *(const USite_Field &s1, const USite_Field &s2){
int sites=0;
Lattice_Vector x;
Umatrix dum;

dum=Umatrix();
while(loop_over_lattice(x,sites)){
dum=dum+s1.get(x)*s2.get(x);
}
return(dum);
}

USite_Field Adj(const USite_Field &l){
int sites;
Lattice_Vector x;
USite_Field dum;

sites=0;
while(loop_over_lattice(x,sites)){
dum.set(x,Adj(l.get(x)));}

return(dum);
}


UPlaq_Field::UPlaq_Field(void){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
for(int nu=0;nu<NUMLINK;nu++){
square[i][mu][nu]=Umatrix();}}}
return;
}

UPlaq_Field::UPlaq_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
square[i][mu][mu]=Umatrix();
for(int nu=mu+1;nu<NUMLINK;nu++){
square[i][mu][nu]=gaussU();
square[i][nu][mu]=-1.0*square[i][mu][nu];
}
}
}
return;
}
cout << "error in link constructor\n" << "\n" << flush;

return;
}

Umatrix UPlaq_Field::get(const Lattice_Vector &x, const int mu, const int nu) const{
int site=0,i;

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(square[site][mu][nu]);
}

void UPlaq_Field::set(const Lattice_Vector &x, const int mu, const int nu, const Umatrix &u){
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

square[site][mu][nu]=u;
return;
}

void UPlaq_Field::print(void){
cout << "link field values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
cout << "j,k= " << j << k << ":" << square[i][j][k] << "\n" << flush;}
cout << "\n" << flush;}}
return;
}

UPlaq_Field operator +(const UPlaq_Field &s1, const UPlaq_Field &s2){
int sites=0;
Lattice_Vector x;
UPlaq_Field dum=UPlaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,s1.get(x,j,k)+s2.get(x,j,k));
}}
}
return(dum);
}

UPlaq_Field operator -(const UPlaq_Field &s1, const UPlaq_Field &s2){
int sites=0;
Lattice_Vector x;
UPlaq_Field dum=UPlaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,s1.get(x,j,k)-s2.get(x,j,k));
}}}

return(dum);
}

UPlaq_Field Plaq(const Gauge_Field &U){
Lattice_Vector x,e_mu,e_nu;
int site=0,mu,nu;
UPlaq_Field P;

while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=mu+1;nu<NUMLINK;nu++){
e_mu=Lattice_Vector(mu);
e_nu=Lattice_Vector(nu);
P.set(x,mu,nu,
U.get(x,mu)*U.get(x+e_mu,nu)*Adj(U.get(x+e_nu,mu))*Adj(U.get(x,nu)));
P.set(x,nu,mu,Adj(P.get(x,mu,nu)));
}}}
return(P);
}

UPlaq_Field operator *(const double o, const UPlaq_Field &s){
int sites=0;
Lattice_Vector x;
UPlaq_Field dum=UPlaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,o*s.get(x,j,k));
}}}

return(dum);
}

UPlaq_Field operator *(const Complex &o, const UPlaq_Field &s){
int sites=0;
Lattice_Vector x;
UPlaq_Field dum=UPlaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,o*s.get(x,j,k));
}}}

return(dum);
}

Umatrix operator *(const UPlaq_Field &s1, const UPlaq_Field &s2){
int sites=0;
Lattice_Vector x;
Umatrix dum;

dum=Umatrix();
while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=j+1;k<NUMLINK;k++){
dum=dum+s1.get(x,j,k)*s2.get(x,j,k);
}}}
return(dum);
}

UPlaq_Field Adj(const UPlaq_Field &l){
int sites=0;
Lattice_Vector x;
UPlaq_Field dum;

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,Adj(l.get(x,j,k)));}
}}

return(dum);
}


Site_Field::Site_Field(void){
for(int i=0;i<SITES;i++){
points[i]=Afield();}
return;
}

Site_Field::Site_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
points[i]=gaussA();
}
return;}

cout << "error in site constructor\n" << flush;

}

Afield Site_Field::get(const Lattice_Vector &x) const{
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(points[site]);
}

void Site_Field::set(const Lattice_Vector &x, const Afield &u){
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

points[site]=u;
return;
}

void Site_Field::print(void){
cout << "site field values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;
cout << points[i] << "\n" << flush;}
return;
}


Site_Field operator +(const Site_Field &s1, const Site_Field &s2){
int sites=0;
Lattice_Vector x;
Site_Field dum=Site_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,s1.get(x)+s2.get(x));
}

return(dum);
}


Site_Field operator -(const Site_Field &s1, const Site_Field &s2){
int sites=0;
Lattice_Vector x;
Site_Field dum=Site_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,s1.get(x)-s2.get(x));
}

return(dum);
}

Site_Field operator *(const double o, const Site_Field &s){
int sites=0;
Lattice_Vector x;
Site_Field dum=Site_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,o*s.get(x));
}

return(dum);
}


Site_Field operator *(const Complex &o, const Site_Field &s){
int sites=0;
Lattice_Vector x;
Site_Field dum=Site_Field();

while(loop_over_lattice(x,sites)){
dum.set(x,o*s.get(x));
}

return(dum);
}

Complex operator *(const Site_Field &s1, const Site_Field &s2){
int sites=0;
Lattice_Vector x;
Complex dum;

dum=Complex();
while(loop_over_lattice(x,sites)){
dum=dum+s1.get(x)*s2.get(x);
}
return(dum);
}

Site_Field Cjg(const Site_Field &l){
int sites;
Lattice_Vector x;
Site_Field dum;

sites=0;
while(loop_over_lattice(x,sites)){
dum.set(x,Cjg(l.get(x)));}

return(dum);
}


Link_Field::Link_Field(void){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
links[i][mu]=Afield();}}
return;
}

Link_Field::Link_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
links[i][mu]=gaussA();
}
}
return;
}
cout << "error in link constructor\n" << "\n" << flush;

return;
}

Afield Link_Field::get(const Lattice_Vector &x, const int mu) const{
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(links[site][mu]);
}

void Link_Field::set(const Lattice_Vector &x, const int mu, const Afield &u){
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

links[site][mu]=u;
return;
}

void Link_Field::print(void){
cout << "link field values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;
for(int j=0;j<NUMLINK;j++){
cout << "j= " << j << ":" << links[i][j] << "\n" << flush;}
cout << "\n" << flush;}
return;
}

Link_Field operator +(const Link_Field &s1, const Link_Field &s2){
int sites=0;
Lattice_Vector x;
Link_Field dum=Link_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum.set(x,j,s1.get(x,j)+s2.get(x,j));
}}

return(dum);
}

Link_Field operator -(const Link_Field &s1, const Link_Field &s2){
int sites=0;
Lattice_Vector x;
Link_Field dum=Link_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum.set(x,j,s1.get(x,j)-s2.get(x,j));
}}

return(dum);
}


Link_Field operator *(const double o, const Link_Field &s){
int sites=0;
Lattice_Vector x;
Link_Field dum=Link_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum.set(x,j,o*s.get(x,j));
}}

return(dum);
}

Link_Field operator *(const Complex &o, const Link_Field &s){
int sites=0;
Lattice_Vector x;
Link_Field dum=Link_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum.set(x,j,o*s.get(x,j));
}}

return(dum);
}

Complex operator *(const Link_Field &s1, const Link_Field &s2){
int sites=0;
Lattice_Vector x;
Complex dum;

dum=Complex();
while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum=dum+s1.get(x,j)*s2.get(x,j);
}}
return(dum);
}

Link_Field Cjg(const Link_Field &l){
int sites,mu;
Lattice_Vector x;
Link_Field dum;

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
dum.set(x,mu,Cjg(l.get(x,mu)));}
}

return(dum);
}



Plaq_Field::Plaq_Field(void){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
for(int nu=0;nu<NUMLINK;nu++){
square[i][mu][nu]=Afield();}}}
return;
}

Plaq_Field::Plaq_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
square[i][mu][mu]=Afield();
for(int nu=mu+1;nu<NUMLINK;nu++){
square[i][mu][nu]=gaussA();
square[i][nu][mu]=-1.0*square[i][mu][nu];
}
}
}
return;
}
cout << "error in link constructor\n" << "\n" << flush;

return;
}

Afield Plaq_Field::get(const Lattice_Vector &x, const int mu, const int nu) const{
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(square[site][mu][nu]);
}

void Plaq_Field::set(const Lattice_Vector &x, const int mu, const int nu,
const Afield &u){
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

square[site][mu][nu]=u;
return;
}

void Plaq_Field::print(void){
cout << "plaq field values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
cout << "j,k= " << j << k << ":" << square[i][j][k] << "\n" << flush;}
cout << "\n" << flush;}}
return;
}

Plaq_Field operator +(const Plaq_Field &s1, const Plaq_Field &s2){
int sites=0;
Lattice_Vector x;
Plaq_Field dum=Plaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,s1.get(x,j,k)+s2.get(x,j,k));
}}
}
return(dum);
}

Plaq_Field operator -(const Plaq_Field &s1, const Plaq_Field &s2){
int sites=0;
Lattice_Vector x;
Plaq_Field dum=Plaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,s1.get(x,j,k)-s2.get(x,j,k));
}}}

return(dum);
}


Plaq_Field operator *(const double o, const Plaq_Field &s){
int sites=0;
Lattice_Vector x;
Plaq_Field dum=Plaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,o*s.get(x,j,k));
}}}

return(dum);
}

Plaq_Field operator *(const Complex &o, const Plaq_Field &s){
int sites=0;
Lattice_Vector x;
Plaq_Field dum=Plaq_Field();

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,o*s.get(x,j,k));
}}}

return(dum);
}

Complex operator *(const Plaq_Field &s1, const Plaq_Field &s2){
int sites=0;
Lattice_Vector x;
Complex dum;

dum=Complex();
while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=j+1;k<NUMLINK;k++){
dum=dum+s1.get(x,j,k)*s2.get(x,j,k);
}}}
return(dum);
}

Plaq_Field Cjg(const Plaq_Field &l){
int sites=0;
Lattice_Vector x;
Plaq_Field dum;

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,Cjg(l.get(x,j,k)));}
}}

return(dum);
}



Twist_Fermion::Twist_Fermion(void){
L=Link_Field();
S=Site_Field();
C=Plaq_Field();
return;
}

Twist_Fermion::Twist_Fermion(int c){
if(c==1){
L=Link_Field(1);
S=Site_Field(1);
C=Plaq_Field(1);
return;}
cout << "error in Twist_Fermion constructor\n" << flush;

}

Twist_Fermion Cjg(const Twist_Fermion &K){
Twist_Fermion dum;

dum.setS(Cjg(K.getS()));
dum.setL(Cjg(K.getL()));
dum.setC(Cjg(K.getC()));

return(dum);
}


const Link_Field& Twist_Fermion::getL(void) const{
return(L);
}

const Site_Field& Twist_Fermion::getS(void) const{
return(S);
}

const Plaq_Field& Twist_Fermion::getC(void) const{
return(C);
}

void Twist_Fermion::setL(const Link_Field &l){
L=l;
}

void Twist_Fermion::setS(const Site_Field &s){
S=s;
}

void Twist_Fermion::setC(const Plaq_Field &c){
C=c;
}

void Twist_Fermion::print(void){
cout << "Twist_Fermion: \n" << flush;

S.print();
L.print();
C.print();
}

Twist_Fermion operator +(const Twist_Fermion &k1, const Twist_Fermion &k2){
Twist_Fermion dum=Twist_Fermion();

dum.setS(k1.getS()+k2.getS());
dum.setL(k1.getL()+k2.getL());
dum.setC(k1.getC()+k2.getC());

return(dum);
}

Twist_Fermion operator -(const Twist_Fermion &k1, const Twist_Fermion &k2){
Twist_Fermion dum=Twist_Fermion();

dum.setS(k1.getS()-k2.getS());
dum.setL(k1.getL()-k2.getL());
dum.setC(k1.getC()-k2.getC());

return(dum);
}


Twist_Fermion operator *(const double o, const Twist_Fermion &k){
Twist_Fermion dum=Twist_Fermion();

dum.setS(o*k.getS());
dum.setL(o*k.getL());
dum.setC(o*k.getC());

return(dum);
}

Complex operator *(const Twist_Fermion &k1, const Twist_Fermion &k2){
Complex tmp;

tmp=k1.getS()*k2.getS()+k1.getL()*k2.getL()+k1.getC()*k2.getC();
return(tmp);
}

Adjoint_Matrix::Adjoint_Matrix(void){
for(int i=0;i<NUMGEN;i++)
for(int j=0;j<NUMGEN;j++){
amat[i][j]=Complex();}
return;
}

Complex Adjoint_Matrix::get(const int i, const int j) const{
return(amat[i][j]);}

void Adjoint_Matrix::set(const int i, const int j, const Complex &c){
amat[i][j]=c;
return;
}

Adjoint_Matrix::Adjoint_Matrix(Complex m[NUMGEN][NUMGEN]){
for(int i=0;i<NUMGEN;i++){
for(int j=0;j<NUMGEN;j++){
amat[i][j]= m[i][j];
}}
return;
}

ostream& operator<<(ostream& out, Adjoint_Matrix s){
for(int i=0;i<NUMGEN;i++){
for(int j=0;j<NUMGEN;j++){
  out<<s.get(i,j)<<'\t';}}
  return out;}

istream& operator>>(istream& in, Adjoint_Matrix & s){
  Complex v[NUMGEN][NUMGEN];
  for(int j=0;j<NUMGEN;j++){
  for(int i=0;i<NUMGEN;i++){
  in>>v[j][i];}}
  s=Adjoint_Matrix(v);
  return in;
  }


Adjoint_Links::Adjoint_Links(void){
for(int i=0;i<SITES;i++)
for(int j=0;j<NUMLINK;j++){
alinks[i][j]=Adjoint_Matrix();
}
return;
}

Adjoint_Matrix Adjoint_Links::get(const Lattice_Vector &x, const int mu) const{
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(alinks[site][mu]);
}


void Adjoint_Links::set(const Lattice_Vector &x, const int mu,
const Adjoint_Matrix &u){
int site=0,i;

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

alinks[site][mu]=u;
return;
}

void Adjoint_Links::print(void){
cout << "adjoint links values\n" << flush;
for(int i=0;i<SITES;i++){
cout << "site= " << i << "\n" << flush;
for(int j=0;j<D;j++){
cout << "j= " << j << ":" << alinks[i][j]<< "\n" << flush;}}
return;
}

void compute_Adjoint_Links(const Gauge_Field &U, Adjoint_Links &V){
int mu,a,b,sites;
Lattice_Vector x;
Adjoint_Matrix tmp;

//cout << "in compute_Adjoint_Links\n" << flush;

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){

for(a=0;a<NUMGEN;a++){
for(b=0;b<NUMGEN;b++){
tmp.set(a,b,Tr(Lambda[a]*U.get(x,mu)*Lambda[b]));
}
}
V.set(x,mu,tmp);

}
}

return;
}



Plaq_Field Dplus(const Adjoint_Links &V, const Link_Field &L){
Lattice_Vector x,e_mu,e_nu;
int mu,nu,sites,a,b;
Complex tmp;
Afield atmp=Afield();
Plaq_Field dum=Plaq_Field();


sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=mu+1;nu<NUMLINK;nu++){

e_mu=Lattice_Vector(mu);
e_nu=Lattice_Vector(nu);

for(a=0;a<NUMGEN;a++){
tmp=Complex();
for(b=0;b<NUMGEN;b++){
tmp=tmp+
    V.get(x,mu).get(a,b)*L.get(x+e_mu,nu).get(b)*BC(x,e_mu)-
    V.get(x+e_nu,mu).get(b,a)*L.get(x,nu).get(b);
tmp=tmp-
    V.get(x,nu).get(a,b)*L.get(x+e_nu,mu).get(b)*BC(x,e_nu)+
    V.get(x+e_mu,nu).get(b,a)*L.get(x,mu).get(b);
    }
atmp.set(a,tmp);
}

dum.set(x,mu,nu,atmp);
dum.set(x,nu,mu,-1.0*atmp);

}
}
}

return(dum);
}

Link_Field Dminus(const Adjoint_Links &V, const Plaq_Field &P){
Lattice_Vector x,e_mu,e_nu;
int sites,mu,nu,a,b;
Complex tmp;
Afield atmp;
Link_Field dum=Link_Field();

sites=0;
while(loop_over_lattice(x,sites)){
for(nu=0;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
atmp=Afield();

for(mu=0;mu<NUMLINK;mu++){
if(mu==nu) continue;

e_mu=Lattice_Vector(mu);

for(a=0;a<NUMGEN;a++){
tmp=Complex();
for(b=0;b<NUMGEN;b++){
tmp=tmp+
V.get(x+e_nu,mu).get(a,b)*P.get(x,mu,nu).get(b)-
V.get(x-e_mu,mu).get(b,a)*P.get(x-e_mu,mu,nu).get(b)*BC(x,-e_mu);
}
atmp.set(a,atmp.get(a)+tmp);
}

}
dum.set(x,nu,atmp);

}
}
return(dum);
}


double order(const int i, const int j, const int k, const int l, const int m){
int seq[5];
seq[0]=i;seq[1]=j;seq[2]=k;seq[3]=l;seq[4]=m;
int swap,tmp,p;
double permutation=1.0;
do{
swap=0;
for(p=0;p<4;p++){
if(seq[p]>seq[p+1]) {tmp=seq[p];seq[p]=seq[p+1];seq[p+1]=tmp;
                     swap++;permutation*=(-1.0);}
        }
}
while(swap>0);
return(permutation);
}

void epsilon(void){
int i,j,k,l,m;
for(i=0;i<NUMLINK;i++)
for(j=0;j<NUMLINK;j++)
for(k=0;k<NUMLINK;k++)
for(l=0;l<NUMLINK;l++)
for(m=0;m<NUMLINK;m++)
perm[i][j][k][l][m]=0.0;

for(i=0;i<NUMLINK;i++){
for(j=0;j<NUMLINK;j++){
if(j==i) continue;
for(k=0;k<NUMLINK;k++){
if((k==j)||(k==i)) continue;
for(l=0;l<NUMLINK;l++){
if((l==k)||(l==j)||(l==i)) continue;
for(m=0;m<NUMLINK;m++){
if((m==l)||(m==k)||(m==j)||(m==i)) continue;
perm[i][j][k][l][m]=order(i,j,k,l,m);
//cout << "eps(" << i << "\t" << j << "\t" << k << "\t" << l << "\t" << m << ")" <<
//perm[i][j][k][l][m] << "\n";
}}}}}

return;
}

Plaq_Field Dbminus(const Adjoint_Links &V, const Plaq_Field &p){
Lattice_Vector x,e_a,e_b,e_c;
int a,b,c,d,e,l,m,sites;
Plaq_Field dum=Plaq_Field();
Complex tmp;
Afield atmp;


sites=0;
while(loop_over_lattice(x,sites)){
for(d=0;d<NUMLINK;d++){
for(e=d+1;e<NUMLINK;e++){

for(l=0;l<NUMGEN;l++){
tmp=Complex();

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
tmp=tmp+perm[d][e][c][a][b]*(
conjug(V.get(x-e_a-e_b-e_c,c).get(l,m))*BC(x,-e_a,-e_b)*
p.get(x-e_a-e_b,a,b).get(m)-
conjug(V.get(x-e_c,c).get(m,l))*BC(x,-e_a,-e_b,-e_c)*
p.get(x-e_a-e_b-e_c,a,b).get(m));}

}}
}

atmp.set(l,tmp);
}

dum.set(x,d,e,atmp);
dum.set(x,e,d,-1.0*atmp);
}}}

return(dum);
}



Plaq_Field Dbplus(const Adjoint_Links &V, const Plaq_Field &p){
Lattice_Vector x,e_a,e_b,e_c;
int a,b,c,d,e,l,m,sites;
Plaq_Field dum=Plaq_Field();
Complex tmp;
Afield atmp;

sites=0;
while(loop_over_lattice(x,sites)){
for(a=0;a<NUMLINK;a++){
e_a=Lattice_Vector(a);
for(b=a+1;b<NUMLINK;b++){
e_b=Lattice_Vector(b);

for(l=0;l<NUMGEN;l++){
tmp=Complex();

for(c=0;c<NUMLINK;c++){
if((c==a)||(c==b)) continue;

e_c=Lattice_Vector(c);

for(d=0;d<NUMLINK;d++){
if((d==c)||(d==a)||(d==b)) continue;
for(e=d+1;e<NUMLINK;e++){
if((e==c)||(e==a)||(e==b)) continue;

for(m=0;m<NUMGEN;m++){

tmp=tmp+perm[a][b][c][d][e]*(
conjug(V.get(x+e_a+e_b,c).get(l,m))*BC(x,e_a,e_b,e_c)*
p.get(x+e_a+e_b+e_c,d,e).get(m)-
conjug(V.get(x-e_c,c).get(m,l))*BC(x,e_a,e_b)*
p.get(x+e_a+e_b,d,e).get(m));}

}}
}
atmp.set(l,tmp);
}

dum.set(x,a,b,atmp);
dum.set(x,b,a,-1.0*atmp);
}}}

return(dum);
}


Link_Field Dbplus(const Adjoint_Links &V, const Site_Field &S){
Lattice_Vector x,e_mu;
int mu,sites,a,b;
Complex tmp;
Afield atmp;
Link_Field dum=Link_Field();


sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){

e_mu=Lattice_Vector(mu);

for(a=0;a<NUMGEN;a++){
tmp=Complex();
for(b=0;b<NUMGEN;b++){
tmp=tmp+
    conjug(V.get(x,mu).get(a,b))*S.get(x+e_mu).get(b)*BC(x,e_mu)-
    conjug(V.get(x,mu).get(b,a))*S.get(x).get(b);
    }
atmp.set(a,tmp);
}

dum.set(x,mu,atmp);

}
}

return(dum);
}

Site_Field Dbminus(const Adjoint_Links &V, const Link_Field &L){
Lattice_Vector x,e_mu;
int mu,sites,a,b;
Complex tmp;
Afield atmp;
Site_Field dum=Site_Field();

sites=0;
while(loop_over_lattice(x,sites)){

atmp=Afield();
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);

for(a=0;a<NUMGEN;a++){
tmp=Complex();
for(b=0;b<NUMGEN;b++){
tmp=tmp+
conjug(V.get(x,mu).get(a,b))*L.get(x,mu).get(b)-
conjug(V.get(x-e_mu,mu).get(b,a))*L.get(x-e_mu,mu).get(b)*BC(x,-e_mu);
}
atmp.set(a,atmp.get(a)+tmp);
}

}

dum.set(x,atmp);
}
return(dum);
}


Twist_Fermion betaterm(const Adjoint_Links &V, const Twist_Fermion &F){
  int a,b,sites,mu;
  Lattice_Vector x;
  Twist_Fermion dum=Twist_Fermion();
  Site_Field s=Site_Field();
  Link_Field l=Link_Field();
  Complex tmp;
  Afield tmp2;

  // only relevant if U(N)

  if(NUMGEN==(NCOLOR*NCOLOR-1)){return(dum);}
  //SIMON: beta term modification

  sites=0;
  while(loop_over_lattice(x,sites)){
    tmp2=Afield();

    for(mu=0;mu<NUMLINK;mu++){

      for(a=0;a<NUMGEN-1;a++){
        tmp=Complex();
        for(b=0;b<NUMGEN;b++){
          tmp=tmp+
            conjug(V.get(x,mu).get(a,b))*F.getL().get(x,mu).get(b);}

        tmp2.set(a,tmp);}

      s.set(x,C1*0.5*tmp2);}

    dum.setS(s);
  }

  sites=0;
  while(loop_over_lattice(x,sites)){
    for(mu=0;mu<NUMLINK;mu++){
      tmp2=Afield();
      for(b=0;b<NUMGEN;b++){
        tmp=Complex();
        for(a=0;a<NUMGEN-1;a++){
          tmp=tmp+
            conjug(V.get(x,mu).get(a,b))*F.getS().get(x).get(a);}
        tmp2.set(b,tmp);
      }

      l.set(x,mu,-0.5*C1*tmp2);
    }
  }
  dum.setL(l);
  return dum;
}


Twist_Fermion Fermion_op(const Adjoint_Links &V, const Twist_Fermion &F){
  Twist_Fermion F2=Twist_Fermion(),F3;

  F2.setC(Dplus(V,F.getL()));
  F2.setL(Dminus(V,F.getC()));
  F2.setL(F2.getL()+0.5*Dbplus(V,F.getS()));
  F2.setS(0.5*Dbminus(V,F.getL()));

  F3=betaterm(V,F);
  F2.setS(F2.getS()+F3.getS());
  F2.setL(F2.getL()+F3.getL());

  // Q-closed piece
  if(NUMLINK==5){
    F2.setC(F2.getC()+0.5*Dbminus(V,F.getC())+0.5*Dbplus(V,F.getC()));}

  return F2;
}

Twist_Fermion Adj_Fermion_op(const Adjoint_Links &V, const Twist_Fermion &F){
  Twist_Fermion F2=Twist_Fermion(),F3,F4;

  F3=Cjg(F);

  F2.setC(Dplus(V,F3.getL()));
  F2.setL(Dminus(V,F3.getC()));
  F2.setL(F2.getL()+0.5*Dbplus(V,F3.getS()));
  F2.setS(0.5*Dbminus(V,F3.getL()));

  F4=betaterm(V,F3);
  F2.setS(F2.getS()+F4.getS());
  F2.setL(F2.getL()+F4.getL());

  // Q-closed piece
  if(NUMLINK==5){
    F2.setC(F2.getC()+0.5*Dbminus(V,F3.getC())+0.5*Dbplus(V,F3.getC()));}

  F2=-1.0*Cjg(F2);

  return F2;
}


UPlaq_Field Field_Strength(const Gauge_Field &U){
int site,mu,nu;
Lattice_Vector x,e_mu,e_nu;
Umatrix Fmunu;
UPlaq_Field F=UPlaq_Field();

site=0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<NUMLINK;mu++){
e_mu=Lattice_Vector(mu);
for(nu=mu+1;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
Fmunu=U.get(x,mu)*U.get(x+e_mu,nu)-U.get(x,nu)*U.get(x+e_nu,mu);
F.set(x,mu,nu,Fmunu);
F.set(x,nu,mu,-1.0*Fmunu);
}
}
}
return(F);
}


UPlaq_Field Bianchi(const Gauge_Field &U){
UPlaq_Field F,B;
int a,b,c,d,e,site;
Lattice_Vector x,e_a,e_b,e_c,e_d,e_e;
Umatrix tmp;

F=Field_Strength(U);
B=UPlaq_Field();

site=0;
while(loop_over_lattice(x,site)){
for(a=0;a<NUMLINK;a++){
e_a=Lattice_Vector(a);
for(b=a+1;b<NUMLINK;b++){
e_b=Lattice_Vector(b);

tmp=Umatrix();
for(c=0;c<NUMLINK;c++){
if((c==a)||(c==b)) continue;

e_c=Lattice_Vector(c);

for(d=0;d<NUMLINK;d++){
if((d==c)||(d==a)||(d==b)) continue;
e_d=Lattice_Vector(d);
for(e=d+1;e<NUMLINK;e++){
if((e==c)||(e==a)||(e==b)) continue;
e_e=Lattice_Vector(e);

tmp=tmp+perm[a][b][c][d][e]*
(
U.get(x,c)*F.get(x+e_c,d,e)-
F.get(x,d,e)*U.get(x+e_d+e_e,c)
);

}}
}

//cout << "B is " << Tr(tmp*Adj(tmp)).real() << "\n" << flush;

B.set(x,a,b,tmp);
B.set(x,b,a,-1.0*tmp);
}}


}

return(B);
}



void eigsrt(double d[], int n)
{
        int k,j,i;
        double p;

        for (i=1;i<n;i++) {
                p=d[k=i];
                for (j=i+1;j<=n;j++)
                        if (d[j] >= p) p=d[k=j];
                if (k != i) {
                        d[k]=d[i];
                        d[i]=p;
//                      for (j=1;j<=n;j++) {
//                              p=v[j][i];
//                              v[j][i]=v[j][k];
//                              v[j][k]=p;
//                      }
                }
        }
}

void ceigsrt(Complex d[], int n)
{
        int k,j,i;
        Complex p;

        for (i=1;i<n;i++) {
                p=d[k=i];
                for (j=i+1;j<=n;j++)
                        if (d[j].norm() >= p.norm()) p=d[k=j];
                if (k != i) {
                        d[k]=d[i];
                        d[i]=p;
//                      for (j=1;j<=n;j++) {
//                              p=v[j][i];
//                              v[j][i]=v[j][k];
//                              v[j][k]=p;
//                      }
                }
        }
}

double gasdev(void){
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  if(iset==0){
    do{
      v1=2.0*rand()/(double)RAND_MAX-1.0;
      v2=2.0*rand()/(double)RAND_MAX-1.0;
      rsq=v1*v1+v2*v2;}
    while(rsq>=1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return(v2*fac);}
  else {
    iset=0;
    return(gset);}}

//hardcoded for U(2) for now

Complex det(const Umatrix &u){
        return(u.get(0,0)*u.get(1,1)-u.get(0,1)*u.get(1,0));
}

Umatrix adjugate(const Umatrix &u){
        Umatrix dum;
        dum.set(0,0,u.get(1,1));
        dum.set(0,1,-1.0*u.get(0,1));
        dum.set(1,0,-1.0*u.get(1,0));
        dum.set(1,1,u.get(0,0));
        return(dum);
}

Umatrix inverse(const Umatrix &u){
Umatrix dum;
        dum.set(0,0,u.get(1,1));
        dum.set(0,1,-1.0*u.get(0,1));
        dum.set(1,0,-1.0*u.get(1,0));
        dum.set(1,1,u.get(0,0));
        dum=(Complex(1.0,0.0)/det(u))*dum;
        return(dum);
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(double data[],int nn[],int ndim,int isign)
{
  int i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  int ibit,idim,k1,k2,n,nprev,nrem,ntot;
  double tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;

  ntot=1;
  for (idim=1;idim<=ndim;idim++)
    ntot *= nn[idim];
        nprev=1;
        for (idim=ndim;idim>=1;idim--) {
            n=nn[idim];
            nrem=ntot/(n*nprev);
            ip1=nprev << 1;
            ip2=ip1*n;
            ip3=ip2*nrem;
            i2rev=1;
            for (i2=1;i2<=ip2;i2+=ip1) {
                if (i2 < i2rev) {
                    for (i1=i2;i1<=i2+ip1-2;i1+=2) {
                        for (i3=i1;i3<=ip3;i3+=ip2) {
                            i3rev=i2rev+i3-i2;
                            SWAP(data[i3],data[i3rev]);
                            SWAP(data[i3+1],data[i3rev+1]);
                        }
                    }
                }
                ibit=ip2 >> 1;
                while (ibit >= ip1 && i2rev > ibit) {
                    i2rev -= ibit;
                    ibit >>= 1;
                }
                i2rev += ibit;
            }
            ifp1=ip1;
            while (ifp1 < ip2) {
                ifp2=ifp1 << 1;
                theta=isign*6.28318530717959/(ifp2/ip1);
                wtemp=sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi=sin(theta);
                wr=1.0;
                wi=0.0;
                for (i3=1;i3<=ifp1;i3+=ip1) {
                    for (i1=i3;i1<=i3+ip1-2;i1+=2) {
                        for (i2=i1;i2<=ip3;i2+=ifp2) {
                            k1=i2;
                            k2=k1+ifp1;
                            tempr=wr*data[k2]-wi*data[k2+1];
                            tempi=wr*data[k2+1]+wi*data[k2];
                            data[k2]=data[k1]-tempr;
                            data[k2+1]=data[k1+1]-tempi;
                            data[k1] += tempr;
                            data[k1+1] += tempi;
                        }
                    }
                    wr=(wtemp=wr)*wpr-wi*wpi+wr;
                    wi=wi*wpr+wtemp*wpi+wi;
                }
                ifp1=ifp2;
            }
            nprev *= n;
        }
}

#undef SWAP

