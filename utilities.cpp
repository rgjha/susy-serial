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

Umatrix::Umatrix(double d)
{
  for (int i=0;i<NCOLOR;i++)
  for (int j=0;j<NCOLOR;j++) {
    if (i == j) mat[i][j]=Complex(d,0.0);
    else mat[i][j]=Complex(0.0,0.0);
  }
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
        if(counter==10000){
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
	for(i=0;i<(NCOLOR*NCOLOR-1);i++){
	tmp=tmp+1.0/sqrt(2.0)*Complex(gasdev(),gasdev())*Lambda[i];}
	return(tmp);
	}

Umatrix gaussUtr(void){
   int i;
   Umatrix tmp=Umatrix();
   for(i=0;i<NCOLOR*NCOLOR;i++){
    tmp=tmp+1.0/sqrt(2.0)*(Complex(gasdev(),gasdev())*Lambda[i]);}
    return(tmp);
    }

 Umatrix traceless(const Umatrix &U){
 Umatrix dum;
 dum=U-(1.0/NCOLOR)*Tr(U)*Umatrix(1);
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

Umatrix twist(const Lattice_Vector &x, const Lattice_Vector &v){
    int coord,sum,i;
    Umatrix dum=Umatrix(1);
    coord=x.get(D-1)+v.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ dum=PBC*dum; }
    
    if(!TWIST) return dum;
 
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    if(sum<0) return (dum);
    if(sum>0) return (Adj(dum));

    return(dum);
}

Umatrix twist(const Lattice_Vector &x, const Lattice_Vector &v1 ,const Lattice_Vector &v2){
    int coord,sum,i;
    Umatrix dum=Umatrix(1);
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ dum=PBC*dum; }

    if(!TWIST) return dum;
    
    sum=0;
    for(i=0;i<(D-1);i++){
    coord=x.get(i)+v1.get(i)+v2.get(i);
    if ((coord >= side[i]) || (coord <0)) {
    dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));sum+=coord;}
    }
    if(sum<0) return (dum);
    if(sum>0) return (Adj(dum));

    return(dum);
}

Umatrix twist(const Lattice_Vector &x, const Lattice_Vector &v1 ,const Lattice_Vector &v2, const Lattice_Vector &v3){
    int coord,sum,i;
    Umatrix dum=Umatrix(1);
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1)+v3.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ dum=PBC*dum; }
    
    if(!TWIST) return dum;
    
    sum=0;
    for(i=0;i<(D-1);i++){
    coord=x.get(i)+v1.get(i)+v2.get(i)+v3.get(i);
    if ((coord >= side[i]) || (coord <0)) {
    dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));sum+=coord;}
    }
    
    if(sum<0) return (dum);
    if(sum>0) return (Adj(dum));
    return(dum);
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
link[i][j]=exp(0.1*gaussU());}
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



//J  Return U_\mu(x+v), including the effect of twisted BCs.
Umatrix Gauge_Field::get(const Lattice_Vector &x, const Lattice_Vector &v, const int mu) const{
  Lattice_Vector xp;
  int coord,sum;
  Umatrix dum;
  xp=x+v;
  int site=0,i;
  for(i=0;i<D;i++){
  site=site+xp.get(i)*Lattice_Map[i];}
  
  Umatrix Uloc=link[site][mu];
    if(!TWIST) return (Uloc);
    
  //J  Not twisted in the temporal direction.
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
    coord=x.get(i)+v.get(i);
    if ((coord >= side[i]) || (coord <0)) {
    dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
    sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum); 
    return(Uloc);
}

Umatrix Gauge_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2, const int mu) const{
  Lattice_Vector xp;
  int coord,sum;
  Umatrix dum;
  xp=x+v1+v2;
  int site=0,i;
  for(i=0;i<D;i++){
  site=site+xp.get(i)*Lattice_Map[i];}
  
  Umatrix Uloc=link[site][mu];
    if(!TWIST) return (Uloc);
  
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
    coord=x.get(i)+v1.get(i)+v2.get(i);
    if ((coord >= side[i]) || (coord <0)) {
    dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
    sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);  
  return(Uloc);
}


Umatrix Gauge_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2, const Lattice_Vector &v3, const int mu) const{
  Lattice_Vector xp;
  int coord,sum;
  Umatrix dum;

  xp=x+v1+v2+v3;
  int site=0,i;
  for(i=0;i<D;i++){
  site=site+xp.get(i)*Lattice_Map[i];}
  
  Umatrix Uloc=link[site][mu];
    if(!TWIST) return (Uloc);

    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
    coord=x.get(i)+v1.get(i)+v2.get(i)+v3.get(i);
    if ((coord >= side[i]) || (coord <0)) {
    dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
    sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);  
  return(Uloc);
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

Umatrix USite_Field::get(const Lattice_Vector &x, const Lattice_Vector &v) const{
  Lattice_Vector xp;
  int coord,sum;
  Umatrix dum;
  xp=x+v;
  int site=0,i;
  for(i=0;i<D;i++){
  site=site+xp.get(i)*Lattice_Map[i];}

  Umatrix Uloc=points[site];
    if(!TWIST) return (Uloc);

  //J  Not twisted in the temporal direction.
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
    coord=x.get(i)+v.get(i);
    if ((coord >= side[i]) || (coord <0)) {
    dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
    sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);
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

Umatrix UPlaq_Field::get(const Lattice_Vector &x, 
const Lattice_Vector &v, const int mu, const int nu) const{
  Lattice_Vector xp;
  int coord,sum;
  Umatrix dum;
  xp=x+v;
  int site=0,i;
  for(i=0;i<D;i++){
  site=site+xp.get(i)*Lattice_Map[i];}

  Umatrix Uloc=square[site][mu][nu];
  if(!TWIST) return (Uloc);

    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
    coord=x.get(i)+v.get(i);
    if ((coord >= side[i]) || (coord <0)) {
    dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
    sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);
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
points[i]=Umatrix();}
return;
}

Site_Field::Site_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
points[i]=gaussUtr();
}
return;}

if(c==0){
for(int i=0;i<SITES;i++){
points[i]=Umatrix(1);}

return;}

cout << "error in site constructor\n" << flush;

}



Umatrix Site_Field::get(const Lattice_Vector &x) const{
int site=0,i;

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(points[site]);
}

Umatrix Site_Field::get(const Lattice_Vector &x, const Lattice_Vector &v) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v;
    int site=0,i;
    for(i=0;i<D;i++){
        site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=points[site];
    coord=x.get(D-1)+v.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc; }
    
    if(!TWIST) return (Uloc);
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);
    }


Umatrix Site_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v1+v2;
    int site=0,i;
    for(i=0;i<D;i++){
        site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=points[site];
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc; }
    
    if(!TWIST) return (Uloc);
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v1.get(i)+v2.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);

    }


Umatrix Site_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2, const Lattice_Vector &v3) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v1+v2+v3;
    int site=0,i;
    for(i=0;i<D;i++){
        site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=points[site];
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1)+v3.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc; }
    
    if(!TWIST) return (Uloc);
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v1.get(i)+v2.get(i)+v3.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;
        }
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);
    
}

void Site_Field::set(const Lattice_Vector &x, const Umatrix &u){
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

Umatrix operator *(const Site_Field &s1, const Site_Field &s2){
int sites=0;
Lattice_Vector x;
Umatrix dum;

dum=Umatrix();
while(loop_over_lattice(x,sites)){
dum=dum+s1.get(x)*s2.get(x);
}
return(dum);
}

Site_Field Adj(const Site_Field &l){
int sites;
Lattice_Vector x;
Site_Field dum;

sites=0;
while(loop_over_lattice(x,sites)){
dum.set(x,Adj(l.get(x)));}

return(dum);
}

Link_Field::Link_Field(void){
for(int i=0;i<SITES;i++)
for(int j=0;j<NUMLINK;j++){
flink[i][j]=Umatrix();}
}

Link_Field::Link_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
for(int j=0;j<NUMLINK;j++){
flink[i][j]=gaussUtr();}}
return;
}
if(c==0){
for(int i=0;i<SITES;i++){
for(int j=0;j<NUMLINK;j++){
flink[i][j]=Umatrix(1);
}
}
return;
}
    
cout << "error in Link field constructor " << "\n" << flush;
    
return;
}


Umatrix Link_Field::get(const Lattice_Vector &x, const int mu) const{
    int site=0,i;
    
    for(i=0;i<D;i++)
    site=site+x.get(i)*Lattice_Map[i];
    
    return(flink[site][mu]);
}



Umatrix Link_Field::get(const Lattice_Vector &x, const Lattice_Vector &v, const int mu) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v;
    int site=0,i;
    for(i=0;i<D;i++){
    site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=flink[site][mu];
    coord=x.get(D-1)+v.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc;}
    
    if(!TWIST) return (Uloc);
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);

}

Umatrix Link_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2, const int mu) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v1+v2;
    int site=0,i;
    for(i=0;i<D;i++){
    site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=flink[site][mu];
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc;}
    
    if(!TWIST) return (Uloc);
    
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v1.get(i)+v2.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);

    }


Umatrix Link_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2, const Lattice_Vector &v3, const int mu) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v1+v2+v3;
    int site=0,i;
    for(i=0;i<D;i++){
    site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=flink[site][mu];
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1)+v3.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc;}
    
    if(!TWIST) return (Uloc);
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v1.get(i)+v2.get(i)+v3.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);

}


void Link_Field::set(const Lattice_Vector &x, const int mu, const Umatrix &u){
    int site=0,i;
    
    for(i=0;i<D;i++)
    site=site+x.get(i)*Lattice_Map[i];
    
    flink[site][mu]=u;
    return;
}

void Link_Field::print(void){
cout << "link field values\n" << flush;
for(int i=0;i<SITES;i++){
for(int j=0;j<NUMLINK;j++){
cout << "site= " << i << "\n" << flush;
cout << flink[i][j] << "\n" << flush;}}
return;
}

Link_Field operator *(const double o, const Link_Field &s){
int sites=0;
Lattice_Vector x;
Link_Field dum=Link_Field();

while(loop_over_lattice(x,sites)){
for(int mu=0;mu<NUMLINK;mu++){
dum.set(x,mu,o*s.get(x,mu));
}}

return(dum);
}


Link_Field operator *(const Complex &o, const Link_Field &s){
int sites=0;
Lattice_Vector x;
Link_Field dum=Link_Field();

while(loop_over_lattice(x,sites)){
for(int mu=0;mu<NUMLINK;mu++){
dum.set(x,mu,o*s.get(x,mu));
}}

return(dum);
}

Umatrix operator *(const Link_Field &s1, const Link_Field &s2){
int sites=0;
Lattice_Vector x;
Umatrix dum;

dum=Umatrix();
while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum=dum+s1.get(x,j)*s2.get(x,j);
}}
return(dum);
}

Link_Field operator +(const Link_Field &s1, const Link_Field &s2){
int sites=0;
Lattice_Vector x;
Link_Field dum;
while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum.set(x,j,s1.get(x,j)+s2.get(x,j));
}}
return(dum);
}


Link_Field operator -(const Link_Field &s1, const Link_Field &s2){
int sites=0;
Lattice_Vector x;
Link_Field dum;
while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
dum.set(x,j,s1.get(x,j)-s2.get(x,j));
}}
return(dum);
}

Link_Field Adj(const Link_Field & U){
int site;
Link_Field W;
Lattice_Vector x;
site=0;
while(loop_over_lattice(x,site)){
for(int mu=0;mu<NUMLINK;mu++){
W.set(x,mu,Adj(U.get(x,mu)));
}}
return W;
}



Plaq_Field::Plaq_Field(void){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
for(int nu=0;nu<NUMLINK;nu++){
square[i][mu][nu]=Umatrix();}}}
return;
}

Plaq_Field::Plaq_Field(int c){
if(c==1){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
square[i][mu][mu]=Umatrix();
for(int nu=mu+1;nu<NUMLINK;nu++){
square[i][mu][nu]=gaussUtr();
square[i][nu][mu]=-1.0*square[i][mu][nu];
}
}
}
return;
}
    
if(c==0){
for(int i=0;i<SITES;i++){
for(int mu=0;mu<NUMLINK;mu++){
square[i][mu][mu]=Umatrix();
for(int nu=mu+1;nu<NUMLINK;nu++){
square[i][mu][nu]=Umatrix(1);
square[i][nu][mu]=-1.*Umatrix(1);
}
}
}

return;
}
cout << "error in plaq constructor\n" << "\n" << flush;

return;
}

Umatrix Plaq_Field::get(const Lattice_Vector &x, const int mu, const int nu) const{
int site=0,i;

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(square[site][mu][nu]);
}

Umatrix Plaq_Field::get(const Lattice_Vector &x, const Lattice_Vector &v, const int mu, const int nu) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v;
    int site=0,i;
    for(i=0;i<D;i++){
    site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=square[site][mu][nu];
    coord=x.get(D-1)+v.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc;}
    
    if(!TWIST) return (Uloc);
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    //dum=Adj(dum);
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);

    }


Umatrix Plaq_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1,
const Lattice_Vector &v2, const int mu, const int nu) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v1+v2;
    int site=0,i;
    for(i=0;i<D;i++){
    site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=square[site][mu][nu];
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc;}
    
    if(!TWIST) return (Uloc);
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v1.get(i)+v2.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    //dum=Adj(dum);
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);
    
}



Umatrix Plaq_Field::get(const Lattice_Vector &x, const Lattice_Vector &v1,
const Lattice_Vector &v2, const Lattice_Vector &v3, const int mu, const int nu) const{
    Lattice_Vector xp;
    int coord,sum;
    Umatrix dum;
    xp=x+v1+v2+v3;
    int site=0,i;
    for(i=0;i<D;i++){
    site=site+xp.get(i)*Lattice_Map[i];}
    
    Umatrix Uloc=square[site][mu][nu];
    coord=x.get(D-1)+v1.get(D-1)+v2.get(D-1)+v3.get(D-1);
    if((coord >= side[D-1]) || (coord <0)){ Uloc=PBC*Uloc;}
    
    if(!TWIST) return (Uloc);
    
    dum=Umatrix(1);
    sum=0;
    for(i=0;i<(D-1);i++){
        coord=x.get(i)+v1.get(i)+v2.get(i)+v3.get(i);
        if ((coord >= side[i]) || (coord <0)) {
            dum=dum*Lambda[i]*Complex(0.0,sqrt(2.0));
            sum+=coord;}
    }
    //dum=Adj(dum);
    if(sum<0) return (dum*Uloc*Adj(dum));
    if(sum>0) return (Adj(dum)*Uloc*dum);
    return(Uloc);

    }

void Plaq_Field::set(const Lattice_Vector &x, const int mu, const int nu, const Umatrix &u){
int site=0,i;


for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

square[site][mu][nu]=u;
return;
}

void Plaq_Field::print(void){
cout << "link field values\n" << flush;
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

Umatrix operator *(const Plaq_Field &s1, const Plaq_Field &s2){
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

Plaq_Field Adj(const Plaq_Field &l){
int sites=0;
Lattice_Vector x;
Plaq_Field dum;

while(loop_over_lattice(x,sites)){
for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum.set(x,j,k,Adj(l.get(x,j,k)));}
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

Twist_Fermion Adj(const Twist_Fermion &K){
Twist_Fermion dum;

dum.setS(Adj(K.getS()));
dum.setL(Adj(K.getL()));
dum.setC(Adj(K.getC()));

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

void check_trace(Twist_Fermion &sol){
Site_Field s;
Link_Field l;
Plaq_Field p;
Umatrix dum;
Lattice_Vector x;
int sites;

s=sol.getS();
l=sol.getL();
p=sol.getC();

sites=0;
while(loop_over_lattice(x,sites)){
dum=s.get(x);
if(Tr(dum).norm()>0.00001){cout << "site field not tracelss" << endl;}

for(int j=0;j<NUMLINK;j++){
dum=l.get(x,j);
if(Tr(dum).norm()>0.00001){cout << "link field not tracelss" << endl;}
}

for(int j=0;j<NUMLINK;j++){
for(int k=0;k<NUMLINK;k++){
dum=p.get(x,j,k);
if(Tr(dum).norm()>0.00001){cout << "plaquette field not tracelss" << endl;}
}}
}


return;
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

Umatrix operator *(const Twist_Fermion &k1, const Twist_Fermion &k2){
Umatrix tmp;

tmp=k1.getS()*k2.getS()+k1.getL()*k2.getL()+k1.getC()*k2.getC();
return(tmp);
}



Plaq_Field Dplus(const Gauge_Field &U, const Link_Field &L){
Lattice_Vector x,e_mu,e_nu;
int mu,nu,sites;
Umatrix tmp;
Plaq_Field dum=Plaq_Field();


sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<NUMLINK;mu++){
for(nu=mu+1;nu<NUMLINK;nu++){

e_mu=Lattice_Vector(mu);
e_nu=Lattice_Vector(nu);

tmp=Umatrix();
tmp=tmp+
    U.get(x,mu)*L.get(x,e_mu,nu)-
    L.get(x,nu)*U.get(x,e_nu,mu);
tmp=tmp-
    U.get(x,nu)*L.get(x,e_nu,mu)+
    L.get(x,mu)*U.get(x,e_mu,nu);

dum.set(x,mu,nu, (tmp));
dum.set(x,nu,mu,-1.0*(tmp));
}}
}

return(dum);
}

Link_Field Dminus(const Gauge_Field &U, const Plaq_Field &P){
Lattice_Vector x,e_mu,e_nu;
int sites,mu,nu;
Umatrix tmp;
Link_Field dum=Link_Field();

sites=0;
while(loop_over_lattice(x,sites)){
for(nu=0;nu<NUMLINK;nu++){
e_nu=Lattice_Vector(nu);
tmp=Umatrix();

for(mu=0;mu<NUMLINK;mu++){
if(mu==nu) continue;
e_mu=Lattice_Vector(mu);

tmp=tmp+
U.get(x,e_nu,mu)*P.get(x,mu,nu)-
P.get(x,-e_mu,mu,nu)*U.get(x,-e_mu,mu);
}

dum.set(x,nu,(tmp));
}}
return(dum);
}


Link_Field Dbplus(const Gauge_Field &U, const Site_Field &S){
    Lattice_Vector x,e_mu;
    int mu,sites;
    Umatrix tmp;
    Link_Field dum=Link_Field();
    Gauge_Field Udag;
    
    Udag=Adj(U);
    
    sites=0;
    while(loop_over_lattice(x,sites)){
    for(mu=0;mu<NUMLINK;mu++){
    e_mu=Lattice_Vector(mu);
    tmp=S.get(x,e_mu)*Udag.get(x,mu)-
        Udag.get(x,mu)*S.get(x);
            
        dum.set(x,mu,(tmp));
        }
    }
    return(dum);
}

Site_Field Dbminus(const Gauge_Field &U, const Link_Field &L){
    Lattice_Vector x,e_mu;
    int mu,sites;
    Umatrix tmp;
    Site_Field dum=Site_Field();
    Gauge_Field Udag;
    
    Udag=Adj(U);
    sites=0;
    while(loop_over_lattice(x,sites)){
        
        tmp=Umatrix();
        for(mu=0;mu<NUMLINK;mu++){
        e_mu=Lattice_Vector(mu);
        tmp=tmp+
        L.get(x,mu)*Udag.get(x,mu)-
        Udag.get(x,-e_mu,mu)*L.get(x,-e_mu,mu);
        }
        dum.set(x,(tmp));
    }
    return(dum);
}

Plaq_Field Dbminus(const Gauge_Field &U, const Plaq_Field &p){
    Lattice_Vector x,e_a,e_b,e_c;
    int a,b,c,d,e,sites;
    Plaq_Field dum=Plaq_Field();
    Umatrix tmp;
    Gauge_Field Udag;
    Udag=Adj(U);

    sites=0;
    while(loop_over_lattice(x,sites)){
    for(d=0;d<NUMLINK;d++){
    for(e=d+1;e<NUMLINK;e++){
    tmp=Umatrix();
                    
    for(c=0;c<NUMLINK;c++){
    if((c==d)||(c==e)) continue;
    e_c=Lattice_Vector(c);
                        
    for(a=0;a<NUMLINK;a++){
    if((a==c)||(a==d)||(a==e)) continue;
    e_a=Lattice_Vector(a);
    for(b=a+1;b<NUMLINK;b++){
    if((b==c)||(b==d)||(b==e)) continue;
    e_b=Lattice_Vector(b);
    tmp=tmp+perm[d][e][c][a][b]*(
    p.get(x,-e_a,-e_b,a,b)*Udag.get(x,-e_a,-e_b,-e_c,c)-
    Udag.get(x,-e_c,c)*p.get(x,-e_a,-e_b,-e_c,a,b)
    );
    }}}
        dum.set(x,d,e,(tmp));
        dum.set(x,e,d,-1.0*(tmp));
    }}}
    
    return(dum);
}



Plaq_Field Dbplus(const Gauge_Field &U, const Plaq_Field &p){
    Lattice_Vector x,e_a,e_b,e_c;
    int a,b,c,d,e,sites;
    Plaq_Field dum=Plaq_Field();
    Umatrix tmp;
    Gauge_Field Udag;
    Udag=Adj(U);
    
    sites=0;
    while(loop_over_lattice(x,sites)){
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
    for(e=d+1;e<NUMLINK;e++){
    if((e==c)||(e==a)||(e==b)) continue;
                                    
    tmp=tmp+perm[a][b][c][d][e]*(
    p.get(x,e_a,e_b,e_c,d,e)*Udag.get(x,e_a,e_b,c)-
    Udag.get(x,-e_c,c)*p.get(x,e_a,e_b,d,e));
                                
    }}}
        dum.set(x,a,b,(tmp));
        dum.set(x,b,a,-1.0*(tmp));
    }}}
    
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


Twist_Fermion Fermion_op(const Gauge_Field &U,
const Twist_Fermion &F){
Twist_Fermion F2=Twist_Fermion();

F2.setC(Dplus(U,F.getL()));
F2.setL(Dminus(U,F.getC()));
F2.setL(F2.getL()+0.5*Dbplus(U,F.getS()));
F2.setS(0.5*Dbminus(U,F.getL()));

// Q-closed piece
if(NUMLINK==5){
F2.setC(F2.getC()+0.5*Dbminus(U,F.getC())+0.5*Dbplus(U,F.getC()));}

return(F2);
}

Twist_Fermion Adj_Fermion_op(const Gauge_Field &U,
const Twist_Fermion &F){
Twist_Fermion F2=Twist_Fermion(),F3;

F3=Adj(F);

F2.setC(Dplus(U,F3.getL()));
F2.setL(Dminus(U,F3.getC()));
F2.setL(F2.getL()+0.5*Dbplus(U,F3.getS()));
F2.setS(0.5*Dbminus(U,F3.getL()));

// Q-closed piece
if(NUMLINK==5){
F2.setC(F2.getC()+0.5*Dbminus(U,F3.getC())+0.5*Dbplus(U,F3.getC()));}

F2=-1.0*Adj(F2);

return(F2);
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
Fmunu=U.get(x,mu)*U.get(x,e_mu,nu)-U.get(x,nu)*U.get(x,e_nu,mu);
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
U.get(x,c)*F.get(x,e_c,d,e)-
F.get(x,d,e)*U.get(x,e_d,e_e,c)
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

