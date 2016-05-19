#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <stdlib.h>
#include <iomanip>

#define PREC 10
#define GO 1

#define FMASS 0.0

const int FERMIONS = 1;
const int T = 8;
const int NCOLOR = 3;
const int NUMGEN = (NCOLOR*NCOLOR);
//const int NUMGEN = (NCOLOR*NCOLOR)-1;  // Forcing SU(N) theory

#define GPU
#define Q16 

// Q=16 parameters
#ifdef Q16
const int D = 4;
const int NUMLINK = 5;
const int LX = 1;
const int LY = 1;
const int LZ = 32;
const int SITES = (LX*LY*LZ*T);
const unsigned int LEN = (16*NUMGEN*SITES);       // This is the size of fermion matrix // 
const int supercharges =16;
#else
// Q=4 parameters
const int D = 2;
const int NUMLINK = 2;
const int LX = 2;
const int LY = 1;
const int LZ = 1;
const int SITES = (LX*T);
const int LEN = (4*NUMGEN*SITES);
const int supercharges = 4;
#endif

//const int NONZEROES = LEN/SITES/2;
const unsigned int NONZEROES = (NUMGEN*NUMLINK*8);
 
const double GAUGETOL = 0.00000000000001;
const int DEGREE = 10;      // Number of terms in  Remez approximation
const double NORM = (1.0/sqrt(4.0*D));
const double PBC = -1; // Set to 1 for periodic b.c.
const int SMALLEIG=1;  // Set to 1 for constructing polar projected 'U'



// Omelyan
const double INT_LAMBDA = 0.193;
const double INT_LAMBDA_CONT = 0.386;
const double INT_LAMBDA_MID = 0.614;

//#define FULLMATRIX


extern double ampdeg,amp[DEGREE],shift[DEGREE];
extern int num_in_row[LEN],SIMULATING,SWEEPNO,TOTALNONZEROES;
extern double KAPPA,DT,ALPHA,G,BMASS,C2,MASS;
extern int SWEEPS,GAP,THERM,READIN,SEED, TRAJECTORY_LENGTH,OLDSWEEPNO;
extern double SMALLCUT, LARGECUT;
extern double TIME;
extern double perm[NUMLINK][NUMLINK][NUMLINK][NUMLINK][NUMLINK];
extern int Lattice_Map[D];
extern int side[D], BLOCK_MEASURE,BLOCKING;

class Complex{
	private:
		double re,im;
	public:
		Complex();
		Complex(double, double);
		double real(void) const;
		double imag(void) const;
		double norm(void);
		void print(void) const;
		friend ostream& operator<<(ostream&,Complex);
		friend istream& operator>>(istream&,Complex &);};

inline Complex conjug(const Complex &o1){return(Complex(o1.real(),-o1.imag()));}
inline Complex operator +(const Complex &o1, const Complex &o2){
	return(Complex(o1.real()+o2.real(),o1.imag()+o2.imag()));}
inline Complex operator -(const Complex &o1, const Complex &o2){
	return(Complex(o1.real()-o2.real(),o1.imag()-o2.imag()));}
inline Complex operator *(const Complex &o1, const Complex &o2){
	return(Complex(o1.real()*o2.real()-o1.imag()*o2.imag(),
		o1.real()*o2.imag()+o1.imag()*o2.real()));}
inline Complex operator *(const Complex &o1, const double o2){
	return(Complex(o1.real()*o2,o1.imag()*o2));}
inline Complex operator *(const double o1, const Complex &o2){
	return(Complex(o2.real()*o1,o2.imag()*o1));}

	
Complex operator /(const Complex &, const Complex &);
Complex pow(const Complex &, const int);

class Umatrix{
	private:
		Complex mat[NCOLOR][NCOLOR];
	public:
		Umatrix();
		Umatrix(int);
		Umatrix(Complex [NCOLOR][NCOLOR]);
		Complex get(int,int) const;
		void set(int,int,const Complex);
		void print(void);
		friend ostream& operator<<(ostream &, Umatrix);
		friend istream& operator>>(istream &, Umatrix &);};
		
Umatrix operator +(const Umatrix &o1, const Umatrix &o2);
Umatrix operator -(const Umatrix &o1, const Umatrix &o2);
Umatrix operator *(const Umatrix &, const Umatrix &);
Umatrix operator *(const Umatrix &, const Complex &);
Umatrix operator *(const Complex &, const Umatrix &);
Umatrix operator *(const Umatrix &, const double);
Umatrix operator *(const double, const Umatrix &);
Umatrix gaussU(void);
Umatrix exp(const Umatrix &);
Umatrix Adj(const Umatrix &);
Complex Tr(const Umatrix &);
Umatrix Trans(const Umatrix &);

class Lattice_Vector{
private:
	int coords[D];
public:
	Lattice_Vector(void);
	Lattice_Vector(int);
	void set(int, int);
	int get(int) const;
	void print(void) const;
	};

Lattice_Vector operator +(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x, const Lattice_Vector &y);
Lattice_Vector operator -(const Lattice_Vector &x);
double BC(const Lattice_Vector &x, const Lattice_Vector &y);
double BC(const Lattice_Vector &x, const Lattice_Vector &y, const Lattice_Vector
&z);
double BC(const Lattice_Vector &x, const Lattice_Vector &y, const Lattice_Vector
&z, const Lattice_Vector &w);
int loop_over_lattice(Lattice_Vector &, int &);
int blocklatsite(Lattice_Vector &);

int length(const Lattice_Vector &x);


class Gauge_Field{
private:
	Umatrix link[SITES][NUMLINK];
	
public:
        Gauge_Field(void);
	Gauge_Field(int);
	Umatrix get(const Lattice_Vector &, const int) const;
	void set(const Lattice_Vector &, const int, const Umatrix &);
	};
	
Gauge_Field Adj(const Gauge_Field &);
Gauge_Field Transpose(const Gauge_Field &);

class Afield{
	private:
		Complex afield[NUMGEN];
	public:
		Afield(void);
		Afield(int);
		Afield(Complex [NUMGEN]);
		Complex get(int) const;
		void set(int,const Complex);
		void print(void);
		friend ostream& operator<<(ostream &, Afield);
		friend istream& operator>>(istream &, Afield &);};
		
Afield operator +(const Afield &o1, const Afield &o2);
Afield operator -(const Afield &o1, const Afield &o2);
Afield operator *(const Afield &, const Complex &);
Afield operator *(const Complex &, const Afield &);
Afield operator *(const Afield &, const double);
Afield operator *(const double, const Afield &);
Complex operator *(const Afield &, const Afield &);
Afield Cjg(const Afield &u);
Afield gaussA(void);

class Scalar_Plaquette{
private:
	Complex data[SITES][NUMLINK][NUMLINK];
public:
	Scalar_Plaquette(void);
   	Complex get(const Lattice_Vector &, const int, const int) const;
	void set(const Lattice_Vector &, const int, const int, const Complex &);
	};


Scalar_Plaquette operator +(const Scalar_Plaquette &, const Scalar_Plaquette &);
Scalar_Plaquette operator -(const Scalar_Plaquette &, const Scalar_Plaquette &);
Scalar_Plaquette operator *(const double, const Scalar_Plaquette &);
//Scalar_Plaquette operator *(const Complex &, const Scalar_Plaquette &);
Scalar_Plaquette mydiff(const Scalar_Plaquette &, const Scalar_Plaquette &);
Scalar_Plaquette mydiff2(const Scalar_Plaquette &, const Scalar_Plaquette &);

class USite_Field{
private:
	Umatrix points[SITES];
public:
	USite_Field(void);
	USite_Field(int);
   	Umatrix get(const Lattice_Vector &) const;
	void set(const Lattice_Vector &, const Umatrix &);
	void print(void);
	};

USite_Field Adj(const USite_Field &);	

USite_Field operator +(const USite_Field &, const USite_Field &);
USite_Field operator -(const USite_Field &, const USite_Field &);
USite_Field operator *(const double, const USite_Field &);
USite_Field operator *(const Complex &, const USite_Field &);
Umatrix operator *(const USite_Field &, const USite_Field &);

class UPlaq_Field{
private:
	Umatrix square[SITES][NUMLINK][NUMLINK];
public:
	UPlaq_Field(void);
	UPlaq_Field(int);
	Umatrix get(const Lattice_Vector &, const int, const int) const;
	void set(const Lattice_Vector &, const int, const int,  const Umatrix &);
	void print(void);
	};


UPlaq_Field Adj(const UPlaq_Field &);

UPlaq_Field operator +(const UPlaq_Field &, const UPlaq_Field &);
UPlaq_Field operator -(const UPlaq_Field &, const UPlaq_Field &);
UPlaq_Field operator *(const double, const UPlaq_Field &);
UPlaq_Field operator *(const Complex &, const UPlaq_Field &);
Umatrix operator *(const UPlaq_Field &, const UPlaq_Field &);
UPlaq_Field Plaq(const Gauge_Field &);

class Site_Field{
private:
	Afield points[SITES];
public:
	Site_Field(void);
	Site_Field(int);
	Afield get(const Lattice_Vector &) const;
	void set(const Lattice_Vector &, const Afield &);
	void print(void);
	};

Site_Field Cjg(const Site_Field &);	

Site_Field operator +(const Site_Field &, const Site_Field &);
Site_Field operator -(const Site_Field &, const Site_Field &);
Site_Field operator *(const double, const Site_Field &);
Site_Field operator *(const Complex &, const Site_Field &);
Complex operator *(const Site_Field &, const Site_Field &);


class Link_Field{
private:
	Afield links[SITES][NUMLINK];
public:
	Link_Field(void);
	Link_Field(int);
	Afield get(const Lattice_Vector &, const int) const;
	void set(const Lattice_Vector &, const int, const Afield &);
	void print(void);
	};


Link_Field Cjg(const Link_Field &);

Link_Field operator +(const Link_Field &, const Link_Field &);
Link_Field operator -(const Link_Field &, const Link_Field &);
Link_Field operator *(const double, const Link_Field &);
Link_Field operator *(const Complex &, const Link_Field &);
Complex operator *(const Link_Field &, const Link_Field &);


class Plaq_Field{
private:
	Afield square[SITES][NUMLINK][NUMLINK];
public:
	Plaq_Field(void);
	Plaq_Field(int);
	Afield get(const Lattice_Vector &, const int, const int) const;
	void set(const Lattice_Vector &, const int, const int,  const Afield &);
	void print(void);
	};


Plaq_Field Cjg(const Plaq_Field &);

Plaq_Field operator +(const Plaq_Field &, const Plaq_Field &);
Plaq_Field operator -(const Plaq_Field &, const Plaq_Field &);
Plaq_Field operator *(const double, const Plaq_Field &);
Plaq_Field operator *(const Complex &, const Plaq_Field &);
Complex operator *(const Plaq_Field &, const Plaq_Field &);


class Adjoint_Matrix{
private:
       Complex amat[NUMGEN][NUMGEN];
public:
       Adjoint_Matrix(void);
       Adjoint_Matrix(Complex m[NUMGEN][NUMGEN]);
       Complex get(const int, const int) const;
       void set(const int, const int, const Complex &);
       friend ostream& operator<<(ostream &, Adjoint_Matrix);
       friend istream& operator>>(istream &, Adjoint_Matrix &);
       };
       

class Adjoint_Links{
private:
       Adjoint_Matrix alinks[SITES][NUMLINK];
       
public:
       Adjoint_Links(void);
       Adjoint_Matrix get(const Lattice_Vector &, const int) const;
       void set(const Lattice_Vector &, const int, const Adjoint_Matrix &);
       void print(void);
       };

void compute_Adjoint_Links(const Gauge_Field &U, Adjoint_Links &V);

Link_Field Dbplus(const Adjoint_Links &, const Site_Field &);
Plaq_Field Dplus(const Adjoint_Links &, const Link_Field &);
Site_Field Dbminus(const Adjoint_Links &, const Link_Field &);
Link_Field Dminus(const Adjoint_Links &, const Plaq_Field &);
Plaq_Field Dbminus(const Adjoint_Links &, const Plaq_Field &);
Plaq_Field Dbplus(const Adjoint_Links &, const Plaq_Field &);

class Twist_Fermion{
private:
	Site_Field S;
	Link_Field L;
	Plaq_Field C;
	
public:
	Twist_Fermion(void);
	Twist_Fermion(int);
	const Site_Field& getS(void) const;
	const Link_Field& getL(void) const;
	const Plaq_Field& getC(void) const;
	void setS(const Site_Field &);
	void setL(const Link_Field &);
	void setC(const Plaq_Field &);
	void print(void);
	};


Twist_Fermion Cjg(const Twist_Fermion &);

Twist_Fermion operator +(const Twist_Fermion &, const Twist_Fermion &);
Twist_Fermion operator -(const Twist_Fermion &, const Twist_Fermion &);
Twist_Fermion operator *(const double, const Twist_Fermion &);
Twist_Fermion operator *(const Complex &, const Twist_Fermion &);
Complex operator *(const Twist_Fermion &, const Twist_Fermion &);

Twist_Fermion Fermion_op(const Adjoint_Links &, const Gauge_Field &,const Twist_Fermion &);
Twist_Fermion Adj_Fermion_op(const Adjoint_Links &, const Gauge_Field &,
const Twist_Fermion &);
UPlaq_Field Field_Strength(const Gauge_Field &U);
UPlaq_Field Bianchi(const Gauge_Field &U);

void eigsrt(double d[], int n);
void ceigsrt(Complex d[], int n);
void epsilon(void);
	
extern Umatrix Lambda[NUMGEN];

double gasdev(void);

Complex det(const Umatrix &u);

Umatrix adjugate(const Umatrix &u);
Umatrix inverse(const Umatrix &);
void fourn(double data[], int nn[], int ndum, int isign);
#endif
