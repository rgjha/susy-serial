#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <fstream>
using namespace std;

#include <math.h>
#include <stdlib.h>
#include <iomanip>

#define PREC 10
#define FF 1.0

const int FERMIONS = 1;
const int T = 2;
const int NCOLOR = 2;
const int NUMGEN = (NCOLOR*NCOLOR);

#define Q16 

// Q=16 parameters
#ifdef Q16
const int D = 4;
const int NUMLINK = 5;
const int LX = 2;
const int LY = 2;
const int LZ = 2;
const int SITES = (LX*LY*LZ*T);
const unsigned int LEN = (16*NUMGEN*SITES);
#else
// Q=4 parameters
const int D = 2;
const int NUMLINK = 2;
const int LX = 2;
const int LY = 1;
const int LZ = 1;
const int SITES = (LX*T);
const int LEN = (4*NUMGEN*SITES);
#endif
 
const double GAUGETOL = 0.00000000000001;
const int DEGREE = 10;
const double PBC = 1; // set to 1 for periodic b.c.
const int SMALLEIG=0;
const int TWIST = 0;


// Omelyan
const double INT_LAMBDA = 0.193;
const double INT_LAMBDA_CONT = 0.386;
const double INT_LAMBDA_MID = 0.614;

#define FULLMATRIX


extern double ampdeg,amp[DEGREE],shift[DEGREE];
extern int SWEEPNO;
extern double KAPPA,DT,G,BMASS;
extern int SWEEPS,GAP,THERM,READIN,SEED, TRAJECTORY_LENGTH;
extern double SMALLCUT, LARGECUT;
extern double perm[NUMLINK][NUMLINK][NUMLINK][NUMLINK][NUMLINK];
extern int Lattice_Map[D];
extern int side[D];

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
    Umatrix(double d);
		Umatrix(Complex [NCOLOR][NCOLOR]);
		Complex get(int,int) const;
		void set(int,int,const Complex);
		void print(void);
		friend ostream& operator<<(ostream &, Umatrix);
		friend istream& operator>>(istream &, Umatrix &);

};
		
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
Umatrix traceless(const Umatrix &);

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
Umatrix twist(const Lattice_Vector &x, const Lattice_Vector &v);
Umatrix twist(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2);
Umatrix twist(const Lattice_Vector &x, const Lattice_Vector &v1, const Lattice_Vector &v2, const Lattice_Vector &v3);

class Gauge_Field{
private:
	Umatrix link[SITES][NUMLINK];
	
public:
    Gauge_Field(void);
	Gauge_Field(int);
	Umatrix get(const Lattice_Vector &, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const int) const;
  	void set(const Lattice_Vector &, const int, const Umatrix &);
	};
	
Gauge_Field Adj(const Gauge_Field &);
Gauge_Field Transpose(const Gauge_Field &);

class USite_Field{
private:
	Umatrix points[SITES];
public:
	USite_Field(void);
	USite_Field(int);
   	Umatrix get(const Lattice_Vector &) const;
        Umatrix get(const Lattice_Vector &, const Lattice_Vector &) const;
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
        Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const int, const int) const;
	void set(const Lattice_Vector &, const int, const int,  const Umatrix &);
	void print(void);
	};


UPlaq_Field Adj(const UPlaq_Field &);

UPlaq_Field operator +(const UPlaq_Field &, const UPlaq_Field &);
UPlaq_Field operator -(const UPlaq_Field &, const UPlaq_Field &);
UPlaq_Field operator *(const double, const UPlaq_Field &);
UPlaq_Field operator *(const Complex &, const UPlaq_Field &);
Umatrix operator *(const UPlaq_Field &, const UPlaq_Field &);


class Site_Field{
private:
	Umatrix points[SITES];
public:
	Site_Field(void);
	Site_Field(int);
	Umatrix get(const Lattice_Vector &) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &) const;
	void set(const Lattice_Vector &, const Umatrix &);
	void print(void);
	};

Site_Field Adj(const Site_Field &);

Site_Field operator +(const Site_Field &, const Site_Field &);
Site_Field operator -(const Site_Field &, const Site_Field &);
Site_Field operator *(const double, const Site_Field &);
Site_Field operator *(const Complex &, const Site_Field &);
Umatrix operator *(const Site_Field &, const Site_Field &);



class Link_Field{
private:
	Umatrix flink[SITES][NUMLINK];
public:
	Link_Field(void);
	Link_Field(int);
	Umatrix get(const Lattice_Vector &, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const int) const;
	void set(const Lattice_Vector &, const int, const Umatrix &);
	void print(void);
	};


Link_Field Adj(const Link_Field &);

Link_Field operator +(const Link_Field &, const Link_Field &);
Link_Field operator -(const Link_Field &, const Link_Field &);
Link_Field operator *(const double, const Link_Field &);
Link_Field operator *(const Complex &, const Link_Field &);
Umatrix operator *(const Link_Field &, const Link_Field &);


class Plaq_Field{
private:
	Umatrix square[SITES][NUMLINK][NUMLINK];
public:
	Plaq_Field(void);
	Plaq_Field(int);
	Umatrix get(const Lattice_Vector &, const int, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const int, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const int, const int) const;
    Umatrix get(const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &, const Lattice_Vector &,const int, const int) const;
	void set(const Lattice_Vector &, const int, const int,  const Umatrix &);
	void print(void);
	};


Plaq_Field Adj(const Plaq_Field &);

Plaq_Field operator +(const Plaq_Field &, const Plaq_Field &);
Plaq_Field operator -(const Plaq_Field &, const Plaq_Field &);
Plaq_Field operator *(const double, const Plaq_Field &);
Plaq_Field operator *(const Complex &, const Plaq_Field &);
Umatrix operator *(const Plaq_Field &, const Plaq_Field &);

UPlaq_Field Field_Strength(const Gauge_Field &U);
UPlaq_Field Bianchi(const Gauge_Field &U);

Link_Field Dbplus(const Gauge_Field &, const Site_Field &);
Plaq_Field Dplus(const Gauge_Field &, const Link_Field &);
Site_Field Dbminus(const Gauge_Field &, const Link_Field &);
Link_Field Dminus(const Gauge_Field &, const Plaq_Field &);
Plaq_Field Dbminus(const Gauge_Field &, const Plaq_Field &);
Plaq_Field Dbplus(const Gauge_Field &, const Plaq_Field &);

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


Twist_Fermion Adj(const Twist_Fermion &);

Twist_Fermion operator +(const Twist_Fermion &, const Twist_Fermion &);
Twist_Fermion operator -(const Twist_Fermion &, const Twist_Fermion &);
Twist_Fermion operator *(const double, const Twist_Fermion &);
Twist_Fermion operator *(const Complex &, const Twist_Fermion &);
Umatrix operator *(const Twist_Fermion &, const Twist_Fermion &);

Twist_Fermion Fermion_op(const Gauge_Field &U, const Twist_Fermion &);
Twist_Fermion Adj_Fermion_op(const Gauge_Field &U, const Twist_Fermion &);

void check_trace(Twist_Fermion &);

void eigsrt(double d[], int n);
void ceigsrt(Complex d[], int n);
void epsilon(void);
	
extern Umatrix Lambda[NUMGEN];

double gasdev(void);

void fourn(double data[], int nn[], int ndum, int isign);
#endif
