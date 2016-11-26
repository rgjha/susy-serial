#include "setup.h"

void setup(void){
int i,j,k;

// note: Tr(TaTb)=-delta_ab

cout << "Computing generators for SU(N)\n" << flush;
if(NUMGEN>1){(void)my_gen();}


if(NUMGEN==NCOLOR*NCOLOR){
Lambda[NUMGEN-1]=Umatrix();

for(i=0;i<NCOLOR;i++){
Lambda[NUMGEN-1].set(i,i,(1.0/sqrt(NCOLOR))*Complex(0.0,1.0));
}
}

// test orthogonality
Complex trace,trace2;
for(i=0;i<NUMGEN;i++){
for(j=0;j<NUMGEN;j++){
trace=Tr(Lambda[i]*Lambda[j]);
//trace2=Tr(Lambda[i]*Trans(Lambda[j]));
if(trace.norm()>0.000001){
cout << "TrT_"<< i << "T_" << j << "= " <<  trace << "\n" <<
flush;}
//if(trace2.norm()>0.00001){
//cout << "TrT_"<< i << "T^T_" << j << "= " << trace2 << "\n" << flush;}
}}

epsilon();


/*SMALLCUT=1.0e-7;
LARGECUT=1000.0;
double ERR=1.0e-5;
cout << "min and max eigenvalue are " << SMALLCUT << "\t" << LARGECUT << "\n";
cout << "relative error is " << ERR << " in (15,15) approx\n";


// 0.0000001-->1000 1e-5 error

ampdeg = 9.9997112279957390e-02;
amp[0] = 3.6832229992796258e-07; shift[0] = 3.7549480881878877e-09;
amp[1] = 1.3567666284582589e-06; shift[1] = 4.3373206800920752e-08;
amp[2] = 5.1757466437096689e-06; shift[2] = 2.8970668616661478e-07;
amp[3] = 2.0060578172377753e-05; shift[3] = 1.7970117665113235e-06;
amp[4] = 7.7976055655092961e-05; shift[4] = 1.1016220840281374e-05;
amp[5] = 3.0323983201324125e-04; shift[5] = 6.7403510935204510e-05;
amp[6] = 1.1793570136758038e-03; shift[6] = 4.1228407663619111e-04;
amp[7] = 4.5868079395172696e-03; shift[7] = 2.5216729791207432e-03;
amp[8] = 1.7839421514438226e-02; shift[8] = 1.5423383071080004e-02;
amp[9] = 6.9386638849859295e-02; shift[9] = 9.4337434071853923e-02;
amp[10] = 2.6997414529708952e-01; shift[10] = 5.7713151658913675e-01;
amp[11] = 1.0526731536884490e+00; shift[11] = 3.5350396633271388e+00;
amp[12] = 4.1584233028628317e+00; shift[12] = 2.1815101171813343e+01;
amp[13] = 1.7800823020581991e+01; shift[13] = 1.4102992696626504e+02;
amp[14] = 1.2795681699057995e+02; shift[14] = 1.2544425313051306e+03;
*/

// 10 term approx:
SMALLCUT=1.0e-4;
LARGECUT=1000.0;
double ERR=4.0e-5;
cout << "min and max eigenvalue are " << SMALLCUT << "\t" << LARGECUT << "\n";
cout << "relative error is " << ERR << " in (10,10) approx\n";
ampdeg = 9.9586400598319705e-02;
amp[0] = 3.6274357276778920e-04; shift[0] = 3.6870476914394509e-05;
amp[1] = 1.3221487344738422e-03; shift[1] = 4.2221229509987286e-04;
amp[2] = 4.9863095947783572e-03; shift[2] = 2.7843257882578039e-03;
amp[3] = 1.9111315256425405e-02; shift[3] = 1.7026645143951285e-02;
amp[4] = 7.3471253678035947e-02; shift[4] = 1.0286504285489308e-01;
amp[5] = 2.8269836634904605e-01; shift[5] = 6.2033453997467913e-01;
amp[6] = 1.0902839420690116e+00; shift[6] = 3.7445215727458883e+00;
amp[7] = 4.2622626238054480e+00; shift[7] = 2.2777671191312244e+01;
amp[8] = 1.8089088885699283e+01; shift[8] = 1.4530492478812033e+02;
amp[9] = 1.2952606180372499e+02; shift[9] = 1.2795625017612963e+03;


if(D==2){
side[0]=LX;
side[1]=T;
Lattice_Map[0]=1;
Lattice_Map[1]=LX;}

if(D==4){
side[0]=LX;
side[1]=LY;
side[2]=LZ;
side[3]=T;
Lattice_Map[0]=1;
Lattice_Map[1]=LX;
Lattice_Map[2]=LY*LX;
Lattice_Map[3]=LZ*LY*LX;
}

return;
}

