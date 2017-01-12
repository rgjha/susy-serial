#include "action.h"

double action(const Gauge_Field &U, const Twist_Fermion F) {
  Lattice_Vector x,e_mu,e_nu;
  Gauge_Field Udag;
  Complex dum,d,trace;
  Umatrix p,dummy;
  USite_Field DmuUmu;
  UPlaq_Field P;
  double act_s = 0.0, act_F = 0.0;
  int mu,nu,site;
  Twist_Fermion sol[DEGREE],psol[DEGREE];
  Umatrix Fmunu;

  Udag = Adj(U);


  P = Plaq(U);
  site = 0;
  while(loop_over_lattice(x,site)) {
    dum = Complex();
    for (mu = 0;mu<NUMLINK;mu++) {
      for (nu = 0;nu<NUMLINK;nu++) {
        if (mu == nu) continue;
        dum = dum+(det(P.get(x,mu,nu))-Complex(1.0,0.0));
      }
    }
    DmuUmu.set(x,G*dum*Umatrix(1));
  }

  // usual piece
  site = 0;
  while(loop_over_lattice(x,site)) {
    for (mu = 0;mu<NUMLINK;mu++) {
      e_mu = Lattice_Vector(mu);
      DmuUmu.set(x,DmuUmu.get(x)+U.get(x,mu)*Udag.get(x,mu)-
          Udag.get(x-e_mu,mu)*U.get(x-e_mu,mu));}

      act_s = act_s+0.5*C2*Tr(DmuUmu.get(x)*DmuUmu.get(x)).real();
  }

  // Konishi mass term - single trace operator - 'eig'  //
  site = 0;
  dum = Complex();
  while(loop_over_lattice(x,site)) {
    for (mu = 0;mu<NUMLINK;mu++) {
      dummy = U.get(x,mu)*Udag.get(x,mu)-Umatrix(1);
      dum = dum+Tr(dummy*dummy);
    }}

  act_s = act_s+(BMASS*BMASS)*dum.real();

  // Below : Double trace mass term - "trace" //
  /* site = 0;
     while(loop_over_lattice(x,site)) {
     for (mu = 0;mu<NUMLINK;mu++) {
     dum =(1.0/NCOLOR)*Tr(Udag.get(x,mu)*U.get(x,mu)).real()-1.0;
     act_s = act_s+(BMASS*BMASS)*dum*dum;
     }}  */

  site = 0;
  while(loop_over_lattice(x,site)) {
    for (mu = 0;mu<NUMLINK;mu++) {
      e_mu = Lattice_Vector(mu);
      for (nu = mu+1;nu<NUMLINK;nu++) {
        e_nu = Lattice_Vector(nu);
        Fmunu = U.get(x,mu)*U.get(x+e_mu,nu)-U.get(x,nu)*U.get(x+e_nu,mu);
        act_s = act_s+2.0*Tr(Fmunu*Adj(Fmunu)).real();

      }
    }
  }


  act_s = KAPPA*act_s;


  // Now, the fermions //

  // Pseudofermion contribution


  if (FERMIONS) {
    act_F = ampdeg*(Cjg(F)*F).real();

    MCG_solver(U,F,shift,sol,psol);

    for (int n = 0;n<DEGREE;n++) {
      act_F = act_F+amp[n]*(Cjg(F)*sol[n]).real();}
  }

  return (act_s + act_F);
}
