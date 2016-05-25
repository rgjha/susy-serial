/* Omelyan or 2MN integrator; Takaishi & de Forcrand, Phys Rev E 73, 036706 (2006).
T. D. 25 July 2012

*/


#include "evolve_fields.h"

void evolve_fields(Gauge_Field &U, Gauge_Field &p_U, Gauge_Field &f_U, 
		   Twist_Fermion &F, Twist_Fermion &p_F, Twist_Fermion &f_F){
  
  Lattice_Vector x,y;
  Gauge_Field nf_U;
  Twist_Fermion nf_F;
  int mu,i,site;

  
  void evolve_gauge_step(Gauge_Field &U, Gauge_Field &p_U, Gauge_Field &f_U, 
			 Twist_Fermion &F, Twist_Fermion &p_F, Twist_Fermion &f_F,double eps);
  
  // multi timestep leapfrog update: note convention that TRAJECTORY_LENGTH*DT is hardwired:
  // read_param.cpp:TRAJECTORY_LENGTH=(int)(0.5/DT);
  
  
  
  update_gauge_momenta(p_U,f_U,INT_LAMBDA*DT);
  p_F=p_F+INT_LAMBDA*DT*f_F;
  
  for(int i=1;i<=TRAJECTORY_LENGTH;i++){
    
    
    evolve_gauge_step(U, p_U,f_U, F, p_F, f_F, 0.5*DT);
    
    // pure S_F contribution to force
    force(U,f_U,F,f_F,1);
    
    // intermediate step for momenta
    update_gauge_momenta(p_U,f_U,INT_LAMBDA_MID*DT);
    p_F=p_F+INT_LAMBDA_MID*DT*f_F;
    
    evolve_gauge_step(U, p_U,f_U, F, p_F, f_F, 0.5*DT);
    
    // pure S_F contribution to force
    force(U,f_U,F,f_F,1);
    
    if(i<TRAJECTORY_LENGTH){
      update_gauge_momenta(p_U,f_U,INT_LAMBDA_CONT*DT);
      p_F=p_F+INT_LAMBDA_CONT*DT*f_F;
    }
    else{
      update_gauge_momenta(p_U,f_U,INT_LAMBDA*DT);
      p_F=p_F+INT_LAMBDA*DT*f_F;
    }
    
    
  } // outer loop over steps
  
  return;
}

void evolve_gauge_step(Gauge_Field &U, Gauge_Field &p_U, Gauge_Field &f_U, 
		       Twist_Fermion &F, Twist_Fermion &p_F, Twist_Fermion &f_F, double eps)
{
  int STEPS=5;
  int j;
  
  // inner S_B force update carries forward a time interval `eps'
  
  // pure S_B contribution to force
  force(U,f_U,F,f_F,0);
  update_gauge_momenta(p_U,f_U,INT_LAMBDA*eps/STEPS);
  
  for(j=1;j<=STEPS;j++){
    //cout << "boson step " << endl;
    
    
    update_gauge_field(U,p_U,0.5*eps/STEPS);
    F=F+0.5*eps/STEPS*p_F;
    force(U,f_U,F,f_F,0);
    
    update_gauge_momenta(p_U,f_U,INT_LAMBDA_MID*eps/STEPS);
    
    update_gauge_field(U,p_U,0.5*eps/STEPS);
    F=F+0.5*eps/STEPS*p_F;
    force(U,f_U,F,f_F,0);
    
    if(j<STEPS){
      update_gauge_momenta(p_U,f_U,INT_LAMBDA_CONT*eps/STEPS);
   }
    else{
      update_gauge_momenta(p_U,f_U,INT_LAMBDA*eps/STEPS);
    }
    
    
  }
  
  return;
}
