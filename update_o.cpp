
#include "update.h"

void update(Gauge_Field &U, Twist_Fermion &F){

  Gauge_Field p_U,old_U,old_f_U;
  Twist_Fermion p_F,old_F,old_f_F;
  
  double H_old,K_old,K_new,S_new,H_new,hmc_test;
  
  static Gauge_Field f_U;
  static Twist_Fermion f_F;
  
  static double S_old=0.0;
  static int first_time=1,accept=0,no_calls=0;
  static ofstream f_hmc;
  
  int i;
  
  no_calls++;
  
  // refresh momenta
  p_U=Gauge_Field(2);
  
  p_F=Twist_Fermion(1);
  
  K_old=kinetic_energy(p_U,p_F);
  
  
  if(first_time){
    f_hmc.open("hmc_test");
    if(f_hmc.bad()){
      cout << "failed to open hmc_test file\n" << flush;}
    
    S_old=action(U,F);
    force(U,f_U,F,f_F,1);
    first_time=0;
  }	
  
  if((no_calls%100==0)&&(!first_time)){
    cout << "acceptance rate " << (double)accept/(double)no_calls << "\n" <<
      flush;
    no_calls=0;
    accept=0;
  }
  
  
  
  
  //save copies of fields
  old_U=U;
  old_F=F;
  
  H_old=S_old+K_old;
  
  old_f_U=f_U;
  old_f_F=f_F;
  
  

  
  
  /* 
   classical evolution: this is different from evolve_field in that the loop is in
   evolve_fields_o.cpp
  */
  
    evolve_fields(U,p_U,f_U,F,p_F,f_F);
    
  
  K_new=kinetic_energy(p_U,p_F);

  //cout << "  MOM " << mid1 << "\n";
  S_new=action(U,F);
  
  
  H_new=S_new+K_new;
  cout << "  ACTION BEGIN TRAJ " << H_old << "\n" << flush;
  cout << "  ACTION END TRAJ " << H_new << "\n" << flush;
  
  hmc_test=exp(-H_new+H_old);
  f_hmc << hmc_test << "\n" << flush;
  
  //metropolis test
  if((rand()/(double)RAND_MAX)<exp(H_old-H_new)){
    cout << "ACCEPT  " << H_old << "\t" << H_new <<"\t" << H_old-H_new << "\n" << flush;
    cout << "++++++++++++++++++++++++++++++++++++" << endl ; 
    S_old=S_new;
    accept++;
    return;
  }
  else{
    cout << "hmc_test " << hmc_test << " failed\n" << flush;
    cout << "REJECT  " << H_old << "\t" << H_new <<"\t" << H_old-H_new << "\n" << flush;
    cout << "++++++++++++++++++++++++++++++++++++" << endl ; 
	// if fails copy back fields
    U=old_U;
    F=old_F;
    f_F=old_f_F;
    f_U=old_f_U;
    return;
  }	   
  
  
  
  return;	    
}	    
	    
	    
	    
	    
	    
