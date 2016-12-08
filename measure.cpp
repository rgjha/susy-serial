#include <cstdio>
#include "measure.h"  

void measure(const Gauge_Field &U, const Twist_Fermion &F,int &num){

	static int first_time=1;
	static ofstream f_data,f_av,f_scalars;
	static ofstream f_line,f_kon1,f_kon2;
        static ofstream f_loop,f_loop2,f_det;
        double act_s,act_F,mass,t1,t2;
	double eigenvals[SITES][NUMLINK][NCOLOR];
	int i,j,k,sites,mu,r,m;
	Lattice_Vector x;
	Gauge_Field Utmp,U2,U3,U4;
        Complex av2,d;
	UPlaq_Field B;
        double wilson[LX][T],wilson2[LX][T];
#ifdef FULLMATRIX
Complex M[LEN][LEN];
#endif 
        cout << " in measure \n" << flush;	
	if(first_time){
	f_data.open("data",ios::app);
	if(f_data.bad()){ 
	cout << "failed to open data file\n" << flush ;}

        f_scalars.open("scalars",ios::app);
        if(f_scalars.bad()){
        cout << "failed to open scalars file\n" << flush;}

	f_line.open("lines",ios::app);
	if(f_line.bad()){
	cerr << "failed to open lines_t file" << "\n";exit(1);}

        f_loop.open("loops",ios::app);
        if(f_loop.bad()){
        cerr << "failed to open loops file" << "\n";exit(1);}
       
	first_time=0;
	}

	// check Bianchi
	if(NUMLINK==5){
	B=Bianchi(U);
	cout << "Bianchi is " << Tr(B*Adj(B)).real() << flush << "\n";
	}
	
	f_line  << line(U,D-1) << "\n" << flush;
    
  
        loop(U,wilson);  

        for(r=1;r<=(LX/2);r++){
        for(m=1;m<=(T/2);m++){
	f_loop << r << "\t" << m <<  "\t" << wilson[r][m]/((D-1)*SITES) << "\t" << flush;}
        }
        f_loop  << "\n" << flush;


        obs(U,F,act_s,mass,d,act_F,eigenvals);

        
#ifdef FULLMATRIX
        if((num%(10*GAP)==0)&&(FERMIONS==1)){
	full_fermion_op(U,M);
        eigenvalues(M);
	(void)Pfaffian(M);
}
#endif

	f_data << act_s << "\t" << act_F << "\n" << flush;
    
        scalars(U,t1,t2);
	
        f_scalars << t1 << "\t" << t2 <<  "\n" << flush;

        if(SMALLEIG==1){	
	sites=0;
	while(loop_over_lattice(x,sites)){
	for(j=0;j<NUMLINK;j++){
	for(k=0;k<NCOLOR;k++){
        f_av << eigenvals[sites-1][j][k] << "\t" ;}
	f_av << "\n";}}
	f_av << flush;}

	return;
}
