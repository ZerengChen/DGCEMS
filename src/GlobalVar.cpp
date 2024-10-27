#include"GlobalVar.h"
//const MeshUnion *meshunion = &mesh;
int *dg_demo_comp_id=new int;//for coupler
int *decomp_id=new int;
int *grid_h2d_id=new int;//for coupler
int *time_step=new int;//for coupler
int *coupling_freq=new int;//for coupler
int First_CouplingStep;//for coupler
int Calculating_CouplingSteps;//for coupler

double ftime;
double dt;
double tidalinterval;
int  NOut;
double  gra;
double  Hcrit;
int  Nvar;
int  Nfield3d; //{'hu'  'hv'  'omega'  'h'  'nv'  'z'  'eta'  'zx'  'zy'  'w'  'hw'  'hc'}
int  Nfield2d;//{'h'  'hu'  'hv'  'z'  'zx'  'zy'}

//Barotropic
//double  rho;
double  rho0;
double T0;//Init temperature(oC)
double S0;//Init salinity(g/kg)
double alphaT;
double betaS;
//Vertical diffusion
//char buf;//GotmFile
double z0b;//Value of the physical bottom roughness length:0.0015
double z0s;//Value of the physical surface roughness length:0.02
double WindTauxC;//constant wind stress at x (tau/rho),Np2d*K2d
double WindTauyC;//constant wind stress at y (tau/rho),Np2d*K2d
int Switch_Limiter3D;
int Switch_GOTM;
double cf;//Drag force coefficient
double nv_con;//constant vertical diffusion
double Latitude;
//for tecplot output
int StepNumber;//How many steps to output

