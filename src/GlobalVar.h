#pragma once
//#include"NdgPhysMat.h"

extern int *dg_demo_comp_id;//for coupler
extern int *decomp_id;
extern int *grid_h2d_id;//for coupler
extern int *time_step;//for coupler
extern int *coupling_freq;//for coupler
extern int First_CouplingStep;//for coupler
extern int Calculating_CouplingSteps;//for coupler

extern double ftime;
extern double dt;
extern double tidalinterval;
extern int  NOut;
extern double  gra;
extern double  Hcrit;
extern int  Nvar;
extern int  Nfield3d;
extern int  Nfield2d;

//extern double  rho;
extern double  rho0;
extern double  T0;
extern double  S0;
extern double  alphaT;
extern double  betaS;

//extern char buf;
extern double  z0b;
extern double  z0s;
extern double  WindTauxC;
extern double  WindTauyC;
extern int  Switch_Limiter3D;
extern int  Switch_GOTM;
extern double  cf;
extern double  nv_con;
extern double  Latitude;
//for tecplot output
extern int StepNumber;