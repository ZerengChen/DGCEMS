#pragma once

#ifndef Turbulence
#define Turbulence

#include"GlobalVar.h"
#include <omp.h>

//extern MeshUnion mesh;
//extern const MeshUnion *meshunion;

//#define max(a, b) ((a > b) ? a : b)
//#define min(a, b) ((a < b) ? a : b)

#define TURBULENCE_mp_INIT_TURBULENCE __turbulence_MOD_init_turbulence
#define MTRIDIAGONAL_mp_INIT_TRIDIAGONAL __mtridiagonal_MOD_init_tridiagonal
#define TURBULENCE_mp_DO_TURBULENCE __turbulence_MOD_do_turbulence
#define MTRIDIAGONAL_mp_CLEAN_TRIDIAGONAL __mtridiagonal_MOD_clean_tridiagonal
#define TURBULENCE_mp_CLEAN_TURBULENCE __turbulence_MOD_clean_turbulence
#define TURBULENCE_mp_TKE __turbulence_MOD_tke
#define TURBULENCE_mp_EPS __turbulence_MOD_eps
#define TURBULENCE_mp_L __turbulence_MOD_l
#define TURBULENCE_mp_NUM __turbulence_MOD_num
#define TURBULENCE_mp_NUH __turbulence_MOD_nuh

struct F90Array {
	double * dataPtr; //base address typically the first descriptor element
	int metadata[10]; //this should be only big enough to hold the remainder of the descriptor
};

	void UpdateEddyViscosity(double *fphys_, double ImplicitA_, double *fphys2d_, int Np2d, int K2d, int Np3d, int K3d, long long int nlev, double* VCV, double* , double* , double *Hhuv2d, int);
	void UpdateEddyViscosity_CW(double *fphys_, double ImplicitA_, double *fphys2d_, int Np2d, int K2d, int Np3d, int K3d, long long int nlev, double* VCV, double*, double*, double *Hhuv2d, int, double *UBOT, double *TMBOT);

	/*This function is used to interpolate the physical value from interpolation point to central point in vertical direction*/
	void InterpolationToCentralPoint(double *, double *, int *, int *, int *, double *);
	/*This function is used to map the date located at the central point to the vertical line that the GOTM adapted*/
	void mapCentralPointDateToVerticalDate(double *, double *, int, long long int, int);
	/*This function is used to calculate the shear production term*/
	void CalculateShearFrequencyDate(double *, int, int, double, long long int);

	//void CalculateBuoyanceFrequencyDate(int , int , long long int );
	void CalculateBuoyanceFrequencyDate(double *, int, int, double, long long int, double, double);
	// CalculateBuoyanceProductionTerm to be added 
	//Here, z0b is the bottom roughness, utaub is the friction velocity, z0s is the surface roughness
	void CalculateLengthScaleAndShearVelocity(double, double, double *, double, double *, double *, double *, int, int, long long int);//void DGDoTurbulence(double *dt, double *H2d, double *ShearProductionDate, double *buoyanceProductionDate, double *utaus, double *utaub, double *z0s, double *z0b);

	void DGDoTurbulence(double *, double *, double, double *, int, int, long long int);
	/*This function is used to map the date calculated by GOTM to the output matrix*/
	void mapVedgeDateToDof(double *, double *, int, int, int, long long int);
	/*This function is used to initialize the GOTM module*/
	void InitTurbulenceModelGOTM(long long int *, char *, long long int, long long int, int, int);
	/*This function is used to calculate the water depth for each layer*/
	void CalculateWaterDepth(double *, int, int, double, long long int);

	void getGotmDate(int, long long int);

	void setGotmDate(int, long long int);

	void MatrixMultiplyForGOTM(double *, double *, double *, const int, const int, const int, const double);

	/*The following is the GOTM part*/
extern "C" {
	void TURBULENCE_mp_INIT_TURBULENCE(long long int *, char *, long long int *, long long int);

	void MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(long long int *);

	void TURBULENCE_mp_DO_TURBULENCE(long long int *, double *, double *, double *,
		double *, double *, double *, double *, double *, double *, double *);
	//void do_turbulence_mp_test_do_turbulence_(int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
	void MTRIDIAGONAL_mp_CLEAN_TRIDIAGONAL();

	void TURBULENCE_mp_CLEAN_TURBULENCE();

	void CleanGOTM();

	extern F90Array TURBULENCE_mp_TKE;
	extern F90Array TURBULENCE_mp_EPS;
	extern F90Array TURBULENCE_mp_L;
	extern F90Array TURBULENCE_mp_NUM;
	extern F90Array TURBULENCE_mp_NUH;
}
#endif	
//};
