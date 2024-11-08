#include "NdgMemory.h"
#include "NdgMath.h"
#include "eqstate.h"
#include "CalculateDensityField.h"
#include "VertLimit3d.h"
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

extern double *BaroclinicT, *BaroclinicS, *BaroclinicDTS;
extern double gra;
extern double Hcrit;
extern int Nvar;
extern int Nfield3d;
extern signed char *Status3d;
extern double  rho0;
extern double  T0;
extern double  S0;
extern double  alphaT;
extern double  betaS;


void CalculateDensityField(double *fphys) {
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	
	double *height = fphys + 3 * Np * K;
    double *hT = fphys + 14 * Np * K;
    double *hS = fphys + 15 * Np * K;
	double *z = meshunion->z;
    double hcrit = Hcrit;
    
    double *rho = fphys + 13 * Np * K;
    
	/***  Wet and Dry  ***/
	signed char *Status2d = meshunion->mesh2d_p->status;
	typedef enum {
		NdgRegionNormal = 1,
		NdgRegionRefine = 2,
		NdgRegionSponge = 3,
		NdgRegionWet = 4,
		NdgRegionDry = 5,
		NdgRegionPartialWet = 6,
		NdgRegionPartialWetFlood = 7,
		NdgRegionPartialWetDamBreak = 8
	} NdgRegionType;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		if ((NdgRegionType)Status3d[i] != NdgRegionWet) {
			continue;
		}
		else {
			DotCriticalDivide(BaroclinicT + i * Np, hT + i * Np, &hcrit, height + i * Np, Np);
			DotCriticalDivide(BaroclinicS + i * Np, hS + i * Np, &hcrit, height + i * Np, Np);
		}
	}

		Limiter3d(BaroclinicT);
		Limiter3d(BaroclinicS);
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K; k++) {
			for (int n = 0; n < Np; n++) {
				BaroclinicT[k * Np + n] = fmax(0.0, BaroclinicT[k * Np + n]);
				BaroclinicS[k * Np + n] = fmax(0.0, BaroclinicS[k * Np + n]);
			}
		}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		if ((NdgRegionType)Status3d[i] != NdgRegionWet) {
			for (int p = 0; p < Np; p++) {
				rho[i * Np + p] = rho0;
				hT[i * Np + p] = 0.0;
				hS[i * Np + p] = 0.0;
			}
		}
		else {
			//		if (!strcmp(EosType, "Jackett05")) {
						//for (int p = 0; p < Np; p++) {
						//	EosByFeistel(rho + i*Np + p, fmax(*(BaroclinicT + i*Np + p), 0.0), fmax(*(BaroclinicS + i*Np + p), 0.0) );
						//    hT[i * Np + p] = fmax(0, BaroclinicT[i * Np + p]) * height[i * Np + p];
						//    hS[i * Np + p] = fmax(0, BaroclinicS[i * Np + p]) * height[i * Np + p];
						//}
			//		}
					//else if (!strcmp(EosType, "UNESCO83")) {
					//	for (int p = 0; p < Np; p++) {
					//		EosByUNESCO(rho + i*Np + p, fmax(*(BaroclinicT + i*Np + p), 0.0), fmax(*(BaroclinicS + i*Np + p), 0.0));
					//      hT[i * Np + p] = fmax(0, BaroclinicT[i * Np + p]) * height[i * Np + p];
					//      hS[i * Np + p] = fmax(0, BaroclinicS[i * Np + p]) * height[i * Np + p];
					//	}
					//}
					//else if (!strcmp(EosType, "Linear")) {
			for (int p = 0; p < Np; p++) {
				EosByLinear(rho + i * Np + p, fmax(*(BaroclinicT + i * Np + p), 0.0), fmax(*(BaroclinicS + i * Np + p), 0.0), rho0, T0, S0, alphaT, betaS);
				BaroclinicT[i * Np + p] = fmin(BaroclinicT[i * Np + p], 20);
				BaroclinicS[i * Np + p] = fmin(BaroclinicS[i * Np + p], 35);
				hT[i * Np + p] = fmax(0, BaroclinicT[i * Np + p]) * height[i * Np + p];
				hS[i * Np + p] = fmax(0, BaroclinicS[i * Np + p]) * height[i * Np + p];
			}
			//}
			//else {
			//	printf("Equation of state(EOS) needs to be pointed for this part!\n");
			//	break;
			//}
		}
    }
}