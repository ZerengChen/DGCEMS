#include "NdgSWEVertGOTMDiffSolver.h"
#include"UpdateEddyViscosity.h"
#include "NdgMath.h"

NdgSWEVertGOTMDiffSolver::NdgSWEVertGOTMDiffSolver()
{
}

NdgSWEVertGOTMDiffSolver::~NdgSWEVertGOTMDiffSolver()
{
}
extern double *BBE, *SBE, *Hhuv2d;
extern double Hcrit;
extern int Nvar;
extern double nv_con;
extern double cf;
extern double WindTauxC;
extern double WindTauyC;
extern int Switch_GOTM;
extern signed char *Status3d;

void NdgSWEVertGOTMDiffSolver::EvaluateVertDiffRHS(double *fphys_, double *frhs_, double *time_, double *fphys2d, double ImplicitA_, int *varIndex)
{
	int K2d = *(meshunion->mesh2d_p->K2d);
	int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	int Np3d = *(meshunion->cell_p->Np);
	int K3d = *(meshunion->K);
	int NLayer = *(meshunion->Nlayer);
	int Nz = *meshunion->cell_p->Nz;
	long long int nlev = *(meshunion->Nlayer);
    double *VCV = meshunion->cell_p->VCV;
	int oneI = 1;
	double *h = fphys_ + 3 * K3d * Np3d;
	double *nv_v = fphys_ + 4 * K3d * Np3d;

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
	/***  Wet and Dry  ***/
	memset(Hhuv2d, 0, sizeof(double)*(Np2d * K2d * 3));

	if ((int)Switch_GOTM == 1) {
		UpdateEddyViscosity(fphys_, ImplicitA_, fphys2d, Np2d, K2d, Np3d, K3d, nlev, VCV, BBE, SBE, Hhuv2d, NLayer);

#ifdef _BAROCLINIC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			for (int n = 0; n < Np2d; n++) {
				BBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				BBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;

				SBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				SBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;
			}
		}
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K3d; k++) {
			if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = 0.0;
				}
			}
			else {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = nv_v[k * Np3d + n] / h[k * Np3d + n] / h[k * Np3d + n];
					//nv_v[k * Np3d + n] = 0.0001 / h[k * Np3d + n] / h[k * Np3d + n];// For SaltyWater Case
				}
			}
		}
	}
	else if ((int)Switch_GOTM == 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			if ((NdgRegionType)Status2d[k] == NdgRegionDry) {
				continue;
			}
			else {
				MatrixMultiply(VCV, fphys_ + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hu
				MatrixMultiply(VCV, fphys_ + K3d * Np3d + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hv
				MatrixMultiply(VCV, fphys_ + K3d * Np3d * 3 + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d * 2 + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom h
			}
		}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			for (int n = 0; n < Np2d; n++) {
				if (Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] > Hcrit) {
					BBE[k * Np2d + n] = cf * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					BBE[K2d * Np2d + k * Np2d + n] = cf * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);

					////--------------For wind driven flow
					//BBE[k * Np2d + n] = cf * Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					//BBE[K2d * Np2d + k * Np2d + n] = 0.0;//cf * Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					////--------------For wind driven flow

					SBE[k * Np2d + n] = WindTauxC;
					SBE[K2d * Np2d + k * Np2d + n] = WindTauyC;
				}
				else {
					BBE[k * Np2d + n] = 0.0;
					BBE[K2d * Np2d + k * Np2d + n] = 0.0;
					SBE[k * Np2d + n] = 0.0;
					SBE[K2d * Np2d + k * Np2d + n] = 0.0;
				}

			}
		}

#ifdef _BAROCLINIC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			for (int n = 0; n < Np2d; n++) {
				BBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				BBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;

				SBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				SBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;
			}
		}
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K3d; k++) {
			if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = 0.0;
				}
			}
			else {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = nv_con / h[k * Np3d + n] / h[k * Np3d + n];
				}
			}
		}

	}
	else {
		std::cout << "Warning! The turbulent model is undefined, please check." << std::endl;
	}

	updateimplicitrhs.EvaluateupdateimplicitRHS(fphys_, nv_v, frhs_, ImplicitA_, BBE, SBE, varIndex);

}

void NdgSWEVertGOTMDiffSolver::EvaluateVertDiffRHS_CW(double *fphys_, double *frhs_, double *time_, double *fphys2d, double ImplicitA_, double *UBOT, double *TMBOT, int *varIndex)
{
	int K2d = *(meshunion->mesh2d_p->K2d);
	int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	int Np3d = *(meshunion->cell_p->Np);
	int K3d = *(meshunion->K);
	int NLayer = *(meshunion->Nlayer);
	int Nz = *meshunion->cell_p->Nz;
	long long int nlev = *(meshunion->Nlayer);
	double *VCV = meshunion->cell_p->VCV;
	int oneI = 1;
	double *h = fphys_ + 3 * K3d * Np3d;
	double *nv_v = fphys_ + 4 * K3d * Np3d;

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
	/***  Wet and Dry  ***/
	memset(Hhuv2d, 0, sizeof(double)*(Np2d * K2d * 3));

	if ((int)Switch_GOTM == 1) {
		UpdateEddyViscosity_CW(fphys_, ImplicitA_, fphys2d, Np2d, K2d, Np3d, K3d, nlev, VCV, BBE, SBE, Hhuv2d, NLayer, UBOT, TMBOT);

#ifdef _BAROCLINIC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			for (int n = 0; n < Np2d; n++) {
				BBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				BBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;

				SBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				SBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;
			}
		}
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K3d; k++) {
			if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = 0.0;
				}
			}
			else {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = nv_v[k * Np3d + n] / h[k * Np3d + n] / h[k * Np3d + n];
				}
			}
		}
	}
	else if ((int)Switch_GOTM == 0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			if ((NdgRegionType)Status2d[k] == NdgRegionDry) {
				continue;
			}
			else {
				MatrixMultiply(VCV, fphys_ + K3d * Np3d * 9 + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hu
				MatrixMultiply(VCV, fphys_ + K3d * Np3d * 10+ k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hv
				MatrixMultiply(VCV, fphys_ + K3d * Np3d * 3 + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d * 2 + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom h
			}
		}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			for (int n = 0; n < Np2d; n++) {
				//if (Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] > Hcrit) {
					BBE[k * Np2d + n] = cf * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					BBE[K2d * Np2d + k * Np2d + n] = cf * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);

					////--------------For wind driven flow
					//BBE[k * Np2d + n] = cf * Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					//BBE[K2d * Np2d + k * Np2d + n] = 0.0;//cf * Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					////--------------For wind driven flow

					SBE[k * Np2d + n] = WindTauxC;
					SBE[K2d * Np2d + k * Np2d + n] = WindTauyC;
				//}
				//else {
				//	BBE[k * Np2d + n] = 0.0;
				//	BBE[K2d * Np2d + k * Np2d + n] = 0.0;
				//	SBE[k * Np2d + n] = 0.0;
				//	SBE[K2d * Np2d + k * Np2d + n] = 0.0;
				//}

			}
		}

#ifdef _BAROCLINIC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K2d; k++) {
			for (int n = 0; n < Np2d; n++) {
				BBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				BBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;

				SBE[K2d * Np2d * 2 + k * Np2d + n] = 0.0;
				SBE[K2d * Np2d * 3 + k * Np2d + n] = 0.0;
			}
		}
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < K3d; k++) {
			if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = 0.0;
				}
			}
			else {
				for (int n = 0; n < Np3d; n++) {
					nv_v[k * Np3d + n] = nv_con / h[k * Np3d + n] / h[k * Np3d + n];
				}
			}
		}


	}
	else {
		std::cout << "Warning! The turbulent model is undefined, please check." << std::endl;
	}

	updateimplicitrhs.EvaluateupdateimplicitRHS(fphys_, nv_v, frhs_, ImplicitA_, BBE, SBE, varIndex);

}