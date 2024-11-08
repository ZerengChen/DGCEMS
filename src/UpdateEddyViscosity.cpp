#include<string.h>
//#include "cblas.h"
#include "NdgMemory.h"
#include"UpdateEddyViscosity.h"
#include "NdgMath.h"
#include <cmath>
#include <iostream>

extern double gra;
extern double Hcrit;
extern int  Nvar;
extern double z0b;
extern double dt;
extern double z0s;
//extern char buf;
extern double rho0;
extern double WindTauxC;
extern double WindTauyC;

extern double *tkeGOTM, *epsGOTM, *LGOTM, *nuhGOTM, *numGOTM, *layerHeight, *huCentralDate, \
*hvCentralDate, *huVerticalLine, *hvVerticalLine, *shearFrequencyDate, *buoyanceFrequencyDate, *BottomFrictionLength, \
*BottomFrictionVelocity, *SurfaceFrictionLength, *SurfaceFrictionVelocity, *eddyViscosityDate, \
*rhoCentralDate, *rhoVerticalLine, *eddyDiffusionDate, *eddyTKEDate, *eddyLengthDate, *eddyEPSDate;

/*the Von kamma constant*/
double kappa = 0.41;
const char *GOTMInitialized = "False";

void getGotmDate(int index, long long int nlev) {
	for (int i = 0; i < (int)nlev + 1; i++) {
		tkeGOTM[index*((int)nlev + 1) + i] = *(TURBULENCE_mp_TKE.dataPtr + i);
		epsGOTM[index*((int)nlev + 1) + i] = *(TURBULENCE_mp_EPS.dataPtr + i);
		LGOTM[index*((int)nlev + 1) + i] = *(TURBULENCE_mp_L.dataPtr + i);
		nuhGOTM[index*((int)nlev + 1) + i] = *(TURBULENCE_mp_NUH.dataPtr + i);
		numGOTM[index*((int)nlev + 1) + i] = *(TURBULENCE_mp_NUM.dataPtr + i);
	}
}

void setGotmDate(int index, long long int nlev) {
	for (int i = 0; i < (int)nlev + 1; i++) {
		*(TURBULENCE_mp_TKE.dataPtr + i) = tkeGOTM[i + index * ((int)nlev + 1)];
		*(TURBULENCE_mp_EPS.dataPtr + i) = epsGOTM[i + index * ((int)nlev + 1)];
		*(TURBULENCE_mp_L.dataPtr + i) = LGOTM[i + index * ((int)nlev + 1)];
		*(TURBULENCE_mp_NUH.dataPtr + i) = nuhGOTM[i + index * ((int)nlev + 1)];
		*(TURBULENCE_mp_NUM.dataPtr + i) = numGOTM[i + index * ((int)nlev + 1)];
	}
}


void InitTurbulenceModelGOTM(long long int *NameList, char * buf, long long int buflen, long long int nlev, int Np2d, int K2d) {

	TURBULENCE_mp_INIT_TURBULENCE(NameList, buf, &nlev, buflen);

	MTRIDIAGONAL_mp_INIT_TRIDIAGONAL(&nlev);

	for (int i = 0; i < Np2d*K2d; i++) {
		getGotmDate(i, nlev);
	}
}

void InterpolationToCentralPoint(double *fphys, double *dest, int *Np2d, int *K3d, int *Np3d, double *VCV,int*pE3d,int MyID) {
	double alpha = 1;
	int Col = 1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < (int)(*K3d); i++) {
		if (MyID == pE3d[i]) {
			MatrixMultiplyForGOTM(VCV, fphys + i * (*Np3d), dest + i * (*Np2d), *Np2d, Col, *Np3d, alpha);
		}
	}
}

void mapCentralPointDateToVerticalDate(double *centralDate, double *verticalLineDate, int K2d, long long int nlev, int Np2d, int*pE2d, int MyID) {
	//This has been verified by tests
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			for (int L = 1; L < nlev + 1; L++) {
				for (int p = 0; p < Np2d; p++) {
					verticalLineDate[k*(nlev + 1)*Np2d + p * (nlev + 1) + L] = \
						centralDate[k*nlev*Np2d + (nlev - L)*Np2d + p];
				}
			}
		}
	}
}

void CalculateWaterDepth(double *H2d, int Np2d, int K2d, double hcrit, long long int nlev, int*pE2d, int MyID) {
	/*if the water depth is larger than the threshold value, the layer height is calculated by the ratio of water depth to number of lyers.
	*For the current version, we only consider the equalspace division in vertical. When the water depth is less than the threshold, the water
	*depth is set to be zero
	*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int p = 0; p < Np2d * K2d; p++) {
		if (MyID == pE2d[p/Np2d]) {
			if (H2d[p] >= hcrit) {
				for (int L = 1; L < nlev + 1; L++) {
					layerHeight[p * (nlev + 1)] = 0.0;
					layerHeight[p * (nlev + 1) + L] = H2d[p] / nlev;
				}
			}
		}
	}

}


void CalculateShearFrequencyDate(double *H2d, int Np2d, int K2d, double hcrit, long long int nlev, int*pE2d, int MyID) {
	//SS = $(\frac{\partial u}{\partial x})^2+(\frac{\partial v}{\partial y})^2$
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int p = 0; p < Np2d*K2d; p++) {
		if (MyID == pE2d[p/Np2d]) {
			if (H2d[p] >= hcrit) {
				for (int L = 1; L < nlev; L++) {
					shearFrequencyDate[p*(nlev + 1) + L] = pow((huVerticalLine[p*(nlev + 1) + L + 1] - huVerticalLine[p*(nlev + 1) + L]) / H2d[p] / (0.5*(layerHeight[p*(nlev + 1) + L + 1] + layerHeight[p*(nlev + 1) + L])), 2) \
						+ pow((hvVerticalLine[p*(nlev + 1) + L + 1] - hvVerticalLine[p*(nlev + 1) + L]) / H2d[p] / (0.5*(layerHeight[p*(nlev + 1) + L + 1] + layerHeight[p*(nlev + 1) + L])), 2);
				}
				//For each vertical segment, we have SS(0) = SS(1), SS(nlev) = SS(nlev - 1)
				shearFrequencyDate[p*(nlev + 1)] = shearFrequencyDate[p*(nlev + 1) + 1];
				shearFrequencyDate[p*(nlev + 1) + nlev] = shearFrequencyDate[p*(nlev + 1) + nlev - 1];
			}
		}
	}

}

void CalculateLengthScaleAndShearVelocity(double z0b, double z0s, double *H2d, double hcrit, double *DragCoefficient, \
	double *Taux, double *Tauy, int Np2d, int K2d, long long int nlev, int*pE2d, int MyID) {
	/*for surface friction length, another way is the charnock method*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int p = 0; p < Np2d * K2d; p++) {
		if (MyID == pE2d[p/Np2d]) {
			SurfaceFrictionLength[p] = z0s;
			SurfaceFrictionVelocity[p] = 0;
			DragCoefficient[p] = 0;
			double rr = 0;
			if (H2d[p] >= hcrit) {
				rr = kappa / log((z0b + layerHeight[p*(nlev + 1) + 1] / 2) / z0b);
				BottomFrictionVelocity[p] = rr * sqrt(pow((huVerticalLine[p*(nlev + 1) + 1] / H2d[p]), 2) + pow((hvVerticalLine[p*(nlev + 1) + 1] / H2d[p]), 2));
				BottomFrictionLength[p] = z0b;

				SurfaceFrictionVelocity[p] = sqrt(pow(Taux[p], 2) + pow(Tauy[p], 2));

				DragCoefficient[p] = pow(rr, 2);   //Cd = rr^2;

			}
		}
	}
}

void DGDoTurbulence(double *TimeStep, double *H2d, double hcrit, double *Grass, int Np2d, int K2d, long long int nlev, int*pE2d, int MyID) {
	//For the current version, grass is not considered
	//int _nlev = nlev - 1;
	for (int p = 0; p < Np2d * K2d; p++) {
		if (MyID == pE2d[p / Np2d]) {
			if (H2d[p] >= hcrit) {
				setGotmDate(p, nlev);
				TURBULENCE_mp_DO_TURBULENCE(&nlev, TimeStep, H2d + p, SurfaceFrictionVelocity + p, BottomFrictionVelocity + p, SurfaceFrictionLength + p, \
					BottomFrictionLength + p, layerHeight + p * (nlev + 1), buoyanceFrequencyDate + p * (nlev + 1), shearFrequencyDate + p * (nlev + 1), Grass);
				getGotmDate(p, nlev);

				for (int L = 0; L < nlev + 1; L++) {
					eddyTKEDate[p*(nlev + 1) + L] = *(TURBULENCE_mp_TKE.dataPtr + L);
					eddyEPSDate[p*(nlev + 1) + L] = *(TURBULENCE_mp_EPS.dataPtr + L);
					eddyLengthDate[p*(nlev + 1) + L] = *(TURBULENCE_mp_L.dataPtr + L);
					eddyViscosityDate[p*(nlev + 1) + L] = *(TURBULENCE_mp_NUM.dataPtr + L);
					eddyDiffusionDate[p*(nlev + 1) + L] = *(TURBULENCE_mp_NUH.dataPtr + L);
				}
			}
		}
	}

}

void mapVedgeDateToDof(double *SourceDate, double *DestinationDate, int Np2d, int K2d, int Np3d, long long int nlev, int*pE2d, int MyID) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			for (int p = 0; p < Np2d; p++) {
				DestinationDate[k*nlev*Np3d + (nlev - 1)*Np3d + p] = SourceDate[k*Np2d*(nlev + 1) + p * (nlev + 1)];//the down face of the bottommost cell for each column
				DestinationDate[k*nlev*Np3d + p + Np2d] = SourceDate[k*Np2d*(nlev + 1) + p * (nlev + 1) + nlev];//the upper face of the topmost cell for each column
				for (int L = 1; L < nlev; L++) {
					DestinationDate[k*nlev*Np3d + (nlev - L)*Np3d + p + Np2d] = SourceDate[k*Np2d*(nlev + 1) + p * (nlev + 1) + L];  //The top layer of the down cell
					DestinationDate[k*nlev*Np3d + (nlev - L - 1)*Np3d + p] = SourceDate[k*Np2d*(nlev + 1) + p * (nlev + 1) + L];//The bottom layer of the up cell
				}
			}
		}
	}
}

void CalculateBuoyanceFrequencyDate(double *H2d, int Np2d, int K2d, double hcrit, long long int nlev, \
	double gra, double rho0, int*pE2d, int MyID) {

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int p = 0; p < Np2d*K2d; p++) {
		if (MyID == pE2d[p/Np2d]) {
			if (H2d[p] >= hcrit) {
				for (int L = 1; L < nlev; L++) {
					buoyanceFrequencyDate[p*(nlev + 1) + L] = -1 * gra / rho0 * (rhoVerticalLine[p*(nlev + 1) + L + 1] - rhoVerticalLine[p*(nlev + 1) + L]) / (0.5*(layerHeight[p*(nlev + 1) + L + 1] + layerHeight[p*(nlev + 1) + L]));
					buoyanceFrequencyDate[p*(nlev + 1) + L] = fmax(buoyanceFrequencyDate[p*(nlev + 1) + L], 0.0);
				}
				//For each vertical segment, we have NN(0) = NN(1), NN(nlev) = NN(nlev - 1)
				buoyanceFrequencyDate[p*(nlev + 1)] = buoyanceFrequencyDate[p*(nlev + 1) + 1];
				buoyanceFrequencyDate[p*(nlev + 1) + nlev] = buoyanceFrequencyDate[p*(nlev + 1) + nlev - 1];
			}
			else
			{
				for (int L = 0; L < nlev; L++) {
					//buoyanceFrequencyDate[p*(nlev + 1)] = 0.0;
					buoyanceFrequencyDate[p*(nlev + 1) + L] = 0;
				}
			}
		}
	}
}

void MatrixMultiplyForGOTM(double *matrix1, double *matrix2, double *matrix3, const int M_, const int N_, const int K_, const double Alpha)
{
	int i, j, k;
	double *A1 = (double *)malloc((int)M_ * (int)K_ * sizeof(double));
	double *B1 = (double *)malloc((int)K_ * (int)N_ * sizeof(double));
	double *C1 = (double *)malloc((int)M_ * (int)N_ * sizeof(double));
	memset(A1, 0, (int)M_ * (int)K_ * sizeof(double));
	memset(B1, 0, (int)K_ * (int)N_ * sizeof(double));
	memset(C1, 0, (int)M_ * (int)N_ * sizeof(double));

	for (i = 0; i < (int)M_; i++) {
		for (k = 0; k < (int)K_; k++) {
			A1[i * K_ + k] = matrix1[k * M_ + i];
		}
	}

	for (k = 0; k < (int)K_; k++) {
		for (j = 0; j < (int)N_; j++) {
			B1[k * N_ + j] = Alpha * matrix2[j * K_ + k];
		}
	}

	for (i = 0; i < (int)M_; i++) {
		for (j = 0; j < (int)N_; j++) {
			for (k = 0; k < (int)K_; k++) {
				C1[j * M_ + i] = C1[j * M_ + i] + A1[i * K_ + k] * B1[k * N_ + j];
			}
		}
	}

	for (i = 0; i < (int)M_; i++) {
		for (j = 0; j < (int)N_; j++) {
			matrix3[i * N_ + j] = C1[i * N_ + j];
		}
	}

	free(A1); A1 = NULL;
	free(B1); B1 = NULL;
	free(C1); C1 = NULL;
}

void CleanGOTM() {
	TURBULENCE_mp_CLEAN_TURBULENCE();
    MTRIDIAGONAL_mp_CLEAN_TRIDIAGONAL();
}

void UpdateEddyViscosity(double *fphys_, double ImplicitA_, double *fphys2d_, int Np2d, int K2d, int Np3d, int K3d, long long int nlev, double* VCV, double*BBE,double*SBE,double *Hhuv2d,int NLayer,int*pE2d,int*pE3d,int MyID)
{
	double* h = fphys2d_;
	double* hu = fphys_;
	double* hv = fphys_ + Np3d * K3d;
	double* WindTaux = (double*)malloc(Np2d*K2d * sizeof(double));
	double* WindTauy = (double*)malloc(Np2d*K2d * sizeof(double));
	double* rho = fphys_ + 13 * Np3d * K3d;
	int Num2d = Np2d * K2d;
	long long int Interface = nlev + 1;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			for (int n = 0; n < Np2d; n++) {
				WindTaux[k * Np2d + n] = WindTauxC;
				WindTauy[k * Np2d + n] = WindTauyC;
			}
		}
	}
#ifndef _BAROCLINIC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		if (MyID == pE3d[k]) {
			for (int n = 0; n < Np3d; n++) {
				rho[k * Np3d + n] = rho0;
			}
		}
	}
#endif
	
	char buf0[] = "/data/home/DG/czr/opt/GOTM/lib/gotmturb.nml";
	//char buf0[] = "/data2/DG/CZR/opt/GOTM/lib/gotmturb.nml";
	long long int buflen = sizeof(buf0);
	char *buf = (char *)malloc(buflen);
	buf = buf0;
	long long int _nNamelist = 2;
	int oneI = 1;

	double *PtrOutEddyViscosity = (double *)malloc(Np3d*K3d * sizeof(double));
	double *PtrOutDragCoefficient = (double *)malloc(Np2d*K2d * sizeof(double));

	double *PtrOutTKE = (double *)malloc(Np3d*K3d * sizeof(double));

	double *PtrOutEPS = (double *)malloc(Np3d*K3d * sizeof(double));

	if (GOTMInitialized == "False") {
		InitTurbulenceModelGOTM(&_nNamelist, buf, buflen, nlev, Np2d, K2d);

		GOTMInitialized = "true";
	}
	//else {
		//std::cout<< "GOTM has already initialized! " << std::endl;
	//}

	int TempNp2d = Np2d;
	int TempNp3d = Np3d;
	int TempK3d = K3d;

	InterpolationToCentralPoint(hu, huCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV, pE3d, MyID);

	InterpolationToCentralPoint(hv, hvCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV, pE3d, MyID);
	/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S.
	Details about the latter manner can be found in Tuomas and Vincent(2012, Ocean modelling)
	*/
	InterpolationToCentralPoint(rho, rhoCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV, pE3d, MyID);

	//Tc to be continued
	//Sc to be continued
	mapCentralPointDateToVerticalDate(huCentralDate, huVerticalLine, K2d, nlev, Np2d, pE2d, MyID);

	mapCentralPointDateToVerticalDate(hvCentralDate, hvVerticalLine, K2d, nlev, Np2d, pE2d, MyID);
	/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S*/
	mapCentralPointDateToVerticalDate(rhoCentralDate, rhoVerticalLine, K2d, nlev, Np2d, pE2d, MyID);
	//Tvl to be continued
	//Svl to be continued
	CalculateWaterDepth(h, Np2d, K2d, Hcrit, nlev, pE2d, MyID);

	CalculateShearFrequencyDate(h, Np2d, K2d, Hcrit, nlev, pE2d, MyID);

	CalculateBuoyanceFrequencyDate(h, Np2d, K2d, Hcrit, nlev, gra, rho0, pE2d, MyID);

	CalculateLengthScaleAndShearVelocity(z0b, z0s, h, Hcrit, PtrOutDragCoefficient, WindTaux, WindTauy, Np2d, K2d, nlev, pE2d, MyID);

//-----------------------------------------------------------------------------------Neumann BC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			MatrixMultiply(VCV, fphys_ + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hu
			MatrixMultiply(VCV, fphys_ + K3d * Np3d + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hv
			MatrixMultiply(VCV, fphys_ + K3d * Np3d * 3 + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d * 2 + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom h	
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			for (int n = 0; n < Np2d; n++) {
				if (Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] > Hcrit) {
					BBE[k * Np2d + n] = PtrOutDragCoefficient[k * Np2d + n] * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					BBE[K2d * Np2d + k * Np2d + n] = PtrOutDragCoefficient[k * Np2d + n] * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);

					SBE[k * Np2d + n] = WindTaux[k * Np2d + n];
					SBE[K2d * Np2d + k * Np2d + n] = WindTauy[k * Np2d + n];
				}
				else {
					BBE[k * Np2d + n] = 0.0;
					BBE[K2d * Np2d + k * Np2d + n] = 0.0;
					SBE[k * Np2d + n] = 0.0;
					SBE[K2d * Np2d + k * Np2d + n] = 0.0;
				}

			}
		}
	}
//-----------------------------------------------------------------------------------END Neumann BC

	DGDoTurbulence(&dt, h, Hcrit, NULL, Np2d, K2d, nlev, pE2d, MyID);

	mapVedgeDateToDof(eddyViscosityDate, PtrOutEddyViscosity, Np2d, K2d, Np3d, nlev, pE2d, MyID);
	/* If the following parts are to be exported, the corresponding parts in DGDoTurbulence in file mxGOTM.c need to be activated*/


	mapVedgeDateToDof(eddyTKEDate, PtrOutTKE, Np2d, K2d, Np3d, nlev, pE2d, MyID);


	mapVedgeDateToDof(eddyEPSDate, PtrOutEPS, Np2d, K2d, Np3d, nlev, pE2d, MyID);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		if (MyID == pE3d[k]) {
			for (int n = 0; n < Np3d; n++) {
				fphys_[4 * K3d*Np3d + k * Np3d + n] = PtrOutEddyViscosity[k * Np3d + n];

			}
		}
	}
	free(PtrOutEddyViscosity);
	free(PtrOutDragCoefficient);
	free(PtrOutTKE);
	free(PtrOutEPS);
	free(WindTaux);
	free(WindTauy);
	//free(rho);
	//free(buf);
}

void UpdateEddyViscosity_CW(double *fphys_, double ImplicitA_, double *fphys2d_, int Np2d, int K2d, int Np3d, int K3d, long long int nlev, double* VCV, double*BBE, double*SBE, double *Hhuv2d, int NLayer, double *UBOT, double *TMBOT,int*pE2d,int*pE3d,int MyID)
{
	double* h = fphys2d_;
	double* hu = fphys_ + 9 * Np3d * K3d;
	double* hv = fphys_ + 10 * Np3d * K3d;

	double* WindTaux = (double*)malloc(Np2d*K2d * sizeof(double));
	double* WindTauy = (double*)malloc(Np2d*K2d * sizeof(double));
	double* rho = fphys_ + 13 * Np3d * K3d;
	int Num2d = Np2d * K2d;
	long long int Interface = nlev + 1;
	const double eps1 = 1.0E-14;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			for (int n = 0; n < Np2d; n++) {
				WindTaux[k * Np2d + n] = WindTauxC;
				WindTauy[k * Np2d + n] = WindTauyC;
			}
		}
	}
#ifndef _BAROCLINIC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		if (MyID == pE3d[k]) {
			for (int n = 0; n < Np3d; n++) {
				rho[k * Np3d + n] = rho0;
			}
		}
	}
#endif


	char buf0[] = "/data/home/DG/czr/opt/GOTM/lib/gotmturb.nml";
	//char buf0[] = "/data2/DG/CZR/opt/GOTM/lib/gotmturb.nml";
	long long int buflen = sizeof(buf0);
	char *buf = (char *)malloc(buflen);
	buf = buf0;
	long long int _nNamelist = 2;
	int oneI = 1;

	double *PtrOutEddyViscosity = (double *)malloc(Np3d*K3d * sizeof(double));
	double *PtrOutDragCoefficient = (double *)malloc(Np2d*K2d * sizeof(double));
	double *Tau_w = (double *)malloc(Np2d*K2d * sizeof(double));
	double *Tau_CW = (double *)malloc(Np2d*K2d * sizeof(double));
	double *fw = (double *)malloc(Np2d*K2d * sizeof(double));
	double *Ab = (double *)malloc(Np2d*K2d * sizeof(double));
	double *O_Ab = (double *)malloc(Np2d*K2d * sizeof(double));
	double *midTau = (double *)malloc(Np2d*K2d * sizeof(double));
	memset(Tau_w, 0, Np2d*K2d * sizeof(double));
	memset(Tau_CW, 0, Np2d*K2d * sizeof(double));
	memset(fw, 0, Np2d*K2d * sizeof(double));
	memset(Ab, 0, Np2d*K2d * sizeof(double));
	memset(O_Ab, 0, Np2d*K2d * sizeof(double));
	memset(midTau, 0, Np2d*K2d * sizeof(double));

	double *PtrOutTKE = (double *)malloc(Np3d*K3d * sizeof(double));

	double *PtrOutEPS = (double *)malloc(Np3d*K3d * sizeof(double));

	if (GOTMInitialized == "False") {
		InitTurbulenceModelGOTM(&_nNamelist, buf, buflen, nlev, Np2d, K2d);

		GOTMInitialized = "true";
	}
	//else {
		//std::cout<< "GOTM has already initialized! " << std::endl;
	//}

	int TempNp2d = Np2d;
	int TempNp3d = Np3d;
	int TempK3d = K3d;

	InterpolationToCentralPoint(hu, huCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV, pE3d, MyID);

	InterpolationToCentralPoint(hv, hvCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV, pE3d, MyID);
	/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S.
	Details about the latter manner can be found in Tuomas and Vincent(2012, Ocean modelling)
	*/
	InterpolationToCentralPoint(rho, rhoCentralDate, &TempNp2d, &TempK3d, &TempNp3d, VCV, pE3d, MyID);

	//Tc to be continued
	//Sc to be continued
	mapCentralPointDateToVerticalDate(huCentralDate, huVerticalLine, K2d, nlev, Np2d, pE2d, MyID);

	mapCentralPointDateToVerticalDate(hvCentralDate, hvVerticalLine, K2d, nlev, Np2d, pE2d, MyID);
	/*The gradient about rho in vertical direction is calculated according to rho directly, not T and S*/
	mapCentralPointDateToVerticalDate(rhoCentralDate, rhoVerticalLine, K2d, nlev, Np2d, pE2d, MyID);
	//Tvl to be continued
	//Svl to be continued
	CalculateWaterDepth(h, Np2d, K2d, Hcrit, nlev, pE2d, MyID);

	CalculateShearFrequencyDate(h, Np2d, K2d, Hcrit, nlev, pE2d, MyID);

	CalculateBuoyanceFrequencyDate(h, Np2d, K2d, Hcrit, nlev, gra, rho0, pE2d, MyID);

	CalculateLengthScaleAndShearVelocity(z0b, z0s, h, Hcrit, PtrOutDragCoefficient, WindTaux, WindTauy, Np2d, K2d, nlev, pE2d, MyID);

#ifdef COUPLING_SWAN
	//wave ShearVelocity and Tau_CW
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			for (int n = 0; n < Np2d; n++) {
				if (h[k * Np2d + n] <= Hcrit) {
					continue;
				}
				else {
					Ab[k * Np2d + n] = UBOT[k * NLayer * Np2d + (NLayer - 1) * Np3d + n] * TMBOT[k * NLayer * Np2d + (NLayer - 1) * Np3d + n] / 2.0 / 3.1415926;
					O_Ab[k * Np2d + n] = z0b / (Ab[k * Np2d + n] + eps1);
					fw[k * Np2d + n] = 1.39 * pow(O_Ab[k * Np2d + n], 0.52);
					Tau_w[k * Np2d + n] = 0.5 * fw[k * Np2d + n] * UBOT[k * NLayer * Np2d + (NLayer - 1) * Np3d + n] * UBOT[k * NLayer * Np2d + (NLayer - 1) * Np3d + n];
					midTau[k * Np2d + n] = Tau_w[k * Np2d + n] / (Tau_w[k * Np2d + n] + PtrOutDragCoefficient[k * Np2d + n] + eps1);
					Tau_CW[k * Np2d + n] = PtrOutDragCoefficient[k * Np2d + n] * (1.0 + 1.2 * pow(midTau[k * Np2d + n], 3.2));
				}
			}
		}
	}
#endif

	//-----------------------------------------------------------------------------------Neumann BC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			MatrixMultiply(VCV, hu + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hu
			MatrixMultiply(VCV, hv + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom hv
			MatrixMultiply(VCV, fphys_ + K3d * Np3d * 3 + k * NLayer*Np3d + (NLayer - 1) * Np3d, Hhuv2d + K2d * Np2d * 2 + k * Np2d, Np2d, oneI, Np3d, 1.0);//bottom h		
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if (MyID == pE2d[k]) {
			for (int n = 0; n < Np2d; n++) {
				//if (Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] > Hcrit) {
					//BBE[k * Np2d + n] = Tau_CW[k * Np2d + n] * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
					//	pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
					//	Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					//BBE[K2d * Np2d + k * Np2d + n] = Tau_CW[k * Np2d + n] * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
					//	pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
					//	Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);

					BBE[k * Np2d + n] = PtrOutDragCoefficient[k * Np2d + n] * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);
					BBE[K2d * Np2d + k * Np2d + n] = PtrOutDragCoefficient[k * Np2d + n] * sqrt(pow(Hhuv2d[k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2) + \
						pow(Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n], 2)) * \
						Hhuv2d[K2d * Np2d + k * Np2d + n] / Hhuv2d[K2d * Np2d * 2 + k * Np2d + n] * (-1);

					SBE[k * Np2d + n] = WindTaux[k * Np2d + n];
					SBE[K2d * Np2d + k * Np2d + n] = WindTauy[k * Np2d + n];
				//}
				//else {
				//	BBE[k * Np2d + n] = 0.0;
				//	BBE[K2d * Np2d + k * Np2d + n] = 0.0;
				//	SBE[k * Np2d + n] = 0.0;
				//	SBE[K2d * Np2d + k * Np2d + n] = 0.0;
				//}

			}
		}
	}
	//-----------------------------------------------------------------------------------END Neumann BC

	DGDoTurbulence(&dt, h, Hcrit, NULL, Np2d, K2d, nlev, pE2d, MyID);

	mapVedgeDateToDof(eddyViscosityDate, PtrOutEddyViscosity, Np2d, K2d, Np3d, nlev, pE2d, MyID);

	mapVedgeDateToDof(eddyTKEDate, PtrOutTKE, Np2d, K2d, Np3d, nlev, pE2d, MyID);

	mapVedgeDateToDof(eddyEPSDate, PtrOutEPS, Np2d, K2d, Np3d, nlev, pE2d, MyID);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		if (MyID == pE3d[k]) {
			for (int n = 0; n < Np3d; n++) {
				fphys_[4 * K3d*Np3d + k * Np3d + n] = PtrOutEddyViscosity[k * Np3d + n];

			}
		}
	}
	free(PtrOutEddyViscosity);
	free(PtrOutDragCoefficient);
	free(PtrOutTKE);
	free(PtrOutEPS);
	free(WindTaux);
	free(WindTauy);
	//free(rho);
	free(Tau_w); Tau_w = NULL;
	free(Tau_CW); Tau_CW = NULL;
	free(fw); fw = NULL;
	free(Ab); Ab = NULL;
	free(O_Ab); O_Ab = NULL;
	free(midTau); midTau = NULL;
	//free(buf);
}