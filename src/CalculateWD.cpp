#include "NdgMath.h"
#include "CalculateWD.h"
#include <stdbool.h>
#include <stdlib.h>
#include <omp.h>

extern int Nvar;
extern double Hcrit;
extern signed char *Status3d;
extern double  rho0;
extern double  T0;
extern double  S0;
extern double *Eta2d;
#define EPSILON 1.0e-12
#define INF 10.0e5
//extern double *h_mean, *hu_mean, *hv_mean;

void JudgingNodeAndElementWD(double *source_h2d, double *source_hu2d, double *source_hv2d, double *source_z, \
	double *hu3d, double *hv3d, double *h3d, double *Limited_huhv2D,signed char *dest, int Np2d, int K2d, int Np3d, int NLayer, int k) {
	//Define some initial value
	int _Nz = *(meshunion->cell_p->Nz) + 1;
	int wetNodeNum = 0; // The total number of wet nodes in an element
	NdgRegionType type = NdgRegionWet; // initialize the element statu to wet element
	double b_max = source_z[0]; // the initial maximal bottom elevation
	double zeta_max = source_z[0] + source_h2d[0]; // the initial maximal water elevation

	//Judge WD of Nodes in an element
	for (int n = 0; n < Np2d; n++) {
		if (source_h2d[n] > Hcrit) { // wet nodes
			wetNodeNum += 1;
		}
		else {
			continue;
		}
	}

	//Judge WD of the element and update the WD statu in cellType
	if (wetNodeNum == 0) { // dry element
		type = NdgRegionDry;
		//Correct the h, hu2d, hv2d to zero
		for (int n1 = 0; n1 < Np2d; n1++) {
			if (source_h2d[n1] < 0.0) {
				source_h2d[n1] = 0.0;
			}
			source_hu2d[n1] = 0.0;
			source_hv2d[n1] = 0.0;
		}
		//Correct the h3d, hu3d, hv3d to zero
		for (int i = 0; i < NLayer; i++) {
			for (int j = 0; j < _Nz; j++) {
				for (int n = 0; n < Np2d; n++) {
					h3d[Np3d * i + Np2d * j + n] = source_h2d[n];
					hu3d[Np3d * i + Np2d * j + n] = source_hu2d[n];//0
					hv3d[Np3d * i + Np2d * j + n] = source_hv2d[n];//0
#ifdef _BAROCLINIC
					hu3d[14 * Np3d * K2d * NLayer + Np3d * i + Np2d * j + n] = 0.0;//hT==0
					hu3d[15 * Np3d * K2d * NLayer + Np3d * i + Np2d * j + n] = 0.0;//hS==0
#endif
				}
			}
		}
	}
	else if (wetNodeNum < Np2d) { // partial wet element
		/***  Restructing the element  ***/
		/* 1. Calculate h_mean, hu_mean, hv_mean of the partial WD element.*/
		/* 2. Calculate theta.*/
		/* 3. Restruct h2d,hu2d and hv2d.*/
		WDelementRestructing(source_h2d, source_hu2d, source_hv2d, Np3d, NLayer, Np2d, K2d, k, type);
		type = NdgRegionPartialWetDamBreak;
		//if (type == NdgRegionPartialWet) {
		//	/****** A refined version to judge the flood and dambreak ******/
		//	for (int n = 0; n < Np2d; n++) {
		//		b_max = fmax(source_z[n], b_max); // find the maximal bottom elevation
		//		zeta_max = fmax((source_z[n] + source_h2d[n]), zeta_max);// find the maximal water elevation
		//	}
		//	//For Flooding, the maximal bottom elevation is higher than the maximal water elevation
		//	if (b_max + Hcrit > zeta_max) {
		//		type = NdgRegionPartialWetFlood;
		//		//type = NdgRegionDry;
		//	}
		//	//For dambreak, the maximal bottom elevation is lower than the maximal water elevation
		//	else {
		//		type = NdgRegionPartialWetDamBreak;
		//		//type = NdgRegionDry;
		//	}
		//}

		/* 4. Extend to 3d variables.*/
		for (int i = 0; i < NLayer; i++) {
			for (int j = 0; j < _Nz; j++) {
				for (int n = 0; n < Np2d; n++) {
					h3d[Np3d * i + Np2d * j + n] = source_h2d[n];
#ifdef _BAROCLINIC
					hu3d[14 * Np3d * K2d * NLayer + Np3d * i + Np2d * j + n] = 0.0; //hT==0
					hu3d[15 * Np3d * K2d * NLayer + Np3d * i + Np2d * j + n] = 0.0; //hS==0
#endif
					//if (source_h2d[n] > Hcrit) {
					//	hu3d[Np3d * i + Np2d * j + n] = hu3d[Np3d * i + Np2d * j + n] + source_hu2d[n] - Limited_huhv2D[n];
					//	hv3d[Np3d * i + Np2d * j + n] = hv3d[Np3d * i + Np2d * j + n] + source_hv2d[n] - Limited_huhv2D[K2d * Np2d + n];
					//}
					//else {
						hu3d[Np3d * i + Np2d * j + n] = source_hu2d[n];//0
						hv3d[Np3d * i + Np2d * j + n] = source_hv2d[n];//0
					//}
				}
			}
		}

        /****** ********************** END ********************** ******/		
	}

	*dest = (signed char)type;
}

/* The elementRestructing is just for partial WD elements. */
void WDelementRestructing(double *h2d, double *hu2d, double *hv2d, int Np3d, int NLayer, int Np2d, int K2d, int k, NdgRegionType type2d_) {
	double *wq2d = meshunion->mesh2d_p->mesh2dcell_p->wq2d;
	double *Jacobian2d = meshunion->mesh2d_p->J2d;
	double *Vq2d = meshunion->mesh2d_p->mesh2dcell_p->Vq2d;
	int *Nq2d = meshunion->mesh2d_p->mesh2dcell_p->Nq2d;
	double *LAV2d = meshunion->mesh2d_p->LAV2d;

	int oneI = 1;
	double _h_mean = 0.0;
	double _hu_mean = 0.0;
	double _hv_mean = 0.0;
	double h_min = h2d[0];
	//double h_max = h2d[0];
	double theta = 1.0; // A parameter to keep positive

	// Get the average h2d of an element 
	GetMeshAverageValue(&_h_mean, LAV2d + k, Nq2d, &oneI, &Np2d, Vq2d, h2d + k * Np2d, Jacobian2d + k * Np2d, Nq2d, wq2d);
	GetMeshAverageValue(&_hu_mean, LAV2d + k, Nq2d, &oneI, &Np2d, Vq2d, hu2d + k * Np2d, Jacobian2d + k * Np2d, Nq2d, wq2d);
    GetMeshAverageValue(&_hv_mean, LAV2d + k, Nq2d, &oneI, &Np2d, Vq2d, hv2d + k * Np2d, Jacobian2d + k * Np2d, Nq2d, wq2d);

	// Calculate theta and restruct the h2d, hu2d, hv2d
	for (int i = 1; i < Np2d; i++) {
		h_min = fmin(h_min, h2d[i]);
		//h_max = max(h_max, h2d[i]);
	}

	if (_h_mean > Hcrit) {
		type2d_ = NdgRegionPartialWet;
		theta = fmin(_h_mean / (_h_mean - h_min + EPSILON), 1.0);
		theta = fmax(0.0, theta);
		
		for (int n = 0; n < Np2d; n++) {
			h2d[n] = _h_mean + theta * (h2d[n] - _h_mean);
			//hu2d[n] = _hu_mean / _h_mean * h2d[n];
			//hv2d[n] = _hv_mean / _h_mean * h2d[n];
			hu2d[n] = _hu_mean + theta * (hu2d[n] - _hu_mean);
			hv2d[n] = _hv_mean + theta * (hv2d[n] - _hv_mean);

			if (h2d[n] <= Hcrit) {
				hu2d[n] = 0.0;
				hv2d[n] = 0.0;
			}
		}
	}
	else {
		type2d_ = NdgRegionDry;
		//for (int n = 0; n < Np2d; n++) {
		//	if (h2d[n] < 0.0) {
		//		h2d[n] = 0.0;
		//	}

		//	hu2d[n] = 0.0;
		//	hv2d[n] = 0.0;
		//}
	}

}

void UpdateWetDryState(double *fphys_, double *fphys2d_,double *Limited_huhv2D, int *NLayer_, signed char *status_, int Np3d, int K3d, int Np2d, int K2d) {
	double *h2d = fphys2d_;
	double *hu3d = fphys_;
	double *hv3d = fphys_ + Np3d * K3d;
	double *h3d = fphys_ + Np3d * K3d * 3;
	int NLayer = *NLayer_;
	signed char *cellType = status_;

	/*Judge the WD status in each element. Firstly, confirm wet, dry and patial status.*/
	/*** For partial WD elements, we must restruct the elements. ***/
	/*Finally, confirm the patial status Flooding or Dambreaking.*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
        JudgingNodeAndElementWD(h2d + k * Np2d, h2d + K2d * Np2d + k * Np2d, h2d + K2d * Np2d * 2 + k * Np2d, h2d + K2d * Np2d * 3 + k * Np2d, \
			hu3d + k * Np3d * NLayer, hv3d + k * Np3d * NLayer, h3d + k * Np3d * NLayer, Limited_huhv2D + k * Np2d,cellType + k, Np2d, K2d, Np3d, NLayer,k);

		for (int i = 0; i < NLayer; i++) {
			Status3d[k * NLayer + i] = cellType[k];
		}
	}
}

void Limiter2d(double *fphys2d_, int fieldID, double *Limited_huhv2D) {

	/* get inputs */
	int Np = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	int K = *(meshunion->mesh2d_p->K2d);
	double *fphys = fphys2d_;
	double *avar = (double *)malloc(K * sizeof(double));
	memset(avar, 0, K * sizeof(double));
	double *EToE = meshunion->mesh2d_p->EToE2d;
	int Nface = *(meshunion->mesh2d_p->mesh2dcell_p->Nface2d); // number of faces of the computation cell

	double one = 1;
	double zero = 0;
	int one_ptrdiff = 1;
	double *wq = meshunion->mesh2d_p->mesh2dcell_p->wq2d;
	double *Jacobian = meshunion->mesh2d_p->J2d;
	double *Vq = meshunion->mesh2d_p->mesh2dcell_p->Vq2d;
	int *Nq = meshunion->mesh2d_p->mesh2dcell_p->Nq2d;
	double *LAV = meshunion->mesh2d_p->LAV2d;

	signed char *status = meshunion->mesh2d_p->status;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int n = 0; n < Np; n++) {
			Eta2d[k * Np + n] = fphys[k * Np + n] + fphys[3 * Np * K + k * Np + n];//use fphys_5 to save zeta2d
	}
}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif 
	for (int k = 0; k < K; k++) {
		GetMeshAverageValue(avar + k, LAV + k, Nq, &one_ptrdiff, &Np, Vq, Eta2d + k * Np, Jacobian + k * Np, Nq, wq);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif    
	for (int k = 0; k < K; k++) {
		/*This part is used to determine the maxmum and minimum allowable value of the studied cell with index k */
		double amax = -1 * pow(10, 10), amin = pow(10, 10);
		for (int f = 0; f < Nface; f++) {
			/*Only the adjacent cell is considered, and the average value of the studied cell is not included*/
			if ((int)EToE[k*Nface + f] - 1 != k) {
				amax = fmax(amax, avar[(int)EToE[k*Nface + f] - 1]);
				amin = fmin(amin, avar[(int)EToE[k*Nface + f] - 1]);
			}
		}
		/*This part is used to decide whether the studied cell is problematic following Cockburn and Shu*/
		int flag = 0;
		for (int i = 0; i < Np; i++) {
			if (Eta2d[k * Np + i] > amax || Eta2d[k * Np + i] < amin) {
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			continue;
		else {
			double Lambda = 1;
			for (int i = 0; i < Np; i++) {
				if (Eta2d[k * Np + i] > amax)
					Lambda = fmin(Lambda, (amax - avar[k]) / (Eta2d[k*Np + i] - avar[k] + pow(10, -10)));
				//Lambda = fmin(Lambda, (amax - avar[k]) / (fphys[k * Np + i] - avar[k]));
				else if (Eta2d[k * Np + i] < amin)
					Lambda = fmin(Lambda, (avar[k] - amin) / (avar[k] - Eta2d[k*Np + i] + pow(10, -10)));
				//Lambda = fmin(Lambda, (avar[k] - amin) / (avar[k] - fphys[k * Np + i]));
			}
			Lambda = fmax(0.0, Lambda);
			for (int i = 0; i < Np; i++) {
				Eta2d[k * Np + i] = Lambda * Eta2d[k * Np + i] + (1 - Lambda)*avar[k];
			}
		}

	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {	
		if (status[k] == NdgRegionWet) {
			for (int j = 0; j < Np; j++) {			
				fphys[k * Np + j] = Eta2d[k * Np + j] - fphys[3 * Np * K + k * Np + j];
			}
		}
	}

	free(avar);
}