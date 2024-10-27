#include "VertLimit3d.h"
#include "NdgMath.h"
#include <stdlib.h>

extern int Nvar;
extern signed char *Status3d;
/*
 * Purpose: This function is used to limit the physical value to satisfy the maxmium principle
 * This function is programmed according to ( Philippe Delandmeter, 2017 ).
 *
 * Input:
 *      double[Np x K] fphys the physical field to be limited.
 * 		double[1 x K]  avar the average value of each cell for the physical field
 * 		double[Nface x K]  EToE the topological relation of the studied mesh
 * 		double[Npz x Npz]  V1d the inversed one-dimensional Vandmonde matrix, used to calculate the average value in vertical direction for each line
 *      int[1] Npz number of interpolation points in vertical direction
 * 		int[1] Nph number of interpolation points in horizontal direction
 *      double[Npz x Npz]  OV1d the one-dimensional Vandmonde matrix, the first value of this matrix is used to calculate the average value
 * Output:
 * 		double[Np x K] limfphys the limited physical field.
 */

void Limiter3d(double *fphys_) {
    
    /* get inputs */
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
    double *fphys = fphys_;
    double *avar = (double *)malloc(K*sizeof(double));
	memset(avar, 0, K * sizeof(double));
    double *EToE = meshunion->EToE;
	int Nface = *(meshunion->cell_p->Nface); // number of faces of the computation cell

    /*This is the inversed one-dimensional Vandmande matrix*/
	int Nz = *(meshunion->cell_p->Nz);
    double *V1d = (double*)malloc((Nz + 1)*(Nz + 1) * sizeof(double));
    const int Npz = *(meshunion->cell_p->Npz);
    const int Nph = *(meshunion->cell_p->Nph);
    /*This is the one-dimensional VandMande matrix*/
 //   double *OV1d = meshunion->cell_p->V1d;   
	//cblas_dcopy((Nz + 1)*(Nz + 1), OV1d, 1, V1d, 1);
	//MatrixInverse(V1d, Nz + 1);


    double one = 1;
    double zero = 0;
    int one_ptrdiff = 1;
    int Np_ptrdiff = Npz;
    int Nq_ptrdiff = Npz;
	double *wq = meshunion->cell_p->wq;
	double *Jacobian = meshunion->J;
	double *Vq = meshunion->cell_p->Vq;
	int *Nq = meshunion->cell_p->Nq;
	double *LAV = meshunion->LAV;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif 
	for (int k = 0; k < K; k++) {
		GetMeshAverageValue(avar + k, LAV + k, Nq, &one_ptrdiff, &Np, Vq, fphys + k * Np, Jacobian + k * Np, Nq, wq);
	}

	//double *avL = (double*)malloc(sizeof(double) * Nph), \
	//	*tempValue = (double*)malloc(sizeof(double) * Npz), \
	//	*fmod = (double*)malloc(sizeof(double)*Npz);
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
/*---------------------------------------------------------------------------------------------------------------------------------*/
		//for (int i = 0; i < Np; i++) {
		//	if (fphys[k * Np + i] > amax || fphys[k * Np + i] < amin) {
		//		flag = 1;
		//		break;
		//	}
		//}

		//if (flag == 0) //the studied cell satisfy the maximum principle
		//	continue;
		//else {// the studied cell is problematic
		//	// Calculate the average value of the vertical line
		//	memset(avL, 0, sizeof(double) * Nph);
		//	memset(tempValue, 0, sizeof(double) * Nph);
		//	memset(fmod, 0, sizeof(double) * Nph);
		//	for (int i = 0; i < Nph; i++) {
		//		for (int j = 0; j < Npz; j++) {
		//			//Fetch the original value in each line and store them in tempValue
		//			*(tempValue + j) = fphys[k * Np + Nph * (Npz - 1) + i - j * Nph];
		//		}
		//		// Calculate the corresponding mode coefficients for each vertical line
		//		MatrixMultiplyTT(V1d, tempValue, fmod, Nq_ptrdiff, one_ptrdiff, Nq_ptrdiff, 1.0);
		//		//MatrixMultiply(V1d, tempValue, fmod, Nq_ptrdiff, one_ptrdiff, Nq_ptrdiff, 1.0);
		//		// the final average data for each vertical line
		//		*(avL + i) = fmod[0] * OV1d[0];
		//		//Next limit the corresponding date
		//		/*Calculate the slope parameter first*/
		//		double Lambda = 1.0;
		//		if (fphys[k*Np + Nph * (Npz - 1) + i] > amax)
		//			Lambda = fmin(Lambda, (amax - avL[i]) / ( fphys[k*Np + Nph*(Npz - 1) + i] - avL[i] + pow(10,-10.0) ));
		//			//Lambda = fmin(Lambda, (amax - avL[i]) / (fphys[k*Np + Nph * (Npz - 1) + i] - avL[i]));
		//		else if (fphys[k*Np + Nph * (Npz - 1) + i] < amin)
		//			Lambda = fmin(Lambda, (avL[i] - amin) / (avL[i] - fphys[k*Np + Nph*(Npz - 1) + i] + pow(10, -10.0)) );
		//			//Lambda = fmin(Lambda, (avL[i] - amin) / (avL[i] - fphys[k*Np + Nph * (Npz - 1) + i]));
		//		else
		//			Lambda = 1.0;
		//		/*Limit the value in vertical direction*/
		//		Lambda = fmax(0, Lambda);
		//		for (int j = 0; j < Npz; j++) {
		//			fphys[k*Np + Nph * (Npz - 1) + i - j * Nph] = Lambda * fphys[k*Np + Nph * (Npz - 1) + i - j * Nph]\
		//				+ (1 - Lambda) * avL[i];
		//		}
		//	}
		//}
/*---------------------------------------------------------------------------------------------------------------------------------*/
			/*Next to limit the whole computation cell follows Cockburn and Shu*/
			/*This part is used to decide whether the studied cell is problematic following Cockburn and Shu*/
			flag = 0;
			for (int i = 0; i < Np; i++) {
				if (fphys[k * Np + i] > amax || fphys[k * Np + i] < amin) {
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				continue;
			else {
				double Lambda = 1;
				for (int i = 0; i < Np; i++) {
					if (fphys[k * Np + i] > amax)
						Lambda = fmin(Lambda, (amax - avar[k]) / (fphys[k*Np + i] - avar[k] + pow(10, -10)));
						//Lambda = fmin(Lambda, (amax - avar[k]) / (fphys[k * Np + i] - avar[k]));
					else if (fphys[k * Np + i] < amin)
						Lambda = fmin(Lambda, (avar[k] - amin) / (avar[k] - fphys[k*Np + i] + pow(10, -10)));
						//Lambda = fmin(Lambda, (avar[k] - amin) / (avar[k] - fphys[k * Np + i]));
				}
				Lambda = fmax(0.0, Lambda);
				for (int i = 0; i < Np; i++) {
					fphys[k * Np + i] = Lambda * fphys[k * Np + i] + (1 - Lambda)*avar[k];
				}
			}

    }

	free(V1d);
	//free(avL);
	//free(tempValue);
	//free(fmod);
	free(avar);
}
