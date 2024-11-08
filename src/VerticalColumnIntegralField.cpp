//#include <mex.h>
#include "NdgSWE3D.h"
#include "NdgMath.h"
#include "VerticalColumnIntegralField.h"

VerticalColumnIntegralField::VerticalColumnIntegralField()
{
}

VerticalColumnIntegralField::~VerticalColumnIntegralField()
{
}

void VerticalColumnIntegralField::EvaluateVerticalIntegral(double *fphys2d, double *fphys,int*pE2d,int MyID){
	int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	int K2d = *(meshunion->mesh2d_p->K2d);
	double *V2d = meshunion->mesh2d_p->mesh2dcell_p->V2d;
	double *V3d = meshunion->cell_p->V;
	double *Jz = meshunion->Jz;
	double *field3d = fphys;
	int NLayer = *(meshunion->Nlayer);
	int Np3d = *(meshunion->cell_p->Np);
	int K3d = *(meshunion->K);

	double *field2d = fphys2d;
	double *Tempfield2d = (double *)malloc(Np2d*K2d*sizeof(double));
	memset(Tempfield2d, 0, Np2d*K2d*sizeof(double));
	double *Tempfield3d = (double *)malloc(Np3d*K3d*sizeof(double));
	double *fmod = (double *)malloc(Np3d*K3d*sizeof(double));
	double *InvV3d = (double *)malloc(Np3d*Np3d*sizeof(double));
	memcpy(InvV3d, V3d, Np3d*Np3d*sizeof(double));
	MatrixInverse(InvV3d, Np3d);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++){
		if (MyID == pE2d[i]) {
			VerticalColumnIntegralField3d(field2d + i * Np2d, Np2d, V2d, Tempfield2d + i * Np2d, \
				Tempfield3d + i * Np3d*NLayer, field3d + i * Np3d*NLayer, Jz + i * Np3d*NLayer, \
				fmod + i * Np3d*NLayer, InvV3d, Np3d, NLayer);
		}
	}
	free(Tempfield3d);
	free(Tempfield2d);
	free(fmod);
	free(InvV3d);
}