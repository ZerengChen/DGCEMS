//#include "mex.h"
#include "NdgSWE3D.h"
#include "NdgMath.h"
#include "VerticalRepmatFacialValue.h"

VerticalRepmatFacialValue::VerticalRepmatFacialValue()
{
}

VerticalRepmatFacialValue::~VerticalRepmatFacialValue()
{
}

void VerticalRepmatFacialValue::EvaluateRepmatFacialValue(double *field2d_, double *field3d_) {

	int Nfp = *(meshunion->boundaryedge_p->Nfp);
	int Ne = (*meshunion->boundaryedge_p->Ne);
	int NLayer = (*meshunion->boundaryedge_p->Nlayer);
	double *FToF = meshunion->boundaryedge_p->FToF;

	int Npz = (*meshunion->cell_p->Npz);
	int LNfp = Nfp / Npz;

	double *field2d = field2d_;
	int Ne2d = (*meshunion->mesh2d_p->mesh2dboundaryedge_p->Ne2d);

	double *DepthExtendValue = field3d_;
	//void VerticalRepmatFacialValue(double *dest, double *source, int NLayer, int Nfp, int LNfp, int Npz, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < Ne2d; e++){
		VerticalRepmatFacialValue_(DepthExtendValue + e*NLayer*Nfp, field2d + e*LNfp, NLayer, Nfp, LNfp, Npz, (int)(*(FToF + 2 * e*NLayer)));
	}
}