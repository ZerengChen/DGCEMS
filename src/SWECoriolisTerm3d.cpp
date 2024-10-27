#include "SWECoriolisTerm3d.h"
#include "NdgMath.h"

extern double Latitude;
extern signed char *Status3d;

SWECoriolisTerm3d::SWECoriolisTerm3d()
{
}

SWECoriolisTerm3d::~SWECoriolisTerm3d()
{
}
extern double *BBE, *BBE3d;

void SWECoriolisTerm3d::EvaluateCoriolisTermRHS(double *fphys_, double *frhs_)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	int K2d = *(meshunion->mesh2d_p->K2d);
	int NLayer = *(meshunion->Nlayer);
	int Nz = *meshunion->cell_p->Nz;
	double *OutputRHS = frhs_;
	double *hu = fphys_;
	double *hv = fphys_ + K * Np;
	double *h = fphys_ + K * Np * 3;
	double *MRHS = (double *)malloc(Np*K*2 * sizeof(double));
	memset(MRHS, 0, Np * K * 2 * sizeof(double));

	const double Deg2Rad = 0.0174532925;
	const double omega2 = 0.0001454441043;//4pi/86400

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {		
		NdgExtend2dField(BBE3d, BBE, Np2d, i, Np, NLayer, Nz);
		NdgExtend2dField(BBE3d + Np * K, BBE + Np2d * K2d, Np2d, i, Np, NLayer, Nz);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionDry) {
			for (int n = 0; n < Np; n++) {
				if (h[k * Np + n] >= 2 * Hcrit) {
					MRHS[k * Np + n] = hv[k * Np + n] * sin(Latitude * Deg2Rad) * omega2;
					MRHS[K * Np + k * Np + n] = -hu[k * Np + n] * sin(Latitude * Deg2Rad) * omega2;
				}
				else {
					MRHS[k * Np + n] = hv[k * Np + n] * sin(Latitude * Deg2Rad) * omega2 - BBE3d[k * Np + n];
					MRHS[K * Np + k * Np + n] = -hu[k * Np + n] * sin(Latitude * Deg2Rad) * omega2 - BBE3d[K * Np + k * Np + n];
				}
			}

			Add(OutputRHS + k * Np, MRHS + k * Np, OutputRHS + k * Np, Np);
			Add(OutputRHS + Np * K + k * Np, MRHS + Np * K + k * Np, OutputRHS + Np * K + k * Np, Np);
		}

	}

	free(MRHS);
}
