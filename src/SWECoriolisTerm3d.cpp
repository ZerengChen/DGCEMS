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


void SWECoriolisTerm3d::EvaluateCoriolisTermRHS(double *fphys_, double *frhs_, int*pE3d, int MyID)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	double *OutputRHS = frhs_;
	double *hu = fphys_;
	double *hv = fphys_ + K * Np;
	double *MRHS = (double *)malloc(Np*K*2 * sizeof(double));
	memset(MRHS, 0, Np * K * 2 * sizeof(double));

	const double Deg2Rad = 0.0174532925;
	const double omega2 = 0.0001454441043;//4pi/86400

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if (MyID == pE3d[k]) {
			if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
				for (int n = 0; n < Np; n++) {
					MRHS[k * Np + n] = hv[k * Np + n] * sin(Latitude * Deg2Rad) * omega2;
					MRHS[K * Np + k * Np + n] = -hu[k * Np + n] * sin(Latitude * Deg2Rad) * omega2;
				}

				Add(OutputRHS + k * Np, MRHS + k * Np, OutputRHS + k * Np, Np);
				Add(OutputRHS + Np * K + k * Np, MRHS + Np * K + k * Np, OutputRHS + Np * K + k * Np, Np);
			}
		}
	}

	free(MRHS);
}
