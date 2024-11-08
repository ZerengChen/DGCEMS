#include "NdgSourceTermSolver3d.h"
#include "NdgMath.h"

extern double gra;
extern double Hcrit;
extern int Nvar;

SWETopographySourceTerm3d::SWETopographySourceTerm3d()
{
}

SWETopographySourceTerm3d::~SWETopographySourceTerm3d()
{
}

extern signed char *Status3d;

void SWETopographySourceTerm3d::EvaluateTopographySourceTerm(double *fphys_, double *frhs_, int*pE3d, int MyID)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	//int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	int K2d = *(meshunion->mesh2d_p->K2d);
	double *OutputRHS = frhs_;
	double *MRHS = (double *)malloc(Np*K*2 * sizeof(double));
	memset(MRHS, 0.0, Np * K * 2 * sizeof(double));
	double *fphys = fphys_;
	double *bx = fphys + 7 * K * Np;
	double *by = fphys + 8 * K * Np;
	double *eta_ = fphys + 6 * K * Np;
	signed char *status = meshunion->mesh2d_p->status;
	int NLayer = *(meshunion->Nlayer);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if (MyID == pE3d[k]) {
			if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
				for (int n = 0; n < Np; n++) {
					MRHS[k * Np + n] = -gra * eta_[k * Np + n] * bx[k * Np + n];
					MRHS[K * Np + k * Np + n] = -gra * eta_[k * Np + n] * by[k * Np + n];
				}
				Add(OutputRHS + k * Np, MRHS + k * Np, OutputRHS + k * Np, Np);

				Add(OutputRHS + Np * K + k * Np, MRHS + Np * K + k * Np, OutputRHS + Np * K + k * Np, Np);
			}

			else if ((NdgRegionType)Status3d[k] == NdgRegionPartialWetDamBreak) {
				for (int n = 0; n < Np; n++) {
					MRHS[k * Np + n] = -gra * eta_[k * Np + n] * bx[k * Np + n];
					MRHS[K * Np + k * Np + n] = -gra * eta_[k * Np + n] * by[k * Np + n];
				}
				Add(OutputRHS + k * Np, MRHS + k * Np, OutputRHS + k * Np, Np);

				Add(OutputRHS + Np * K + k * Np, MRHS + Np * K + k * Np, OutputRHS + Np * K + k * Np, Np);
			}

			else {
				continue;
			}
		}
	}

	free(MRHS);
}
