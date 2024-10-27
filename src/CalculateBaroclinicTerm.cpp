#include "NdgMath.h"
#include "NdgSWE3D.h"
#include "NdgMemory.h"
#include "CalculateBaroclinicTerm.h"
#include "eqstate.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double *BaroclinicPDPX, *BaroclinicPDPY, *BaroclinicPRHOPX, *BaroclinicPRHOPY, \
*BaroclinicPRHOPS,*BaroclinicInXPartOne, *BaroclinicInXPartTwo, *BaroclinicInYPartOne, \
*BaroclinicInYPartTwo, *BaroclinicInXTempRHS,*BaroclinicInYTempRHS, *Baroclinicfmod, \
*BaroclinicBotEfm, *BaroclinicBotEfp, *BaroclinicBotEFluxM,*BaroclinicBotEFluxP, \
*BaroclinicBotEFluxS, *BaroclinicIEfm, *BaroclinicIEfp,*BaroclinicIEfluxM, \
*BaroclinicIEfluxP, *BaroclinicIEfluxS, *BaroclinicERHS, *BaroclinicTempFacialIntegral,\
*BaroclinicTempVolumeIntegral, *BaroclinicBEfm, *BaroclinicBEfp, *BaroclinicBEfluxM,\
*BaroclinicBEfluxS;

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

void EvaluateBaroclinicTerm(double *fphys, double *frhs, double *fext_) {
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);

	double *hurhs = frhs;
	double *hvrhs = frhs + Np * K;
	double *h = fphys + 3 * Np * K;
	double *rho = fphys + 13 * Np * K;

	double *fext3d = fext_;	

	int NLayer = *(meshunion->Nlayer);
	double *z = meshunion->z;
	double *Jz = meshunion ->Jz;

	int Nface = *(meshunion->cell_p->Nface);

	int IENe = *(meshunion->inneredge_p->Ne);
	int IENfp = *(meshunion->inneredge_p->Nfp);

	int BENe = *(meshunion->boundaryedge_p->Ne);
	int BENfp = *(meshunion->boundaryedge_p->Nfp);

	int Npz = *(meshunion->cell_p ->Npz);
	int Nph = *(meshunion->cell_p->Nph);

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

#ifdef _BAROCLINIC
	double *VintU = meshunion->cell_p->VintU;//垂向从顶至底的积分矩阵
#endif
	double *V = meshunion->cell_p ->V;

	int K2d = K/ NLayer;
	int BotENe = *(meshunion->bottomedge_p->Ne);

	GetFirstOrderPartialDerivativeInHorizontalDirection(BaroclinicPDPX, BaroclinicPDPY, BaroclinicPRHOPX, BaroclinicPRHOPY, h, rho, fext3d);
	GetFirstOrderPartialDerivativeInVerticalDirection(BaroclinicPRHOPS, rho);

	double *InvV = (double*)malloc(Np*Np * sizeof(double));
	memcpy(InvV, V, Np*Np * sizeof(double));
	MatrixInverse(InvV, Np);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			for (int p = 0; p < Np; p++) {
				DotProduct(BaroclinicInXPartOne + k * Np, h + k * Np, BaroclinicPRHOPX + k * Np, Np);
				DotProduct(BaroclinicInXPartTwo + k * Np, z + k * Np, BaroclinicPDPX + k * Np, Np);
				DotProduct(BaroclinicInXPartTwo + k * Np, BaroclinicInXPartTwo + k * Np, BaroclinicPRHOPS + k * Np, Np);
				Minus(BaroclinicInXPartOne + k * Np, BaroclinicInXPartOne + k * Np, BaroclinicInXPartTwo + k * Np, Np);
				MultiplyByConstant(BaroclinicInXPartOne + k * Np, BaroclinicInXPartOne + k * Np, -1 * gra / rho0, Np);

				DotProduct(BaroclinicInYPartOne + k * Np, h + k * Np, BaroclinicPRHOPY + k * Np, Np);
				DotProduct(BaroclinicInYPartTwo + k * Np, z + k * Np, BaroclinicPDPY + k * Np, Np);
				DotProduct(BaroclinicInYPartTwo + k * Np, BaroclinicInYPartTwo + k * Np, BaroclinicPRHOPS + k * Np, Np);
				Minus(BaroclinicInYPartOne + k * Np, BaroclinicInYPartOne + k * Np, BaroclinicInYPartTwo + k * Np, Np);
				MultiplyByConstant(BaroclinicInYPartOne + k * Np, BaroclinicInYPartOne + k * Np, -1 * gra / rho0, Np);
			}
		}
	}

#ifdef _BAROCLINIC
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if ((NdgRegionType)Status2d[k] != NdgRegionWet) {
			continue;
		}
		else {
			VerticalIntegralFromSurface(BaroclinicInXTempRHS + k * NLayer*Np, BaroclinicInXPartOne + k * NLayer*Np, Jz + k * NLayer*Np, Baroclinicfmod + k * Np, NLayer, Np, InvV, Nph, Npz, VintU);
			VerticalIntegralFromSurface(BaroclinicInYTempRHS + k * NLayer*Np, BaroclinicInYPartOne + k * NLayer*Np, Jz + k * NLayer*Np, Baroclinicfmod + k * Np, NLayer, Np, InvV, Nph, Npz, VintU);
		}
	}
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			for (int p = 0; p < Np; p++)
			{
				BaroclinicInXTempRHS[k*Np + p] = h[k*Np + p] * BaroclinicInXTempRHS[k*Np + p];
				hurhs[k*Np + p] = hurhs[k*Np + p] + BaroclinicInXTempRHS[k*Np + p];
				BaroclinicInYTempRHS[k*Np + p] = h[k*Np + p] * BaroclinicInYTempRHS[k*Np + p];
				hvrhs[k*Np + p] = hvrhs[k*Np + p] + BaroclinicInYTempRHS[k*Np + p];
			}
		}
	}
	free(InvV);
}

void GetFirstOrderPartialDerivativeInVerticalDirection(double *BaroclinicPRHOPS, double *rho) {
	double *tz = meshunion->tz;
	double *J = meshunion->J;
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	double *Dt = meshunion->cell_p->Dt;
	int Nface = *(meshunion->cell_p->Nface);
	double *invM = meshunion->cell_p->invM;

	/*For bottom edge object*/
	int BotENe = *(meshunion->bottomedge_p->Ne);
	int BotENfp = *(meshunion->bottomedge_p->Nfp);
	double *BotEMb = meshunion->bottomedge_p->M;
	double *BotEJs = meshunion->bottomedge_p->Js;
	double *BotEnz = meshunion->bottomedge_p->nz;
	double *BotEFToE = meshunion->bottomedge_p->FToE;
	double *BotEFToF = meshunion->bottomedge_p->FToF;
	double *BotEFToN1 = meshunion->bottomedge_p->FToN1;
	double *BotEFToN2 = meshunion->bottomedge_p->FToN2;

	memset(BaroclinicERHS, 0, Np*K * Nface * sizeof(double));

	double *rhoM = BaroclinicBotEfm, *rhoP = BaroclinicBotEfp, *rhoFluxM = BaroclinicBotEFluxM, \
		*rhoFluxP = BaroclinicBotEFluxP, *rhoFluxS = BaroclinicBotEFluxS;

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
	for (int face = 0; face < BotENe; face++) {
		/*Fetch variable BotEfm and BotEfp first*/
		int adjacentE = (int)BotEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			FetchInnerEdgeFacialValue(rhoM + face * BotENfp, rhoP + face * BotENfp, rho, BotEFToE + 2 * face, \
				BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);

			EvaluateFaceSurfFlux(rhoFluxM + face * BotENfp, rhoM + face * BotENfp, BotEnz + face * BotENfp, BotENfp);

			EvaluateFaceSurfFlux(rhoFluxP + face * BotENfp, rhoP + face * BotENfp, BotEnz + face * BotENfp, BotENfp);

			EvaluateFaceNumFlux_Central(rhoFluxS + face * BotENfp, rhoM + face * BotENfp, rhoP + face * BotENfp, BotEnz + face * BotENfp, BotENfp);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++) {
		int adjacentE = (int)BotEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			for (int field = 0; field < 1; field++) {
				StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, BaroclinicBotEFluxM + field * BotENe*BotENfp, \
					BaroclinicBotEFluxP + field * BotENe*BotENfp, BaroclinicBotEFluxS + field * BotENe*BotENfp, BotEJs, BotEMb, BaroclinicERHS + field * Np*K*Nface);
			}
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			for (int field = 0; field < 1; field++) {
				for (int face = 1; face < Nface; face++) {
					Add(BaroclinicERHS + field * Np*K*Nface + k * Np, BaroclinicERHS + field * Np*K*Nface + k * Np, BaroclinicERHS + field * Np*K*Nface + face * Np*K + k * Np, Np);
				}
			}
		}
		else {
			continue;
		}
	}

	int np = Np;
	int oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			for (int field = 0; field < 1; field++) {
				MultiEdgeContributionByLiftOperator(BaroclinicERHS + field * Np*K*Nface + k * Np, BaroclinicTempFacialIntegral + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
			}
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/*$\bold{t_z}\cdot (Dt*rho)$*/
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			GetVolumnIntegral1d(BaroclinicPRHOPS + k * Np, &np, &oneI, &np, &one, \
				Dt, &np, rho + k * Np, &np, &zero, &np, tz + k * Np);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			Minus(BaroclinicPRHOPS + k * Np, BaroclinicPRHOPS + k * Np, BaroclinicERHS + k * Np, Np);
		}
		else {
			continue;
		}
	}

}

void GetFirstOrderPartialDerivativeInHorizontalDirection(double *BaroclinicPDPX, double *BaroclinicPDPY, \
	double *BaroclinicPRHOPX, double *BaroclinicPRHOPY, double *h, double *rho, double *fext)
{
	double *rx = meshunion->rx;
	double *sx = meshunion->sx;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;
	double *J = meshunion->J;
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);

	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	double *Dt = meshunion->cell_p->Dt;
	int Nface = *(meshunion->cell_p->Nface);
	double *invM = meshunion->cell_p->invM;

	int IENe = *(meshunion->inneredge_p->Ne);
	int IENfp = *(meshunion->inneredge_p->Nfp);
	double *IEMb = meshunion->inneredge_p->M;
	double *IEJs = meshunion->inneredge_p->Js;
	double *IEnx = meshunion->inneredge_p->nx;
	double *IEny = meshunion->inneredge_p->ny;
	double *IEFToE = meshunion->inneredge_p->FToE;
	double *IEFToF = meshunion->inneredge_p->FToF;
	double *IEFToN1 = meshunion->inneredge_p->FToN1;
	double *IEFToN2 = meshunion->inneredge_p->FToN2;
	double *IELAV = meshunion->inneredge_p ->LAV ;

	int BENe = *(meshunion->boundaryedge_p->Ne);
	int BENfp = *(meshunion->boundaryedge_p->Nfp);
	double *BEMb = meshunion->boundaryedge_p->M;
	double *BEJs = meshunion->boundaryedge_p->Js;
	double *BEnx = meshunion->boundaryedge_p->nx;
	double *BEny = meshunion->boundaryedge_p->ny;
	double *BEFToE = meshunion->boundaryedge_p->FToE;
	double *BEFToF = meshunion->boundaryedge_p->FToF;
	double *BEFToN1 = meshunion->boundaryedge_p->FToN1;
	double *BELAV = meshunion->boundaryedge_p ->LAV ;

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

	double *ftype = meshunion->boundaryedge_p->ftype;

	double *hM = BaroclinicIEfm, *rhoM = BaroclinicIEfm + IENe*IENfp;
	double *hP = BaroclinicIEfp, *rhoP = BaroclinicIEfp + IENe*IENfp;
	double *HIEfluxMx = BaroclinicIEfluxM, *rhoIEfluxMx = BaroclinicIEfluxM + IENe*IENfp, \
		*HIEfluxMy = BaroclinicIEfluxM + 2*IENe*IENfp, *rhoIEfluxMy = BaroclinicIEfluxM + 3*IENe*IENfp;
	double *HIEfluxPx = BaroclinicIEfluxP, *rhoIEfluxPx = BaroclinicIEfluxP + IENe*IENfp, \
		*HIEfluxPy = BaroclinicIEfluxP + 2 * IENe*IENfp, *rhoIEfluxPy = BaroclinicIEfluxP + 3 * IENe*IENfp;
	double *HIEfluxSx = BaroclinicIEfluxS, *rhoIEfluxSx = BaroclinicIEfluxS + IENe*IENfp, \
	    *HIEfluxSy = BaroclinicIEfluxS + 2 * IENe*IENfp, *rhoIEfluxSy = BaroclinicIEfluxS + 3 * IENe*IENfp;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		/*Fetch variable IEfm and IEfp first*/
		int adjacentE = (int)IEFToE[2 * face];
		int adjacentE2 = (int)IEFToE[2 * face + 1];
		if ((NdgRegionType)Status3d[adjacentE - 1] != NdgRegionWet || (NdgRegionType)Status3d[adjacentE2 - 1] != NdgRegionWet) {
			continue;
		}
		else {
			FetchInnerEdgeFacialValue(hM + face * IENfp, hP + face * IENfp, h, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
			FetchInnerEdgeFacialValue(rhoM + face * IENfp, rhoP + face * IENfp, rho, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);

			EvaluateFaceSurfFlux(HIEfluxMx + face * IENfp, hM + face * IENfp, IEnx + face * IENfp, IENfp);
			EvaluateFaceSurfFlux(HIEfluxPx + face * IENfp, hP + face * IENfp, IEnx + face * IENfp, IENfp);

			EvaluateFaceSurfFlux(HIEfluxMy + face * IENfp, hM + face * IENfp, IEny + face * IENfp, IENfp);
			EvaluateFaceSurfFlux(HIEfluxPy + face * IENfp, hP + face * IENfp, IEny + face * IENfp, IENfp);

			EvaluateFaceSurfFlux(rhoIEfluxMx + face * IENfp, rhoM + face * IENfp, IEnx + face * IENfp, IENfp);
			EvaluateFaceSurfFlux(rhoIEfluxPx + face * IENfp, rhoP + face * IENfp, IEnx + face * IENfp, IENfp);

			EvaluateFaceSurfFlux(rhoIEfluxMy + face * IENfp, rhoM + face * IENfp, IEny + face * IENfp, IENfp);
			EvaluateFaceSurfFlux(rhoIEfluxPy + face * IENfp, rhoP + face * IENfp, IEny + face * IENfp, IENfp);

			EvaluateFaceNumFlux_Central(HIEfluxSx + face * IENfp, hM + face * IENfp, hP + face * IENfp, IEnx + face * IENfp, IENfp);
			EvaluateFaceNumFlux_Central(HIEfluxSy + face * IENfp, hM + face * IENfp, hP + face * IENfp, IEny + face * IENfp, IENfp);
			EvaluateFaceNumFlux_Central(rhoIEfluxSx + face * IENfp, rhoM + face * IENfp, rhoP + face * IENfp, IEnx + face * IENfp, IENfp);
			EvaluateFaceNumFlux_Central(rhoIEfluxSy + face * IENfp, rhoM + face * IENfp, rhoP + face * IENfp, IEny + face * IENfp, IENfp);
		}
	}

	hM = BaroclinicBEfm, rhoM = BaroclinicBEfm + BENe*BENfp, hP = BaroclinicBEfp;
	// The second part space store hT first, then the density rho next.
	double *hTrhoP = BaroclinicBEfp + BENe*BENfp, *hSp = BaroclinicBEfp + 2*BENe*BENfp;
	double *HBEfluxMx = BaroclinicBEfluxM, *rhoBEfluxMx = BaroclinicBEfluxM + BENe*BENfp, \
		*HBEfluxMy = BaroclinicBEfluxM + 2 * BENe*BENfp, *rhoBEfluxMy = BaroclinicBEfluxM + 3 * BENe*BENfp,\
		*HBEfluxSx = BaroclinicBEfluxS, *rhoBEfluxSx = BaroclinicBEfluxS + BENe*BENfp, \
		*HBEfluxSy = BaroclinicBEfluxS + 2 * BENe*BENfp, *rhoBEfluxSy = BaroclinicBEfluxS + 3 * BENe*BENfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++) {
		NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			FetchBoundaryEdgeFacialValue(hM + face * BENfp, h, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(rhoM + face * BENfp, rho, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
			ImposeBcsForRhoAndH(hTrhoP + face * BENfp, hSp + face * BENfp, hM + face * BENfp, hP + face * BENfp, rhoM + face * BENfp, \
				fext + face * BENfp, BENfp, BENe, type);
			/*Fetch variable IEfm and IEfp first*/
			EvaluateFaceSurfFlux(HBEfluxMx + face * BENfp, hM + face * BENfp, BEnx + face * BENfp, BENfp);
			EvaluateFaceSurfFlux(HBEfluxMy + face * BENfp, hM + face * BENfp, BEny + face * BENfp, BENfp);
			EvaluateFaceSurfFlux(rhoBEfluxMx + face * BENfp, rhoM + face * BENfp, BEnx + face * BENfp, BENfp);
			EvaluateFaceSurfFlux(rhoBEfluxMy + face * BENfp, rhoM + face * BENfp, BEny + face * BENfp, BENfp);
			EvaluateFaceNumFlux_Central(HBEfluxSx + face * BENfp, hM + face * BENfp, hP + face * BENfp, BEnx + face * BENfp, BENfp);
			EvaluateFaceNumFlux_Central(HBEfluxSy + face * BENfp, hM + face * BENfp, hP + face * BENfp, BEny + face * BENfp, BENfp);
			EvaluateFaceNumFlux_Central(rhoBEfluxSx + face * BENfp, rhoM + face * BENfp, hTrhoP + face * BENfp, BEnx + face * BENfp, BENfp);
			EvaluateFaceNumFlux_Central(rhoBEfluxSy + face * BENfp, rhoM + face * BENfp, hTrhoP + face * BENfp, BEny + face * BENfp, BENfp);
		}
	}

	memset(BaroclinicERHS, 0, Np * K * 4 * (Nface - 2) * sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		int adjacentE = (int)IEFToE[2 * face];
		int adjacentE2 = (int)IEFToE[2 * face + 1];
		if ((NdgRegionType)Status3d[adjacentE - 1] != NdgRegionWet || (NdgRegionType)Status3d[adjacentE2 - 1] != NdgRegionWet) {
			continue;
		}
		else {
			for (int field = 0; field < 4; field++) {
				StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, BaroclinicIEfluxM + field * IENe*IENfp, \
					BaroclinicIEfluxP + field * IENe*IENfp, BaroclinicIEfluxS + field * IENe*IENfp, IEJs, IEMb, BaroclinicERHS + field * Np*K*(Nface - 2));
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++) {
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < 4; field++) {
				StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, BENfp, BEFToN1, BaroclinicBEfluxM + field * BENe*BENfp, \
					BaroclinicBEfluxS + field * BENe*BENfp, BEJs, BEMb, BaroclinicERHS + field * Np*K*(Nface - 2));
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			for (int field = 0; field < 4; field++) {
				for (int face = 1; face < Nface - 2; face++) {
					Add(BaroclinicERHS + field * Np*K*(Nface - 2) + k * Np, BaroclinicERHS + field * Np*K*(Nface - 2) + k * Np, BaroclinicERHS + field * Np*K*(Nface - 2) + face * Np*K + k * Np, Np);
				}
			}
		}
		else {
			continue;
		}
	}

	int np = Np;
	int oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			for (int field = 0; field < 4; field++) {
				MultiEdgeContributionByLiftOperator(BaroclinicERHS + field * Np*K*(Nface - 2) + k * Np, BaroclinicTempFacialIntegral + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
			}
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/*$\bold{r_x}\cdot (Dr*h)+\bold{s_x}\cdot (Ds*h)$*/
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			GetVolumnIntegral2d(BaroclinicPDPX + k * Np, BaroclinicTempVolumeIntegral + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, h + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);

			GetVolumnIntegral2d(BaroclinicPDPY + k * Np, BaroclinicTempVolumeIntegral + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, h + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);

			GetVolumnIntegral2d(BaroclinicPRHOPX + k * Np, BaroclinicTempVolumeIntegral + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, rho + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);

			GetVolumnIntegral2d(BaroclinicPRHOPY + k * Np, BaroclinicTempVolumeIntegral + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, rho + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			Minus(BaroclinicPDPX + k * Np, BaroclinicPDPX + k * Np, BaroclinicERHS + k * Np, Np);

			Minus(BaroclinicPRHOPX + k * Np, BaroclinicPRHOPX + k * Np, BaroclinicERHS + Np * K*(Nface - 2) + k * Np, Np);

			Minus(BaroclinicPDPY + k * Np, BaroclinicPDPY + k * Np, BaroclinicERHS + 2 * Np*K*(Nface - 2) + k * Np, Np);

			Minus(BaroclinicPRHOPY + k * Np, BaroclinicPRHOPY + k * Np, BaroclinicERHS + 3 * Np*K*(Nface - 2) + k * Np, Np);
		}
		else {
			continue;
		}
	}

}

/*Note: Only part of the boundary type is implemented at present*/
void ImposeBcsForRhoAndH(double *srcAndDest, double *hS, double *hM, double *hP, double *rhoM, double *fext, \
	int BENfp, int BENe, NdgEdgeType type) 
{
	if (type == NdgEdgeZeroGrad) {
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
	}
	else if (type == NdgEdgeClampedVel) {
		//T
		//DotCriticalDivide(srcAndDest, fext + 3*BENe*BENfp, &Hcrit, hM, BENfp);
		//S
		//DotCriticalDivide(hS, fext + 4 * BENe*BENfp, &Hcrit, hM, BENfp);
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
		//if (!strcmp(EosType, "Jackett05")) {
		//	for (int i = 0; i < BENfp; i++) {
		//		hP[i] = hM[i];
		//		EosByFeistel(srcAndDest + i, fmax(*(srcAndDest + i), 0.0), fmax(*(hS + i), 0.0));
		//	}
		//}
		//else if (!strcmp(EosType, "UNESCO83")) {
		//	for (int i = 0; i < BENfp; i++) {
		//		hP[i] = hM[i];
		//		EosByUNESCO(srcAndDest + i, fmax(*(srcAndDest + i), 0.0), fmax(*(hS + i), 0.0));
		//	}
		//}
		//else if (!strcmp(EosType, "Linear")) {
			//for (int i = 0; i < BENfp; i++) {
			//	hP[i] = hM[i];
			//	EosByLinear(srcAndDest + i, fmax(*(srcAndDest + i), 0.0), fmax(*(hS + i), 0.0), rho0, T0, S0, alphaT, betaS);
			//}
		//}
		//else {
		//	printf("Equation of state(EOS) needs to be pointed for this part!\n");
		//}
	}
	else if (type == NdgEdgeClampedDepth) {
		//T. The exterior value is stored as hu, hv, h, ht, and hs
		//DotCriticalDivide(srcAndDest, fext + 3 * BENe*BENfp, &Hcrit, fext + 2 * BENe*BENfp, BENfp);
		//S
		//DotCriticalDivide(hS, fext + 4 * BENe*BENfp, &Hcrit, fext + 2 * BENe*BENfp, BENfp);
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
		//if (!strcmp(EosType, "Jackett05")) {
		//	for (int i = 0; i < BENfp; i++) {
		//		hP[i] = fext[2 * BENe*BENfp +i];
		//		EosByFeistel(srcAndDest + i, fmax(*(srcAndDest + i), 0.0), fmax(*(hS + i), 0.0));
		//	}
		//}
		//else if (!strcmp(EosType, "UNESCO83")) {
		//	for (int i = 0; i < BENfp; i++) {
		//		hP[i] = fext[2 * BENe*BENfp + i];
		//		EosByUNESCO(srcAndDest + i, fmax(*(srcAndDest + i), 0.0), fmax(*(hS + i), 0.0));
		//	}
		//}
		//else if (!strcmp(EosType, "Linear")) {
			//for (int i = 0; i < BENfp; i++) {
			//	hP[i] = fext[2 * BENe*BENfp + i];
			//	EosByLinear(srcAndDest + i, fmax(*(srcAndDest + i), 0.0), fmax(*(hS + i), 0.0), rho0, T0, S0, alphaT, betaS);
			//}
		//}
		//else {
		//	printf("Equation of state(EOS) needs to be pointed for this part!\n");
		//}
	}
	else if (type == NdgEdgeClamped) {
		//T
		//DotCriticalDivide(srcAndDest, fext + 3 * BENe*BENfp, &Hcrit, fext + 2*BENfp*BENe, BENfp);
		//S
		//DotCriticalDivide(hS, fext + 4 * BENe*BENfp, &Hcrit, fext + 2 * BENfp*BENe, BENfp);
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
		//if (!strcmp(EosType, "Jackett05")) {
		//	for (int i = 0; i < BENfp; i++) {
		//		hP[i] = fext[2 * BENfp*BENe + i];
		//		EosByFeistel(srcAndDest + i, fmax(*(srcAndDest + i), 0.0), fmax(*(hS + i), 0.0));
		//	}
		//}
		//else if (!strcmp(EosType, "UNESCO83")) {
		//	for (int i = 0; i < BENfp; i++) {
		//		hP[i] = fext[2 * BENfp*BENe + i];
		//		EosByUNESCO(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0));
		//	}
		//}
		//else if (!strcmp(EosType, "Linear")) {
		//	for (int i = 0; i < BENfp; i++) {
		//		hP[i] = fext[2 * BENfp*BENe + i];
		//		EosByLinear(srcAndDest + i, max(*(srcAndDest + i), 0.0), max(*(hS + i), 0.0), rho0, T0, S0, alphaT, betaS);
		//	}
		//}
		//else {
		//	printf("Equation of state(EOS) needs to be pointed for this part!\n");
		//}
	}
	else if (type == NdgEdgeSlipWall) {
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
	}
	else if (type == NdgEdgeNonSlipWall) {
		for (int i = 0; i < BENfp; i++) {
			srcAndDest[i] = rhoM[i];
			hP[i] = hM[i];
		}
	}
	else if (type == NdgEdgeFlather) {
	    for (int i = 0; i < BENfp; i++) {
		    srcAndDest[i] = rhoM[i];
		    hP[i] = hM[i];
	    }
	}
	else if (type == NdgEdgeNonLinearFlather) {
	    for (int i = 0; i < BENfp; i++) {
		    srcAndDest[i] = rhoM[i];
		    hP[i] = hM[i];
	    }
	}
	else if (type == NdgEdgeNonLinearFlatherFlow) {
	    for (int i = 0; i < BENfp; i++) {
		    srcAndDest[i] = rhoM[i];
		    hP[i] = hM[i];
	    }
	}
}

void EvaluateFaceSurfFlux(double *dest, double *source, double *vector, int size) {
	for (int i = 0; i < size; i++)
	{
		dest[i] = source[i] * vector[i];
	}
}

void EvaluateFaceNumFlux_Central(double *dest, double *sourcefm, double *sourcefp, double *vector, int size) {
	for (int i = 0; i < size; i++)
	{
		dest[i] = (sourcefm[i] + sourcefp[i]) / 2.0 * vector[i];
	}
}