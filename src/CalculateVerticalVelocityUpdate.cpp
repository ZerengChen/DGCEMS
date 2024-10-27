
#include "NdgMath.h"
#include "NdgSWE.h"
#include "NdgSWE3D.h"
//#include "NdgMemory.h"
#include "CalculateVerticalVelocityUpdate.h"

CalculateVerticalVelocity::CalculateVerticalVelocity()
{
}

CalculateVerticalVelocity::~CalculateVerticalVelocity()
{
}

extern double *UpdatedVSrhs2d, *UpdatedVSIEfm2d, *UpdatedVSIEfp2d, *UpdatedVSIEFluxM2d, \
*UpdatedVSIEFluxP2d, *UpdatedVSIEFluxS2d, *UpdatedVSERHS2d, *UpdatedVSVolumeIntegralX, \
*UpdatedVSTempVolumeIntegralX, *UpdatedVSVolumeIntegralY, *UpdatedVSTempVolumeIntegralY, \
*UpdatedVSBEfm2d, *UpdatedVSBEzM2d, *UpdatedVSBEfp2d, *UpdatedVSBEzP2d, *UpdatedVSBEFluxS2d, \
*UpdatedVSBEFluxM2d, *UpdatedVSTempFacialIntegral, *UpdatedVSfield2d, *UpdatedVSrhs3d, \
*UpdatedVSIEfm3d, *UpdatedVSIEfp3d, *UpdatedVSIEFluxM3d, *UpdatedVSIEFluxP3d, *UpdatedVSIEFluxS3d, \
*UpdatedVSERHS3d, *UpdatedVSVolumeIntegralX3d, *UpdatedVSTempVolumeIntegralX3d, \
*UpdatedVSVolumeIntegralY3d, *UpdatedVSTempVolumeIntegralY3d, *UpdatedVSBEfm3d, \
*UpdatedVSBEzM3d, *UpdatedVSBEfp3d, *UpdatedVSBEzP3d, *UpdatedVSBEFluxS3d, *UpdatedVSBEFluxM3d, \
*UpdatedVSTempFacialIntegral3d, *UpdatedVSIEfmod, *UpdatedVSBEfmod, *Updatedfmod;

//extern char *UpdatedVertVelocityInitialized;
extern double Hcrit;
extern double gra;
extern signed char *Status3d;
//void MyExit()
//{
//	if (!strcmp("True", UpdatedVertVelocityInitialized)){
//		UpdatedVertVelocitySolverMemoryDeAllocation();
//		UpdatedVertVelocityInitialized = "False";
//	}
//	return;
//}


void CalculateVerticalVelocity::EvaluateVerticalVelocity(double *fphys2d_, double *fphys_, double *fext2d_, double *fext3d_)
{
	//mexAtExit(&MyExit);
	/*Properties contained in mesh2d*/
	double *rx2d = meshunion->mesh2d_p->rx2d;
	int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	int K2d = *(meshunion->mesh2d_p->K2d);
	double *sx2d = meshunion->mesh2d_p->sx2d;
	double *ry2d = meshunion->mesh2d_p->ry2d;
	double *sy2d = meshunion->mesh2d_p->sy2d;
	double *J2d = meshunion->mesh2d_p->J2d;

	/*Properties contained in mesh3d*/
	double *rx3d = meshunion->rx;
	int Np3d = *(meshunion->cell_p->Np);
	int K3d = *(meshunion->K);
	double *sx3d = meshunion->sx;
	double *ry3d = meshunion->ry;
	double *sy3d = meshunion->sy;
	double *J3d = meshunion->J;
	int NLayer = *(meshunion->Nlayer);
	double *Jz = meshunion->Jz;

	/*Properties contained in two dimensional inner edge*/
	int IENe2d = *(meshunion->mesh2d_p->mesh2dinneredge_p->Ne2d);
	int IENfp2d = *(meshunion->mesh2d_p->mesh2dinneredge_p->Nfp2d);
	double *IEFToE2d = meshunion->mesh2d_p->mesh2dinneredge_p->FToE2d;
	double *IEFToF2d = meshunion->mesh2d_p->mesh2dinneredge_p->FToF2d;
	double *IEFToN12d = meshunion->mesh2d_p->mesh2dinneredge_p->FToN12d;
	double *IEFToN22d = meshunion->mesh2d_p->mesh2dinneredge_p->FToN22d;
	double *IEnx2d = meshunion->mesh2d_p->mesh2dinneredge_p->nx2d;
	double *IEny2d = meshunion->mesh2d_p->mesh2dinneredge_p->ny2d;
	double *IEJs2d = meshunion->mesh2d_p->mesh2dinneredge_p->Js2d;
	double *IEMb2d = meshunion->mesh2d_p->mesh2dinneredge_p->M2d;

	/*Properties contained in two dimensional boundary edge*/
	int BENe2d = *(meshunion->mesh2d_p->mesh2dboundaryedge_p->Ne2d);
	int BENfp2d = *(meshunion->mesh2d_p->mesh2dboundaryedge_p->Nfp2d);
	double *BEFToE2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->FToE2d;
	double *BEFToF2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->FToF2d;
	double *BEFToN12d = meshunion->mesh2d_p->mesh2dboundaryedge_p->FToN12d;
	double *BEnx2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->nx2d;
	double *BEny2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->ny2d;
	double *BEJs2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->Js2d;
	double *BEMb2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->M2d;

	/*Properties contained in three dimensional inner edge*/
	int IENe3d = *(meshunion->inneredge_p->Ne);
	int IENfp3d = *(meshunion->inneredge_p->Nfp);
	double *IEFToE3d = meshunion->inneredge_p->FToE;
	double *IEFToF3d = meshunion->inneredge_p->FToF;
	double *IEFToN13d = meshunion->inneredge_p->FToN1;
	double *IEFToN23d = meshunion->inneredge_p->FToN2;
	double *IEnx3d = meshunion->inneredge_p->nx;
	double *IEny3d = meshunion->inneredge_p->ny;
	double *IEMb3d = meshunion->inneredge_p->M;
	double *IEJs3d = meshunion->inneredge_p->Js;
	double *IEJz3d = meshunion->inneredge_p->Jz;
	double *V1d = meshunion->inneredge_p->V1d;
	double *V2d = meshunion->inneredge_p->V2d;

	double *InvV2d = (double *)malloc(IENfp3d*IENfp3d * sizeof(double));
	memcpy(InvV2d, V2d, IENfp3d*IENfp3d * sizeof(double));
	MatrixInverse(InvV2d, IENfp3d);

	/*Properties contained in three dimensional boundary edge*/
	int BENe3d = *(meshunion->boundaryedge_p->Ne);
	int BENfp3d = *(meshunion->boundaryedge_p->Nfp);
	double *BEFToE3d = meshunion->boundaryedge_p->FToE;
	double *BEFToF3d = meshunion->boundaryedge_p->FToF;
	double *BEFToN13d = meshunion->boundaryedge_p->FToN1;
	double *BEnx3d = meshunion->boundaryedge_p->nx;
	double *BEny3d = meshunion->boundaryedge_p->ny;
	double *BEMb3d = meshunion->boundaryedge_p->M;
    double *BEJs3d = meshunion->boundaryedge_p->Js;
	double *BEJz3d = meshunion->boundaryedge_p->Jz;

	/*Data contained in two dimensional physical field and three dimensional physical field*/
	double *fphys2d = fphys2d_;
	double *h2d = fphys2d;
	double *hu2d = fphys2d + Np2d*K2d;
	double *hv2d = fphys2d + 2 * Np2d*K2d;
	double *z2d = fphys2d + 3 * Np2d*K2d;
	double *fphys3d = fphys_;
	double *hu3d = fphys3d;
	double *hv3d = fphys3d + Np3d * K3d;
	double *h3d = fphys3d + 3 * Np3d*K3d;
	double *z3d = fphys3d + 5 * Np3d*K3d;
	
	/*Data contained in two-dimensional standard cell*/
	double *Dr2d = meshunion->mesh2d_p->mesh2dcell_p->Dr2d;
	double *Ds2d = meshunion->mesh2d_p->mesh2dcell_p->Ds2d;
	int Nface = *(meshunion->mesh2d_p->mesh2dcell_p->Nface2d);
	double *invM2d = meshunion->mesh2d_p->mesh2dcell_p->invM2d;

	/*Data contained in three-dimensional standard cell*/
	double *Dr3d = meshunion->cell_p->Dr;
	double *Ds3d = meshunion->cell_p->Ds;
	double *invM3d = meshunion->cell_p->invM;
	/*Approximation order in vertical direction*/
	int Nz = *(meshunion->cell_p->Nz);
	int Npz = *(meshunion->cell_p->Npz);
	double *Vint = meshunion->cell_p->Vint;
	double *V3d = meshunion->cell_p->V;

	double *InvV3d = (double *)malloc(Np3d*Np3d*sizeof(double));
	memcpy(InvV3d, V3d, Np3d*Np3d*sizeof(double));
	MatrixInverse(InvV3d, Np3d);


	/*The two-dimensional external value*/
	double *fext2d = fext2d_;
//	double *BEhE2d = fext2d, *BEhuE2d = fext2d + BENe2d * BENfp2d, *BEhvE2d = fext2d + 2 * BENe2d * BENfp2d;
	
	/*The three-dimensional external value*/	
	double *fext3d = fext3d_;

	double *ftype2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->ftype2d;
	double *ftype3d = meshunion->boundaryedge_p->ftype;

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
	//int adjacentE = 0; //store the element number of the current face
	/***  Wet and Dry  ***/

	double *VerticalVelocity = fphys3d + K3d*Np3d*2;

	/**********************************************************  Three Dimensional Part  *******************************************************************************/

	memset(UpdatedVSrhs3d, 0, Np3d*K3d*sizeof(double));
	double *IEhuM3d = UpdatedVSIEfm3d, *IEhvM3d = UpdatedVSIEfm3d + IENfp3d * IENe3d, \
		*IEhM3d = UpdatedVSIEfm3d + 2 * IENfp3d*IENe3d;
	double *IEhuP3d = UpdatedVSIEfp3d, *IEhvP3d = UpdatedVSIEfp3d + IENfp3d * IENe3d, \
		*IEhP3d = UpdatedVSIEfp3d + 2 * IENfp3d*IENe3d;
	memset(UpdatedVSIEFluxM3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSIEFluxP3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSIEFluxS3d, 0, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSERHS3d, 0, Np3d*K3d*Nface*sizeof(double));

	int np = Np3d;
	int oneI = 1;
	double one = 1.0, zero = 0.0;
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
			GetVolumnIntegral2d(UpdatedVSVolumeIntegralX3d + k * Np3d, UpdatedVSTempVolumeIntegralX3d + k * Np3d, &np, &oneI, &np, &one, \
				Dr3d, Ds3d, &np, hu3d + k * Np3d, &np, &zero, &np, rx3d + k * Np3d, sx3d + k * Np3d);
			/*$\bold{r_y}\cdot (Dr*hv)+\bold{s_y}\cdot (Ds*hv)$*/
			GetVolumnIntegral2d(UpdatedVSVolumeIntegralY3d + k * Np3d, UpdatedVSTempVolumeIntegralY3d + k * Np3d, &np, &oneI, &np, &one, \
				Dr3d, Ds3d, &np, hv3d + k * Np3d, &np, &zero, &np, ry3d + k * Np3d, sy3d + k * Np3d);

			Add(UpdatedVSrhs3d + k * Np3d, UpdatedVSVolumeIntegralX3d + k * Np3d, UpdatedVSVolumeIntegralY3d + k * Np3d, Np3d);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		int adjacentE = (int)IEFToE3d[2 * e];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			FetchInnerEdgeFacialValue(IEhM3d + e * IENfp3d, IEhP3d + e * IENfp3d, h3d, IEFToE3d + 2 * e, IEFToN13d + e * IENfp3d, IEFToN23d + e * IENfp3d, Np3d, IENfp3d);
			FetchInnerEdgeFacialValue(IEhuM3d + e * IENfp3d, IEhuP3d + e * IENfp3d, hu3d, IEFToE3d + 2 * e, IEFToN13d + e * IENfp3d, IEFToN23d + e * IENfp3d, Np3d, IENfp3d);
			FetchInnerEdgeFacialValue(IEhvM3d + e * IENfp3d, IEhvP3d + e * IENfp3d, hv3d, IEFToE3d + 2 * e, IEFToN13d + e * IENfp3d, IEFToN23d + e * IENfp3d, Np3d, IENfp3d);
			GetFacialFluxTerm2d(UpdatedVSIEFluxM3d + e * IENfp3d, IEhuM3d + e * IENfp3d, IEhvM3d + e * IENfp3d, IEnx3d + e * IENfp3d, IEny3d + e * IENfp3d, IENfp3d);
			GetFacialFluxTerm2d(UpdatedVSIEFluxP3d + e * IENfp3d, IEhuP3d + e * IENfp3d, IEhvP3d + e * IENfp3d, IEnx3d + e * IENfp3d, IEny3d + e * IENfp3d, IENfp3d);
			GetPCENumericalFluxTerm_HLLC_LAI(UpdatedVSIEFluxS3d + e * IENfp3d, UpdatedVSIEfm3d + e * IENfp3d, UpdatedVSIEfp3d + e * IENfp3d, IEnx3d + e * IENfp3d, IEny3d + e * IENfp3d, &gra, Hcrit, IENfp3d, IENe3d);
		}
		else {
			continue;
		}
	}

	memset(UpdatedVSIEfmod, 0, IENe2d*IENfp3d*sizeof(double));

	//void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		int adjacentE = (int)IEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionWet) {
			VerticalFaceColumnIntegral(UpdatedVSIEFluxS2d + e * IENfp2d, UpdatedVSIEFluxS3d + e * NLayer*IENfp3d, UpdatedVSIEfmod + e * IENfp3d, InvV2d, IENfp3d, IEJz3d + e * NLayer*IENfp3d, NLayer, V1d, IENfp2d, (int)(*(IEFToF3d + e * NLayer * 2)));
		}
		else {
			continue;
		}
	}


	double *BEhuM3d = UpdatedVSBEfm3d, *BEhvM3d = UpdatedVSBEfm3d + BENe3d * BENfp3d, \
		*BEhM3d = UpdatedVSBEfm3d + 2 * BENe3d * BENfp3d;

	int Nfield = 3;

	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		int adjacentE = (int)BEFToE3d[2 * e];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			NdgEdgeType type = (NdgEdgeType)ftype3d[e];  // boundary condition
			FetchBoundaryEdgeFacialValue(BEhuM3d + e * BENfp3d, hu3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);
			FetchBoundaryEdgeFacialValue(BEhvM3d + e * BENfp3d, hv3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);
			FetchBoundaryEdgeFacialValue(BEhM3d + e * BENfp3d, h3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);
			FetchBoundaryEdgeFacialValue(UpdatedVSBEzM3d + e * BENfp3d, z3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);

			ImposeBoundaryCondition(&gra, type, BEnx3d + e * BENfp3d, BEny3d + e * BENfp3d, UpdatedVSBEfm3d + e * BENfp3d, UpdatedVSBEfp3d + e * BENfp3d, \
				UpdatedVSBEzM3d + e * BENfp3d, UpdatedVSBEzP3d + e * BENfp3d, fext3d + e * BENfp3d, BENfp3d, Nfield, BENe3d);
			EvaluateHydroStaticReconstructValue(Hcrit, UpdatedVSBEfm3d + e * BENfp3d, UpdatedVSBEfp3d + e * BENfp3d, UpdatedVSBEzM3d + e * BENfp3d, UpdatedVSBEzP3d + e * BENfp3d, BENfp3d, Nfield, BENe3d);
			GetFacialFluxTerm2d(UpdatedVSBEFluxM3d + e * BENfp3d, BEhuM3d + e * BENfp3d, BEhvM3d + e * BENfp3d, BEnx3d + e * BENfp3d, BEny3d + e * BENfp3d, BENfp3d);
			GetPCENumericalFluxTerm_HLLC_LAI(UpdatedVSBEFluxS3d + e * BENfp3d, UpdatedVSBEfm3d + e * BENfp3d, UpdatedVSBEfp3d + e * BENfp3d, BEnx3d + e * BENfp3d, BEny3d + e * BENfp3d, &gra, Hcrit, BENfp3d, BENe3d);
		}
		else {
			continue;
		}
	}

	memset(UpdatedVSBEfmod, 0, BENe2d*BENfp3d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		int adjacentE = (int)BEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionWet) {
			VerticalFaceColumnIntegral(UpdatedVSBEFluxS2d + e * BENfp2d, UpdatedVSBEFluxS3d + e * NLayer*BENfp3d, UpdatedVSBEfmod + e * BENfp3d, InvV2d, BENfp3d, BEJz3d + e * NLayer*BENfp3d, NLayer, V1d, BENfp2d, (int)(*(BEFToF3d + e * NLayer * 2)));
		}
		else {
			continue;
		}
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		int adjacentE = (int)IEFToE3d[2 * e];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			StrongFormInnerEdgeRHS(e, IEFToE3d, IEFToF3d, Np3d, K3d, IENfp3d, IEFToN13d, IEFToN23d, UpdatedVSIEFluxM3d, UpdatedVSIEFluxP3d, UpdatedVSIEFluxS3d, IEJs3d, IEMb3d, UpdatedVSERHS3d);
		}
		else {
			continue;
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		int adjacentE = (int)BEFToE3d[2 * e];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			StrongFormBoundaryEdgeRHS(e, BEFToE3d, BEFToF3d, Np3d, K3d, BENfp3d, BEFToN13d, UpdatedVSBEFluxM3d, UpdatedVSBEFluxS3d, BEJs3d, BEMb3d, UpdatedVSERHS3d);
		}
		else {
			continue;
		}
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k = 0; k < K3d; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			for (int face = 1; face < Nface; face++) {
				Add(UpdatedVSERHS3d + k * Np3d, UpdatedVSERHS3d + k * Np3d, UpdatedVSERHS3d + face * Np3d*K3d + k * Np3d, Np3d);
			}
		}
		else {
			continue;
		}
    }
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			MultiEdgeContributionByLiftOperator(UpdatedVSERHS3d + k * Np3d, UpdatedVSTempFacialIntegral3d + k * Np3d, &np, &oneI, &np, \
				&one, invM3d, &np, &np, &zero, &np, J3d + k * Np3d, Np3d);
		}
		else {
			continue;
		}
	}

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			Minus(UpdatedVSrhs3d + k * Np3d, UpdatedVSERHS3d + k * Np3d, UpdatedVSrhs3d + k * Np3d, Np3d);
		}
		else {
			continue;
		}
	}

	/**********************************************************  Three Dimensional Part Finished  *******************************************************************************/

	/**********************************************************  Two Dimensional Part  *******************************************************************************/
	memset(UpdatedVSrhs2d, 0, Np2d*K2d*sizeof(double));
	double *IEhuM2d = UpdatedVSIEfm2d, *IEhvM2d = UpdatedVSIEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = UpdatedVSIEfm2d + 2 * IENfp2d*IENe2d;
	double *IEhuP2d = UpdatedVSIEfp2d, *IEhvP2d = UpdatedVSIEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = UpdatedVSIEfp2d + 2 * IENfp2d*IENe2d;
	memset(UpdatedVSIEFluxM2d, 0, IENfp2d*IENe2d*sizeof(double));
	memset(UpdatedVSIEFluxP2d, 0, IENfp2d*IENe2d*sizeof(double));
	memset(UpdatedVSERHS2d, 0, Np2d*K2d*Nface*sizeof(double));

	np = Np2d;
	oneI = 1;
	one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		if ((NdgRegionType)Status2d[k] == NdgRegionWet) {
			/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
			GetVolumnIntegral2d(UpdatedVSVolumeIntegralX + k * Np2d, UpdatedVSTempVolumeIntegralX + k * Np2d, &np, &oneI, &np, &one, \
				Dr2d, Ds2d, &np, hu2d + k * Np2d, &np, &zero, &np, rx2d + k * Np2d, sx2d + k * Np2d);
			/*$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$*/
			GetVolumnIntegral2d(UpdatedVSVolumeIntegralY + k * Np2d, UpdatedVSTempVolumeIntegralY + k * Np2d, &np, &oneI, &np, &one, \
				Dr2d, Ds2d, &np, hv2d + k * Np2d, &np, &zero, &np, ry2d + k * Np2d, sy2d + k * Np2d);

			Add(UpdatedVSrhs2d + k * Np2d, UpdatedVSVolumeIntegralX + k * Np2d, UpdatedVSVolumeIntegralY + k * Np2d, Np2d);
		}
		else {
			continue;
		}
	}
	/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		int adjacentE = (int)IEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionWet) {
			FetchInnerEdgeFacialValue(IEhM2d + e * IENfp2d, IEhP2d + e * IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e * IENfp2d, IEFToN22d + e * IENfp2d, Np2d, IENfp2d);
			FetchInnerEdgeFacialValue(IEhuM2d + e * IENfp2d, IEhuP2d + e * IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e * IENfp2d, IEFToN22d + e * IENfp2d, Np2d, IENfp2d);
			FetchInnerEdgeFacialValue(IEhvM2d + e * IENfp2d, IEhvP2d + e * IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e * IENfp2d, IEFToN22d + e * IENfp2d, Np2d, IENfp2d);
			GetFacialFluxTerm2d(UpdatedVSIEFluxM2d + e * IENfp2d, IEhuM2d + e * IENfp2d, IEhvM2d + e * IENfp2d, IEnx2d + e * IENfp2d, IEny2d + e * IENfp2d, IENfp2d);
			GetFacialFluxTerm2d(UpdatedVSIEFluxP2d + e * IENfp2d, IEhuP2d + e * IENfp2d, IEhvP2d + e * IENfp2d, IEnx2d + e * IENfp2d, IEny2d + e * IENfp2d, IENfp2d);
		}
		else {
			continue;
		}
	}

	double *BEhuM2d = UpdatedVSBEfm2d, *BEhvM2d = UpdatedVSBEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = UpdatedVSBEfm2d + 2 * BENe2d * BENfp2d;
	memset(UpdatedVSBEFluxM2d, 0, BENe2d*BENfp2d*sizeof(double));

	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		int adjacentE = (int)BEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionWet) {
			NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
			FetchBoundaryEdgeFacialValue(BEhM2d + e * BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			FetchBoundaryEdgeFacialValue(BEhuM2d + e * BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			FetchBoundaryEdgeFacialValue(BEhvM2d + e * BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			FetchBoundaryEdgeFacialValue(UpdatedVSBEzM2d + e * BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			ImposeBoundaryCondition(&gra, type, BEnx2d + e * BENfp2d, BEny2d + e * BENfp2d, UpdatedVSBEfm2d + e * BENfp2d, UpdatedVSBEfp2d + e * BENfp2d, \
				UpdatedVSBEzM2d + e * BENfp2d, UpdatedVSBEzP2d + e * BENfp2d, fext2d + e * BENfp2d, BENfp2d, Nfield, BENe2d);
			EvaluateHydroStaticReconstructValue(Hcrit, UpdatedVSBEfm2d + e * BENfp2d, UpdatedVSBEfp2d + e * BENfp2d, UpdatedVSBEzM2d + e * BENfp2d, UpdatedVSBEzP2d + e * BENfp2d, BENfp2d, Nfield, BENe2d);
			GetFacialFluxTerm2d(UpdatedVSBEFluxM2d + e * BENfp2d, BEhuM2d + e * BENfp2d, BEhvM2d + e * BENfp2d, BEnx2d + e * BENfp2d, BEny2d + e * BENfp2d, BENfp2d);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		int adjacentE = (int)IEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionWet) {
			StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, UpdatedVSIEFluxM2d, UpdatedVSIEFluxP2d, UpdatedVSIEFluxS2d, IEJs2d, IEMb2d, UpdatedVSERHS2d);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		int adjacentE = (int)BEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionWet) {
			StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, UpdatedVSBEFluxM2d, UpdatedVSBEFluxS2d, BEJs2d, BEMb2d, UpdatedVSERHS2d);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		if ((NdgRegionType)Status2d[k] == NdgRegionWet) {
			for (int face = 1; face < Nface; face++) {
				Add(UpdatedVSERHS2d + k * Np2d, UpdatedVSERHS2d + k * Np2d, UpdatedVSERHS2d + face * Np2d*K2d + k * Np2d, Np2d);
			}
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if ((NdgRegionType)Status2d[k] == NdgRegionWet) {
			MultiEdgeContributionByLiftOperator(UpdatedVSERHS2d + k * Np2d, UpdatedVSTempFacialIntegral + k * Np2d, &np, &oneI, &np, \
				&one, invM2d, &np, &np, &zero, &np, J2d + k * Np2d, Np2d);
		}
		else {
			continue;
		}
	}

	/*Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		if ((NdgRegionType)Status2d[k] == NdgRegionWet) {
			Minus(UpdatedVSrhs2d + k * Np2d, UpdatedVSERHS2d + k * Np2d, UpdatedVSrhs2d + k * Np2d, Np2d);
		}
		else {
			continue;
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		if ((NdgRegionType)Status2d[k] == NdgRegionWet) {
			NdgExtend2dField(UpdatedVSfield2d, UpdatedVSrhs2d, Np2d, k, Np3d, NLayer, Nz);
		}
		else {
			continue;
		}
	}

	/**********************************************************  Two Dimensional Part Finished  *******************************************************************************/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			/*Substract the VSfield2d from VSrhs3d to assemble the final right hand side*/
			Minus(UpdatedVSrhs3d + k * Np3d, UpdatedVSrhs3d + k * Np3d, UpdatedVSfield2d + k * Np3d, Np3d);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		if ((NdgRegionType)Status2d[k] == NdgRegionWet) {
			VerticalIntegralFromBottom(VerticalVelocity + k * NLayer*Np3d, UpdatedVSrhs3d + k * NLayer*Np3d, Jz + k * NLayer*Np3d, Updatedfmod + k * Np3d, NLayer, Np3d, InvV3d, Np2d, Npz, Vint);
		}
		//For partial dry and dry element, the VerticalVelocity in each layer is set to 0.
		else {
			for (int n = 0; n < NLayer; n++) {
				for (int i = 0; i < Np3d; i++) {
					VerticalVelocity[k * NLayer * Np3d + n * Np3d + i] = 0.0;
				}
			}
		}
	}

	free(InvV2d);
	free(InvV3d);
}