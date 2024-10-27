#include "NdgMath.h"
#include "NdgSWE.h"
#include "NdgSWE3D.h"
#include "NdgMemory.h"
#include "NdgQuadFreeStrongFormPECSolver2d.h"
#include <stdio.h>

extern double *PCEUpdatedIEfm2d, *PCEUpdatedIEfp2d, *PCEUpdatedIEFluxM2d, *PCEUpdatedIEFluxP2d, *PCEUpdatedIEFluxS2d, \
*PCEUpdatedERHS2d, *PCEUpdatedVolumeIntegralX, *PCEUpdatedTempVolumeIntegralX, *PCEUpdatedVolumeIntegralY, \
*PCEUpdatedTempVolumeIntegralY, *PCEUpdatedBEfm2d, *PCEUpdatedBEzM2d, *PCEUpdatedBEfp2d, *PCEUpdatedBEzP2d, \
*PCEUpdatedBEFluxS2d, *PCEUpdatedBEFluxM2d, *PCEUpdatedPCETempFacialIntegral, *PCEUpdatedIEfmod, *PCEUpdatedBEfmod,\
*PCEUpdatedIEfm3d, *PCEUpdatedIEfp3d, *PCEUpdatedIEFluxS3d, *PCEUpdatedBEFluxS3d, *PCEUpdatedBEfm3d, *PCEUpdatedBEfp3d, \
*PCEUpdatedBEzM3d, *PCEUpdatedBEzP3d;

//extern char *PCEUpdatedInitialized;
extern double gra;
extern double Hcrit;
extern int Nfield3d;
extern signed char *Status3d;
//void MyExit()
//{
//	if (!strcmp("True", PCEUpdatedInitialized)){
//		PCEUpdatedMemoryDeAllocation();
//		PCEUpdatedInitialized = "False";
//	}
//	return;
//}
//PCEUpdatedMemoryDeAllocation();

NdgQuadFreeStrongFormPECSolver2d::NdgQuadFreeStrongFormPECSolver2d()
{
}

NdgQuadFreeStrongFormPECSolver2d::~NdgQuadFreeStrongFormPECSolver2d()
{
};

void NdgQuadFreeStrongFormPECSolver2d::evaluatePCERHSUpdated(double *fphys_, double *frhs_, double *fext_, int *varFieldIndex_, double *fphys2d_, double *fext2d_)
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

	int K3d = *(meshunion->K);
	int NLayer = *(meshunion->Nlayer);

	int Np3d = *(meshunion->cell_p->Np);

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

	/*Data contained in two-dimensional standard cell*/
	double *Dr2d = meshunion->mesh2d_p->mesh2dcell_p->Dr2d;
	double *Ds2d = meshunion->mesh2d_p->mesh2dcell_p->Ds2d;
	int Nface = *(meshunion->mesh2d_p->mesh2dcell_p->Nface2d);
	double *invM2d = meshunion->mesh2d_p->mesh2dcell_p->invM2d;

	double *ftype2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->ftype2d;
	double *ftype3d = meshunion->boundaryedge_p->ftype;

	/*Data contained in two dimensional physical field*/
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

	/*The two-dimensional external value*/
	double *fext2d = fext2d_;
	double *fext3d = fext_;
//	double *BEhE2d = fext2d, *BEhuE2d = fext2d + BENe2d * BENfp2d, *BEhvE2d = fext2d + 2 * BENe2d * BENfp2d;

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
	double *IEJz3d = meshunion->inneredge_p->Jz;
	double *V1d = meshunion->inneredge_p->V1d;
	double *V2d = meshunion->inneredge_p->V2d;

	double *InvV2d = (double *)malloc(IENfp3d*IENfp3d*sizeof(double));
	memcpy(InvV2d, V2d, IENfp3d*IENfp3d*sizeof(double));
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
	double *BEJz3d = meshunion->boundaryedge_p->Jz;

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
	//int adjacentE2 = 0;//store the adjacent element number of the current face
	/***  Wet and Dry  ***/

	double *RHS = frhs_;

	/*if (!strcmp("False", PCEUpdatedInitialized)){
		PCEUpdatedMemoryAllocation(IENfp2d, IENe2d, Np2d, K2d, Nface, BENe2d, BENfp2d,\
			IENfp3d, IENe3d, BENe3d, BENfp3d);
	}*/
	/*PCEUpdatedMemoryAllocation(IENfp2d, IENe2d, Np2d, K2d, Nface, BENe2d, BENfp2d, \
		IENfp3d, IENe3d, BENe3d, BENfp3d);*/

	double *IEhuM2d = PCEUpdatedIEfm2d, *IEhvM2d = PCEUpdatedIEfm2d + IENfp2d * IENe2d, \
		*IEhM2d = PCEUpdatedIEfm2d + 2 * IENfp2d*IENe2d;
	
	double *IEhuP2d = PCEUpdatedIEfp2d, *IEhvP2d = PCEUpdatedIEfp2d + IENfp2d * IENe2d, \
		*IEhP2d = PCEUpdatedIEfp2d + 2 * IENfp2d*IENe2d;
	
	memset(PCEUpdatedIEFluxM2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(PCEUpdatedIEFluxP2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(PCEUpdatedIEFluxS2d, 0, IENfp2d*IENe2d*sizeof(double));
	
	memset(PCEUpdatedERHS2d, 0, Np2d*K2d*Nface*sizeof(double));


	int np = Np2d;
	int oneI = 1;
	double one = 1.0, zero = 0.0;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){                /*  2D Volume Integral Part  */
		if ((NdgRegionType)Status2d[k] == NdgRegionDry) {
			continue;
		}
		else  {
			/*$\bold{r_x}\cdot (Dr*hu2d)+\bold{s_x}\cdot (Ds*hu2d)$*/
			GetVolumnIntegral2d(PCEUpdatedVolumeIntegralX + k * Np2d, PCEUpdatedTempVolumeIntegralX + k * Np2d, &np, &oneI, &np, &one, \
				Dr2d, Ds2d, &np, hu2d + k * Np2d, &np, &zero, &np, rx2d + k * Np2d, sx2d + k * Np2d);
			/*$\bold{r_y}\cdot (Dr*hv2d)+\bold{s_y}\cdot (Ds*hv2d)$*/
			GetVolumnIntegral2d(PCEUpdatedVolumeIntegralY + k * Np2d, PCEUpdatedTempVolumeIntegralY + k * Np2d, &np, &oneI, &np, &one, \
				Dr2d, Ds2d, &np, hv2d + k * Np2d, &np, &zero, &np, ry2d + k * Np2d, sy2d + k * Np2d);

			Add(RHS + k * Np2d, PCEUpdatedVolumeIntegralX + k * Np2d, PCEUpdatedVolumeIntegralY + k * Np2d, Np2d);
		}
	}
	/*Two dimensional inner edge flux part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		int adjacentE = (int)IEFToE2d[2 * e];
		int adjacentE2 = (int)IEFToE2d[2 * e + 1];
		if ((NdgRegionType)Status3d[adjacentE - 1] != NdgRegionWet || (NdgRegionType)Status3d[adjacentE2 - 1] != NdgRegionWet) {
			continue; //no wet element 
		}
		else {
			FetchInnerEdgeFacialValue(IEhM2d + e * IENfp2d, IEhP2d + e * IENfp2d, h2d, IEFToE2d + 2 * e, IEFToN12d + e * IENfp2d, IEFToN22d + e * IENfp2d, Np2d, IENfp2d);
			FetchInnerEdgeFacialValue(IEhuM2d + e * IENfp2d, IEhuP2d + e * IENfp2d, hu2d, IEFToE2d + 2 * e, IEFToN12d + e * IENfp2d, IEFToN22d + e * IENfp2d, Np2d, IENfp2d);
			FetchInnerEdgeFacialValue(IEhvM2d + e * IENfp2d, IEhvP2d + e * IENfp2d, hv2d, IEFToE2d + 2 * e, IEFToN12d + e * IENfp2d, IEFToN22d + e * IENfp2d, Np2d, IENfp2d);
			GetFacialFluxTerm2d(PCEUpdatedIEFluxM2d + e * IENfp2d, IEhuM2d + e * IENfp2d, IEhvM2d + e * IENfp2d, IEnx2d + e * IENfp2d, IEny2d + e * IENfp2d, IENfp2d);
			GetFacialFluxTerm2d(PCEUpdatedIEFluxP2d + e * IENfp2d, IEhuP2d + e * IENfp2d, IEhvP2d + e * IENfp2d, IEnx2d + e * IENfp2d, IEny2d + e * IENfp2d, IENfp2d);
		}
	}

	
	double *BEhuM2d = PCEUpdatedBEfm2d, *BEhvM2d = PCEUpdatedBEfm2d + BENe2d * BENfp2d, \
		*BEhM2d = PCEUpdatedBEfm2d + 2 * BENe2d * BENfp2d;

	memset(PCEUpdatedBEFluxS2d, 0, BENe2d*BENfp2d*sizeof(double));
	
	memset(PCEUpdatedBEFluxM2d, 0, BENe2d*BENfp2d*sizeof(double));

	int Nfield = 3;
	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		int adjacentE = (int)BEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionDry) {
			continue; //not wet 
		}
		else {
			NdgEdgeType type = (NdgEdgeType)ftype2d[e];  // boundary condition
			FetchBoundaryEdgeFacialValue(BEhM2d + e * BENfp2d, h2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			FetchBoundaryEdgeFacialValue(BEhuM2d + e * BENfp2d, hu2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			FetchBoundaryEdgeFacialValue(BEhvM2d + e * BENfp2d, hv2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			FetchBoundaryEdgeFacialValue(PCEUpdatedBEzM2d + e * BENfp2d, z2d, BEFToE2d + 2 * e, BEFToN12d + e * BENfp2d, Np2d, BENfp2d);
			ImposeBoundaryCondition(&gra, type, BEnx2d + e * BENfp2d, BEny2d + e * BENfp2d, PCEUpdatedBEfm2d + e * BENfp2d, PCEUpdatedBEfp2d + e * BENfp2d, \
				PCEUpdatedBEzM2d + e * BENfp2d, PCEUpdatedBEzP2d + e * BENfp2d, fext2d + e * BENfp2d, BENfp2d, Nfield, BENe2d);
			EvaluateHydroStaticReconstructValue(Hcrit, PCEUpdatedBEfm2d + e * BENfp2d, PCEUpdatedBEfp2d + e * BENfp2d, PCEUpdatedBEzM2d + e * BENfp2d, PCEUpdatedBEzP2d + e * BENfp2d, BENfp2d, Nfield, BENe2d);
			GetFacialFluxTerm2d(PCEUpdatedBEFluxM2d + e * BENfp2d, BEhuM2d + e * BENfp2d, BEhvM2d + e * BENfp2d, BEnx2d + e * BENfp2d, BEny2d + e * BENfp2d, BENfp2d);
		}
	}

	double *IEhuM3d = PCEUpdatedIEfm3d, *IEhvM3d = PCEUpdatedIEfm3d + IENfp3d * IENe3d, \
		*IEhM3d = PCEUpdatedIEfm3d + 2 * IENfp3d*IENe3d;
	double *IEhuP3d = PCEUpdatedIEfp3d, *IEhvP3d = PCEUpdatedIEfp3d + IENfp3d * IENe3d, \
		*IEhP3d = PCEUpdatedIEfp3d + 2 * IENfp3d*IENe3d;
	memset(PCEUpdatedIEFluxS3d, 0, IENfp3d*IENe3d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe3d; e++){
		int adjacentE = (int)IEFToE3d[2 * e];
		int adjacentE2 = (int)IEFToE3d[2 * e + 1];
		if ((NdgRegionType)Status3d[adjacentE - 1] != NdgRegionWet || (NdgRegionType)Status3d[adjacentE2 - 1] != NdgRegionWet) {
			continue; //不是湿单元旁边单元为湿单元才计算面积分
		}
		else {
			FetchInnerEdgeFacialValue(IEhM3d + e * IENfp3d, IEhP3d + e * IENfp3d, h3d, IEFToE3d + 2 * e, IEFToN13d + e * IENfp3d, IEFToN23d + e * IENfp3d, Np3d, IENfp3d);
			FetchInnerEdgeFacialValue(IEhuM3d + e * IENfp3d, IEhuP3d + e * IENfp3d, hu3d, IEFToE3d + 2 * e, IEFToN13d + e * IENfp3d, IEFToN23d + e * IENfp3d, Np3d, IENfp3d);
			FetchInnerEdgeFacialValue(IEhvM3d + e * IENfp3d, IEhvP3d + e * IENfp3d, hv3d, IEFToE3d + 2 * e, IEFToN13d + e * IENfp3d, IEFToN23d + e * IENfp3d, Np3d, IENfp3d);
			GetPCENumericalFluxTerm_HLLC_LAI(PCEUpdatedIEFluxS3d + e * IENfp3d, PCEUpdatedIEfm3d + e * IENfp3d, PCEUpdatedIEfp3d + e * IENfp3d, IEnx3d + e * IENfp3d, IEny3d + e * IENfp3d, &gra, Hcrit, IENfp3d, IENe3d);
		}
	}

	memset(PCEUpdatedIEfmod, 0, IENe2d*IENfp3d*sizeof(double));

	//void VerticalFaceColumnIntegral(double *dest, double *source, double *fmod, double *InvV2d, int Nfp2d, double *Jz, int Nlayer, double *V1d, int LNfp2d, int FToF)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		int adjacentE = (int)IEFToE2d[2 * e];
		int adjacentE2 = (int)IEFToE2d[2 * e + 1];
		if ((NdgRegionType)Status3d[adjacentE - 1] != NdgRegionWet || (NdgRegionType)Status3d[adjacentE2 - 1] != NdgRegionWet) {
			continue; 
		}
		else {
			VerticalFaceColumnIntegral(PCEUpdatedIEFluxS2d + e * IENfp2d, PCEUpdatedIEFluxS3d + e * NLayer*IENfp3d, PCEUpdatedIEfmod + e * IENfp3d, InvV2d, (int)IENfp3d, IEJz3d + e * NLayer*IENfp3d, NLayer, V1d, (int)IENfp2d, (int)(*(IEFToF3d + e * NLayer * 2)));
		}
	}


	double *BEhuM3d = PCEUpdatedBEfm3d, *BEhvM3d = PCEUpdatedBEfm3d + BENe3d * BENfp3d, \
		*BEhM3d = PCEUpdatedBEfm3d + 2 * BENe3d * BENfp3d;

	/*fetch boundary edge value h, hu, hv and z, apply hydrostatic construction at the boundary and compute the numerical flux*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe3d; e++){
		int adjacentE = (int)BEFToE3d[2 * e];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue;  
		}
		else {
			NdgEdgeType type = (NdgEdgeType)ftype3d[e];  // boundary condition
			FetchBoundaryEdgeFacialValue(BEhuM3d + e * BENfp3d, hu3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);
			FetchBoundaryEdgeFacialValue(BEhvM3d + e * BENfp3d, hv3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);
			FetchBoundaryEdgeFacialValue(BEhM3d + e * BENfp3d, h3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);
			FetchBoundaryEdgeFacialValue(PCEUpdatedBEzM3d + e * BENfp3d, z3d, BEFToE3d + 2 * e, BEFToN13d + e * BENfp3d, Np3d, BENfp3d);

			ImposeBoundaryCondition(&gra, type, BEnx3d + e * BENfp3d, BEny3d + e * BENfp3d, PCEUpdatedBEfm3d + e * BENfp3d, PCEUpdatedBEfp3d + e * BENfp3d, \
				PCEUpdatedBEzM3d + e * BENfp3d, PCEUpdatedBEzP3d + e * BENfp3d, fext3d + e * BENfp3d, BENfp3d, Nfield, BENe3d);
			EvaluateHydroStaticReconstructValue(Hcrit, PCEUpdatedBEfm3d + e * BENfp3d, PCEUpdatedBEfp3d + e * BENfp3d, PCEUpdatedBEzM3d + e * BENfp3d, PCEUpdatedBEzP3d + e * BENfp3d, BENfp3d, Nfield, BENe3d);
			GetPCENumericalFluxTerm_HLLC_LAI(PCEUpdatedBEFluxS3d + e * BENfp3d, PCEUpdatedBEfm3d + e * BENfp3d, PCEUpdatedBEfp3d + e * BENfp3d, BEnx3d + e * BENfp3d, BEny3d + e * BENfp3d, &gra, Hcrit, BENfp3d, BENe3d);
		}
	}

	memset(PCEUpdatedBEfmod, 0, BENe2d*BENfp3d*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < BENe2d; e++){
		int adjacentE = (int)BEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionDry) {
			continue;  
		}
		else {
			VerticalFaceColumnIntegral(PCEUpdatedBEFluxS2d + e * BENfp2d, PCEUpdatedBEFluxS3d + e * NLayer*BENfp3d, PCEUpdatedBEfmod + e * BENfp3d, InvV2d, (int)BENfp3d, BEJz3d + e * NLayer*BENfp3d, NLayer, V1d, (int)BENfp2d, (int)(*(BEFToF3d + e * NLayer * 2)));
		}
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int e = 0; e < IENe2d; e++){
		int adjacentE = (int)IEFToE2d[2 * e];
		int adjacentE2 = (int)IEFToE2d[2 * e + 1];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionDry && (NdgRegionType)Status2d[adjacentE2 - 1] == NdgRegionDry) {
			continue; 
		}
		else {
			StrongFormInnerEdgeRHS(e, IEFToE2d, IEFToF2d, Np2d, K2d, IENfp2d, IEFToN12d, IEFToN22d, PCEUpdatedIEFluxM2d, PCEUpdatedIEFluxP2d, PCEUpdatedIEFluxS2d, IEJs2d, IEMb2d, PCEUpdatedERHS2d);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif   
	for (int e = 0; e < BENe2d; e++){
		int adjacentE = (int)BEFToE2d[2 * e];
		if ((NdgRegionType)Status2d[adjacentE - 1] == NdgRegionDry) {
			continue; //Both dry 
		}
		else {
			StrongFormBoundaryEdgeRHS(e, BEFToE2d, BEFToF2d, Np2d, K2d, BENfp2d, BEFToN12d, PCEUpdatedBEFluxM2d, PCEUpdatedBEFluxS2d, BEJs2d, BEMb2d, PCEUpdatedERHS2d);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K2d; k++){
		for (int face = 1; face<Nface; face++){
			Add(PCEUpdatedERHS2d + k*Np2d, PCEUpdatedERHS2d + k*Np2d, PCEUpdatedERHS2d + face*Np2d*K2d + k*Np2d, Np2d);
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if ((NdgRegionType)Status2d[k] == NdgRegionDry) {
			continue;
		}
		else {
			MultiEdgeContributionByLiftOperator(PCEUpdatedERHS2d + k * Np2d, PCEUpdatedPCETempFacialIntegral + k * Np2d, &np, &oneI, &np, \
				&one, invM2d, &np, &np, &zero, &np, J2d + k * Np2d, Np2d);
		}
	}

	/*Add face integral and volume integral up to form the right hand side corresponding to the discretization of the depth-averaged part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++){
		if ((NdgRegionType)Status2d[k] == NdgRegionDry) {
			continue;
		}
		else {
			Minus(RHS + k * Np2d, PCEUpdatedERHS2d + k * Np2d, RHS + k * Np2d, Np2d);
		}
	}

	free(InvV2d);

}