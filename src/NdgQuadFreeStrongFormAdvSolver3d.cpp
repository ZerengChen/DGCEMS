#include "NdgMath.h"
#include "NdgSWE.h"
#include "NdgSWE3D.h"
#include "NdgQuadFreeStrongFormAdvSolver3d.h"
#include <omp.h>
#include <iostream>
using namespace std;

NdgQuadFreeStrongFormAdvSolver3d::NdgQuadFreeStrongFormAdvSolver3d()
{
}

NdgQuadFreeStrongFormAdvSolver3d::~NdgQuadFreeStrongFormAdvSolver3d()
{
}



extern double *TempFacialIntegral, *IEfm, *IEfp, *IEFluxM, *IEFluxP, \
*IEFluxS, *ERHS, *BEfm, *BEfp, *AdvzM, *AdvzP, *BEFluxM, \
*BEFluxS, *BotEfm, *BotEfp, *BotEFluxM, *BotEFluxP, *BotEFluxS, \
*BotBEfm, *BotBEFluxM, *BotBEFluxS, *SurfBEfm, *SurfBEFluxM, \
*SurfBEFluxS, *E, *G, *H, *TempVolumeIntegral;

//extern char *AdvInitialized;
extern double gra;
extern double Hcrit;
extern int Nvar;
extern int Nfield3d;
extern signed char *Status3d;

//void MyExit()
//{
	//if (!strcmp("True", AdvInitialized)){
		//AdvMemoryDeAllocation();
		//AdvInitialized = "False";
	//}
	//return;
//}

void NdgQuadFreeStrongFormAdvSolver3d::evaluateAdvectionRHS(double *fphys_, double *frhs_, double *fext, int *varFieldIndex_, int*pE3d, int MyID)
{
	double *rx = meshunion->rx;
	double *sx = meshunion->sx;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;
	double *tz = meshunion->tz;
	double *J = meshunion->J;
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	/*For cell object*/
	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	double *Dt = meshunion->cell_p->Dt;
	int Nface = *(meshunion->cell_p->Nface);
	double *invM = meshunion->cell_p->invM;
	/*For inner edge object*/
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
	/*For boundary edge object*/
	int BENe = *(meshunion->boundaryedge_p->Ne);
	int BENfp = *(meshunion->boundaryedge_p->Nfp);
	double *BEMb = meshunion->boundaryedge_p->M;
	double *BEJs = meshunion->boundaryedge_p->Js;
	double *BEnx = meshunion->boundaryedge_p->nx;
	double *BEny = meshunion->boundaryedge_p->ny;
	double *BEFToE = meshunion->boundaryedge_p->FToE;
	double *BEFToF = meshunion->boundaryedge_p->FToF;
	double *BEFToN1 = meshunion->boundaryedge_p->FToN1;
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
	/*For bottom boundary edge object*/
	int BotBENe = *(meshunion->bottomboundaryedge_p->Ne);
	int BotBENfp = *(meshunion->bottomboundaryedge_p->Nfp);
	double *BotBEMb = meshunion->bottomboundaryedge_p->M;
	double *BotBEJs = meshunion->bottomboundaryedge_p->Js;
	double *BotBEnz = meshunion->bottomboundaryedge_p->nz;
	double *BotBEFToE = meshunion->bottomboundaryedge_p->FToE;
	double *BotBEFToF = meshunion->bottomboundaryedge_p->FToF;
	double *BotBEFToN1 = meshunion->bottomboundaryedge_p->FToN1;
    /*For surface boundary edge object*/
	int SurfBENe = *(meshunion->surfaceboundaryedge_p->Ne);
	int SurfBENfp = *(meshunion->surfaceboundaryedge_p->Nfp);
	double *SurfBEMb = meshunion->surfaceboundaryedge_p->M;
	double *SurfBEJs = meshunion->surfaceboundaryedge_p->Js;
	double *SurfBEnz = meshunion->surfaceboundaryedge_p->nz;
	double *SurfBEFToE = meshunion->surfaceboundaryedge_p->FToE;
	double *SurfBEFToF = meshunion->surfaceboundaryedge_p->FToF;
	double *SurfBEFToN1 = meshunion->surfaceboundaryedge_p->FToN1;

	int *varFieldIndex = varFieldIndex_;

	double *fphys = fphys_;
	double *hu = fphys_, *hv = fphys_ + Np*K, *omega = fphys_ + 2 * Np*K,  \
		*h = fphys_ + 3 * Np*K, *z = fphys_ + 5 * Np*K;
	//double *fext = fext;
	//double gra = gra_;
	//double Hcrit = Hcrit_;

	double *ftype = meshunion->boundaryedge_p->ftype;

	int np = Np;
	int oneI = 1;
	double one = 1.0, zero = 0.0;
 //   size_t NdimOut = 3;
	//mwSize dimOut[3] = { Np, K, Nvar };
	//plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
	double *OutputRHS = frhs_;

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
	//int adjacentE2 = 0; //store the adjacent element number of the current face
	double no_gra = 0.0; //For Flooding elements
	/***  Wet and Dry  ***/
                                                                  /************************  Face Integral Part  ****************************/
/*************************************************************************************************************************************/
	/**************************************Inner Edge Part*******************************************************/
	            /********************************************************************/
	/*Allocate memory for fm and fp defined over inner edges. Here, variables correspond to hu, hv, hT, hS, sediment and other passive transport material, vertical velocity omega is not included*/
	
	double *huM = IEfm, *hvM = IEfm + IENfp * IENe, *hM = IEfm + 2 * IENfp*IENe;
	double *huP = IEfp, *hvP = IEfp + IENfp * IENe, *hP = IEfp + 2 * IENfp*IENe;
	/*Allocate memory for IEFluxM, IEFluxP and IEFluxS, and calculate these flux term*/
	memset(IEFluxM, 0, IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxP, 0, IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxS, 0, IENfp*IENe*Nvar*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		int adjacentE = (int)IEFToE[2 * face];
		int adjacentE2 = (int)IEFToE[2 * face + 1];
		if ((MyID == pE3d[adjacentE - 1]) || (MyID == pE3d[adjacentE2 - 1])) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry && (NdgRegionType)Status3d[adjacentE2 - 1] == NdgRegionDry) {
				continue; //不是湿单元旁边单元为湿单元才计算面积分//换成干单元旁边是干单元才不计算
			}
			else {
				//printf("The order of thread is:%d\n",omp_get_thread_num());  // This is not thread-safe
				/*Fetch variable IEfm and IEfp first*/
				FetchInnerEdgeFacialValue(hM + face * IENfp, hP + face * IENfp, h, IEFToE + 2 * face, \
					IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
				FetchInnerEdgeFacialValue(huM + face * IENfp, huP + face * IENfp, hu, IEFToE + 2 * face, \
					IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
				FetchInnerEdgeFacialValue(hvM + face * IENfp, hvP + face * IENfp, hv, IEFToE + 2 * face, \
					IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
				/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
				here 1 stands for the space occupied by water depth h.*/
				for (int field = 2; field < Nvar; field++) {
					FetchInnerEdgeFacialValue(IEfm + (field + 1)*IENe*IENfp + face * IENfp, \
						IEfp + (field + 1)*IENe*IENfp + face * IENfp, fphys + (varFieldIndex[field] - 1)*Np*K, \
						IEFToE + 2 * face, IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
				}

				EvaluateVerticalFaceSurfFlux(IEFluxM + face * IENfp, IEfm + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, &gra, Hcrit, IENfp, Nvar, IENe);
				EvaluateVerticalFaceSurfFlux(IEFluxP + face * IENfp, IEfp + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, &gra, Hcrit, IENfp, Nvar, IENe);

				EvaluateVerticalFaceNumFlux_HLLC_LAI(IEFluxS + face * IENfp, IEfm + face * IENfp, IEfp + face * IENfp, \
					IEnx + face * IENfp, IEny + face * IENfp, &gra, Hcrit, IENfp, Nvar, IENe);

			}
		}
	}

/*Allocate memory for contribution to RHS due to inner edge facial integral, and
 calculate contribution to RHS due to inner edge facial integral in strong form manner*/
	memset(ERHS, 0, Np*K*Nvar*Nface*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < IENe; face++){
		int adjacentE = (int)IEFToE[2 * face];
		int adjacentE2 = (int)IEFToE[2 * face + 1];
		if ((MyID == pE3d[adjacentE - 1]) || (MyID == pE3d[adjacentE2 - 1])) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry && (NdgRegionType)Status3d[adjacentE2 - 1] == NdgRegionDry) {
				continue; //不是湿单元旁边单元为湿单元才计算面积分
			}
			else {
				for (int field = 0; field < Nvar; field++) {
					StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, IEFluxM + field * IENe*IENfp, \
						IEFluxP + field * IENe*IENfp, IEFluxS + field * IENe*IENfp, IEJs, IEMb, ERHS + field * Np*K*Nface);
				}
			}
		}
	}
	/*************************************************************************************************************************************/
	/*************************************************************************************************************************************/
	            /**************************************Boundary Edge Part*******************************************************/
	                           /********************************************************************/
	/*Allocate memory for fm and fp defined over boundary edges. Here, variables correspond to hu, hv, h, hT, hS, 
	sediment and other passive transport material, vertical velocity omega is not included for boundary edges*/
	huM = BEfm, hvM = BEfm + BENfp*BENe, hM = BEfm + 2 * BENfp*BENe;
	huP = BEfp, hvP = BEfp + BENfp*BENe, hP = BEfp + 2 * BENfp*BENe;
	memset(BEFluxM, 0, BENfp*BENe*Nvar*sizeof(double));
	memset(BEFluxS, 0, BENfp*BENe*Nvar*sizeof(double));

	/*Fetch variable BEfm and BEfp first, then impose boundary condition and conduct hydrostatic reconstruction.
	Finally, calculate local flux term, adjacent flux term and numerical flux term*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		int adjacentE = (int)BEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
				continue; //BEFToE is just itself
			}
			else {
				NdgEdgeType type = (NdgEdgeType)ftype[face];  // boundary condition
				FetchBoundaryEdgeFacialValue(huM + face * BENfp, hu, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				FetchBoundaryEdgeFacialValue(hvM + face * BENfp, hv, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				FetchBoundaryEdgeFacialValue(hM + face * BENfp, h, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				FetchBoundaryEdgeFacialValue(AdvzM + face * BENfp, z, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
				here 1 stands for the memory occupied by water depth h*/
				for (int field = 2; field < Nvar; field++) {
					FetchBoundaryEdgeFacialValue(BEfm + (field + 1)*BENe*BENfp + face * BENfp, \
						fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
						BEFToE + 2 * face, BEFToN1 + BENfp * face, Np, BENfp);
				}
				ImposeBoundaryCondition(&gra, type, BEnx + face * BENfp, BEny + face * BENfp, BEfm + face * BENfp, BEfp + face * BENfp, \
					AdvzM + face * BENfp, AdvzP + face * BENfp, fext + face * BENfp, BENfp, Nvar + 1, BENe);
				EvaluateHydroStaticReconstructValue(Hcrit, BEfm + face * BENfp, BEfp + face * BENfp, AdvzM + face * BENfp, AdvzP + face * BENfp, BENfp, Nvar + 1, BENe);

				EvaluateVerticalFaceSurfFlux(BEFluxM + face * BENfp, BEfm + face * BENfp, BEnx + face * BENfp, BEny + face * BENfp, &gra, Hcrit, BENfp, Nvar, BENe);

				EvaluateVerticalFaceNumFlux_HLLC_LAI(BEFluxS + face * BENfp, BEfm + face * BENfp, BEfp + face * BENfp, \
					BEnx + face * BENfp, BEny + face * BENfp, &gra, Hcrit, BENfp, Nvar, BENe);
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		int adjacentE = (int)BEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
				continue;
			}
			else {
				for (int field = 0; field < Nvar; field++) {
					StrongFormBoundaryEdgeRHS(face, BEFToE, BEFToF, Np, K, BENfp, BEFToN1, BEFluxM + field * BENe*BENfp, BEFluxS + field * BENe*BENfp, BEJs, BEMb, ERHS + field * Np*K*Nface);
				}
			}
		}
	}
	/*************************************************************************************************************************************/
	/*************************************************************************************************************************************/
	           /**************************************Bottom Edge Part*******************************************************/
	                              /********************************************************************/
	double *omegaM = BotEfm + 2 * BotENfp*BotENe;
	huM = BotEfm, hvM = BotEfm + BotENfp*BotENe, hM = BotEfm + 3 * BotENfp*BotENe;
	double *omegaP = BotEfp + 2 * BotENfp*BotENe;
	huP = BotEfp, hvP = BotEfp + BotENfp*BotENe, hP = BotEfp + 3 * BotENfp*BotENe;
	/*Allocate memory for BotEFluxM, BotEFluxP and BotEFluxS, and calculate these flux term*/
	memset(BotEFluxM, 0, BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxP, 0, BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxS, 0, BotENfp*BotENe*Nvar*sizeof(double));
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++){
		int adjacentE = (int)BotEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
				/*Fetch variable BotEfm and BotEfp first*/
				FetchInnerEdgeFacialValue(hM + face * BotENfp, hP + face * BotENfp, h, BotEFToE + 2 * face, \
					BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);
				FetchInnerEdgeFacialValue(huM + face * BotENfp, huP + face * BotENfp, hu, BotEFToE + 2 * face, \
					BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);
				FetchInnerEdgeFacialValue(hvM + face * BotENfp, hvP + face * BotENfp, hv, BotEFToE + 2 * face, \
					BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);
				FetchInnerEdgeFacialValue(omegaM + face * BotENfp, omegaP + face * BotENfp, omega, BotEFToE + 2 * face, \
					BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);
				/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
				here 2 stands for the memory occupied by water depth h and omega*/
				for (int field = 2; field < Nvar; field++) {
					FetchInnerEdgeFacialValue(BotEfm + (field + 2)*BotENe*BotENfp + face * BotENfp, \
						BotEfp + (field + 2)*BotENe*BotENfp + face * BotENfp, fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
						BotEFToE + 2 * face, BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);
				}
				EvaluateHorizontalFaceSurfFlux(BotEFluxM + face * BotENfp, BotEfm + face * BotENfp, BotEnz + face * BotENfp, Hcrit, BotENfp, Nvar, BotENe);
				EvaluateHorizontalFaceSurfFlux(BotEFluxP + face * BotENfp, BotEfp + face * BotENfp, BotEnz + face * BotENfp, Hcrit, BotENfp, Nvar, BotENe);
				EvaluateHorizontalFaceNumFlux(BotEFluxS + face * BotENfp, BotEfm + face * BotENfp, BotEfp + face * BotENfp, \
					BotEnz + face * BotENfp, Hcrit, BotENfp, Nvar, BotENe);
			}
			else {
				continue;
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < BotENe; face++){
		int adjacentE = (int)BotEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
				for (int field = 0; field < Nvar; field++) {
					StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, BotEFluxM + field * BotENe*BotENfp, \
						BotEFluxP + field * BotENe*BotENfp, BotEFluxS + field * BotENe*BotENfp, BotEJs, BotEMb, ERHS + field * Np*K*Nface);
				}
			}
			else {
				continue;
			}
		}
	}
	/*************************************************************************************************************************************/
	/*************************************************************************************************************************************/
	             /**************************************Bottom Boundary Edge Part*******************************************************/
	                            /********************************************************************/
	huM = BotBEfm;
	hvM = BotBEfm + BotBENfp*BotBENe;
	omegaM = BotBEfm + 2 * BotBENfp*BotBENe;
	hM = BotBEfm + 3 * BotBENfp*BotBENe;

	/*Allocate memory for BotBEFluxM, fluxP and BotBEFluxS, and calculate these flux term*/
	memset(BotBEFluxM, 0, BotBENfp*BotBENe*Nvar*sizeof(double));
	memset(BotBEFluxS, 0, BotBENfp*BotBENe*Nvar*sizeof(double));
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		int adjacentE = (int)BotBEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
				FetchBoundaryEdgeFacialValue(huM + face * BotBENfp, hu, BotBEFToE + 2 * face, BotBEFToN1 + face * BotBENfp, Np, BotBENfp);
				FetchBoundaryEdgeFacialValue(hvM + face * BotBENfp, hv, BotBEFToE + 2 * face, BotBEFToN1 + face * BotBENfp, Np, BotBENfp);
				FetchBoundaryEdgeFacialValue(omegaM + face * BotBENfp, omega, BotBEFToE + 2 * face, BotBEFToN1 + face * BotBENfp, Np, BotBENfp);
				FetchBoundaryEdgeFacialValue(hM + face * BotBENfp, h, BotBEFToE + 2 * face, BotBEFToN1 + face * BotBENfp, Np, BotBENfp);

				/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
				2 stands for the memory occupied by water depth h and omega*/
				for (int field = 2; field < Nvar; field++) {
					FetchBoundaryEdgeFacialValue(BotBEfm + (field + 2)*BotBENe*BotBENfp + face * BotBENfp, \
						fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
						BotBEFToE + 2 * face, BotBEFToN1 + BotBENfp * face, Np, BotBENfp);
				}

				EvaluateHorizontalFaceSurfFlux(BotBEFluxM + face * BotBENfp, BotBEfm + face * BotBENfp, BotBEnz + face * BotBENfp, Hcrit, BotBENfp, Nvar, BotBENe);
			}
			else {
				continue;
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++){
		int adjacentE = (int)BotBEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
				for (int field = 0; field < Nvar; field++) {
					StrongFormBoundaryEdgeRHS(face, BotBEFToE, BotBEFToF, Np, K, BotBENfp, BotBEFToN1, BotBEFluxM + field * BotBENe*BotBENfp, BotBEFluxS + field * BotBENe*BotBENfp, BotBEJs, BotBEMb, ERHS + field * Np*K*Nface);
				}
			}
			else {
				continue;
			}
		}
	}
	/*************************************************************************************************************************************/
	/*************************************************************************************************************************************/
	         /**************************************Surface Boundary Edge Part*******************************************************/
	                               /********************************************************************/
	huM = SurfBEfm;
	hvM = SurfBEfm + SurfBENfp*SurfBENe;
	omegaM = SurfBEfm + 2 * SurfBENfp*SurfBENe;
	hM = SurfBEfm + 3 * SurfBENfp*SurfBENe;

	/*Allocate memory for SurfBEFluxM, fluxP and SurfBEFluxS, and calculate these flux term*/
	memset(SurfBEFluxM, 0, SurfBENfp*SurfBENe*Nvar*sizeof(double));
	/*WE NOTE THAT, FOR THE CONVERGENCE TEST, THIS FLUX IS NOT TAKEN AS ZERO AND SHOULD BE TAKEN FORM THE INPUT*/
	//double *SurfFluxS = mxGetPr(prhs[13]);
	memset(SurfBEFluxS, 0, SurfBENfp*SurfBENe*Nvar*sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++){
		int adjacentE = (int)SurfBEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
				FetchBoundaryEdgeFacialValue(huM + face * SurfBENfp, hu, SurfBEFToE + 2 * face, SurfBEFToN1 + face * SurfBENfp, Np, SurfBENfp);
				FetchBoundaryEdgeFacialValue(hvM + face * SurfBENfp, hv, SurfBEFToE + 2 * face, SurfBEFToN1 + face * SurfBENfp, Np, SurfBENfp);
				FetchBoundaryEdgeFacialValue(omegaM + face * SurfBENfp, omega, SurfBEFToE + 2 * face, SurfBEFToN1 + face * SurfBENfp, Np, SurfBENfp);
				FetchBoundaryEdgeFacialValue(hM + face * SurfBENfp, h, SurfBEFToE + 2 * face, SurfBEFToN1 + face * SurfBENfp, Np, SurfBENfp);

				/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included
				2 stands for the memory occupied by water depth h and omega*/
				for (int field = 2; field < Nvar; field++) {
					FetchBoundaryEdgeFacialValue(SurfBEfm + (field + 2)*SurfBENe*SurfBENfp + face * SurfBENfp, \
						fphys + ((int)varFieldIndex[field] - 1)*Np*K, \
						SurfBEFToE + 2 * face, SurfBEFToN1 + SurfBENfp * face, Np, SurfBENfp);
				}

				EvaluateHorizontalFaceSurfFlux(SurfBEFluxM + face * SurfBENfp, SurfBEfm + face * SurfBENfp, SurfBEnz + face * SurfBENfp, Hcrit, SurfBENfp, Nvar, SurfBENe);
			}
			else {
				continue;
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < SurfBENe; face++){
		int adjacentE = (int)SurfBEFToE[2 * face];
		if (MyID == pE3d[adjacentE - 1]) {
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
				for (int field = 0; field < Nvar; field++) {
					StrongFormBoundaryEdgeRHS(face, SurfBEFToE, SurfBEFToF, Np, K, SurfBENfp, SurfBEFToN1, SurfBEFluxM + field * SurfBENe*SurfBENfp, SurfBEFluxS + field * SurfBENe*SurfBENfp, SurfBEJs, SurfBEMb, ERHS + field * Np*K*Nface);
				}
			}
			else {
				continue;
			}
		}
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k < K; k++){
		if (MyID == pE3d[k]) {
			for (int field = 0; field < Nvar; field++) {
				for (int face = 1; face < Nface; face++) {
					Add(ERHS + field * Np*K*Nface + k * Np, ERHS + field * Np*K*Nface + k * Np, ERHS + field * Np*K*Nface + face * Np*K + k * Np, Np);
				}
			}
		}
    }    
    

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if (MyID == pE3d[k]) {
			for (int field = 0; field < Nvar; field++) {
				MultiEdgeContributionByLiftOperator(ERHS + field * Np*K*Nface + k * Np, TempFacialIntegral + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
			}
		}
	}

	/***************************************************************************************************************/

	                                                                                    /************************  Volume Integral Part  ****************************/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		if (MyID == pE3d[k]) {
			if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
				continue;
			}
			else if ((NdgRegionType)Status3d[k] == NdgRegionPartialWetFlood) {
				EvaluatePrebalanceVolumeTerm(E + k * Np, G + k * Np, H + k * Np, fphys + k * Np, \
					varFieldIndex, Nvar, &no_gra, Np, K, Hcrit);

				GetVolumnIntegral3d(OutputRHS + k * Np, TempVolumeIntegral + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, Dt, E + k * Np, G + k * Np, H + k * Np, &np, &np, &zero, \
					&np, rx + k * Np, sx + k * Np, ry + k * Np, sy + k * Np, tz + k * Np, Nvar, Np, K);
			}
			else { //DamBreak or wet
				EvaluatePrebalanceVolumeTerm(E + k * Np, G + k * Np, H + k * Np, fphys + k * Np, \
					varFieldIndex, Nvar, &gra, Np, K, Hcrit);

				GetVolumnIntegral3d(OutputRHS + k * Np, TempVolumeIntegral + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, Dt, E + k * Np, G + k * Np, H + k * Np, &np, &np, &zero, \
					&np, rx + k * Np, sx + k * Np, ry + k * Np, sy + k * Np, tz + k * Np, Nvar, Np, K);
			}
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		if (MyID == pE3d[k]) {
			for (int n = 0; n < Nvar; n++) {
				Minus(OutputRHS + n * Np*K + k * Np, \
					ERHS + n * Np*K*Nface + k * Np, OutputRHS + n * Np*K + k * Np, Np);
			}
		}
	}

}
