#include "NdgMath.h"
#include "NdgSWE.h"
//#include "NdgMemory.h"
#include "NdgSWEHorizSmagrinskyDiffSolver.h"
#include <cmath>
//#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

//#define max(a, b) ((a > b) ? a : b)
//#define min(a, b) ((a < b) ? a : b)

extern double *HorDiffnv, *HorDiffvariable, *HorDiffBEfp, *HorDiffzM, \
*HorDiffzP, *HorDiffTempBEfp, *HorDiffTempBEfm, *HorDiffAVx, \
*HorDiffAVy, *HorDiffVx, *HorDiffTempVx, *HorDiffVy, *HorDiffTempVy, \
*HorDiffIEfm, *HorDiffAVIEfm, *HorDiffIEfp, *HorDiffAVIEfp, *HorDiffIEFluxM, \
*HorDiffIEFluxP, *HorDiffBEfm, *HorDiffIEFluxS, *HorDiffAVBEfm, *HorDiffBEFluxM, \
*HorDiffBEFluxS, *HorDiffERHSX, *HorDiffERHSY, *HorDiffLocalPrimitiveDiffTermX, \
*HorDiffLocalPrimitiveDiffTermY, *HorDiffLPDTIEfm, *HorDiffLPDTIEfp, *HorDiffLPDTBEfm, \
*HorDiffTempFacialIntegralX, *HorDiffTempFacialIntegralY, *HorDiffInnerEdgeTau, \
*HorDiffBoundaryEdgeTau, *HorDiffIEnvfm, *HorDiffIEnvfp, *HorDiffBEnvfm, *u_u, *v_v, \
*rx_Dr_u, *rx_Dr_v, *sx_Ds_u, *sx_Ds_v, *ry_Dr_u, *ry_Dr_v, *sy_Ds_u, *sy_Ds_v,\
*Hrms, *WaveNumber_, *nv_c, *nv_w, *nv_cw;


int TimeStep = 0;

//extern char *HorDiffInitialized;
extern double gra;
extern double Hcrit;
extern int  Nvar;
extern signed char *Status3d;
//HorizDiffMemoryDeAllocation();

NdgSWEHorizSmagrinskyDiffSolver::NdgSWEHorizSmagrinskyDiffSolver(double c_) :C(c_)
{
}

NdgSWEHorizSmagrinskyDiffSolver::~NdgSWEHorizSmagrinskyDiffSolver()
{
};

void NdgSWEHorizSmagrinskyDiffSolver::EvaluateDiffRHS_Nowave(double *fphys_, double *frhs_, double *fext_, int *varFieldIndex_) {
	//Order of the input variable is 0 hcrit, 1 meshUnion.type, 2 prantl, 3 InnerEdge, 4 BoundaryEdge, 5 nv, 6 frhs, 7 fphys, 8 varIndex, 9 cell, 10 mesh, 11 BoundaryEdgefp
	//mexAtExit(&MyExit);
	//NdgMeshType type = meshunion->type;
	double Prantl = 1.0;//chenzereng
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
	double *IELAV = meshunion->inneredge_p->LAV;

	int BENe = *(meshunion->boundaryedge_p->Ne);
	int BENfp = *(meshunion->boundaryedge_p->Nfp);
	double *BEMb = meshunion->boundaryedge_p->M;
	double *BEJs = meshunion->boundaryedge_p->Js;
	double *BEnx = meshunion->boundaryedge_p->nx;
	double *BEny = meshunion->boundaryedge_p->ny;
	double *BEFToE = meshunion->boundaryedge_p->FToE;
	double *BEFToF = meshunion->boundaryedge_p->FToF;
	double *BEFToN1 = meshunion->boundaryedge_p->FToN1;
	double *BELAV = meshunion->boundaryedge_p->LAV;

	//const size_t *PRHS;
	//PRHS = mxGetDimensions(prhs[6]);
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	double *fphys = fphys_;
	int *varIndex = varFieldIndex_;

	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	int Nface3d = *(meshunion->cell_p->Nface);
	double *invM = meshunion->cell_p->invM;
	int Nface;

	double *rx = meshunion->rx;
	double *sx = meshunion->sx;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;
	double *J = meshunion->J;
	double *LAV = meshunion->LAV;
	double *ftype = meshunion->boundaryedge_p->ftype;
	double *fext = fext_;
    
	//double *Tempnv = new double[Np * K];

	/*Set the output right hand side*/
	double *OutputRHS = frhs_;
	int Nfield;
	//mxArray *TempOrder;
	int MeshType, Order;
	double *hu = NULL, *hv = NULL, *h = NULL, *z = NULL;
	double *huM = NULL, *hvM = NULL, *hM = NULL;

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
	//if (!strcmp("False", HorDiffInitialized)){

		//HorizDiffMemoryAllocation( type, Np, K, Nvar, Nface3d, BENfp, BENe, IENfp, IENe);

	//}

	/*Allocate memory and calcualte variable u, v, $\theta$. Then impose the boudary condition and calculate u, v, $\theta$ defined over the ghost edge*/
//	if (type == Two){
//		/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
//		Nfield = Nvar - 1;
	//if (type == Three){
		Nfield = Nvar;
//#ifdef _BAROCLINIC
//		Nfield = 2;//For cases that set the horizontal diffusion coefficient to zero
//#endif
        /*For 3d shallow water problem, the face number is equal to TempNface - 2, since there 
         * is surface edge and bottom edge is not considered for horizontal diffusion term*/
        Nface = Nface3d - 2;        
		MeshType = 3;
		//TempOrder = mxGetField(prhs[9], 0, "N");
		//Order = meshunion->cell_p->N;
		//TempOrder = mxGetField(prhs[9], 0, "Nz");
		//Order = max(Order, (int)mxGetScalar(TempOrder));
		Order = fmax(*meshunion->cell_p->N,* meshunion->cell_p->Nz);
		
		UpdateViscosity_Nowave(fphys, nv_cw);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		/*  Calculate nv*h  */
		for (int k = 0; k < K; k++){
			DotProduct(HorDiffnv + k*Np, nv_cw + k*Np, fphys + 3 * Np*K + k*Np, Np);
			for (int field = 0; field < Nfield; field++){
				DotCriticalDivide(HorDiffvariable + field*Np*K + k*Np, \
					fphys + (int)(varIndex[field] - 1)*Np*K + k*Np, &Hcrit, \
					fphys + 3*Np*K + k*Np, Np);//For 3d shallow water problem, variable about height is organized as the forth variable
			}
		}
		huM = HorDiffTempBEfm, hvM = HorDiffTempBEfm + BENe*BENfp, hM = HorDiffTempBEfm + 2 * BENe*BENfp;
		hu = fphys;
		hv = fphys + Np*K;
		h = fphys + 3 * Np*K;
		z = fphys + 5 * Np*K;
		/*Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.
		Finally, calculate local flux term, adjacent flux term and numerical flux term*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int face = 0; face < BENe; face++){
			int adjacentE = (int)BEFToE[2 * face];
			if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
				continue; //BEFToE is just itself
			}
			else {
				NdgEdgeType facetype = (NdgEdgeType)ftype[face];  // boundary condition
				FetchBoundaryEdgeFacialValue(huM + face * BENfp, hu, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				FetchBoundaryEdgeFacialValue(hvM + face * BENfp, hv, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				FetchBoundaryEdgeFacialValue(hM + face * BENfp, h, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				FetchBoundaryEdgeFacialValue(HorDiffzM + face * BENfp, z, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
				/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
				here 1 stands for the memory occupied by water depth h*/
				for (int field = 2; field < Nfield; field++) {
					FetchBoundaryEdgeFacialValue(HorDiffTempBEfm + (field + 1)*BENe*BENfp + face * BENfp, \
						fphys + ((int)varIndex[field] - 1)*Np*K, \
						BEFToE + 2 * face, BEFToN1 + BENfp * face, Np, BENfp);
				}
				/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
				ImposeBoundaryCondition(&gra, facetype, BEnx + face * BENfp, BEny + face * BENfp, HorDiffTempBEfm + face * BENfp, HorDiffTempBEfp + face * BENfp, \
					HorDiffzM + face * BENfp, HorDiffzP + face * BENfp, fext + face * BENfp, BENfp, Nfield + 1, BENe);
				/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
				EvaluateHydroStaticReconstructValue(Hcrit, HorDiffTempBEfm + face * BENfp, HorDiffTempBEfp + face * BENfp, HorDiffzM + face * BENfp, HorDiffzP + face * BENfp, BENfp, Nfield + 1, BENe);
				/*We divide the variable by water depth to get the original variable*/
				for (int field = 0; field < 2; field++) {
					DotCriticalDivide(HorDiffBEfp + field * BENfp*BENe + face * BENfp, \
						HorDiffTempBEfp + field * BENfp*BENe + face * BENfp, &Hcrit, \
						HorDiffTempBEfp + 2 * BENfp*BENe + face * BENfp, BENfp);
				}
				/*Water depth is stored as the third variable*/
				for (int field = 2; field < Nfield; field++) {
					DotCriticalDivide(HorDiffBEfp + field * BENfp*BENe + face * BENfp, \
						HorDiffTempBEfp + (field + 1)*BENfp*BENe + face * BENfp, &Hcrit, \
						HorDiffTempBEfp + 2 * BENfp*BENe + face * BENfp, BENfp);
				}
			}

		}
	//}
	//delete [] Tempnv;

	memset(HorDiffERHSX, 0, Np*K*Nfield*Nface*sizeof(double));
	memset(HorDiffERHSY, 0, Np*K*Nfield*Nface*sizeof(double));
   
	int np = Np;
	int oneI = 1;
	double one = 1.0, zero = 0.0;
	/*Volume integral part*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				/*$\bold{r_x}\cdot (Dr*u(v,\theta))+\bold{s_x}\cdot (Ds*u(v,\theta))$*/
				GetVolumnIntegral2d(HorDiffVx + field * Np*K + k * Np, HorDiffTempVx + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffvariable + field * Np*K + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);
				/*$\bold{r_y}\cdot (Dr*u(v,\theta))+\bold{s_y}\cdot (Ds*u(v,\theta))$*/
				GetVolumnIntegral2d(HorDiffVy + field * Np*K + k * Np, HorDiffTempVy + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffvariable + field * Np*K + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);
			}
		}

	}
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				DotProduct(HorDiffLocalPrimitiveDiffTermX + field * Np*K + k * Np, HorDiffVx + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
				DotProduct(HorDiffLocalPrimitiveDiffTermY + field * Np*K + k * Np, HorDiffVy + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
			}
			/*for substance transport, prantl number is considered*/
			for (int field = 2; field < Nfield; field++) {
				DotDivideByConstant(HorDiffLocalPrimitiveDiffTermX + field * Np*K + k * Np, HorDiffLocalPrimitiveDiffTermX + field * Np*K + k * Np, Prantl, Np);
				DotDivideByConstant(HorDiffLocalPrimitiveDiffTermY + field * Np*K + k * Np, HorDiffLocalPrimitiveDiffTermY + field * Np*K + k * Np, Prantl, Np);
			}
		}

	}

/*Inner edge facial integral part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < IENe; face++){
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				FetchInnerEdgeFacialValue(HorDiffIEfm + field * IENe*IENfp + face * IENfp, HorDiffIEfp + field * IENe*IENfp + face * IENfp, HorDiffvariable + field * Np*K, IEFToE + 2 * face, \
					IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
				/*Inner edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$*/
				GetIEContributionToAuxialaryVariable(HorDiffERHSX + field * Np*K*Nface, face, IENe, IENfp, field, HorDiffIEfm, HorDiffIEfp, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, HorDiffIEFluxM, HorDiffIEFluxP, HorDiffIEFluxS, IEnx, IEMb, IEJs);
				/*Inner edge contribution to RHSY of $\frac{\partial u(v,\theta)}{\partial y}$*/
				GetIEContributionToAuxialaryVariable(HorDiffERHSY + field * Np*K*Nface, face, IENe, IENfp, field, HorDiffIEfm, HorDiffIEfp, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, HorDiffIEFluxM, HorDiffIEFluxP, HorDiffIEFluxS, IEny, IEMb, IEJs);
			}
		}
	}
	/*Boundary edge facial integral part*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < BENe; face++){
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				FetchBoundaryEdgeFacialValue(HorDiffBEfm + field * BENe*BENfp + face * BENfp, HorDiffvariable + field * Np*K, BEFToE + 2 * face, BEFToN1 + BENfp * face, Np, BENfp);
				/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$, we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToAuxialaryVariable(HorDiffERHSX + field * Np*K*Nface, face, BENe, BENfp, field, HorDiffBEfm, HorDiffBEfp, BEFToE, BEFToF, BEFToN1, Np, K, HorDiffBEFluxM, HorDiffBEFluxS, BEnx, BEMb, BEJs);
				/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial y}$, we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToAuxialaryVariable(HorDiffERHSY + field * Np*K*Nface, face, BENe, BENfp, field, HorDiffBEfm, HorDiffBEfp, BEFToE, BEFToF, BEFToN1, Np, K, HorDiffBEFluxM, HorDiffBEFluxS, BEny, BEMb, BEJs);
			}
		}

	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k<K; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				for (int face = 1; face < Nface; face++) {
					Add(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + face * Np*K + k * Np, Np);
					Add(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + face * Np*K + k * Np, Np);
				}
			}
		}
    }

	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {

				MultiEdgeContributionByLiftOperator(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralX + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);

				MultiEdgeContributionByLiftOperator(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralY + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
			}
		}
	}
	/*Next, sum contribution from volume integral, inner edge contribution and boundary edge contribution into HorDiffAVx and HorDiffAVy*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				Minus(HorDiffAVx + field * Np*K + k * Np, \
					HorDiffVx + field * Np*K + k * Np, HorDiffERHSX + field * Np*K*Nface + k * Np, Np);
				Minus(HorDiffAVy + field * Np*K + k * Np, \
					HorDiffVy + field * Np*K + k * Np, HorDiffERHSY + field * Np*K*Nface + k * Np, Np);
			}
		}
	}
	/*Finally, multiply each component of HorDiffAVx and HorDiffAVy by its diffusion coefficient to get the final auxialary variable,
	this is conducted according to $(Q,v)=(\nv\hat Q,v)$,where $\hat Q = \frac{\partial u(v,\theta)}{\partial (x,y)}$,
	we note that in this projection precedure, aliasing error is introduced.*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				DotProduct(HorDiffAVx + field * Np*K + k * Np, HorDiffAVx + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
				DotProduct(HorDiffAVy + field * Np*K + k * Np, HorDiffAVy + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
			}
			for (int field = 2; field < Nfield; field++) {
				DotDivideByConstant(HorDiffAVx + field * Np*K + k * Np, HorDiffAVx + field * Np*K + k * Np, Prantl, Np);
				DotDivideByConstant(HorDiffAVy + field * Np*K + k * Np, HorDiffAVy + field * Np*K + k * Np, Prantl, Np);
			}
		}

	}

	/*Calculate the contribution to the right hand side due to the auxialary variable HorDiffAVx and HorDiffAVy with IPDG.*/
	/*Calculate the penalty parameter $\tau$ first, this parameter is calculated as $\tau=\frac{(N+1)(N+d)}{d}\frac{n_0}{2}\frac{A}{V}\nv$*/
	/*Inner edge first*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++){
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			FetchInnerEdgeFacialValue(HorDiffIEnvfm + face * IENfp, HorDiffIEnvfp + face * IENfp, HorDiffnv, \
				IEFToE + face * 2, IEFToN1 + face * IENfp, IEFToN2 + face * IENfp, Np, IENfp);
			double localRatio = IELAV[face] / LAV[(int)IEFToE[face * 2] - 1];
			double adjacentRatio = IELAV[face] / LAV[(int)IEFToE[face * 2 + 1] - 1];
			for (int p = 0; p < IENfp; p++) {
				HorDiffInnerEdgeTau[face*IENfp + p] = fmax(localRatio*(Order + 1)*(Order + MeshType) / MeshType * Nface / 2 * HorDiffIEnvfm[face*IENfp + p], \
					adjacentRatio*(Order + 1)*(Order + MeshType) / MeshType * Nface / 2 * HorDiffIEnvfp[face*IENfp + p]);
			}
		}
	}
	/*Boundary edge next*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++){
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			FetchBoundaryEdgeFacialValue(HorDiffBEnvfm + face * BENfp, HorDiffnv, BEFToE + face * 2, BEFToN1 + face * BENfp, Np, BENfp);
			double localRatio = BELAV[face] / LAV[(int)BEFToE[face * 2] - 1];
			for (int p = 0; p < BENfp; p++) {
				HorDiffBoundaryEdgeTau[face*BENfp + p] = localRatio * (Order + 1)*(Order + MeshType) / MeshType * Nface / 2 * HorDiffBEnvfm[face*IENfp + p];
			}
		}
	}

	/*Volume integral of the second order operator*/
	
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++){
		/*Prantl number is not considered here, since we have considered this parameter when calculating the auxialary variable,
		it is not need anymore in the volumn integral part when considering the contribution of second order operator to the RHS*/
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				/*$\bold{r_x}\cdot (Dr*Q_x)+\bold{s_x}\cdot (Ds*Q_x)$*/
				GetVolumnIntegral2d(HorDiffVx + field * Np*K + k * Np, HorDiffTempVx + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffAVx + field * Np*K + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);
				/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$*/
				GetVolumnIntegral2d(HorDiffVy + field * Np*K + k * Np, HorDiffTempVy + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffAVy + field * Np*K + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);
				//if (type == Two){
				//	/*The water depth field is excluded from this part*/
				//	Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		HorDiffVx + field*Np*K + k*Np, Np);
				//	Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		HorDiffVy + field*Np*K + k*Np, Np);
				//}
				//else if (type == Three){
				Add(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffVx + field * Np*K + k * Np, Np);
				Add(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffVy + field * Np*K + k * Np, Np);
				//}
			}
		}
    }
/************************************************************************************************************************************/

///***********************************************************************************************************************************************/
	/*Reset all the data contained in HorDiffERHSX and HorDiffERHSY to zero, becaused these space contains the data left when calculating the auxialary variable*/
	memset(HorDiffERHSX, 0, Np*K*Nfield*Nface*sizeof(double));
	memset(HorDiffERHSY, 0, Np*K*Nfield*Nface*sizeof(double));
	/*Surface integral of the second order operator, this part is calculated seperately for field with index less than 2 and field with index larger than 2 since prantl number should be considered*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < IENe; face++){
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < 2; field++) {
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
				GetIEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEnx, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEnx, HorDiffInnerEdgeTau, 1.0, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial y}$*/
				GetIEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEny, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEny, HorDiffInnerEdgeTau, 1.0, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);

			}
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < IENe; face++){
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 2; field < Nfield; field++) {
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$*/
				GetIEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEnx, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEnx, HorDiffInnerEdgeTau, Prantl, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$*/
				GetIEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEny, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEny, HorDiffInnerEdgeTau, Prantl, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);

			}
		}
	}
	/*Boundary edge facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < BENe; face++){
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < 2; field++) {
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
				GetBEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVx + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial (u,v)}{\partial y}$*/
				GetBEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVy + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
			}
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < BENe; face++){
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 2; field < Nfield; field++) {
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$,
				we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, Prantl, HorDiffAVBEfm, HorDiffAVx + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$,
				we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, Prantl, HorDiffAVBEfm, HorDiffAVy + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
			}
		}
	}
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k<K; k++){
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				for (int face = 1; face < Nface; face++) {
					Add(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + face * Np*K + k * Np, Np);
					Add(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + face * Np*K + k * Np, Np);
				}
			}
		}
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {

				MultiEdgeContributionByLiftOperator(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralX + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);

				MultiEdgeContributionByLiftOperator(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralY + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
			}
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {

				Minus(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffERHSX + field * Np*K*Nface + k * Np, Np);
				Minus(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffERHSX + field * Np*K*Nface + k * Np, Np);

			}
		}
	}
}

void NdgSWEHorizSmagrinskyDiffSolver::UpdateViscosity_Nowave(double *fphys, double *nv_) {
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	double *nv_cw = nv_;
	double *hu = fphys;
	double *hv = fphys + K*Np;
	double *h = fphys + K * Np*3;
	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	double *rx = meshunion->rx;
	double *sx = meshunion->sx;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;
	double *LAV = meshunion->LAV;
	int *Nlayer = meshunion->Nlayer;
	int oneI = 1;

	memset(u_u, 0, Np * K * sizeof(double));
	memset(v_v, 0, Np * K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {	
		for (int j = 0; j < Np; j++) {
			/*  Calculate u,v in a cell   */
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i * Np + j, hu + i * Np + j, u_u + i * Np + j);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i * Np + j, hv + i * Np + j, v_v + i * Np + j);
		}
	}

	memset(rx_Dr_u, 0, Np * K * sizeof(double));
	memset(rx_Dr_v, 0, Np * K * sizeof(double));
	memset(sx_Ds_u, 0, Np * K * sizeof(double));
	memset(sx_Ds_v, 0, Np * K * sizeof(double));
	memset(ry_Dr_u, 0, Np * K * sizeof(double));
	memset(ry_Dr_v, 0, Np * K * sizeof(double));
	memset(sy_Ds_u, 0, Np * K * sizeof(double));
	memset(sy_Ds_v, 0, Np * K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		/*  Calculate MatrixMultiply of Dr*u INTO ry_Dr_u, Dr*v INTO ry_Dr_v, Ds*u INTO sy_Ds_u, Ds*v INTO sy_Ds_v in a cell   */
		MatrixMultiply(Dr, u_u + i * Np, ry_Dr_u + i * Np, Np, oneI, Np, 1.0);
		MatrixMultiply(Dr, v_v + i * Np, ry_Dr_v + i * Np, Np, oneI, Np, 1.0);
		MatrixMultiply(Ds, u_u + i * Np, sy_Ds_u + i * Np, Np, oneI, Np, 1.0);
		MatrixMultiply(Ds, v_v + i * Np, sy_Ds_v + i * Np, Np, oneI, Np, 1.0);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		/*  Calculate Dotproduct of rx_Dr*u to sy_Ds*v in a cell   */
		DotProduct(rx_Dr_u + i * Np, rx, ry_Dr_u + i * Np, Np);
		DotProduct(rx_Dr_v + i * Np, rx, ry_Dr_v + i * Np, Np);
		DotProduct(sx_Ds_u + i * Np, sx, sy_Ds_u + i * Np, Np);
		DotProduct(sx_Ds_v + i * Np, sx, sy_Ds_v + i * Np, Np);
		DotProduct(ry_Dr_u + i * Np, ry, ry_Dr_u + i * Np, Np);
		DotProduct(ry_Dr_v + i * Np, ry, ry_Dr_v + i * Np, Np);
		DotProduct(sy_Ds_u + i * Np, sy, sy_Ds_u + i * Np, Np);
		DotProduct(sy_Ds_v + i * Np, sy, sy_Ds_v + i * Np, Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < Np; j++) {
			nv_cw[i*Np + j] = C * LAV[i] * (*Nlayer) * sqrt(pow((rx_Dr_u[i*Np + j] + sx_Ds_u[i*Np + j]), 2) + 0.5* pow((rx_Dr_v[i*Np + j]+ sx_Ds_v[i*Np + j]+ ry_Dr_u[i*Np + j] + sy_Ds_u[i*Np + j]), 2) + pow((ry_Dr_v[i*Np + j]+ sy_Ds_v[i*Np + j]), 2));
			//nv_cw[i*Np + j] = 200.0;//For SaltyWater Case,1,10,200
		}

	}
	//std::cout << "nv_c = " << nv[0] << " and " << nv[7] << " and " << nv[10] << std::endl;
}


void NdgSWEHorizSmagrinskyDiffSolver::EvaluateDiffRHS(double *fphys_, double *frhs_, double *fext_, int *varFieldIndex_, double *time,  double *HS, double *WLEN, double *UBOT) {
	//Order of the input variable is 0 hcrit, 1 meshUnion.type, 2 prantl, 3 InnerEdge, 4 BoundaryEdge, 5 nv, 6 frhs, 7 fphys, 8 varIndex, 9 cell, 10 mesh, 11 BoundaryEdgefp
	//mexAtExit(&MyExit);
	//NdgMeshType type = meshunion->type;
	double Prantl = 1.0;//chenzereng
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
	double *IELAV = meshunion->inneredge_p->LAV;

	int BENe = *(meshunion->boundaryedge_p->Ne);
	int BENfp = *(meshunion->boundaryedge_p->Nfp);
	double *BEMb = meshunion->boundaryedge_p->M;
	double *BEJs = meshunion->boundaryedge_p->Js;
	double *BEnx = meshunion->boundaryedge_p->nx;
	double *BEny = meshunion->boundaryedge_p->ny;
	double *BEFToE = meshunion->boundaryedge_p->FToE;
	double *BEFToF = meshunion->boundaryedge_p->FToF;
	double *BEFToN1 = meshunion->boundaryedge_p->FToN1;
	double *BELAV = meshunion->boundaryedge_p->LAV;

	//const size_t *PRHS;
	//PRHS = mxGetDimensions(prhs[6]);
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	double *fphys = fphys_;
	int *varIndex = varFieldIndex_;

	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	int Nface3d = *(meshunion->cell_p->Nface);
	double *invM = meshunion->cell_p->invM;
	int Nface;

	double *rx = meshunion->rx;
	double *sx = meshunion->sx;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;
	double *J = meshunion->J;
	double *LAV = meshunion->LAV;
	double *ftype = meshunion->boundaryedge_p->ftype;
	double *fext = fext_;

	//double *Tempnv = new double[Np * K];

	/*Set the output right hand side*/
	double *OutputRHS = frhs_;
	int Nfield;
	//mxArray *TempOrder;
	int MeshType, Order;
	double *hu = NULL, *hv = NULL, *h = NULL, *z = NULL;
	double *huM = NULL, *hvM = NULL, *hM = NULL;

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


	Nfield = Nvar;
//#ifdef _BAROCLINIC
//	Nfield = 2;//For cases that set the horizontal diffusion coefficient to zero
//#endif
	/*For 3d shallow water problem, the face number is equal to TempNface - 2, since there
	 * is surface edge and bottom edge is not considered for horizontal diffusion term*/
	Nface = Nface3d - 2;
	MeshType = 3;
	//TempOrder = mxGetField(prhs[9], 0, "N");
	//Order = meshunion->cell_p->N;
	//TempOrder = mxGetField(prhs[9], 0, "Nz");
	//Order = max(Order, (int)mxGetScalar(TempOrder));
	Order = fmax(*meshunion->cell_p->N, *meshunion->cell_p->Nz);

	UpdateViscosity(fphys, nv_cw, HS, WLEN, UBOT,time);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	/*  Calculate nv*h  */
	for (int k = 0; k < K; k++) {
		DotProduct(HorDiffnv + k * Np, nv_cw + k * Np, fphys + 3 * Np*K + k * Np, Np);
		for (int field = 0; field < Nvar; field++) {
			DotCriticalDivide(HorDiffvariable + field * Np*K + k * Np, \
				fphys + (int)(varIndex[field] - 1)*Np*K + k * Np, &Hcrit, \
				fphys + 3 * Np*K + k * Np, Np);//For 3d shallow water problem, variable about height is organized as the forth variable
		}
	}
	huM = HorDiffTempBEfm, hvM = HorDiffTempBEfm + BENe * BENfp, hM = HorDiffTempBEfm + 2 * BENe*BENfp;
	hu = fphys;
	hv = fphys + Np * K;
	h = fphys + 3 * Np*K;
	z = fphys + 5 * Np*K;
	/*Fetch variable fm and fp first, then impose boundary condition and conduct hydrostatic reconstruction.
	Finally, calculate local flux term, adjacent flux term and numerical flux term*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++) {
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			NdgEdgeType facetype = (NdgEdgeType)ftype[face];  // boundary condition
			FetchBoundaryEdgeFacialValue(huM + face * BENfp, hu, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hvM + face * BENfp, hv, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(hM + face * BENfp, h, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
			FetchBoundaryEdgeFacialValue(HorDiffzM + face * BENfp, z, BEFToE + 2 * face, BEFToN1 + face * BENfp, Np, BENfp);
			/*The following part is used to fetch the field corresponding to temperature, salinity, and sediment if they are included,
			here 1 stands for the memory occupied by water depth h*/
			for (int field = 2; field < Nvar; field++) {
				FetchBoundaryEdgeFacialValue(HorDiffTempBEfm + (field + 1)*BENe*BENfp + face * BENfp, \
					fphys + ((int)varIndex[field] - 1)*Np*K, \
					BEFToE + 2 * face, BEFToN1 + BENfp * face, Np, BENfp);
			}
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			ImposeBoundaryCondition(&gra, facetype, BEnx + face * BENfp, BEny + face * BENfp, HorDiffTempBEfm + face * BENfp, HorDiffTempBEfp + face * BENfp, \
				HorDiffzM + face * BENfp, HorDiffzP + face * BENfp, fext + face * BENfp, BENfp, Nfield + 1, BENe);
			/*Water depth needs to be considered, so we plus Nfield by one to consider the water depth field*/
			EvaluateHydroStaticReconstructValue(Hcrit, HorDiffTempBEfm + face * BENfp, HorDiffTempBEfp + face * BENfp, HorDiffzM + face * BENfp, HorDiffzP + face * BENfp, BENfp, Nfield + 1, BENe);
			/*We divide the variable by water depth to get the original variable*/
			for (int field = 0; field < 2; field++) {
				DotCriticalDivide(HorDiffBEfp + field * BENfp*BENe + face * BENfp, \
					HorDiffTempBEfp + field * BENfp*BENe + face * BENfp, &Hcrit, \
					HorDiffTempBEfp + 2 * BENfp*BENe + face * BENfp, BENfp);
			}
			/*Water depth is stored as the third variable*/
			for (int field = 2; field < Nfield; field++) {
				DotCriticalDivide(HorDiffBEfp + field * BENfp*BENe + face * BENfp, \
					HorDiffTempBEfp + (field + 1)*BENfp*BENe + face * BENfp, &Hcrit, \
					HorDiffTempBEfp + 2 * BENfp*BENe + face * BENfp, BENfp);
		}
	}

}
	//}
	//delete [] Tempnv;

	memset(HorDiffERHSX, 0, Np*K*Nfield*Nface * sizeof(double));
	memset(HorDiffERHSY, 0, Np*K*Nfield*Nface * sizeof(double));

	int np = Np;
	int oneI = 1;
	double one = 1.0, zero = 0.0;
	/*Volume integral part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				/*$\bold{r_x}\cdot (Dr*u(v,\theta))+\bold{s_x}\cdot (Ds*u(v,\theta))$*/
				GetVolumnIntegral2d(HorDiffVx + field * Np*K + k * Np, HorDiffTempVx + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffvariable + field * Np*K + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);
				/*$\bold{r_y}\cdot (Dr*u(v,\theta))+\bold{s_y}\cdot (Ds*u(v,\theta))$*/
				GetVolumnIntegral2d(HorDiffVy + field * Np*K + k * Np, HorDiffTempVy + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffvariable + field * Np*K + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);
			}
		}

	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				DotProduct(HorDiffLocalPrimitiveDiffTermX + field * Np*K + k * Np, HorDiffVx + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
				DotProduct(HorDiffLocalPrimitiveDiffTermY + field * Np*K + k * Np, HorDiffVy + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
	}
			/*for substance transport, prantl number is considered*/
			for (int field = 2; field < Nfield; field++) {
				DotDivideByConstant(HorDiffLocalPrimitiveDiffTermX + field * Np*K + k * Np, HorDiffLocalPrimitiveDiffTermX + field * Np*K + k * Np, Prantl, Np);
				DotDivideByConstant(HorDiffLocalPrimitiveDiffTermY + field * Np*K + k * Np, HorDiffLocalPrimitiveDiffTermY + field * Np*K + k * Np, Prantl, Np);
			}
		}

	}

	/*Inner edge facial integral part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				FetchInnerEdgeFacialValue(HorDiffIEfm + field * IENe*IENfp + face * IENfp, HorDiffIEfp + field * IENe*IENfp + face * IENfp, HorDiffvariable + field * Np*K, IEFToE + 2 * face, \
					IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
				/*Inner edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$*/
				GetIEContributionToAuxialaryVariable(HorDiffERHSX + field * Np*K*Nface, face, IENe, IENfp, field, HorDiffIEfm, HorDiffIEfp, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, HorDiffIEFluxM, HorDiffIEFluxP, HorDiffIEFluxS, IEnx, IEMb, IEJs);
				/*Inner edge contribution to RHSY of $\frac{\partial u(v,\theta)}{\partial y}$*/
				GetIEContributionToAuxialaryVariable(HorDiffERHSY + field * Np*K*Nface, face, IENe, IENfp, field, HorDiffIEfm, HorDiffIEfp, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, HorDiffIEFluxM, HorDiffIEFluxP, HorDiffIEFluxS, IEny, IEMb, IEJs);
			}
		}
	}
	/*Boundary edge facial integral part*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++) {
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
	}
		else {
			for (int field = 0; field < Nfield; field++) {
				FetchBoundaryEdgeFacialValue(HorDiffBEfm + field * BENe*BENfp + face * BENfp, HorDiffvariable + field * Np*K, BEFToE + 2 * face, BEFToN1 + BENfp * face, Np, BENfp);
				/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial x}$, we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToAuxialaryVariable(HorDiffERHSX + field * Np*K*Nface, face, BENe, BENfp, field, HorDiffBEfm, HorDiffBEfp, BEFToE, BEFToF, BEFToN1, Np, K, HorDiffBEFluxM, HorDiffBEFluxS, BEnx, BEMb, BEJs);
				/*Boundary edge contribution to RHSX of $\frac{\partial u(v,\theta)}{\partial y}$, we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToAuxialaryVariable(HorDiffERHSY + field * Np*K*Nface, face, BENe, BENfp, field, HorDiffBEfm, HorDiffBEfp, BEFToE, BEFToF, BEFToN1, Np, K, HorDiffBEFluxM, HorDiffBEFluxS, BEny, BEMb, BEJs);
			}
		}

	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				for (int face = 1; face < Nface; face++) {
					Add(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + face * Np*K + k * Np, Np);
					Add(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + face * Np*K + k * Np, Np);
				}
			}
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {

				MultiEdgeContributionByLiftOperator(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralX + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);

				MultiEdgeContributionByLiftOperator(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralY + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
			}
		}
	}
	/*Next, sum contribution from volume integral, inner edge contribution and boundary edge contribution into HorDiffAVx and HorDiffAVy*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				Minus(HorDiffAVx + field * Np*K + k * Np, \
					HorDiffVx + field * Np*K + k * Np, HorDiffERHSX + field * Np*K*Nface + k * Np, Np);
				Minus(HorDiffAVy + field * Np*K + k * Np, \
					HorDiffVy + field * Np*K + k * Np, HorDiffERHSY + field * Np*K*Nface + k * Np, Np);
			}
		}
	}
	/*Finally, multiply each component of HorDiffAVx and HorDiffAVy by its diffusion coefficient to get the final auxialary variable,
	this is conducted according to $(Q,v)=(\nv\hat Q,v)$,where $\hat Q = \frac{\partial u(v,\theta)}{\partial (x,y)}$,
	we note that in this projection precedure, aliasing error is introduced.*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				DotProduct(HorDiffAVx + field * Np*K + k * Np, HorDiffAVx + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
				DotProduct(HorDiffAVy + field * Np*K + k * Np, HorDiffAVy + field * Np*K + k * Np, HorDiffnv + k * Np, Np);
			}
			for (int field = 2; field < Nfield; field++) {
				DotDivideByConstant(HorDiffAVx + field * Np*K + k * Np, HorDiffAVx + field * Np*K + k * Np, Prantl, Np);
				DotDivideByConstant(HorDiffAVy + field * Np*K + k * Np, HorDiffAVy + field * Np*K + k * Np, Prantl, Np);
			}
		}

	}

	/*Calculate the contribution to the right hand side due to the auxialary variable HorDiffAVx and HorDiffAVy with IPDG.*/
	/*Calculate the penalty parameter $\tau$ first, this parameter is calculated as $\tau=\frac{(N+1)(N+d)}{d}\frac{n_0}{2}\frac{A}{V}\nv$*/
	/*Inner edge first*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			FetchInnerEdgeFacialValue(HorDiffIEnvfm + face * IENfp, HorDiffIEnvfp + face * IENfp, HorDiffnv, \
				IEFToE + face * 2, IEFToN1 + face * IENfp, IEFToN2 + face * IENfp, Np, IENfp);
			double localRatio = IELAV[face] / LAV[(int)IEFToE[face * 2] - 1];
			double adjacentRatio = IELAV[face] / LAV[(int)IEFToE[face * 2 + 1] - 1];
			for (int p = 0; p < IENfp; p++) {
				HorDiffInnerEdgeTau[face*IENfp + p] = fmax(localRatio*(Order + 1)*(Order + MeshType) / MeshType * Nface / 2 * HorDiffIEnvfm[face*IENfp + p], \
					adjacentRatio*(Order + 1)*(Order + MeshType) / MeshType * Nface / 2 * HorDiffIEnvfp[face*IENfp + p]);
			}
		}
		}
	/*Boundary edge next*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++) {
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			FetchBoundaryEdgeFacialValue(HorDiffBEnvfm + face * BENfp, HorDiffnv, BEFToE + face * 2, BEFToN1 + face * BENfp, Np, BENfp);
			double localRatio = BELAV[face] / LAV[(int)BEFToE[face * 2] - 1];
			for (int p = 0; p < BENfp; p++) {
				HorDiffBoundaryEdgeTau[face*BENfp + p] = localRatio * (Order + 1)*(Order + MeshType) / MeshType * Nface / 2 * HorDiffBEnvfm[face*IENfp + p];
			}
		}
	}

	/*Volume integral of the second order operator*/

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/*Prantl number is not considered here, since we have considered this parameter when calculating the auxialary variable,
		it is not need anymore in the volumn integral part when considering the contribution of second order operator to the RHS*/
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				/*$\bold{r_x}\cdot (Dr*Q_x)+\bold{s_x}\cdot (Ds*Q_x)$*/
				GetVolumnIntegral2d(HorDiffVx + field * Np*K + k * Np, HorDiffTempVx + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffAVx + field * Np*K + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);
				/*$\bold{r_y}\cdot (Dr*Q_y)+\bold{s_y}\cdot (Ds*Q_y)$*/
				GetVolumnIntegral2d(HorDiffVy + field * Np*K + k * Np, HorDiffTempVy + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, HorDiffAVy + field * Np*K + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);
				//if (type == Two){
				//	/*The water depth field is excluded from this part*/
				//	Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		HorDiffVx + field*Np*K + k*Np, Np);
				//	Add(OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		OutputRHS + ((int)varIndex[field + 1] - 1)*Np*K + k*Np, \
				//		HorDiffVy + field*Np*K + k*Np, Np);
				//}
				//else if (type == Three){
				Add(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffVx + field * Np*K + k * Np, Np);
				Add(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffVy + field * Np*K + k * Np, Np);
				//}
			}
		}
	}
	/************************************************************************************************************************************/

	///***********************************************************************************************************************************************/
		/*Reset all the data contained in HorDiffERHSX and HorDiffERHSY to zero, becaused these space contains the data left when calculating the auxialary variable*/
	memset(HorDiffERHSX, 0, Np*K*Nfield*Nface * sizeof(double));
	memset(HorDiffERHSY, 0, Np*K*Nfield*Nface * sizeof(double));
	/*Surface integral of the second order operator, this part is calculated seperately for field with index less than 2 and field with index larger than 2 since prantl number should be considered*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < 2; field++) {
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
				GetIEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEnx, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEnx, HorDiffInnerEdgeTau, 1.0, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial y}$*/
				GetIEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEny, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEny, HorDiffInnerEdgeTau, 1.0, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);

			}
		}
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		int adjacentE = (int)IEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 2; field < Nfield; field++) {
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$*/
				GetIEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEnx, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEnx, HorDiffInnerEdgeTau, Prantl, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVx + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);
				/*Inner edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$*/
				GetIEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTIEfm, HorDiffLPDTIEfp, \
					face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, IENe, IENfp, field, HorDiffIEFluxS, \
					IEny, HorDiffIEfm + field * IENe*IENfp, HorDiffIEfp + field * IENe*IENfp, IEny, HorDiffInnerEdgeTau, Prantl, \
					HorDiffAVIEfm, HorDiffAVIEfp, HorDiffAVy + field * Np*K, HorDiffIEFluxM, \
					HorDiffIEFluxP, IEJs, IEMb);

			}
		}
	}
	/*Boundary edge facial integral part*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BENe; face++) {
		int adjacentE = (int)BEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionDry) {
			continue; //BEFToE is just itself
		}
		else {
			for (int field = 0; field < 2; field++) {
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\nv H\frac{\partial (u,v)}{\partial x}$*/
				GetBEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVx + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\nv H\frac{\partial (u,v)}{\partial y}$*/
				GetBEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, 1.0, HorDiffAVBEfm, HorDiffAVy + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
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
			for (int field = 2; field < Nfield; field++) {
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_x}{\partial x}$, here $Q_x=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial x}$,
				we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToRHS(HorDiffERHSX + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermX + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEnx, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEnx, HorDiffBoundaryEdgeTau, Prantl, HorDiffAVBEfm, HorDiffAVx + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
				/*Boundary edge contribution to right hand side due to term $\frac{\partial Q_y}{\partial y}$, here $Q_y=\frac{\nv}{\sigma} H\frac{\partial \theta}{\partial y}$,
				we note that the numerical flux at the boundary is set directly to the outer value*/
				GetBEContributionToRHS(HorDiffERHSY + field * Np*K*Nface, HorDiffLPDTBEfm, face, HorDiffLocalPrimitiveDiffTermY + field * Np*K, \
					BEFToE, BEFToF, BEFToN1, Np, K, BENe, BENfp, field, HorDiffBEFluxS, BEny, HorDiffBEfm + field * BENe*BENfp, \
					HorDiffBEfp + field * BENe*BENfp, BEny, HorDiffBoundaryEdgeTau, Prantl, HorDiffAVBEfm, HorDiffAVy + field * Np*K, \
					HorDiffBEFluxM, BEJs, BEMb);
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {
				for (int face = 1; face < Nface; face++) {
					Add(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffERHSX + field * Np*K*Nface + face * Np*K + k * Np, Np);
					Add(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffERHSY + field * Np*K*Nface + face * Np*K + k * Np, Np);
				}
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {

				MultiEdgeContributionByLiftOperator(HorDiffERHSX + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralX + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);

				MultiEdgeContributionByLiftOperator(HorDiffERHSY + field * Np*K*Nface + k * Np, HorDiffTempFacialIntegralY + field * Np*K + k * Np, &np, &oneI, &np, \
					&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
			}
		}
	}


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionDry) {
			continue;
		}
		else {
			for (int field = 0; field < Nfield; field++) {

				Minus(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffERHSX + field * Np*K*Nface + k * Np, Np);
				Minus(OutputRHS + field * Np*K + k * Np, \
					OutputRHS + field * Np*K + k * Np, \
					HorDiffERHSX + field * Np*K*Nface + k * Np, Np);

			}
		}
	}

}

void NdgSWEHorizSmagrinskyDiffSolver::UpdateViscosity(double *fphys_, double *nv_,double *HS_, double *WLEN_, double *UBOT_, double *time) {
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	double *WLEN = WLEN_;
	double *HS = HS_;
	double *UBOT = UBOT_;
	double *nv_cw = nv_;
	double *hu = fphys_;
	double *hv = fphys_ + K * Np;
	double *h = fphys_ + K * Np * 3;
	double *zeta = fphys_ + K * Np * 6;
	double *Dr = meshunion->cell_p->Dr;
	double *Ds = meshunion->cell_p->Ds;
	double *rx = meshunion->rx;
	double *sx = meshunion->sx;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;
	double *LAV = meshunion->LAV;
	double *sigma = meshunion->z;
	int *Nlayer = meshunion->Nlayer;
	int oneI = 1;
	const double WLEN_min = 0.05;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int n = 0; n < Np; n++) {
			WLEN[k * Np + n] = fmax(WLEN[k * Np + n], WLEN_min);
			WaveNumber_[k * Np + n] = 2.0 * 3.14159265 / WLEN[k * Np + n];
			Hrms[k * Np + n] = HS[k * Np + n] / sqrt(2);
		}

	}

	memset(u_u, 0, Np * K * sizeof(double));
	memset(v_v, 0, Np * K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < Np; j++) {
			/*  Calculate u,v in a cell   */
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i * Np + j, hu + i * Np + j, u_u + i * Np + j);
			EvaluatePhysicalVariableByDepthThreshold(Hcrit, h + i * Np + j, hv + i * Np + j, v_v + i * Np + j);
		}
	}

	memset(rx_Dr_u, 0, Np * K * sizeof(double));
	memset(rx_Dr_v, 0, Np * K * sizeof(double));
	memset(sx_Ds_u, 0, Np * K * sizeof(double));
	memset(sx_Ds_v, 0, Np * K * sizeof(double));
	memset(ry_Dr_u, 0, Np * K * sizeof(double));
	memset(ry_Dr_v, 0, Np * K * sizeof(double));
	memset(sy_Ds_u, 0, Np * K * sizeof(double));
	memset(sy_Ds_v, 0, Np * K * sizeof(double));


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		/*  Calculate MatrixMultiply of Dr*u INTO ry_Dr_u, Dr*v INTO ry_Dr_v, Ds*u INTO sy_Ds_u, Ds*v INTO sy_Ds_v in a cell   */
		MatrixMultiply(Dr, u_u + i * Np, ry_Dr_u + i * Np, Np, oneI, Np, 1.0);
		MatrixMultiply(Dr, v_v + i * Np, ry_Dr_v + i * Np, Np, oneI, Np, 1.0);
		MatrixMultiply(Ds, u_u + i * Np, sy_Ds_u + i * Np, Np, oneI, Np, 1.0);
		MatrixMultiply(Ds, v_v + i * Np, sy_Ds_v + i * Np, Np, oneI, Np, 1.0);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		/*  Calculate Dotproduct of rx_Dr*u to sy_Ds*v in a cell   */
		DotProduct(rx_Dr_u + i * Np, rx, ry_Dr_u + i * Np, Np);
		DotProduct(rx_Dr_v + i * Np, rx, ry_Dr_v + i * Np, Np);
		DotProduct(sx_Ds_u + i * Np, sx, sy_Ds_u + i * Np, Np);
		DotProduct(sx_Ds_v + i * Np, sx, sy_Ds_v + i * Np, Np);
		DotProduct(ry_Dr_u + i * Np, ry, ry_Dr_u + i * Np, Np);
		DotProduct(ry_Dr_v + i * Np, ry, ry_Dr_v + i * Np, Np);
		DotProduct(sy_Ds_u + i * Np, sy, sy_Ds_u + i * Np, Np);
		DotProduct(sy_Ds_v + i * Np, sy, sy_Ds_v + i * Np, Np);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < Np; j++) {
			if (h[i*Np + j] >= Hcrit) {
				nv_c[i*Np + j] = C * LAV[i] * (*Nlayer) * sqrt(pow((rx_Dr_u[i*Np + j] + sx_Ds_u[i*Np + j]), 2) + 0.5* pow((rx_Dr_v[i*Np + j] + sx_Ds_v[i*Np + j] + ry_Dr_u[i*Np + j] + sy_Ds_u[i*Np + j]), 2) + pow((ry_Dr_v[i*Np + j] + sy_Ds_v[i*Np + j]), 2));

				//*****nv_w = lemda * Hrms * Ubot * cosh[k * (z + h)]; lemda = 0.3; z = sigma * h +zeta; k is wavenumber; h means still water level. 
				//nv_w[i*Np + j] = 0.3*Hrms[i*Np + j] * UBOT[i*Np + j] * cosh(WaveNumber_[i*Np + j] * ((1 + sigma[i*Np + j])*h[i*Np + j]));

			}
			else {
				nv_c[i*Np + j] = 0.0;
				nv_w[i*Np + j] = 0.0;
			}

			//nv_cw[i*Np + j] = sqrt(nv_c[i*Np + j]* nv_c[i*Np + j]+ nv_w[i*Np + j]* nv_w[i*Np + j]);//Why the nv_cw is sometimes unsteady? e.g. In the Inlet case!
			nv_cw[i*Np + j] = sqrt(2) * nv_c[i*Np + j]; 
		};
		/************************************************************************************************************************************/
			/*  Thought about the wave influence, we need Hrms, Ubot to calculate the Am,w. Finally, Am=sqrt(Am,c^2+Am,w^2)    */

		/************************************************************************************************************************************/
	}

	//ofstream fout0("nv_w.dat", ios_base::app);
	//if (!fout0.is_open()) {
	//	cerr << "Can't open file for output, exit here.";
	//	exit(EXIT_FAILURE);
	//}
	//	fout0 << "This is the data for time step " << *time << endl;
	//	for (int k = 0; k < K; k++) {
	//		for (int n = 0; n < Np; n++) {
	//			fout0 << nv_w[k * Np + n] << " ";
	//		}
	//		fout0 << endl;
	//	}
	//fout0.close();
	//ofstream fout1("nv_cw.dat", ios_base::app);
	//if (!fout1.is_open()) {
	//	cerr << "Can't open file for output, exit here.";
	//	exit(EXIT_FAILURE);
	//}
	//	fout1 << "This is the data for time step " << *time << endl;
	//	for (int k = 0; k < K; k++) {
	//		for (int n = 0; n < Np; n++) {
	//			fout1 << nv_cw[k * Np + n] << " ";
	//		}
	//		fout1 << endl;
	//	}
	//fout1.close();

	//ofstream fout2("nv_c.dat", ios_base::app);
	//if (!fout2.is_open()) {
	//	cerr << "Can't open file for output, exit here.";
	//	exit(EXIT_FAILURE);
	//}
	//	fout2 << "This is the data for time step " << *time << endl;
	//	for (int k = 0; k < K; k++) {
	//		for (int n = 0; n < Np; n++) {
	//			fout2 << nv_c[k * Np + n] << " ";
	//		}
	//		fout2 << endl;
	//	}
	//fout2.close();
}
