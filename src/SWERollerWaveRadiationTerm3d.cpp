#include "SWERollerWaveRadiationTerm3d.h"
#include <cmath>
#include "NdgMath.h"
#include "NdgSWE3D.h"
//#include "NdgMemory.h"
#include "HorizontalDiffusion.h"
#include <fstream>
#include <cstdlib>

using namespace std;

//#define max(a, b) ((a > b) ? a : b)
//#define min(a, b) ((a < b) ? a : b)


SWERollerWaveRadiationTerm3d::SWERollerWaveRadiationTerm3d()
{
}

SWERollerWaveRadiationTerm3d::~SWERollerWaveRadiationTerm3d()
{
}

extern double *SIN_DIR, *COS_DIR, *WaveNumber, *RSIEfm, *RSIEfp, *RSIEFluxM, *RSIEFluxP, \
*RSIEFluxS, *RSERHS, *RSBotEfm, *RSBotEfp, *RSBotEFluxM, *RSBotEFluxP, *RSBotEFluxS, \
*RSBEfm, *RSBEfp, *RSBEFluxM, *RSBEFluxS, \
*WaveEnergy, *KD, *FSS, *FCS, *FSC, *FCC, *M2, *N1, *N2, *Ar, *Rz, *R, *CFF1, *CFF2, *CFF3, *CFF4,\
*Roller1, *Roller2, *Roller3, *H_Radiation1, *H_Radiation2, *H_Radiation3, *V_Radiation,\
*RS_E,*RS_G,*Hx, *Hy, *tempRHSx, *tempRHSy, *midtempRHSx, *midtempRHSy, *VtempRHSx, *VtempRHSy,\
*midERHS, *RSBotBEfm, *RSBotBEFluxM, *RSBotBEFluxS, *RSSurfBEfm, *RSSurfBEFluxM, *RSSurfBEFluxS, *CFF5,\
*CFF6, *CFF7, *CFF8, *wave_C, *Rzn, *Rz1;

extern double gra;
extern double Hcrit;
extern int Nvar;
extern signed char *Status3d;
//RollerWaveRadiationMemoryDeAllocation();

void SWERollerWaveRadiationTerm3d::EvaluateWaveRadiationRHS(double *fphys_, double *frhs_, double *time_, double *Hs_, double *T_, double *DIR_, double *QB_, double *WLEN_)
{
	int Np = *(meshunion->cell_p->Np);
	int K = *(meshunion->K);
	int K2d = *(meshunion->mesh2d_p->K2d);
	int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	signed char *status = meshunion->mesh2d_p->status;
	int NLayer = *(meshunion->Nlayer);
	int Npz = *(meshunion->cell_p->Npz);

	double *rx = meshunion->rx;
	double *Dr = meshunion->cell_p->Dr;
	double *sx = meshunion->sx;
	double *Ds = meshunion->cell_p->Ds;
	double *ry = meshunion->ry;
	double *sy = meshunion->sy;
	double *tz = meshunion->tz;
	double *Dt = meshunion->cell_p->Dt;
	double *sigma = meshunion->z;
	double *J = meshunion->J;
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

	double *fphys = fphys_;
	double *h = fphys + 3 * Np * K;
	double *OutputRHS = frhs_;
	double *Hs = Hs_;
	double *Tp = T_;
	double *DIR = DIR_;
	double *QB = QB_;
	double *WLEN = WLEN_;
	
	const double Deg2Rad = 0.0174532925;
	const double WLEN_min = 0.5;
	const double KD_max = 5.0;
	const double eps1 = 1.0E-14;
	const double cff_Ar = 23.570226;// ==sqrt(2)/0.06

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

	/*Prepare vital varible: FSS,FCC,FCS,M1=FSS,M2,Ar,Rz,N1,N2*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int n = 0; n < Np; n++) {
			if (Hs[k * Np + n] > 0.005) {
				SIN_DIR[k * Np + n] = sin(DIR[k * Np + n] * Deg2Rad);
				COS_DIR[k * Np + n] = cos(DIR[k * Np + n] * Deg2Rad);
				WLEN[k * Np + n] = fmax(WLEN[k * Np + n], WLEN_min);
				WaveNumber[k * Np + n] = 2.0 * 3.14159265 / WLEN[k * Np + n];
				WaveEnergy[k * Np + n] = 0.0625 * gra * Hs[k * Np + n] * Hs[k * Np + n];
				KD[k * Np + n] = fmin(WaveNumber[k * Np + n] * h[k * Np + n] + eps1, KD_max);
				KD[k * Np + n] = fmax(0, KD[k * Np + n]);
				CFF1[k * Np + n] = KD[k * Np + n] * (1 + sigma[k * Np + n]);//CFF1 = KD(1+sigma) can be zero
				KD[k * Np + n] = fmax(eps1, KD[k * Np + n]);
				wave_C[k * Np + n] = fmax(0.0, sqrt(gra / WaveNumber[k * Np + n] * tanh(KD[k * Np + n])));
				FSS[k * Np + n] = sinh(CFF1[k * Np + n]) / sinh(KD[k * Np + n]);
				FCS[k * Np + n] = cosh(CFF1[k * Np + n]) / sinh(KD[k * Np + n]);
				FSC[k * Np + n] = sinh(CFF1[k * Np + n]) / cosh(KD[k * Np + n]);
				FCC[k * Np + n] = cosh(CFF1[k * Np + n]) / cosh(KD[k * Np + n]);
				//Rz1[k * Np + n] = cosh(2.0*sqrt(2)*3.14159265 / (Hs[k * Np + n] +eps1) * h[k * Np + n] * (1.0 + sigma[k * Np + n]));
				Rz1[k * Np + n] = 2.0 * sigma[k * Np + n] * h[k * Np + n] / Hs[k * Np + n];
				Rz1[k * Np + n] = 1.0 - tanh(pow(Rz[k * Np + n], 4));
				R[k * Np + n] = 1.0 + 2.0 * KD[k * Np + n] / sinh(KD[k * Np + n] * 2.0);
				N1[k * Np + n] = (CFF1[k * Np + n] * CFF1[k * Np + n] / R[k * Np + n] / sinh(2.0 * KD[k * Np + n]) - CFF1[k * Np + n] + 1.0 / tanh(KD[k * Np + n]) / R[k * Np + n] / R[k * Np + n])* FSS[k * Np + n] + (CFF1[k * Np + n] / tanh(KD[k * Np + n]) / R[k * Np + n] / R[k * Np + n] + 2.0*CFF1[k * Np + n] / R[k * Np + n] / sinh(2.0 * KD[k * Np + n]) - 1.0)*FCS[k * Np + n] - sigma[k * Np + n] * FCS[k * Np + n];
				N2[k * Np + n] = (CFF1[k * Np + n] * CFF1[k * Np + n] / R[k * Np + n] / sinh(2.0 * KD[k * Np + n]) - CFF1[k * Np + n])*FCC[k * Np + n] + CFF1[k * Np + n] / tanh(KD[k * Np + n]) / R[k * Np + n] / R[k * Np + n] * FSC[k * Np + n] - N1[k * Np + n];
			}
			else {
				continue;
			}

		}
	}

	verticalColumnIntegralField.EvaluateVerticalIntegral(Rzn, Rz1);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K2d; k++) {
		if ((NdgRegionType)Status2d[k] == NdgRegionWet) {
			for (int i = 0; i < NLayer; i++) {
				for (int j = 0; j < Npz; j++) {
					for (int n = 0; n < Np2d; n++) {
						Rz[k * NLayer * Np + i * Np + j * Np2d + n] = Rz1[k * NLayer * Np + i * Np + j * Np2d + n] / (Rzn[k * Np2d + n] + eps1);
					}
				}
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
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			Minus(M2 + k * Np, FCC + k * Np, FSS + k * Np, Np);
			//Ar = 0.06/sqrt(2) * Hs * L * Qb
			DotProduct(Ar + k * Np, Hs + k * Np, WLEN + k * Np, Np);
			DotProduct(Ar + k * Np, Ar + k * Np, QB + k * Np, Np);
			DotDivideByConstant(Ar + k * Np, Ar + k * Np, cff_Ar, Np);

			DotProduct(CFF2 + k * Np, M2 + k * Np, N1 + k * Np, Np);
			DotProduct(CFF3 + k * Np, FSS + k * Np, N2 + k * Np, Np);
			Minus(CFF2 + k * Np, CFF2 + k * Np, CFF3 + k * Np, Np);//CFF2 = N2M1-M1N2
		}
	}

	int np = Np;
	int oneI = 1;
	double one = 1.0, zero = 0.0;
	memset(RSIEfm, 0, IENfp*IENe * 3 * sizeof(double));
	memset(RSIEfp, 0, IENfp*IENe * 3 * sizeof(double));
	double *Roller1M = RSIEfm, *Roller2M = RSIEfm + IENfp * IENe, *Roller3M = RSIEfm + 2 * IENfp*IENe;
	double *Roller1P = RSIEfp, *Roller2P = RSIEfp + IENfp * IENe, *Roller3P = RSIEfp + 2 * IENfp*IENe;
	memset(RSBEfm, 0, BENfp*BENe * 3 * sizeof(double));
	memset(RSBEfp, 0, BENfp*BENe * 3 * sizeof(double));
	double *Roller1BM = RSBEfm, *Roller2BM = RSBEfm + BENfp * BENe, *Roller3BM = RSBEfm + 2 * BENfp*BENe;
	double *Roller1BP = RSBEfp, *Roller2BP = RSBEfp + BENfp * BENe, *Roller3BP = RSBEfp + 2 * BENfp*BENe;
	memset(RS_E, 0, Np * K * 2 * sizeof(double));
	memset(RS_G, 0, Np * K * 2 * sizeof(double));
	memset(RSERHS, 0, Np*K * 2 * Nface * sizeof(double));
/***********  ELEMENT TEST ***********/
	//ofstream fout("SIN_DIR.dat", ios_base::app);
	//if (!fout.is_open()) {
	//	cerr << "Can't open file for output, exit here.";
	//	exit(EXIT_FAILURE);
	//}

	//fout << "This is the data for time step " << *time_ << endl;
	//for (int k = 0; k < K; k++) {
	//	for (int n = 0; n < Np; n++) {
	//		fout << SIN_DIR[k * Np + Np * K + n] << " ";
	//	}
	//	fout << endl;
	//}
	//fout.close();
/***********  ELEMENT TEST ***********/

	/*-----------------------------------------------------------------------------------------Calculate surface wave roller source term.*/
	/************************  Volume Integral Part 1  ****************************/
	/****************************  (Roller's Volume Integral)  ****************************/

	memset(Roller1, 0, Np*K * sizeof(double));
	memset(Roller2, 0, Np*K * sizeof(double));
	memset(Roller3, 0, Np*K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			for (int n = 0; n < Np; n++) {
				if (Hs[k * Np + n] > 0.005) {
					Roller1[k*Np + n] = h[k * Np + n] * COS_DIR[k * Np + n] * COS_DIR[k * Np + n] / WLEN[k * Np + n] * wave_C[k * Np + n] * wave_C[k * Np + n] * Ar[k * Np + n] * Rz[k * Np + n];
					Roller2[k*Np + n] = h[k * Np + n] * COS_DIR[k * Np + n] * SIN_DIR[k * Np + n] / WLEN[k * Np + n] * wave_C[k * Np + n] * wave_C[k * Np + n] * Ar[k * Np + n] * Rz[k * Np + n];
					Roller3[k*Np + n] = h[k * Np + n] * SIN_DIR[k * Np + n] * SIN_DIR[k * Np + n] / WLEN[k * Np + n] * wave_C[k * Np + n] * wave_C[k * Np + n] * Ar[k * Np + n] * Rz[k * Np + n];
				}
				else {
					Roller1[k*Np + n] = 0.0;
					Roller2[k*Np + n] = 0.0;
					Roller3[k*Np + n] = 0.0;
				}
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/*hRxx = Roller1,hRxy = Roller2,hRyy = Roller3 saved in E and G*/
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			EvaluatePrebalanceVolumeSourceTerm(RS_E + k * Np, RS_G + k * Np, Roller1 + k * Np, Roller2 + k * Np, Roller3 + k * Np, Np, K);
			//把hRxx存进RS_E的前K*Np空间,hRxy存进RS_E的后K*Np空间;hRxy存进RS_G的前K*Np空间,hRyy存进RS_G的后K*Np空间;
		}
	}


	/*************************  Inner,Boundary Edge Part 1  *********************/
	/****************************  (Roller's Inner,Boundary Edge Integral)  ****************************/
	memset(RSIEFluxM, 0, IENfp*IENe*2 * sizeof(double));
	memset(RSIEFluxP, 0, IENfp*IENe*2 * sizeof(double));
	memset(RSIEFluxS, 0, IENfp*IENe*2 * sizeof(double));
	memset(RSBEFluxM, 0, BENfp*BENe * 2 * sizeof(double));

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
			FetchInnerEdgeFacialValue(Roller1M + face * IENfp, Roller1P + face * IENfp, Roller1, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
			FetchInnerEdgeFacialValue(Roller2M + face * IENfp, Roller2P + face * IENfp, Roller2, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
			FetchInnerEdgeFacialValue(Roller3M + face * IENfp, Roller3P + face * IENfp, Roller3, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);

			EvaluateVerticalFaceSurfFluxRandS(RSIEFluxM + face * IENfp, Roller1M + face * IENfp, Roller2M + face * IENfp, Roller3M + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, IENfp, IENe);
			EvaluateVerticalFaceSurfFluxRandS(RSIEFluxP + face * IENfp, Roller1P + face * IENfp, Roller2P + face * IENfp, Roller3P + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, IENfp, IENe);
			EvaluateHorizontalFaceCentralNumFlux(RSIEFluxS + face * IENfp, Roller1M + face * IENfp, Roller2M + face * IENfp, Roller3M + face * IENfp, Roller1P + face * IENfp, Roller2P + face * IENfp, \
				Roller3P + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, IENfp, IENe);
		}
	}

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
			for (int i = 0; i < 2; i++) {//Not Nvar because only huhv have
				StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, RSIEFluxM + i * IENe*IENfp, \
					RSIEFluxP + i * IENe*IENfp, RSIEFluxS + i * IENe*IENfp, IEJs, IEMb, RSERHS + i * Np*K*Nface);
			}
		}
	}

	/*-------------------------------------------------------------------------*******************  End Surface Roller  ****************/

//	std::cout << " Calculate horizontal wave radiation source term. " << std::endl;
	/************************  Volume Integral Part 2  ****************************/
	//把hSxx存进RS_E的前K*Np空间,hSxy存进RS_E的后K*Np空间;hSxy存进RS_G的前K*Np空间,hSyy存进RS_G的后K*Np空间;

	memset(H_Radiation1, 0, Np*K * sizeof(double));
	memset(H_Radiation2, 0, Np*K * sizeof(double));
	memset(H_Radiation3, 0, Np*K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int n = 0; n < Np; n++) {
			if (Hs[k * Np + n] > 0.005) {
				H_Radiation1[k*Np + n] = h[k * Np + n] * WaveNumber[k * Np + n] * WaveEnergy[k * Np + n] * (COS_DIR[k * Np + n] * COS_DIR[k * Np + n] * FCC[k * Np + n] * FCS[k * Np + n] + FCC[k * Np + n] * FCS[k * Np + n] - FSS[k * Np + n] * FCS[k * Np + n]);
				H_Radiation2[k*Np + n] = h[k * Np + n] * WaveNumber[k * Np + n] * WaveEnergy[k * Np + n] * (COS_DIR[k * Np + n] * SIN_DIR[k * Np + n] * FCC[k * Np + n] * FCS[k * Np + n]);
				H_Radiation3[k*Np + n] = h[k * Np + n] * WaveNumber[k * Np + n] * WaveEnergy[k * Np + n] * (SIN_DIR[k * Np + n] * SIN_DIR[k * Np + n] * FCC[k * Np + n] * FCS[k * Np + n] + FCC[k * Np + n] * FCS[k * Np + n] - FSS[k * Np + n] * FCS[k * Np + n]);
			}
			else {
				H_Radiation1[k*Np + n] = 0.0;
				H_Radiation2[k*Np + n] = 0.0;
				H_Radiation3[k*Np + n] = 0.0;
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			EvaluatePrebalanceVolumeSourceTerm(RS_E + k * Np, RS_G + k * Np, H_Radiation1 + k * Np, H_Radiation2 + k * Np, H_Radiation3 + k * Np, Np, K);
		}
	}

	//将水平辐射应力与波面水滚的体积分一同计算，把离散后的结果存入加(实际上是减，负号在水滚求解中)到OutputRHS中
	memset(Hx, 0, Np*K * 2 * sizeof(double));
	memset(Hy, 0, Np*K * 2 * sizeof(double));
	memset(tempRHSx, 0, Np*K * 2 * sizeof(double));
	memset(tempRHSy, 0, Np*K * 2 * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			for (int field = 0; field < 2; field++) {
				/*$\bold{r_x}\cdot (Dr*E)+\bold{s_x}\cdot (Ds*E)$*/
				GetVolumnIntegral2d(Hx + field * Np*K + k * Np, tempRHSx + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, RS_E + field * Np*K + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);
				/*$\bold{r_y}\cdot (Dr*E)+\bold{s_y}\cdot (Ds*G)$*/
				GetVolumnIntegral2d(Hy + field * Np*K + k * Np, tempRHSy + field * Np*K + k * Np, &np, &oneI, &np, &one, \
					Dr, Ds, &np, RS_G + field * Np*K + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);

				Add(OutputRHS + field * Np*K + k * Np, OutputRHS + field * Np*K + k * Np, Hx + field * Np*K + k * Np, Np);
				Add(OutputRHS + field * Np*K + k * Np, OutputRHS + field * Np*K + k * Np, Hy + field * Np*K + k * Np, Np);
			}
		}
	}

	//Here we output the Horizontal Stress to test the function
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//	for (int k = 0; k < K; k++) {
//		for (int n = 0; n < Np; n++) {
//			//fphys[K * Np * 13 + k * Np + n] = Hx[k * Np + n] + Hy[k * Np + n];
//			//fphys[K * Np * 14 + k * Np + n] = Hx[K * Np + k * Np + n] + Hy[K * Np + k * Np + n];
//			fphys[K * Np * 13 + k * Np + n] = H_Radiation1[k * Np + n] / h[k * Np + n];
//			fphys[K * Np * 14 + k * Np + n] = H_Radiation3[k * Np + n] / h[k * Np + n];
//
//		}
//	}

	/*************************  Inner,Boundary Edge Part 2  *********************/
	/*************************  Inner,Boundary Edge of Horizontal Wave Radiation  *********************/
	memset(RSIEFluxM, 0, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxP, 0, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxS, 0, IENfp*IENe * 2 * sizeof(double));
	memset(RSBEFluxM, 0, BENfp*BENe * 2 * sizeof(double));

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
			FetchInnerEdgeFacialValue(Roller1M + face * IENfp, Roller1P + face * IENfp, H_Radiation1, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
			FetchInnerEdgeFacialValue(Roller2M + face * IENfp, Roller2P + face * IENfp, H_Radiation2, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
			FetchInnerEdgeFacialValue(Roller3M + face * IENfp, Roller3P + face * IENfp, H_Radiation3, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);

			EvaluateVerticalFaceSurfFluxRandS(RSIEFluxM + face * IENfp, Roller1M + face * IENfp, Roller2M + face * IENfp, Roller3M + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, IENfp, IENe);
			EvaluateVerticalFaceSurfFluxRandS(RSIEFluxP + face * IENfp, Roller1P + face * IENfp, Roller2P + face * IENfp, Roller3P + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, IENfp, IENe);
			EvaluateHorizontalFaceCentralNumFlux(RSIEFluxS + face * IENfp, Roller1M + face * IENfp, Roller2M + face * IENfp, Roller3M + face * IENfp, Roller1P + face * IENfp, Roller2P + face * IENfp, \
				Roller3P + face * IENfp, IEnx + face * IENfp, IEny + face * IENfp, IENfp, IENe);
		}
	}

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
			for (int i = 0; i < 2; i++) {
				StrongFormInnerEdgeRHS(face, IEFToE, IEFToF, Np, K, IENfp, IEFToN1, IEFToN2, RSIEFluxM + i * IENe*IENfp, \
					RSIEFluxP + i * IENe*IENfp, RSIEFluxS + i * IENe*IENfp, IEJs, IEMb, RSERHS + i * Np*K*Nface);
			}
		}
	}

	/*************************  End Horizental Wave Radiation  *********************/

//	std::cout << " Calculate vertical wave radiation source term. " << std::endl;

	/************************  Volume Integral Part 3  ****************************/
	/*Calculate CFF2, CFF3 and CFF4 in Spx and Spy.*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		for (int n = 0; n < Np; n++) {
			if (Hs[k * Np + n] > 0.005) {
				CFF2[k*Np + n] = CFF2[k * Np + n] * WaveNumber[k * Np + n] * WaveEnergy[k * Np + n];//CFF2 = kE(N2M1-M1N2)
				CFF3[k*Np + n] = -sqrt(2.0) / 8.0*gra*Hs[k * Np + n] * M2[k * Np + n];//CFF3 = -0.5gaM2
				CFF4[k*Np + n] = sqrt(2.0) / 4.0 * Hs[k * Np + n] * FSS[k * Np + n];//CFF4 = a*M1
				CFF5[k*Np + n] = CFF2[k * Np + n] * COS_DIR[k * Np + n] * SIN_DIR[k * Np + n];//CFF5 = kE(N2M1-M1N2)sincos
				CFF6[k*Np + n] = CFF2[k * Np + n] * SIN_DIR[k * Np + n] * SIN_DIR[k * Np + n];//CFF6 = kE(N2M1-M1N2)sin2
				CFF2[k*Np + n] = CFF2[k * Np + n] * COS_DIR[k * Np + n] * COS_DIR[k * Np + n];//CFF2 = kE(N2M1-M1N2)cos2

				CFF7[k*Np + n] = 2 / wave_C[k * Np + n] * cosh(2 * CFF1[k * Np + n]) / sinh(2 * KD[k * Np + n]);
				CFF8[k*Np + n] = 0.0424 * Hs[k * Np + n] * QB[k * Np + n] * wave_C[k * Np + n] * wave_C[k * Np + n];//Roller energy * 2 = Ar * c*c/L
			}
			else {
				continue;
			}
		}
	}

/*--------------------------------------Update the Stokes drift velocity U=fphys10,V=fphys11.------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] == NdgRegionWet) {
			for (int n = 0; n < Np; n++) {
				//With Roller = 2kx/c *cosh[2kh(1+sigma)]/sinh(2kh)*(E+2Er)
				fphys[K * Np * 11 + k * Np + n] = CFF7[k * Np + n] * (WaveEnergy[k * Np + n] + CFF8[k*Np + n]) * WaveNumber[k * Np + n] * COS_DIR[k * Np + n] * h[k * Np + n];
				fphys[K * Np * 12 + k * Np + n] = CFF7[k * Np + n] * (WaveEnergy[k * Np + n] + CFF8[k*Np + n]) * WaveNumber[k * Np + n] * SIN_DIR[k * Np + n] * h[k * Np + n];

				//Without Roller
				//fphys[K * Np * 11 + k * Np + n] = CFF7[k * Np + n] * WaveEnergy[k * Np + n] * WaveNumber[k * Np + n] * COS_DIR[k * Np + n] * h[k * Np + n];
				//fphys[K * Np * 12 + k * Np + n] = CFF7[k * Np + n] * WaveEnergy[k * Np + n] * WaveNumber[k * Np + n] * SIN_DIR[k * Np + n] * h[k * Np + n];
			}
		}
		else {
			for (int n = 0; n < Np; n++) {
				fphys[K * Np * 11 + k * Np + n] = 0.0;
				fphys[K * Np * 12 + k * Np + n] = 0.0;
			}
		}
	}
/*-------------------------------------------------------------------------------------------------------------------------------------*/

	memset(VtempRHSx, 0, Np*K * sizeof(double));
	memset(VtempRHSy, 0, Np*K * sizeof(double));
	memset(midtempRHSx, 0, Np*K * sizeof(double));
	memset(midtempRHSy, 0, Np*K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			/*******  Save the first term of Spx and Spy  *******/
			/*$\bold{r_x}\cdot (Dr*CFF4)+\bold{s_x}\cdot (Ds*CFF4)$*/
			GetVolumnIntegral2d(VtempRHSx + k * Np, midtempRHSx + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, CFF4 + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);
			/*$\bold{r_y}\cdot (Dr*CFF4)+\bold{s_y}\cdot (Ds*CFF4)$*/
			GetVolumnIntegral2d(VtempRHSy + k * Np, midtempRHSy + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, CFF4 + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);
			DotProduct(VtempRHSx + k * Np, VtempRHSx + k * Np, CFF3 + k * Np, Np);
			DotProduct(VtempRHSy + k * Np, VtempRHSy + k * Np, CFF3 + k * Np, Np);
			/*******  END Save the first term of Spx and Spy  *******/
		}
	}

	memset(tempRHSx, 0, Np*K * 2 * sizeof(double));
	memset(tempRHSy, 0, Np*K * 2 * sizeof(double));
	memset(midtempRHSx, 0, Np*K * sizeof(double));
	memset(midtempRHSy, 0, Np*K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/**********  Save the 2-3 term of Spx and Spy  **********/
		/*$\bold{r_x}\cdot (Dr*h)+\bold{s_x}\cdot (Ds*h)$*/
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			GetVolumnIntegral2d(tempRHSx + k * Np, midtempRHSx + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, h + k * Np, &np, &zero, &np, rx + k * Np, sx + k * Np);
			/*$\bold{r_y}\cdot (Dr*h)+\bold{s_y}\cdot (Ds*h)$*/
			GetVolumnIntegral2d(tempRHSx + K * Np + k * Np, midtempRHSy + k * Np, &np, &oneI, &np, &one, \
				Dr, Ds, &np, h + k * Np, &np, &zero, &np, ry + k * Np, sy + k * Np);

			DotProduct(tempRHSx + k * Np, tempRHSx + k * Np, CFF2 + k * Np, Np);//use tempRHSx again to save space. CFF2*pD/px

			DotProduct(tempRHSx + K * Np + k * Np, tempRHSx + K * Np + k * Np, CFF5 + k * Np, Np);//CFF5*pD/py

			DotProduct(tempRHSy + k * Np, tempRHSx + k * Np, CFF5 + k * Np, Np);//use tempRHSy again to save space. CFF5*pD/px

			DotProduct(tempRHSy + K * Np + k * Np, tempRHSx + K * Np + k * Np, CFF6 + k * Np, Np);//CFF6*pD/py
		}
		/**********  END Save the 2-3 term of Spx and Spy  **********/
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		/*******  Save Spx and Spy in VtempRHSx and VtempRHSy *******/
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			Add(VtempRHSx + k * Np, VtempRHSx + k * Np, tempRHSx + k * Np, Np);
			Add(VtempRHSx + k * Np, VtempRHSx + k * Np, tempRHSx + K * Np + k * Np, Np);

			Add(VtempRHSy + k * Np, VtempRHSy + k * Np, tempRHSy + k * Np, Np);
			Add(VtempRHSy + k * Np, VtempRHSy + k * Np, tempRHSy + K * Np + k * Np, Np);
			/*******  END Save Spx and Spy in VtempRHSx and VtempRHSy *******/
		}
	}
/************************  End Volume Integral Part of inner stage  ****************************/

	//std::cout << " Calculate 2-stage inneredge vertical wave radiation. " << std::endl;
//Then, calculate the V_inneredge of inner stage.
	/*************************  Inner Edge Part 3  *********************/
	memset(RS_E, 0, Np * K * 2 * sizeof(double));
	memset(RS_G, 0, Np * K * 4 * sizeof(double));
	memset(Hx, 0, Np * K * 2 * sizeof(double));
	memset(midERHS, 0, Np*K * 2 * Nface * sizeof(double));

	memset(RSIEfm, 0, IENfp*IENe * 3 * sizeof(double));
	memset(RSIEfp, 0, IENfp*IENe * 3 * sizeof(double));
	memset(RSIEFluxM, 0, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxP, 0, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxS, 0, IENfp*IENe * 2 * sizeof(double));

/*Inner edge facial integral part——twice stage term*/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int face = 0; face < IENe; face++) {
		//calculate the fm,fp of CFF4 
		//再用一次空间Roller1M and Roller1P 存放CFF4 的左右值
		int adjacentE = (int)IEFToE[2 * face];
		int adjacentE2 = (int)IEFToE[2 * face + 1];
		if ((NdgRegionType)Status3d[adjacentE - 1] != NdgRegionWet || (NdgRegionType)Status3d[adjacentE2 - 1] != NdgRegionWet) {
			continue; 
		}
		else {
			FetchInnerEdgeFacialValue(Roller1M + face * IENfp, Roller1P + face * IENfp, CFF4, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
			//calculate the FluxM,FluxP,FluxS and the middle varible ERHS of CFF4 saved in midERHS(2*Np*K*Nface)
			GetIEContributionToAuxialaryVariable(midERHS, face, IENe, IENfp, 0, Roller1M, Roller1P, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, RSIEFluxM, RSIEFluxP, RSIEFluxS, IEnx, IEMb, IEJs);
			GetIEContributionToAuxialaryVariable(midERHS + Np * K*Nface, face, IENe, IENfp, 0, Roller1M, Roller1P, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, RSIEFluxM, RSIEFluxP, RSIEFluxS, IEny, IEMb, IEJs);
		}
    }

	//把边界信息都加到第一个Np*K中,前面一个Np*K存x的,后面一个Np*K存y的
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			for (int face = 1; face < Nface; face++) {
				Add(midERHS + k * Np, midERHS + k * Np, midERHS + face * Np*K + k * Np, Np);
				Add(midERHS + Np * K*Nface + k * Np, midERHS + Np * K*Nface + k * Np, midERHS + Np * K*Nface + face * Np*K + k * Np, Np);
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		//左乘M-1,再除Jacobi系数,finish the discretization of twice stage term(CFF4)
		//use E again to save the results
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			MultiEdgeContributionByLiftOperator(midERHS + k * Np, RS_E + k * Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k * Np, Np);

			MultiEdgeContributionByLiftOperator(midERHS + Np * K*Nface + k * Np, RS_E + Np * K + k * Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
		}
	}

	//get CFF3*pian(CFF4)/pian(x);CFF3*pian(CFF4)/pian(y) into E(Np*K*2)
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			DotProduct(RS_E + k * Np, midERHS + k * Np, CFF3 + k * Np, Np);
			DotProduct(RS_E + Np * K + k * Np, midERHS + Np * K*Nface + k * Np, CFF3 + k * Np, Np);
		}
	}

	memset(midERHS, 0, Np*K * 2 * Nface * sizeof(double));
	memset(RSIEFluxM, 0, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxP, 0, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxS, 0, IENfp*IENe * 2 * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < IENe; face++) {
		//calculate the fm,fp of h 
		int adjacentE = (int)IEFToE[2 * face];
		int adjacentE2 = (int)IEFToE[2 * face + 1];
		if ((NdgRegionType)Status3d[adjacentE - 1] != NdgRegionWet || (NdgRegionType)Status3d[adjacentE2 - 1] != NdgRegionWet) {
			continue; 
		}
		else {
			FetchInnerEdgeFacialValue(Roller2M + face * IENfp, Roller2P + face * IENfp, h, IEFToE + 2 * face, \
				IEFToN1 + IENfp * face, IEFToN2 + IENfp * face, Np, IENfp);
			//calculate the FluxM,FluxP,FluxS and the middle varible ERHS of h saved in midERHS(2*Np*K*Nface)
			GetIEContributionToAuxialaryVariable(midERHS, face, IENe, IENfp, 0, Roller2M, Roller2P, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, RSIEFluxM, RSIEFluxP, RSIEFluxS, IEnx, IEMb, IEJs);
			GetIEContributionToAuxialaryVariable(midERHS + Np * K*Nface, face, IENe, IENfp, 0, Roller2M, Roller2P, IEFToE, IEFToF, IEFToN1, IEFToN2, Np, K, RSIEFluxM, RSIEFluxP, RSIEFluxS, IEny, IEMb, IEJs);
		}
	}

	//把边界信息都加到第一个Np*K中,前面一个Np*K存x的,后面一个Np*K存y的
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			for (int face = 1; face < Nface; face++) {
				Add(midERHS + k * Np, midERHS + k * Np, midERHS + face * Np*K + k * Np, Np);
				Add(midERHS + Np * K*Nface + k * Np, midERHS + Np * K*Nface + k * Np, midERHS + Np * K*Nface + face * Np*K + k * Np, Np);
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		//左乘M-1,再除Jacobi系数,finish the discretization of twice stage term(H)
		//use G again to save the results
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			MultiEdgeContributionByLiftOperator(midERHS + k * Np, RS_G + k * Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k * Np, Np);

			MultiEdgeContributionByLiftOperator(midERHS + Np * K*Nface + k * Np, RS_G + Np * K + k * Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
		}
	}

	//get CFF2*pian(H)/pian(x);CFF5*pian(H)/pian(y) into G(Np*K*2);CFF5*pian(H)/pian(x) into G(Np*K*3);CFF6*pian(H)/pian(y) into G(Np*K*4);
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			DotProduct(RS_G + k * Np, midERHS + k * Np, CFF2 + k * Np, Np);
			DotProduct(RS_G + Np * K + k * Np, midERHS + Np * K*Nface + k * Np, CFF5 + k * Np, Np);
			DotProduct(RS_G + Np * K * 2 + k * Np, midERHS + k * Np, CFF5 + k * Np, Np);
			DotProduct(RS_G + Np * K * 3 + k * Np, midERHS + Np * K*Nface + k * Np, CFF6 + k * Np, Np);
		}
	}

	//合并E and G into E
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			Add(RS_E + k * Np, RS_E + k * Np, RS_G + k * Np, Np);
			Add(RS_E + k * Np, RS_E + k * Np, RS_G + Np * K + k * Np, Np);//facial integral of Spx

			Add(RS_E + Np * K + k * Np, RS_E + Np * K + k * Np, RS_G + Np * K * 2 + k * Np, Np);
			Add(RS_E + Np * K + k * Np, RS_E + Np * K + k * Np, RS_G + Np * K * 3 + k * Np, Np);//facial integral of Spy
		}
	}
		/************************* End  Inner Edge Part 3  *********************/

	//将二阶面积分E加入体积分VtempRHSx与VtempRHSy中,注意用minus而不是add,得到二阶项Spx和Spy的整体离散结果
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			Minus(VtempRHSx + k * Np, VtempRHSx + k * Np, RS_E + k * Np, Np);
			Minus(VtempRHSy + k * Np, VtempRHSy + k * Np, RS_E + K * Np + k * Np, Np);
		}
	}
	//End the calculation of the vertical inner stage.

	//std::cout << " Calculate Volume Integral Part of Outer stage. " << std::endl;
	/************************  Volume Integral Part 4  ****************************/
	/*Volume Integral Part of Outer stage*/

	memset(V_Radiation, 0, 2 * Np*K * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			GetVolumnIntegralOnlyVertical(V_Radiation + k * Np, &np, &oneI, &np, &one, Dt, &np, VtempRHSx + k * Np, &np, &zero, &np, tz);
			GetVolumnIntegralOnlyVertical(V_Radiation + K * Np + k * Np, &np, &oneI, &np, &one, Dt, &np, VtempRHSy + k * Np, &np, &zero, &np, tz);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			Minus(OutputRHS + k * Np, OutputRHS + k * Np, V_Radiation + k * Np, Np);
			Minus(OutputRHS + Np * K + k * Np, OutputRHS + Np * K + k * Np, V_Radiation + K * Np + k * Np, Np);
		}
	}
	/* End Volume Integral Part of Outer stage*/
	/************************ End Volume Integral Part 4  ****************************/


	/************************* Inner Edge Part 4  *********************/
	/*Bottom Inner edge facial integral part——outer one stage term*/

	/*Allocate memory for RSBotEFluxM, RSBotEFluxP and RSBotEFluxS, and calculate these flux term*/
	memset(RSBotEFluxM, 0, BotENfp*BotENe*2 * sizeof(double));
	memset(RSBotEFluxP, 0, BotENfp*BotENe*2 * sizeof(double));
	memset(RSBotEFluxS, 0, BotENfp*BotENe*2 * sizeof(double));
	memset(RSBotEfm, 0, BotENfp*BotENe * 2 * sizeof(double));
	memset(RSBotEfp, 0, BotENfp*BotENe * 2 * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotENe; face++) {
		/*Fetch variable RSBotEfm and RSBotEfp first*/
		//Spx and Spy
		int adjacentE = (int)BotEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			FetchInnerEdgeFacialValue(RSBotEfm + face * BotENfp, RSBotEfp + face * BotENfp, VtempRHSx, BotEFToE + 2 * face, \
				BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);
			FetchInnerEdgeFacialValue(RSBotEfm + BotENfp * BotENe + face * BotENfp, RSBotEfp + BotENfp * BotENe + face * BotENfp, VtempRHSy, BotEFToE + 2 * face, \
				BotEFToN1 + BotENfp * face, BotEFToN2 + BotENfp * face, Np, BotENfp);

			DotProduct(RSBotEFluxM + face * BotENfp, RSBotEfm + face * BotENfp, BotEnz + face * BotENfp, BotENfp);
			DotProduct(RSBotEFluxP + face * BotENfp, RSBotEfp + face * BotENfp, BotEnz + face * BotENfp, BotENfp);
			DotProduct(RSBotEFluxM + BotENfp * BotENe + face * BotENfp, RSBotEfm + BotENfp * BotENe + face * BotENfp, BotEnz + face * BotENfp, BotENfp);
			DotProduct(RSBotEFluxP + BotENfp * BotENe + face * BotENfp, RSBotEfp + BotENfp * BotENe + face * BotENfp, BotEnz + face * BotENfp, BotENfp);

			EvaluateVerticalFaceCentralNumFlux(RSBotEFluxS + face * BotENfp, RSBotEfm + face * BotENfp, RSBotEfp + face * BotENfp, \
				BotEnz + face * BotENfp, BotENfp, BotENe);
			EvaluateVerticalFaceCentralNumFlux(RSBotEFluxS + BotENfp * BotENe + face * BotENfp, RSBotEfm + BotENfp * BotENe + face * BotENfp, RSBotEfp + BotENfp * BotENe + face * BotENfp, \
				BotEnz + face * BotENfp, BotENfp, BotENe);
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
			StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, RSBotEFluxM, \
				RSBotEFluxP, RSBotEFluxS, BotEJs, BotEMb, RSERHS);

			StrongFormInnerEdgeRHS(face, BotEFToE, BotEFToF, Np, K, BotENfp, BotEFToN1, BotEFToN2, RSBotEFluxM + BotENe * BotENfp, \
				RSBotEFluxP + BotENe * BotENfp, RSBotEFluxS + BotENe * BotENfp, BotEJs, BotEMb, RSERHS + Np * K*Nface);
		}
		else{
			continue;
		}
	}

	/*Bottom Boundary edge facial integral part——outer one stage term*/
	//fp=-fm
	/*Allocate memory for RSBotBEFluxM, fluxP and RSBotBEFluxS, and calculate these flux term*/
	memset(RSBotBEFluxM, 0, BotBENfp*BotBENe*2 * sizeof(double));
	memset(RSBotBEFluxS, 0, BotBENfp*BotBENe*2 * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++) {
		int adjacentE = (int)BotBEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			FetchBoundaryEdgeFacialValue(RSBotBEfm + face * BotBENfp, VtempRHSx, BotBEFToE + 2 * face, BotBEFToN1 + face * BotBENfp, Np, BotBENfp);
			FetchBoundaryEdgeFacialValue(RSBotBEfm + BotBENfp * BotBENe + face * BotBENfp, VtempRHSy, BotBEFToE + 2 * face, BotBEFToN1 + face * BotBENfp, Np, BotBENfp);

			EvaluateHorizontalFaceSurfFluxRandS(RSBotBEFluxM + face * BotBENfp, RSBotBEfm + face * BotBENfp, RSBotBEfm + BotBENfp * BotBENe + face * BotBENfp, \
				BotBEnz + face * BotBENfp, BotBENfp, BotBENe);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < BotBENe; face++) {
		int adjacentE = (int)BotBEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			for (int field = 0; field < 2; field++) {
				StrongFormBoundaryEdgeRHS(face, BotBEFToE, BotBEFToF, Np, K, BotBENfp, BotBEFToN1, RSBotBEFluxM + field * BotBENe*BotBENfp, \
					RSBotBEFluxS + field * BotBENe*BotBENfp, BotBEJs, BotBEMb, RSERHS + field * Np*K*Nface);
			}
		}
		else {
			continue;
		}
	}

	/*Surface edge facial integral part——outer one stage term*/

	/*Allocate memory for RSSurfBEFluxM, fluxP and RSSurfBEFluxS, and calculate these flux term*/
	memset(RSSurfBEFluxM, 0, SurfBENfp*SurfBENe*2 * sizeof(double));
	memset(RSSurfBEFluxS, 0, SurfBENfp*SurfBENe*2 * sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++) {
		int adjacentE = (int)SurfBEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			FetchBoundaryEdgeFacialValue(RSSurfBEfm + face * SurfBENfp, VtempRHSx, SurfBEFToE + 2 * face, SurfBEFToN1 + face * SurfBENfp, Np, SurfBENfp);
			FetchBoundaryEdgeFacialValue(RSSurfBEfm + SurfBENfp * SurfBENe + face * SurfBENfp, VtempRHSy, SurfBEFToE + 2 * face, SurfBEFToN1 + face * SurfBENfp, Np, SurfBENfp);

			EvaluateHorizontalFaceSurfFluxRandS(RSSurfBEFluxM + face * SurfBENfp, RSSurfBEfm + face * SurfBENfp, RSSurfBEfm + SurfBENfp * SurfBENe + face * SurfBENfp, \
				SurfBEnz + face * SurfBENfp, SurfBENfp, SurfBENe);
		}
		else {
			continue;
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int face = 0; face < SurfBENe; face++) {
		int adjacentE = (int)SurfBEFToE[2 * face];
		if ((NdgRegionType)Status3d[adjacentE - 1] == NdgRegionWet) {
			for (int field = 0; field < 2; field++) {
				StrongFormBoundaryEdgeRHS(face, SurfBEFToE, SurfBEFToF, Np, K, SurfBENfp, SurfBEFToN1, RSSurfBEFluxM + field * SurfBENe*SurfBENfp, \
					RSSurfBEFluxS + field * SurfBENe*SurfBENfp, SurfBEJs, SurfBEMb, RSERHS + field * Np*K*Nface);
			}
		}
		else {
			continue;
		}
	}

	/************************* End Inner Edge Part 4  *********************/

	//把之前所有的123边界信息都加到第一个Np*K与第Nface*Np*K中,前面一个Np*K存x的,后面一个Np*K存y的
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			for (int face = 1; face < Nface; face++) {
				Add(RSERHS + k * Np, RSERHS + k * Np, RSERHS + face * Np*K + k * Np, Np);
				Add(RSERHS + Np * K*Nface + k * Np, RSERHS + Np * K*Nface + k * Np, RSERHS + Np * K*Nface + face * Np*K + k * Np, Np);
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		//左乘M-1,再除Jacobi系数,finish the discretization of one stage term(Spx,Spy)
		//use Hx again to save the results
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			MultiEdgeContributionByLiftOperator(RSERHS + k * Np, Hx + k * Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k * Np, Np);

			MultiEdgeContributionByLiftOperator(RSERHS + Np * K*Nface + k * Np, Hx + Np * K + k * Np, &np, &oneI, &np, \
				&one, invM, &np, &np, &zero, &np, J + k * Np, Np);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K; k++) {
		if ((NdgRegionType)Status3d[k] != NdgRegionWet) {
			continue;
		}
		else {
			Add(OutputRHS + k * Np, RSERHS + k * Np, OutputRHS + k * Np, Np);
			Add(OutputRHS + Np * K + k * Np, RSERHS + Np * K* Nface + k * Np, OutputRHS + Np * K + k * Np, Np);
		}
	}

	//Here we output the Vertical Stress to test the function
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//	for (int k = 0; k < K; k++) {
//		for (int n = 0; n < Np; n++) {
//			fphys[K * Np * 13 + k * Np + n] = - V_Radiation[k * Np + n] + RSERHS[k * Np + n];
//			fphys[K * Np * 14 + k * Np + n] = - V_Radiation[Np * K + k * Np + n] + RSERHS[Np * K* Nface + k * Np + n];
//			fphys[K * Np * 13 + k * Np + n] = fphys[K * Np * 13 + k * Np + n] + RSERHS[k * Np + n];
//			fphys[K * Np * 14 + k * Np + n] = fphys[K * Np * 14 + k * Np + n] + RSERHS[Np * K* Nface + k * Np + n];
//		}		
//	}

}