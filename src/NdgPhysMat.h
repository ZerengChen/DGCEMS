#pragma once
#include"GlobalVar.h"
#include"AbstractOutputFile.h"
#include"NdgSWEHorizSmagrinskyDiffSolver.h"
#include"NdgQuadFreeStrongFormAdvSolver3d.h"
#include"NdgQuadFreeStrongFormPECSolver2d.h"
#include"NdgSourceTermSolver3d.h"
#include "NdgSWEVertGOTMDiffSolver.h"
#include "VerticalColumnIntegralField.h"
#include "CalculateVerticalVelocityUpdate.h"
#include "NdgMemory.h"
#include "VerticalRepmatFacialValue.h"
#include <fstream>
//#include<chrono>
//#include<vector>

using namespace std;
extern MeshUnion mesh;
extern const MeshUnion *meshunion;
class NdgPhysMat 
{
public:
	NdgPhysMat();
	~NdgPhysMat();

	void matSolver();
	void matEvaluateIMEXRK222();

	void UpdateExternalField(double tloc, double *fphys2d, double *fphys);
	void EvaluateRHS_Nowave(double *fphys, double *frhs, double time, double *fext, int *varFieldIndex, double *fphys2d, double *fext2d, double *frhs2d);
#ifdef COUPLING_SWAN
	void EvaluateRHS(double *fphys, double *frhs, double time, double *fext, int *varFieldIndex, double *fphys2d, double *fext2d, double *frhs2d, \
		double*HS, double*T, double*DIR, double*QB, double*WLEN, double*UBOT, double*TMBOT);
#endif
	void UpdateOutputResult(double &time, double *fphys2d, double *fphys);

	void addTecdata(double *, double*, vector<int>, vector<int>, int);

protected:

	//SWEAbstract3d sweabstract3d;
	//SWEConventional3d sweconventional3d;
	NdgQuadFreeStrongFormAdvSolver3d ndgquadfreestrongformadvsolver3d;
	AbstractOutputFile abstractoutputfile;
	NdgSWEHorizSmagrinskyDiffSolver ndgswehorizsmagrinskydiffsolver;
	NdgQuadFreeStrongFormPECSolver2d ndgquadfreestrongformPECsolver2d;
	NdgSourceTermSolver3d ndgsourcetermsolver3d;
	VerticalColumnIntegralField verticalColumnIntegralField;
	CalculateVerticalVelocity calculateVerticalVelocity;
	NdgSWEVertGOTMDiffSolver ndgswevertgotmdiffsolver;
	NdgMemory AllocateMemory;
	VerticalRepmatFacialValue verticalrepmatfacialvalue;
	string casename;

	double *EXfrhs2d;
	double *IMfrhs;
	double *EXfrhs;
	double *fext;
	double *fext2d;
	//double *fphys0;
	double *fphys;
	double *fphys2d;
	vector<double> tidal;
	vector<int> obeindex;


	vector<int> DG_Swan_Node;
	vector<int> Swan_DG_Node;
	vector<int> sizeof_PerNode;

	//double tidalinterval;

	//double ftime;
	int outputIntervalNum;
	int *Np;
	int *K;
	int *K2d;
	int *Np2d;
	int *Nv;
	int *varFieldIndex;
	int *BENfp;
	int *BENe;
	int *BENfp2d;
	int *BENe2d;
	int *IENfp;
	int *IENe;
	int *IENfp2d;
	int *IENe2d;
	int *BotENfp, *BotENe;
	int *BotBENfp, *BotBENe;
	int *SurfBENfp, *SurfBENe;
	int *Nface, *Nface2d;
	int *Nlayer3d;
	double startTime, finalTime;

#ifdef COUPLING_SWAN
	//7 vars for coupling
	float *HS_from_swan,*T_from_swan,*DIR_from_swan,*QB_from_swan,*WLEN_from_swan,*TMBOT_from_swan,*UBOT_from_swan;
	float *H_to_swan, *U_to_swan, *V_to_swan,*test_to_swan;
	float *HS, *T, *DIR, *QB , *WLEN, *TMBOT, *UBOT;
#endif

	typedef enum {
		NdgEdgeInner = 0,
		NdgEdgeGaussEdge = 1,
		NdgEdgeSlipWall = 2,
		NdgEdgeNonSlipWall = 3,
		NdgEdgeZeroGrad = 4,
		NdgEdgeClamped = 5,
		NdgEdgeClampedDepth = 6,
		NdgEdgeClampedVel = 7,
		NdgEdgeFlather = 8,
		NdgEdgeNonLinearFlather = 9,
		NdgEdgeNonLinearFlatherFlow = 10,
		NdgEdgeNonReflectingFlux = 11,
		NdgEdgeBottomBoundary = 12,
		NdgEdgeUpperSurfaceBoundary = 13,
		Newmann = 14,
		Dirichlet = 15
	} NdgEdgeType;

};

