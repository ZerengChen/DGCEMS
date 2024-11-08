#pragma once
#include "MeshUnion.h"
#include"SWETopographySourceTerm3d.h"
#include"SWECoriolisTerm3d.h"

#ifdef COUPLING_SWAN
#include"SWERollerWaveRadiationTerm3d.h"
#endif

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

#ifndef COUPLING_SWAN
class NdgSourceTermSolver3d
{
public:
	NdgSourceTermSolver3d();
	~NdgSourceTermSolver3d();

	void EvaluateSourceTerm(double *fphys, double *frhs, int*pE3d, int MyID);

	SWETopographySourceTerm3d swetopographysourceterm3d;
	SWECoriolisTerm3d swecoriolisterm3d;
};
#endif


#ifdef COUPLING_SWAN
class NdgSourceTermSolver3d
{
public:
	NdgSourceTermSolver3d();
	~NdgSourceTermSolver3d();

	void EvaluateSourceTerm(double *fphys, double *frhs, double *time_, double *HS, double *T, double *DIR, double *QB, double *WLEN, int*pE2d, int*pE3d, int MyID);

	SWETopographySourceTerm3d swetopographysourceterm3d;
	SWECoriolisTerm3d swecoriolisterm3d;
	SWERollerWaveRadiationTerm3d swerollerwaveradiationterm3d;

	
};
#endif
