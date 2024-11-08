#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include "VerticalColumnIntegralField.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class SWERollerWaveRadiationTerm3d
{
public:
	SWERollerWaveRadiationTerm3d();
	~SWERollerWaveRadiationTerm3d();

	void EvaluateWaveRadiationRHS(double *fphys, double *frhs, double *time_, double *Hs, double *T, double *DIR, double *QB, double *WLEN,int*pE2d, int*pE3d, int MyID);
protected:
	VerticalColumnIntegralField verticalColumnIntegralField;
};

