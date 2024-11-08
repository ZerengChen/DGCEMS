#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class SWECoriolisTerm3d
{
public:
	SWECoriolisTerm3d();
	~SWECoriolisTerm3d();

	void EvaluateCoriolisTermRHS(double *fphys, double *frhs, int*pE3d, int MyID);

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
};

