#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class SWETopographySourceTerm3d
{
public:
	SWETopographySourceTerm3d();
	~SWETopographySourceTerm3d();

	void EvaluateTopographySourceTerm(double *fphys, double *frhs,int*pE3d,int MyID);

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

