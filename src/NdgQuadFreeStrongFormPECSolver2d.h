#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
//#include"SWEPreBlanaced3d.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class NdgQuadFreeStrongFormPECSolver2d
{
public:
	NdgQuadFreeStrongFormPECSolver2d();
	~NdgQuadFreeStrongFormPECSolver2d();

	void evaluatePCERHSUpdated(double *fphys, double *frhs, double *fext, int *varFieldIndex, double *fphys2d, double *fext2d);


};

