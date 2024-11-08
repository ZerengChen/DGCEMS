#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
//#include"SWEPreBlanaced3d.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class NdgQuadFreeStrongFormAdvSolver3d
{
public:
	NdgQuadFreeStrongFormAdvSolver3d();
	~NdgQuadFreeStrongFormAdvSolver3d();

	void evaluateAdvectionRHS(double *fphys, double *frhs, double *fext, int *varFieldIndex,int*pE3d,int MyID);

	//SWEAbstract3d sweabstract3d;
	//SWEPreBlanaced2d swepreblanaced2d;
};

