#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include <omp.h>
#include "UpdateImplicitRHS.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;


class NdgSWEVertGOTMDiffSolver
{
public:
	NdgSWEVertGOTMDiffSolver();
	~NdgSWEVertGOTMDiffSolver();

	void EvaluateVertDiffRHS(double *fphys, double *frhs, double *time, double *fphys2d, double ImplicitA, int *varIndex, int*, int*, int);
	void EvaluateVertDiffRHS_CW(double *fphys, double *frhs, double *time, double *fphys2d, double ImplicitA, double *UBOT, double *TMBOT, int *varIndex, int*, int*, int);

protected:
	UpdateImplicitRHS updateimplicitrhs;

};

