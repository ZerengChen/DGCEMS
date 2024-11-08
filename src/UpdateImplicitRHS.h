#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class UpdateImplicitRHS
{
public:
	UpdateImplicitRHS();
	~UpdateImplicitRHS();

	void EvaluateupdateimplicitRHS(double *fphys_, double *nv_v, double *frhs_, double ImplicitA, double *BBE, double *SBE, int *varIndex,int*,int*,int);


};

