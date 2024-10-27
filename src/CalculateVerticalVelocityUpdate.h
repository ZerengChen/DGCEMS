#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class CalculateVerticalVelocity
{
public:
	CalculateVerticalVelocity();
	~CalculateVerticalVelocity();

	void EvaluateVerticalVelocity(double *fphys2d, double *fphys, double *fext2d, double *fext3d);


};

