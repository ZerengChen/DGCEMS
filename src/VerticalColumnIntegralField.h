#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class VerticalColumnIntegralField
{
public:
	VerticalColumnIntegralField();
	~VerticalColumnIntegralField();

	void EvaluateVerticalIntegral(double *fphys2d, double *fphys);


};

