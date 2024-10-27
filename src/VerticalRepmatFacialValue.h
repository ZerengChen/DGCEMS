#pragma once
#include "MeshUnion.h"
#include"GlobalVar.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class VerticalRepmatFacialValue
{
public:
	VerticalRepmatFacialValue();
	~VerticalRepmatFacialValue();

	void EvaluateRepmatFacialValue(double *fext2d, double *fext3d);


};

