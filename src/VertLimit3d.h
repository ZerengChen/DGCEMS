#pragma once
#ifndef _Limiter3D
#define _Limiter3D

#include "MeshUnion.h"
#include"GlobalVar.h"
//#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

void Limiter3d(double *fphys);

#endif

