#ifndef _Density_H
#define _Density_H

#include "MeshUnion.h"
#include "GlobalVar.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

void CalculateDensityField(double *fphys);

#endif