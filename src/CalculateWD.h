#ifndef _WettingAndDrying_H
#define _WettingAndDrying_H

#include "MeshUnion.h"
#include "GlobalVar.h"
#include <omp.h>

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

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

void JudgingNodeAndElementWD(double *source_h2d, double *source_hu2d, double *source_hv2d, double *source_z, \
	double *hu3d, double *hv3d, double *h3d, signed char *dest, int Np2d, int K2d, int Np3d, int NLayer, int k_);

void WDelementRestructing(double *h2d, double *hu2d, double *hv2d, int Np3d, int NLayer, int Np2d, int K2d, int k_, NdgRegionType type);

void UpdateWetDryState(double *fphys_, double *fphys2d_,  int *NLayer_, signed char *status_, int Np3d, int K3d, int Np2d, int K2d);

//****************  Here a 2-dimensional limiter is added to control the steady of field2d. ********************
void Limiter2d(double *);

#endif