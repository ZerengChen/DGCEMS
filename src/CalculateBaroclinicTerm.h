#ifndef _Baroclinic_H
#define _Baroclinic_H

#include "MeshUnion.h"
#include "GlobalVar.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

typedef enum {
	NdgEdgeInner = 0,
	NdgEdgeGaussEdge = 1,
	NdgEdgeSlipWall = 2,
	NdgEdgeNonSlipWall = 3,
	NdgEdgeZeroGrad = 4,
	NdgEdgeClamped = 5,
	NdgEdgeClampedDepth = 6,
	NdgEdgeClampedVel = 7,
	NdgEdgeFlather = 8,
	NdgEdgeNonLinearFlather = 9,
	NdgEdgeNonLinearFlatherFlow = 10,
	NdgEdgeNonReflectingFlux = 11,
	NdgEdgeBottomBoundary = 12,
	NdgEdgeUpperSurfaceBoundary = 13,
	Newmann = 14,
	Dirichlet = 15
} NdgEdgeType;

void GetFirstOrderPartialDerivativeInVerticalDirection(double *BaroclinicPRHOPS, double *rho, int*pE3d, int MyID);

void GetFirstOrderPartialDerivativeInHorizontalDirection(double *, double *, double *, double *, double *, double *, double *, int*pE3d, int MyID);

void ImposeBcsForRhoAndH(double *, double *, double *, double *, double *, double *, \
	int, int, NdgEdgeType type);

void EvaluateFaceSurfFlux(double *, double *, double *, int);

void EvaluateFaceNumFlux_Central(double *, double *, double *, double *, int);

void EvaluateBaroclinicTerm(double *fphys_, double *frhs_, double *fext_, int*pE2d, int*pE3d, int MyID);

#endif