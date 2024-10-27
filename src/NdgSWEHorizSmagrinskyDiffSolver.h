#pragma once

#include "HorizontalDiffusion.h"
#include"GlobalVar.h"
#include "MeshUnion.h"

extern MeshUnion mesh;
extern const MeshUnion *meshunion;

class NdgSWEHorizSmagrinskyDiffSolver //:public NdgHorizDiffSolver
{

public:
	NdgSWEHorizSmagrinskyDiffSolver(double c_);
	~NdgSWEHorizSmagrinskyDiffSolver();

	void EvaluateDiffRHS_Nowave(double *fphys, double *frhs, double *fext, int *varFieldIndex);
	void EvaluateDiffRHS(double *fphys, double *frhs, double *fext, int *varFieldIndex, double *time, double *HS, double *WLEN, double *UBOT);
	void UpdateViscosity_Nowave(double *fphys, double *nv);
	void UpdateViscosity(double *fphys, double *nv, double *HS, double *WLEN, double *UBOT, double *time);
	//void Evaluate_rdhuv(double *r, double *d, double *rdhuv, double *temp);


	double C;

};

