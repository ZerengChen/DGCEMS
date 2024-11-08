#include "NdgSourceTermSolver3d.h"


NdgSourceTermSolver3d::NdgSourceTermSolver3d()
{
}

NdgSourceTermSolver3d::~NdgSourceTermSolver3d()
{
}

#ifndef COUPLING_SWAN
void NdgSourceTermSolver3d::EvaluateSourceTerm(double *fphys, double *frhs, int*pE3d, int MyID)
{
	swetopographysourceterm3d.EvaluateTopographySourceTerm(fphys, frhs, pE3d, MyID);
	swecoriolisterm3d.EvaluateCoriolisTermRHS(fphys, frhs, pE3d, MyID);
}
#endif


#ifdef COUPLING_SWAN
void NdgSourceTermSolver3d::EvaluateSourceTerm(double *fphys, double *frhs, double *time_,double *HS, double *T, double *DIR, double *QB, double *WLEN, int*pE2d, int*pE3d, int MyID)
{
	swetopographysourceterm3d.EvaluateTopographySourceTerm(fphys, frhs, pE3d, MyID);
	swecoriolisterm3d.EvaluateCoriolisTermRHS(fphys, frhs, pE3d, MyID);
	swerollerwaveradiationterm3d.EvaluateWaveRadiationRHS(fphys, frhs, time_, HS, T, DIR, QB, WLEN, pE2d, pE3d, MyID);
}
#endif