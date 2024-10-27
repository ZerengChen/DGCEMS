#include "NdgSourceTermSolver3d.h"


NdgSourceTermSolver3d::NdgSourceTermSolver3d()
{
}

NdgSourceTermSolver3d::~NdgSourceTermSolver3d()
{
}

#ifndef COUPLING_SWAN
void NdgSourceTermSolver3d::EvaluateSourceTerm(double *fphys, double *frhs)
{
	swetopographysourceterm3d.EvaluateTopographySourceTerm(fphys, frhs);
	swecoriolisterm3d.EvaluateCoriolisTermRHS(fphys, frhs);
}
#endif


#ifdef COUPLING_SWAN
void NdgSourceTermSolver3d::EvaluateSourceTerm(double *fphys, double *frhs, double *time_,double *HS, double *T, double *DIR, double *QB, double *WLEN)
{
	swetopographysourceterm3d.EvaluateTopographySourceTerm(fphys, frhs);
	swecoriolisterm3d.EvaluateCoriolisTermRHS(fphys, frhs);
	swerollerwaveradiationterm3d.EvaluateWaveRadiationRHS(fphys, frhs, time_, HS, T, DIR, QB, WLEN);
}
#endif