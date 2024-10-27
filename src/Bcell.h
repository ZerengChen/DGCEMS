#pragma once
#include"MeshUnion_dim.h"
//#include<algorithm>

class Bcell
{
public:
	Bcell();
	~Bcell();

	static double* Dr2d;
	static double* Ds2d;
	static double* Dt2d;
	static double* faceType2d;
	static double* Fmask2d;
	static double* FToV2d;
	static double* invM2d;
	static double* LAV2d;
	static double* M2d;
	static double* N2d;
	static double* Nface2d;
	static double* Nfp2d;
	static double* Nfv2d;
	static double* Np2d;
	static double* Nq2d;
	static double* Nv2d;
	static double* r2d;
	static double* rq2d;
	static double* s2d;
	static double* sq2d;
	static double* t2d;
	static double* TNfp2d;
	static double* tq2d;
	static double* type2d;
	static double* V2d;
	static double* Vq2d;
	static double* vr2d;
	static double* vs2d;
	static double* vt2d;
	static double* wq2d;
};

