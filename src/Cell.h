#pragma once
#include"MeshUnion_dim.h"
#include <iostream>

class Cell
{
public:
	Cell();
	~Cell();

	

	static signed char* type;
	static int* Nq;
	static int* Nv;

	static double* Dr;
	static double* Ds;
	static double* Dt;
	static double* FToV;
	static signed char* faceType;
	static double* Fmask;
	static double* invM;
	static double* LAV;
	static double* M;
	static int* N;
	static int* Nz;
	static int* Nface;
	static int* Nfp;
	static int* Nph;
	static int* Npz;
	static double* V1d;
	static double* Nfv;
	static int* Np;
	static double* r;
	static double* rq;
	static double* s;
	static double* sq;
	static double* t;
	static int* TNfp;
	static double* tq;
	//static double* type;
	static double* V;
	static double* Vq;
	static double* vr;
	static double* vs;
	static double* vt;
	static double* wq;

	static double* r1;
	static double* s1;
	static double* t1;
	static double* Vh;
	static double* Vint;
	static double* VCV;
#ifdef _BAROCLINIC
	static double* VintU;
#endif
};