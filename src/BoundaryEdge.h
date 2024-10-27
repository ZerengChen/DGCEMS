#pragma once
#include"Cell.h"

class BoundaryEdge
{
public:
	BoundaryEdge();
	~BoundaryEdge();

	static double* FToE;
	static double* FToF;
	static double* FToM;
	static double* FToN1;
	static double* FToN2;
	static double* FToV;
	static double* ftype;
	static double* Js;
	static double* LAV;
	static double* M;
	static int* Ne;
	static int* Nfp;
	static double* nx;
	static double* ny;
	static double* nz;
	static double* Jz;
	static double* xb;
	static double* yb;
	static double* zb;
	static double* V1d;
	static double* V2d;
	static int* Nlayer;

};

