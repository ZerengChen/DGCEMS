#pragma once
#include"Cell.h"

class SurfaceBoundaryEdge
{
public:
	SurfaceBoundaryEdge();
	~SurfaceBoundaryEdge();

	static double* ftype;
	static double* FToE;
	static double* FToF;
	static double* FToM;
	static double* FToN1;
	static double* FToN2;
	static double* FToV;
	static double* Js;
	static double* LAV;
	static double* M;
	static int* Ne;
	static int* Nfp;
	static double* nx;
	static double* ny;
	static double* nz;

};

