#pragma once
#include"Bcell.h"

class Mesh2dBoundaryEdge
{
public:
	Mesh2dBoundaryEdge();
	~Mesh2dBoundaryEdge();

	Bcell bcell;
	static Bcell *bcell_p;

	static double* FToE2d;
	static double* FToF2d;
	static double* FToM2d;
	static double* FToN12d;
	static double* FToN22d;
	static double* FToV2d;
	static double* ftype2d;
	static double* Js2d;
	static double* LAV2d;
	static double* M2d;
	static int* Ne2d;
	static int* Nfp2d;
	static double* nx2d;
	static double* ny2d;
	static double* nz2d;
	static double* xb2d;
	static double* yb2d;
	static double* r2d;

};

