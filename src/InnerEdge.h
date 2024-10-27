#pragma once
#include"Cell.h"

class InnerEdge
{
public:
	InnerEdge();
	~InnerEdge();

	static double* FToE;
	static double* FToF;
	static double* FToM;
	static double* FToN1;
	static double* FToN2;
	static double* FToV;
	static double* Js;
	static double* Jz;
	static double* LAV;
	static double* M;
	static int* Ne;
	static int* Nfp;
	static double* nx;
	static double* ny;
	static double* nz;
	//static double* r;
	static double* V1d;
	static double* V2d;
	static int* Nlayer;

};

