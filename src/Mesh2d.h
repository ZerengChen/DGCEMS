#pragma once
#include"MeshUnion_dim.h"
#include"Mesh2dCell.h"
#include"Mesh2dInnerEdge.h"
#include"Mesh2dBoundaryEdge.h"

//extern "C" {
//	void c_GetMeshIntegralValue(double *nodeVal_, double *wq_, double *J_, double *Vq_, int *Np_, int *K_, int *Nq_, double *integralValue_);
//}

class Mesh2d
{
public:
	Mesh2d();
	~Mesh2d();
	//void GetMeshAverageValue(double *nodeVal, double *averageValue);


	Mesh2dInnerEdge mesh2dinneredge;
	static Mesh2dInnerEdge *mesh2dinneredge_p;

	Mesh2dBoundaryEdge mesh2dboundaryedge;
	static Mesh2dBoundaryEdge *mesh2dboundaryedge_p;

	Mesh2dCell mesh2dcell;
	static Mesh2dCell *mesh2dcell_p;

	static int*K2d;
	static int*Nv2d;
	static int Nfield2d;
	static double*charLength2d;
	static double*EToE2d;
	static double*EToF2d;
	static double*EToM2d;
	static double*EToR2d;
	static double*EToV2d;
	static double*ind2d;
	static double*J2d;
	static double*LAV2d;
	static double*rx2d;
	static double*ry2d;
	static double*rz2d;
	static signed char *status;
	static double*sx2d;
	static double*sy2d;
	static double*sz2d;
	static double*tx2d;
	static double*ty2d;
	static double*type2d;
	static double*tz2d;
	static double*vx2d;
	static double*vy2d;
	static double*vz2d;
	static double*x2d;
	static double*xc2d;
	static double*y2d;
	static double*yc2d;
	static double*z2d;
	static double*zc2d;

};