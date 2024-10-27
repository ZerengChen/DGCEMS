#pragma once

#include"BoundaryEdge.h"
#include"InnerEdge.h"
#include"BottomEdge.h"
#include"BottomBoundaryEdge.h"
#include"SurfaceBoundaryEdge.h"
#include"Cell.h"
#include"Mesh2d.h"


class MeshUnion
{
public:
	MeshUnion();
	~MeshUnion();
	//void GetMeshAverageValue(double *nodeVal, double *averageValue);

    Cell cell;
	static Cell *cell_p;

	InnerEdge inneredge;
	static InnerEdge *inneredge_p;

	BoundaryEdge boundaryedge;
	static BoundaryEdge *boundaryedge_p;

	BottomEdge bottomedge;
	static BottomEdge *bottomedge_p;

	BottomBoundaryEdge bottomboundaryedge;
	static BottomBoundaryEdge *bottomboundaryedge_p;

	SurfaceBoundaryEdge surfaceboundaryedge;
	static SurfaceBoundaryEdge *surfaceboundaryedge_p;

	Mesh2d mesh2d;
    static Mesh2d *mesh2d_p;

	//notice the order of definition,which will cause the of constructor. 

	static int*K;
	static int*Nv;
	static int *Nlayer;
	static int Nfield;

	//static double*charLength;
	static double*EToE;
	static double*EToF;
	static double*EToM;
    static double*EToL;
	static double*EToV;
	static double*ind;
	static double*J;
	static double*Jz;
	//static double*K;
	static double*LAV;
	//static double*Nv;
	static double*rx;
	static double*ry;
	static double*rz;
	//static signed char *status;
	static double*sx;
	static double*sy;
	static double*sz;
	static double*tx;
	static double*ty;
	static double*type;
	static double*tz;
	static double*vx;
	static double*vy;
	static double*vz;
	static double*x;
	static double*xc;
	static double*y;
	static double*yc;
	static double*z;
	static double*zc;
};


