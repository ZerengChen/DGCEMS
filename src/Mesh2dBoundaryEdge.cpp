#include "Mesh2dBoundaryEdge.h"
#include<iostream>

double *Mesh2dBoundaryEdge::FToE2d = NULL;
double *Mesh2dBoundaryEdge::FToF2d = NULL;
double *Mesh2dBoundaryEdge::FToM2d = NULL;
double *Mesh2dBoundaryEdge::FToN12d = NULL;
double *Mesh2dBoundaryEdge::FToN22d = NULL;
double *Mesh2dBoundaryEdge::FToV2d = NULL;
double *Mesh2dBoundaryEdge::ftype2d = NULL;
double *Mesh2dBoundaryEdge::Js2d = NULL;
double *Mesh2dBoundaryEdge::LAV2d = NULL;
double *Mesh2dBoundaryEdge::M2d = NULL;
int *Mesh2dBoundaryEdge::Ne2d = NULL;
int *Mesh2dBoundaryEdge::Nfp2d = NULL;
double *Mesh2dBoundaryEdge::nx2d = NULL;
double *Mesh2dBoundaryEdge::ny2d = NULL;
double *Mesh2dBoundaryEdge::nz2d = NULL;
double *Mesh2dBoundaryEdge::xb2d = NULL;
double *Mesh2dBoundaryEdge::yb2d = NULL;
double *Mesh2dBoundaryEdge::r2d = NULL;
Bcell *Mesh2dBoundaryEdge::bcell_p = NULL;

Mesh2dBoundaryEdge::Mesh2dBoundaryEdge() :bcell()
{
	bcell_p = &bcell;

	MeshUnion_dim::ncvar_read(FToE2d, "BoundaryEdge_FToE2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF2d, "BoundaryEdge_FToF2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM2d, "BoundaryEdge_FToM2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToN12d, "BoundaryEdge_FToN12d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(FToN22d, "BoundaryEdge_FToN22d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(FToV2d, "BoundaryEdge_FToV2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(ftype2d, "BoundaryEdge_ftype2d", MeshUnion_dim::one, MeshUnion_dim::Ne_boundary2d);
	MeshUnion_dim::ncvar_read(Js2d, "BoundaryEdge_Js2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(LAV2d, "BoundaryEdge_LAV2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M2d, "BoundaryEdge_M2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(Ne2d, "BoundaryEdge_Ne2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp2d, "BoundaryEdge_Nfp2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx2d, "BoundaryEdge_nx2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(ny2d, "BoundaryEdge_ny2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(nz2d, "BoundaryEdge_nz2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(xb2d, "BoundaryEdge_xb2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(yb2d, "BoundaryEdge_yb2d", MeshUnion_dim::Ne_boundary2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(r2d, "BoundaryEdge_r2d", MeshUnion_dim::one, MeshUnion_dim::Nfp2d);
}


Mesh2dBoundaryEdge::~Mesh2dBoundaryEdge()
{

}


