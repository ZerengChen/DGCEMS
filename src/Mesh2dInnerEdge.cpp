#include "Mesh2dInnerEdge.h"



double *Mesh2dInnerEdge::FToE2d = NULL;
double *Mesh2dInnerEdge::FToF2d = NULL;
double *Mesh2dInnerEdge::FToM2d = NULL;
double *Mesh2dInnerEdge::FToN12d = NULL;
double *Mesh2dInnerEdge::FToN22d = NULL;
double *Mesh2dInnerEdge::FToV2d = NULL;
double *Mesh2dInnerEdge::Js2d = NULL;
double *Mesh2dInnerEdge::LAV2d = NULL;
double *Mesh2dInnerEdge::M2d = NULL;
int *Mesh2dInnerEdge::Ne2d = NULL;
int *Mesh2dInnerEdge::Nfp2d = NULL;
double *Mesh2dInnerEdge::nx2d = NULL;
double *Mesh2dInnerEdge::ny2d = NULL;
double *Mesh2dInnerEdge::nz2d = NULL;
double *Mesh2dInnerEdge::r2d = NULL;
Icell *Mesh2dInnerEdge::icell_p = NULL;

Mesh2dInnerEdge::Mesh2dInnerEdge():icell()
{
	icell_p = &icell;

	MeshUnion_dim::ncvar_read(FToE2d, "InnerEdge_FToE2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF2d, "InnerEdge_FToF2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM2d, "InnerEdge_FToM2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToN12d, "InnerEdge_FToN12d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(FToN22d, "InnerEdge_FToN22d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(FToV2d, "InnerEdge_FToV2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Js2d, "InnerEdge_Js2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(LAV2d, "InnerEdge_LAV2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M2d, "InnerEdge_M2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(Ne2d, "InnerEdge_Ne2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp2d, "InnerEdge_Nfp2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx2d, "InnerEdge_nx2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(ny2d, "InnerEdge_ny2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(nz2d, "InnerEdge_nz2d", MeshUnion_dim::Ne_inner2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(r2d, "InnerEdge_r2d", MeshUnion_dim::one, MeshUnion_dim::Nfp2d);
}


Mesh2dInnerEdge::~Mesh2dInnerEdge()
{

}

