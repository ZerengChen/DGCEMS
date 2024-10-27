#include "BoundaryEdge.h"
using namespace std;

double *BoundaryEdge::FToE = NULL;
double *BoundaryEdge::FToF = NULL;
double *BoundaryEdge::FToM = NULL;
double *BoundaryEdge::FToN1 = NULL;
double *BoundaryEdge::FToN2 = NULL;
double *BoundaryEdge::FToV = NULL;
double *BoundaryEdge::ftype = NULL;
double *BoundaryEdge::Js = NULL;
double *BoundaryEdge::LAV = NULL;
double *BoundaryEdge::M = NULL;
int *BoundaryEdge::Ne = NULL;
int *BoundaryEdge::Nfp = NULL;
double *BoundaryEdge::nx = NULL;
double *BoundaryEdge::ny = NULL;
double *BoundaryEdge::nz = NULL;
double *BoundaryEdge::Jz = NULL;
double *BoundaryEdge::xb = NULL;
double *BoundaryEdge::yb = NULL;
double *BoundaryEdge::zb = NULL;
double *BoundaryEdge::V1d = NULL;
double *BoundaryEdge::V2d = NULL;
int *BoundaryEdge::Nlayer = NULL;

BoundaryEdge::BoundaryEdge() //:bcell()
{
	//cout << "Start BoundaryEdge" << endl;
	MeshUnion_dim::ncvar_read(FToE, "BoundaryEdge_FToE", MeshUnion_dim::Ne_boundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF, "BoundaryEdge_FToF", MeshUnion_dim::Ne_boundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM, "BoundaryEdge_FToM", MeshUnion_dim::Ne_boundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToN1, "BoundaryEdge_FToN1", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToN2, "BoundaryEdge_FToN2", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToV, "BoundaryEdge_FToV", MeshUnion_dim::Ne_boundary, MeshUnion_dim::four);
	MeshUnion_dim::ncvar_read(ftype, "BoundaryEdge_ftype", MeshUnion_dim::one, MeshUnion_dim::Ne_boundary);
	MeshUnion_dim::ncvar_read(Js, "BoundaryEdge_Js", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(LAV, "BoundaryEdge_LAV", MeshUnion_dim::Ne_boundary, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "BoundaryEdge_M", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Ne, "BoundaryEdge_Ne", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "BoundaryEdge_Nfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx, "BoundaryEdge_nx", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(ny, "BoundaryEdge_ny", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(nz, "BoundaryEdge_nz", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Jz, "BoundaryEdge_Jz", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(xb, "BoundaryEdge_xb", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(yb, "BoundaryEdge_yb", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(zb, "BoundaryEdge_zb", MeshUnion_dim::Ne_boundary, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(V1d, "BoundaryEdge_V1d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(V2d, "BoundaryEdge_V2d", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Nlayer, "BoundaryEdge_Nz", MeshUnion_dim::one);

	//cout << "End BoundaryEdge" << endl;
}


BoundaryEdge::~BoundaryEdge()
{

}
