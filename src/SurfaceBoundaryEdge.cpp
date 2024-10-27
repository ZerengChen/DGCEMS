#include "SurfaceBoundaryEdge.h"



double *SurfaceBoundaryEdge::FToE = NULL;
double *SurfaceBoundaryEdge::FToF = NULL;
double *SurfaceBoundaryEdge::FToM = NULL;
double *SurfaceBoundaryEdge::FToN1 = NULL;
double *SurfaceBoundaryEdge::FToN2 = NULL;
double *SurfaceBoundaryEdge::FToV = NULL;
double *SurfaceBoundaryEdge::Js = NULL;
double *SurfaceBoundaryEdge::LAV = NULL;
double *SurfaceBoundaryEdge::M = NULL;
int *SurfaceBoundaryEdge::Ne = NULL;
int *SurfaceBoundaryEdge::Nfp = NULL;
double *SurfaceBoundaryEdge::nx = NULL;
double *SurfaceBoundaryEdge::ny = NULL;
double *SurfaceBoundaryEdge::nz = NULL;
double *SurfaceBoundaryEdge::ftype = NULL;

SurfaceBoundaryEdge::SurfaceBoundaryEdge() //:icell()
{
	//std::cout << "Start SurfaceBoundaryEdge" << std::endl;
	MeshUnion_dim::ncvar_read(FToE, "SurfaceBoundaryEdge_FToE", MeshUnion_dim::Ne_surface, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF, "SurfaceBoundaryEdge_FToF", MeshUnion_dim::Ne_surface, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM, "SurfaceBoundaryEdge_FToM", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToN1, "SurfaceBoundaryEdge_FToN1", MeshUnion_dim::Ne_surface, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(FToN2, "SurfaceBoundaryEdge_FToN2", MeshUnion_dim::Ne_surface, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(FToV, "SurfaceBoundaryEdge_FToV", MeshUnion_dim::Ne_surface, MeshUnion_dim::three);
	MeshUnion_dim::ncvar_read(Js, "SurfaceBoundaryEdge_Js", MeshUnion_dim::Ne_surface, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(LAV, "SurfaceBoundaryEdge_LAV", MeshUnion_dim::Ne_surface, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "SurfaceBoundaryEdge_M", MeshUnion_dim::NfpBottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(Ne, "SurfaceBoundaryEdge_Ne", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "SurfaceBoundaryEdge_Nfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx, "SurfaceBoundaryEdge_nx", MeshUnion_dim::Ne_surface, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(ny, "SurfaceBoundaryEdge_ny", MeshUnion_dim::Ne_surface, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(nz, "SurfaceBoundaryEdge_nz", MeshUnion_dim::Ne_surface, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(ftype, "SurfaceBoundaryEdge_ftype", MeshUnion_dim::one, MeshUnion_dim::Ne_surface);
	//std::cout << "End SurfaceBoundaryEdge" << std::endl;
}


SurfaceBoundaryEdge::~SurfaceBoundaryEdge()
{

}

