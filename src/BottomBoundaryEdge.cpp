#include "BottomBoundaryEdge.h"



double *BottomBoundaryEdge::FToE = NULL;
double *BottomBoundaryEdge::FToF = NULL;
double *BottomBoundaryEdge::FToM = NULL;
double *BottomBoundaryEdge::FToN1 = NULL;
double *BottomBoundaryEdge::FToN2 = NULL;
double *BottomBoundaryEdge::FToV = NULL;
double *BottomBoundaryEdge::Js = NULL;
double *BottomBoundaryEdge::LAV = NULL;
double *BottomBoundaryEdge::M = NULL;
int *BottomBoundaryEdge::Ne = NULL;
int *BottomBoundaryEdge::Nfp = NULL;
double *BottomBoundaryEdge::nx = NULL;
double *BottomBoundaryEdge::ny = NULL;
double *BottomBoundaryEdge::nz = NULL;
double *BottomBoundaryEdge::ftype = NULL;

BottomBoundaryEdge::BottomBoundaryEdge() //:icell()
{
	//std::cout << "Start BottomBoundaryEdge" << std::endl;
	MeshUnion_dim::ncvar_read(FToE, "BottomBoundaryEdge_FToE", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF, "BottomBoundaryEdge_FToF", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM, "BottomBoundaryEdge_FToM", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToN1, "BottomBoundaryEdge_FToN1", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(FToN2, "BottomBoundaryEdge_FToN2", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(FToV, "BottomBoundaryEdge_FToV", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::three);
	MeshUnion_dim::ncvar_read(Js, "BottomBoundaryEdge_Js", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(LAV, "BottomBoundaryEdge_LAV", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "BottomBoundaryEdge_M", MeshUnion_dim::NfpBottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(Ne, "BottomBoundaryEdge_Ne", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "BottomBoundaryEdge_Nfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx, "BottomBoundaryEdge_nx", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(ny, "BottomBoundaryEdge_ny", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(nz, "BottomBoundaryEdge_nz", MeshUnion_dim::Ne_bottomboundary, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(ftype, "BottomBoundaryEdge_ftype", MeshUnion_dim::one, MeshUnion_dim::Ne_bottomboundary);

	//std::cout << "End BottomBoundaryEdge" << std::endl;
}


BottomBoundaryEdge::~BottomBoundaryEdge()
{

}
