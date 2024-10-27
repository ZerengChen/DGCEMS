#include "BottomEdge.h"



double *BottomEdge::FToE = NULL;
double *BottomEdge::FToF = NULL;
double *BottomEdge::FToM = NULL;
double *BottomEdge::FToN1 = NULL;
double *BottomEdge::FToN2 = NULL;
double *BottomEdge::FToV = NULL;
double *BottomEdge::Js = NULL;
double *BottomEdge::LAV = NULL;
double *BottomEdge::M = NULL;
int *BottomEdge::Ne = NULL;
int *BottomEdge::Nfp = NULL;
double *BottomEdge::nx = NULL;
double *BottomEdge::ny = NULL;
double *BottomEdge::nz = NULL;

BottomEdge::BottomEdge() //:icell()
{
	//std::cout<<"Start BottomEdge"<<std::endl;
	MeshUnion_dim::ncvar_read(FToE, "BottomEdge_FToE", MeshUnion_dim::Ne_bottom, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF, "BottomEdge_FToF", MeshUnion_dim::Ne_bottom, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM, "BottomEdge_FToM", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToN1, "BottomEdge_FToN1", MeshUnion_dim::Ne_bottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(FToN2, "BottomEdge_FToN2", MeshUnion_dim::Ne_bottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(FToV, "BottomEdge_FToV", MeshUnion_dim::K, MeshUnion_dim::three);
	MeshUnion_dim::ncvar_read(Js, "BottomEdge_Js", MeshUnion_dim::Ne_bottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(LAV, "BottomEdge_LAV", MeshUnion_dim::Ne_bottom, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "BottomEdge_M", MeshUnion_dim::NfpBottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(Ne, "BottomEdge_Ne", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "BottomEdge_Nfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx, "BottomEdge_nx", MeshUnion_dim::Ne_bottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(ny, "BottomEdge_ny", MeshUnion_dim::Ne_bottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(nz, "BottomEdge_nz", MeshUnion_dim::Ne_bottom, MeshUnion_dim::NfpBottom);
	//MeshUnion_dim::ncvar_read(r, "InnerEdge_r", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	//std::cout << "End BottomEdge" << std::endl;
}


BottomEdge::~BottomEdge()
{

}
