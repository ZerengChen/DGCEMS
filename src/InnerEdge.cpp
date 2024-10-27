#include "InnerEdge.h"
using namespace std;


double *InnerEdge::FToE = NULL;
double *InnerEdge::FToF = NULL;
double *InnerEdge::FToM = NULL;
double *InnerEdge::FToN1 = NULL;
double *InnerEdge::FToN2 = NULL;
double *InnerEdge::FToV = NULL;
double *InnerEdge::Js = NULL;
double *InnerEdge::Jz = NULL;
double *InnerEdge::LAV = NULL;
double *InnerEdge::M = NULL;
int *InnerEdge::Ne = NULL;
int *InnerEdge::Nfp = NULL;
double *InnerEdge::nx = NULL;
double *InnerEdge::ny = NULL;
double *InnerEdge::nz = NULL;
double *InnerEdge::V1d = NULL;
double *InnerEdge::V2d = NULL;
int *InnerEdge::Nlayer = NULL;
//double *InnerEdge::r = NULL;

InnerEdge::InnerEdge() //:icell()
{
	//cout<<"Start InnerEdge"<<endl;
	MeshUnion_dim::ncvar_read(FToE, "InnerEdge_FToE", MeshUnion_dim::Ne_inner, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToF, "InnerEdge_FToF", MeshUnion_dim::Ne_inner, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(FToM, "InnerEdge_FToM", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToN1, "InnerEdge_FToN1", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToN2, "InnerEdge_FToN2", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(FToV, "InnerEdge_FToV", MeshUnion_dim::Ne_inner, MeshUnion_dim::four);
	MeshUnion_dim::ncvar_read(Js, "InnerEdge_Js", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Jz, "InnerEdge_Jz", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(LAV, "InnerEdge_LAV", MeshUnion_dim::Ne_inner, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "InnerEdge_M", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Ne, "InnerEdge_Ne", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "InnerEdge_Nfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(nx, "InnerEdge_nx", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(ny, "InnerEdge_ny", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(nz, "InnerEdge_nz", MeshUnion_dim::Ne_inner, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(V1d, "InnerEdge_V1d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(V2d, "InnerEdge_V2d", MeshUnion_dim::Nfp, MeshUnion_dim::Nfp);
	MeshUnion_dim::ncvar_read(Nlayer, "InnerEdge_Nlayer", MeshUnion_dim::one);
	//MeshUnion_dim::ncvar_read(r, "InnerEdge_r", MeshUnion_dim::one, MeshUnion_dim::Nfp);
	//cout << "End InnerEdge" << endl;
}


InnerEdge::~InnerEdge()
{

}

