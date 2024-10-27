#include "Cell.h"
//#include<netcdfcpp.h>
#include<algorithm>
using namespace std;

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

signed char *Cell::type = NULL;
int *Cell::Nq = NULL;
int *Cell::Nv = NULL;

double *Cell::Dr = NULL;
double *Cell::Ds = NULL;
double *Cell::Dt = NULL;
double *Cell::FToV = NULL;
signed char *Cell::faceType = NULL;
double *Cell::Fmask = NULL;
double *Cell::invM = NULL;
double *Cell::LAV = NULL;
double *Cell::M = NULL;
int *Cell::N = NULL;
int *Cell::Nz = NULL;
int *Cell::Nface = NULL;
int *Cell::Nfp = NULL;
double *Cell::Nfv = NULL;
int *Cell::Np = NULL;
int *Cell::Nph = NULL;
int *Cell::Npz = NULL;
double *Cell::V1d = NULL;
double *Cell::r = NULL;
double *Cell::rq = NULL;
double *Cell::s = NULL;
double *Cell::sq = NULL;
double *Cell::t = NULL;
int *Cell::TNfp = NULL;
double *Cell::tq = NULL;
//double *Cell::type = NULL;
double *Cell::V = NULL;
double *Cell::Vq = NULL;
double *Cell::vr = NULL;
double *Cell::vs = NULL;
double *Cell::vt = NULL;
double *Cell::wq = NULL;
double *Cell::r1 = NULL;
double *Cell::s1 = NULL;
double *Cell::t1 = NULL;
double *Cell::Vh = NULL;
double *Cell::Vint = NULL;
double *Cell::VCV = NULL;
#ifdef _BAROCLINIC
double *Cell::VintU = NULL;
#endif

Cell::Cell()//:Read_NC_dim()
{
	//cout<<"Start Cell"<<endl;
	MeshUnion_dim::ncvar_read(type, "cell_type", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nq, "cell_Nq", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv, "cell_Nv", MeshUnion_dim::one);

	MeshUnion_dim::ncvar_read(Dr, "cell_Dr", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(Ds, "cell_Ds", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(Dt, "cell_Dt", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(FToV, "cell_FToV", MeshUnion_dim::cell_Nface, MeshUnion_dim::four);
	MeshUnion_dim::ncvar_read(faceType, "cell_faceType", MeshUnion_dim::cell_Nface, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Fmask, "cell_Fmask", MeshUnion_dim::cell_Nface, max(MeshUnion_dim::Nfp , MeshUnion_dim::NfpBottom));
	MeshUnion_dim::ncvar_read(invM, "cell_invM", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(LAV, "cell_LAV", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M, "cell_M", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(N, "cell_N", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nz, "cell_Nz", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nface, "cell_Nface", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp, "cell_Nfp", MeshUnion_dim::one, MeshUnion_dim::cell_Nface);
	MeshUnion_dim::ncvar_read(Nfv, "cell_Nfv", MeshUnion_dim::cell_Nface, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Np, "cell_Np", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nph, "cell_Nph", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Npz, "cell_Npz", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(V1d, "cell_V1d", MeshUnion_dim::Nzadd1, MeshUnion_dim::Nzadd1);//Nz+1
	MeshUnion_dim::ncvar_read(r, "cell_r", MeshUnion_dim::one, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(rq, "cell_rq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(s, "cell_s", MeshUnion_dim::one, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(sq, "cell_sq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(t, "cell_t", MeshUnion_dim::one, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(TNfp, "cell_TNfp", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tq, "cell_tq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
	//MeshUnion_dim::ncvar_read(type, "cell_type", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(V, "cell_V", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(Vq, "cell_Vq", MeshUnion_dim::Np, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(vr, "cell_vr", MeshUnion_dim::one, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(vs, "cell_vs", MeshUnion_dim::one, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(vt, "cell_vt", MeshUnion_dim::one, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(wq, "cell_wq", MeshUnion_dim::one, MeshUnion_dim::cell_Nq);
	MeshUnion_dim::ncvar_read(r1, "cell_r1", MeshUnion_dim::one, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(s1, "cell_s1", MeshUnion_dim::one, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(t1, "cell_t1", MeshUnion_dim::one, MeshUnion_dim::Nzadd1);
	MeshUnion_dim::ncvar_read(Vh, "cell_Vh", MeshUnion_dim::NfpBottom, MeshUnion_dim::NfpBottom);
	MeshUnion_dim::ncvar_read(Vint, "cell_Vint", MeshUnion_dim::Np, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(VCV, "cell_VCV", MeshUnion_dim::Np, MeshUnion_dim::NfpBottom);
#ifdef _BAROCLINIC
	MeshUnion_dim::ncvar_read(VintU, "cell_VintU", MeshUnion_dim::Np, MeshUnion_dim::Np);
#endif

	//cout << "End Cell"<< endl;
}


Cell::~Cell()
{

}
