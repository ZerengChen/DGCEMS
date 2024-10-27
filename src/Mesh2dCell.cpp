#include "Mesh2dCell.h"
#include<netcdfcpp.h>
#include<algorithm>


int *Mesh2dCell::Nq2d = NULL;
int *Mesh2dCell::Nv2d = NULL;

double *Mesh2dCell::Dr2d = NULL;
double *Mesh2dCell::Ds2d = NULL;
double *Mesh2dCell::Dt2d = NULL;
double *Mesh2dCell::FToV2d = NULL;
signed char *Mesh2dCell::faceType2d = NULL;
double *Mesh2dCell::Fmask2d = NULL;
double *Mesh2dCell::invM2d = NULL;
double *Mesh2dCell::LAV2d = NULL;
double *Mesh2dCell::M2d = NULL;
int *Mesh2dCell::N2d = NULL;
int *Mesh2dCell::Nface2d = NULL;
int *Mesh2dCell::Nfp2d = NULL;
double *Mesh2dCell::Nfv2d = NULL;
int *Mesh2dCell::Np2d = NULL;
double *Mesh2dCell::r2d = NULL;
double *Mesh2dCell::rq2d = NULL;
double *Mesh2dCell::s2d = NULL;
double *Mesh2dCell::sq2d = NULL;
double *Mesh2dCell::t2d = NULL;
int *Mesh2dCell::TNfp2d = NULL;
double *Mesh2dCell::tq2d = NULL;
double *Mesh2dCell::type2d = NULL;
double *Mesh2dCell::V2d = NULL;
double *Mesh2dCell::Vq2d = NULL;
double *Mesh2dCell::vr2d = NULL;
double *Mesh2dCell::vs2d = NULL;
double *Mesh2dCell::vt2d = NULL;
double *Mesh2dCell::wq2d = NULL;


Mesh2dCell::Mesh2dCell()
{
	MeshUnion_dim::ncvar_read(Nq2d, "cell_Nq2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv2d, "cell_Nv2d", MeshUnion_dim::one);

	MeshUnion_dim::ncvar_read(Dr2d, "cell_Dr2d", MeshUnion_dim::Np2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(Ds2d, "cell_Ds2d", MeshUnion_dim::Np2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(Dt2d, "cell_Dt2d", MeshUnion_dim::Np2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(FToV2d, "cell_FToV2d", MeshUnion_dim::cell_Nv2d, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(faceType2d, "cell_faceType2d", MeshUnion_dim::cell_Nv2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Fmask2d, "cell_Fmask2d", MeshUnion_dim::cell_Nv2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(invM2d, "cell_invM2d", MeshUnion_dim::Np2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(LAV2d, "cell_LAV2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M2d, "cell_M2d", MeshUnion_dim::Np2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(N2d, "cell_N2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nface2d, "cell_Nface2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp2d, "cell_Nfp2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(Nfv2d, "cell_Nfv2d", MeshUnion_dim::cell_Nv2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Np2d, "cell_Np2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(r2d, "cell_r2d", MeshUnion_dim::one, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(rq2d, "cell_rq2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nq2d);
	MeshUnion_dim::ncvar_read(s2d, "cell_s2d", MeshUnion_dim::one, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(sq2d, "cell_sq2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nq2d);
	MeshUnion_dim::ncvar_read(t2d, "cell_t2d", MeshUnion_dim::one, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(TNfp2d, "cell_TNfp2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tq2d, "cell_tq2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nq2d);
	MeshUnion_dim::ncvar_read(type2d, "cell_type2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(V2d, "cell_V2d", MeshUnion_dim::Np2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(Vq2d, "cell_Vq2d", MeshUnion_dim::Np2d, MeshUnion_dim::cell_Nq2d);
	MeshUnion_dim::ncvar_read(vr2d, "cell_vr2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(vs2d, "cell_vs2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(vt2d, "cell_vt2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(wq2d, "cell_wq2d", MeshUnion_dim::one, MeshUnion_dim::cell_Nq2d);
}


Mesh2dCell::~Mesh2dCell()
{


}
