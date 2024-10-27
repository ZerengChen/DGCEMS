#include "Bcell.h"
#include<iostream>

double* Bcell::Dr2d = NULL;
double* Bcell::Ds2d = NULL;
double* Bcell::Dt2d = NULL;
double* Bcell::faceType2d = NULL;
double* Bcell::Fmask2d = NULL;
double* Bcell::FToV2d = NULL;
double* Bcell::invM2d = NULL;
double* Bcell::LAV2d = NULL;
double* Bcell::M2d = NULL;
double* Bcell::N2d = NULL;
double* Bcell::Nface2d = NULL;
double* Bcell::Nfp2d = NULL;
double* Bcell::Nfv2d = NULL;
double* Bcell::Np2d = NULL;
double* Bcell::Nq2d = NULL;
double* Bcell::Nv2d = NULL;
double* Bcell::r2d = NULL;
double* Bcell::rq2d = NULL;
double* Bcell::s2d = NULL;
double* Bcell::sq2d = NULL;
double* Bcell::t2d = NULL;
double* Bcell::TNfp2d = NULL;
double* Bcell::tq2d = NULL;
double* Bcell::type2d = NULL;
double* Bcell::V2d = NULL;
double* Bcell::Vq2d = NULL;
double* Bcell::vr2d = NULL;
double* Bcell::vs2d = NULL;
double* Bcell::vt2d = NULL;
double* Bcell::wq2d = NULL;


Bcell::Bcell()
{
	MeshUnion_dim::ncvar_read(Dr2d, "InnerEdge_cell_Dr2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(Ds2d, "InnerEdge_cell_Ds2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(Dt2d, "InnerEdge_cell_Dt2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(faceType2d,"InnerEdge_cell_faceType2d", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Fmask2d, "InnerEdge_cell_Fmask2d", MeshUnion_dim::two, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(FToV2d, "InnerEdge_cell_FToV2d", MeshUnion_dim::two, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(invM2d, "InnerEdge_cell_invM2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(LAV2d, "InnerEdge_cell_LAV2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(M2d, "InnerEdge_cell_M2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(N2d, "InnerEdge_cell_N2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nface2d, "InnerEdge_cell_Nface2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nfp2d, "InnerEdge_cell_Nfp2d", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Nfv2d, "InnerEdge_cell_Nfv2d", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(Np2d, "InnerEdge_cell_Np2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nq2d, "InnerEdge_cell_Nq2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv2d, "InnerEdge_cell_Nv2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(r2d, "InnerEdge_cell_r2d", MeshUnion_dim::one, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(rq2d, "InnerEdge_cell_rq2d", MeshUnion_dim::one, MeshUnion_dim::icell_Nq);
	MeshUnion_dim::ncvar_read(s2d, "InnerEdge_cell_s2d", MeshUnion_dim::one, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(sq2d, "InnerEdge_cell_sq2d", MeshUnion_dim::one, MeshUnion_dim::icell_Nq);
	MeshUnion_dim::ncvar_read(t2d, "InnerEdge_cell_t2d", MeshUnion_dim::one, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(TNfp2d, "InnerEdge_cell_TNfp2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tq2d, "InnerEdge_cell_tq2d", MeshUnion_dim::one, MeshUnion_dim::icell_Nq);
	MeshUnion_dim::ncvar_read(type2d, "InnerEdge_cell_type2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(V2d, "InnerEdge_cell_V2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::Nfp2d);
	MeshUnion_dim::ncvar_read(Vq2d, "InnerEdge_cell_Vq2d", MeshUnion_dim::Nfp2d, MeshUnion_dim::icell_Nq);
	MeshUnion_dim::ncvar_read(vr2d, "InnerEdge_cell_vr2d", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(vs2d, "InnerEdge_cell_vs2d", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(vt2d, "InnerEdge_cell_vt2d", MeshUnion_dim::one, MeshUnion_dim::two);
	MeshUnion_dim::ncvar_read(wq2d, "InnerEdge_cell_wq2d", MeshUnion_dim::one, MeshUnion_dim::icell_Nq);

}


Bcell::~Bcell()
{

}
