#include "Icell.h"
#include<iostream>

double* Icell::Dr2d = NULL;
double* Icell::Ds2d = NULL;
double* Icell::Dt2d = NULL;
double* Icell::faceType2d = NULL;
double* Icell::Fmask2d = NULL;
double* Icell::FToV2d = NULL;
double* Icell::invM2d = NULL;
double* Icell::LAV2d = NULL;
double* Icell::M2d = NULL;
double* Icell::N2d = NULL;
double* Icell::Nface2d = NULL;
double* Icell::Nfp2d = NULL;
double* Icell::Nfv2d = NULL;
double* Icell::Np2d = NULL;
double* Icell::Nq2d = NULL;
double* Icell::Nv2d = NULL;
double* Icell::r2d = NULL;
double* Icell::rq2d = NULL;
double* Icell::s2d = NULL;
double* Icell::sq2d = NULL;
double* Icell::t2d = NULL;
double* Icell::TNfp2d = NULL;
double* Icell::tq2d = NULL;
double* Icell::type2d = NULL;
double* Icell::V2d = NULL;
double* Icell::Vq2d = NULL;
double* Icell::vr2d = NULL;
double* Icell::vs2d = NULL;
double* Icell::vt2d = NULL;
double* Icell::wq2d = NULL;

Icell::Icell()//:Read_NC_dim2d()
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


Icell::~Icell()
{

}
