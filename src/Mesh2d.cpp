#include "Mesh2d.h"
#include<netcdfcpp.h>
#include<algorithm>


int *Mesh2d::K2d = NULL;
int *Mesh2d::Nv2d = NULL;
int Mesh2d::Nfield2d = 8;
double *Mesh2d::charLength2d = NULL;
double *Mesh2d::EToE2d = NULL;
double *Mesh2d::EToF2d = NULL;
double *Mesh2d::EToM2d = NULL;
double *Mesh2d::EToR2d = NULL;
double *Mesh2d::EToV2d = NULL;
double *Mesh2d::ind2d = NULL;
double *Mesh2d::J2d = NULL;
//double *MeshUnion::K = NULL;
double *Mesh2d::LAV2d = NULL;
//double *MeshUnion::Nv = NULL;
double *Mesh2d::rx2d = NULL;
double *Mesh2d::ry2d = NULL;
double *Mesh2d::rz2d = NULL;
signed char *Mesh2d::status = NULL;
double *Mesh2d::sx2d = NULL;
double *Mesh2d::sy2d = NULL;
double *Mesh2d::sz2d = NULL;
double *Mesh2d::tx2d = NULL;
double *Mesh2d::ty2d = NULL;
double *Mesh2d::type2d = NULL;
double *Mesh2d::tz2d = NULL;
double *Mesh2d::vx2d = NULL;
double *Mesh2d::vy2d = NULL;
double *Mesh2d::vz2d = NULL;
double *Mesh2d::x2d = NULL;
double *Mesh2d::xc2d = NULL;
double *Mesh2d::y2d = NULL;
double *Mesh2d::yc2d = NULL;
double *Mesh2d::z2d = NULL;
double *Mesh2d::zc2d = NULL;
Mesh2dCell *Mesh2d::mesh2dcell_p = NULL;
Mesh2dBoundaryEdge *Mesh2d::mesh2dboundaryedge_p = NULL;
Mesh2dInnerEdge *Mesh2d::mesh2dinneredge_p = NULL;


Mesh2d::Mesh2d() :mesh2dboundaryedge(), mesh2dinneredge(), mesh2dcell()
{
	mesh2dboundaryedge_p = &mesh2dboundaryedge;
	mesh2dinneredge_p = &mesh2dinneredge;
	mesh2dcell_p = &mesh2dcell;

	//std::cout<<"Start Mesh2d"<<std::endl;
	MeshUnion_dim::ncvar_read(K2d, "K2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv2d, "Nv2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(charLength2d, "charLength2d", MeshUnion_dim::K2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(EToE2d, "EToE2d", MeshUnion_dim::K2d, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(EToF2d, "EToF2d", MeshUnion_dim::K2d, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(EToM2d, "EToM2d", MeshUnion_dim::K2d, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(EToR2d, "EToR2d", MeshUnion_dim::one, MeshUnion_dim::K2d);
	MeshUnion_dim::ncvar_read(EToV2d, "EToV2d", MeshUnion_dim::K2d, MeshUnion_dim::cell_Nv2d);
	MeshUnion_dim::ncvar_read(ind2d, "ind2d", MeshUnion_dim::one);

	MeshUnion_dim::ncvar_read(J2d, "J2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(LAV2d, "LAV2d", MeshUnion_dim::K2d, MeshUnion_dim::one);
	//ncvar_read(Nv, "Nv", one);
	MeshUnion_dim::ncvar_read(rx2d, "rx2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(ry2d, "ry2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(rz2d, "rz2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(status, "status", MeshUnion_dim::K2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(sx2d, "sx2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(sy2d, "sy2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(sz2d, "sz2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(tx2d, "tx2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(ty2d, "ty2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(type2d, "type2d", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tz2d, "tz2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(vx2d, "vx2d", MeshUnion_dim::one, MeshUnion_dim::Nv2d);
	MeshUnion_dim::ncvar_read(vy2d, "vy2d", MeshUnion_dim::one, MeshUnion_dim::Nv2d);
	MeshUnion_dim::ncvar_read(vz2d, "vz2d", MeshUnion_dim::one, MeshUnion_dim::Nv2d);
	MeshUnion_dim::ncvar_read(x2d, "x2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(xc2d, "xc2d", MeshUnion_dim::K2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(y2d, "y2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(yc2d, "yc2d", MeshUnion_dim::K2d, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(z2d, "z2d", MeshUnion_dim::K2d, MeshUnion_dim::Np2d);
	MeshUnion_dim::ncvar_read(zc2d, "zc2d", MeshUnion_dim::K2d, MeshUnion_dim::one);
	//std::cout << "End Mesh2d" << std::endl;
}


Mesh2d::~Mesh2d()
{

}
