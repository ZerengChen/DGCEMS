#include "MeshUnion.h"


int *MeshUnion::K = NULL;
int *MeshUnion::Nv = NULL;
int *MeshUnion::Nlayer = NULL;
int MeshUnion::Nfield = 16;


//double *MeshUnion::charLength = NULL;
double *MeshUnion::EToE = NULL;
double *MeshUnion::EToF = NULL;
double *MeshUnion::EToM = NULL;
double *MeshUnion::EToL = NULL;
double *MeshUnion::EToV = NULL;
double *MeshUnion::ind = NULL;
double *MeshUnion::J = NULL;
double *MeshUnion::Jz = NULL;
//double *MeshUnion::K = NULL;
double *MeshUnion::LAV = NULL;
//double *MeshUnion::Nv = NULL;
double *MeshUnion::rx = NULL;
double *MeshUnion::ry = NULL;
double *MeshUnion::rz = NULL;
//signed char *MeshUnion::status = NULL;
double *MeshUnion::sx = NULL;
double *MeshUnion::sy = NULL;
double *MeshUnion::sz = NULL;
double *MeshUnion::tx = NULL;
double *MeshUnion::ty = NULL;
double *MeshUnion::type = NULL;
double *MeshUnion::tz = NULL;
double *MeshUnion::vx = NULL;
double *MeshUnion::vy = NULL;
double *MeshUnion::vz = NULL;
double *MeshUnion::x = NULL;
double *MeshUnion::xc = NULL;
double *MeshUnion::y = NULL;
double *MeshUnion::yc = NULL;
double *MeshUnion::z = NULL;
double *MeshUnion::zc = NULL;
Cell *MeshUnion::cell_p = NULL;
BoundaryEdge *MeshUnion::boundaryedge_p = NULL;
InnerEdge *MeshUnion::inneredge_p = NULL;
BottomEdge *MeshUnion::bottomedge_p = NULL;
BottomBoundaryEdge *MeshUnion::bottomboundaryedge_p = NULL;
SurfaceBoundaryEdge *MeshUnion::surfaceboundaryedge_p = NULL;
Mesh2d *MeshUnion::mesh2d_p = NULL;


MeshUnion::MeshUnion() :boundaryedge(), inneredge(), cell(), bottomedge(), bottomboundaryedge(), mesh2d(), surfaceboundaryedge()
{
	boundaryedge_p = &boundaryedge;
	inneredge_p = &inneredge;
	cell_p = &cell;
	bottomedge_p = &bottomedge;
	bottomboundaryedge_p = &bottomboundaryedge;
	surfaceboundaryedge_p = &surfaceboundaryedge;
	mesh2d_p = &mesh2d;

	//std::cout << "Have a nice time!\n";

	MeshUnion_dim::ncvar_read(K, "K", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nv, "Nv", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(Nlayer, "Nlayer", MeshUnion_dim::one);
	//MeshUnion_dim::ncvar_read(charLength, "charLength", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(EToE, "EToE", MeshUnion_dim::K, MeshUnion_dim::cell_Nface);
	MeshUnion_dim::ncvar_read(EToF, "EToF", MeshUnion_dim::K, MeshUnion_dim::cell_Nface);
	MeshUnion_dim::ncvar_read(EToM, "EToM", MeshUnion_dim::K, MeshUnion_dim::cell_Nface);
	MeshUnion_dim::ncvar_read(EToL, "EToL", MeshUnion_dim::one, MeshUnion_dim::K);
	MeshUnion_dim::ncvar_read(EToV, "EToV", MeshUnion_dim::K, MeshUnion_dim::cell_Nv);
	MeshUnion_dim::ncvar_read(ind, "ind", MeshUnion_dim::one);

	MeshUnion_dim::ncvar_read(J, "J", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(Jz, "Jz", MeshUnion_dim::K, MeshUnion_dim::Np);
	//ncvar_read(K, "K", one);
	MeshUnion_dim::ncvar_read(LAV, "LAV", MeshUnion_dim::K, MeshUnion_dim::one);
	//ncvar_read(Nv, "Nv", one);
	MeshUnion_dim::ncvar_read(rx, "rx", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(ry, "ry", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(rz, "rz", MeshUnion_dim::K, MeshUnion_dim::Np);
	//MeshUnion_dim::ncvar_read(status, "status", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(sx, "sx", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(sy, "sy", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(sz, "sz", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(tx, "tx", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(ty, "ty", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(type, "type", MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(tz, "tz", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(vx, "vx", MeshUnion_dim::one, MeshUnion_dim::Nv);
	MeshUnion_dim::ncvar_read(vy, "vy", MeshUnion_dim::one, MeshUnion_dim::Nv);
	MeshUnion_dim::ncvar_read(vz, "vz", MeshUnion_dim::one, MeshUnion_dim::Nv);
	MeshUnion_dim::ncvar_read(x, "x", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(xc, "xc", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(y, "y", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(yc, "yc", MeshUnion_dim::K, MeshUnion_dim::one);
	MeshUnion_dim::ncvar_read(z, "z", MeshUnion_dim::K, MeshUnion_dim::Np);
	MeshUnion_dim::ncvar_read(zc, "zc", MeshUnion_dim::K, MeshUnion_dim::one);

	//std::cout << "End reading the meshunion!\n";
}


MeshUnion::~MeshUnion()
{
	std::cout << "Free MeshUnion" << std::endl;
}
