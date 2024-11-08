#include "MeshUnion_dim.h"


int MeshUnion_dim::K = 0;
int MeshUnion_dim::Nv = 0;
int MeshUnion_dim::Nlayer = 0;
int MeshUnion_dim::N = 0;
int MeshUnion_dim::Nz = 0;
int	MeshUnion_dim::Ne_inner = 0;
int MeshUnion_dim::Ne_boundary = 0;
int	MeshUnion_dim::Ne_bottom = 0;
int MeshUnion_dim::Ne_bottomboundary = 0;
int	MeshUnion_dim::Ne_surface = 0;
int MeshUnion_dim::Nfp = 0;
int MeshUnion_dim::NfpBottom = 0;
int MeshUnion_dim::Np = 0;
int MeshUnion_dim::cell_Nv = 0;
int MeshUnion_dim::cell_Nq = 0;
int MeshUnion_dim::cell_Nface = 0;
int MeshUnion_dim::one = 1;
int MeshUnion_dim::two = 2;
int MeshUnion_dim::three = 3;
int MeshUnion_dim::four = 4;
int MeshUnion_dim::Nzadd1 = 0;
//int MeshUnion_dim::Nadd1 = 0;

int MeshUnion_dim::K2d = 0;
int MeshUnion_dim::Nv2d = 0;
int	MeshUnion_dim::Ne_inner2d = 0;
int MeshUnion_dim::Ne_boundary2d = 0;
int MeshUnion_dim::Nfp2d = 2;
int MeshUnion_dim::Np2d = 0;
int MeshUnion_dim::cell_Nv2d = 0;
int MeshUnion_dim::cell_Nq2d = 0;
int MeshUnion_dim::icell_Nq = 0;


MeshUnion_dim::MeshUnion_dim()
{
	ncdim_read();
}

MeshUnion_dim::~MeshUnion_dim()
{
}

void MeshUnion_dim::ncdim_read()
{
	static const int NC_ERR = 2;
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	if (!dataFile.is_valid())
	{
		std::cout << "Couldn't open file!\n";
		// FILE *fp = fopen("123.txt", "a+");
		// fclose(fp);
	}

	NcVar *K_v = dataFile.get_var("K");
	NcVar *Nv_v = dataFile.get_var("Nv");
	NcVar *Nlayer_v = dataFile.get_var("Nlayer");//chenzereng
	NcVar *N_v = dataFile.get_var("N");//chenzereng
	NcVar *Nz_v = dataFile.get_var("Nz");//chenzereng
	NcVar *Ne_inner_v = dataFile.get_var("Ne_inner");
	NcVar *Ne_boundary_v = dataFile.get_var("Ne_boundary");
	NcVar *Ne_bottom_v = dataFile.get_var("Ne_bottom");
	NcVar *Ne_bottomboundary_v = dataFile.get_var("Ne_bottomboundary");
	NcVar *Ne_surface_v = dataFile.get_var("Ne_surface");
	NcVar *Nfp_v = dataFile.get_var("Nfp");
	NcVar *NfpBottom_v = dataFile.get_var("NfpBottom");
	NcVar *Np_v = dataFile.get_var("Np");
	NcVar *cell_Nv_v = dataFile.get_var("cell_Nv");
	NcVar *cell_Nq_v = dataFile.get_var("cell_Nq");
	NcVar *cell_Nface_v = dataFile.get_var("cell_Nface");
	NcVar *Nzadd1_v = dataFile.get_var("Nzadd1");//chenzereng
	//NcVar *Nadd1_v = dataFile.get_var("Nadd1");//chenzereng

	NcVar *K2d_v = dataFile.get_var("K2d");
	NcVar *Nv2d_v = dataFile.get_var("Nv2d");
	NcVar *Ne_inner2d_v = dataFile.get_var("Ne_inner2d");
	NcVar *Ne_boundary2d_v = dataFile.get_var("Ne_boundary2d");
	NcVar *Nfp2d_v = dataFile.get_var("InnerEdge_cell_Nfp2d");
	NcVar *Np2d_v = dataFile.get_var("cell_Np2d");
	NcVar *cell_Nv2d_v = dataFile.get_var("cell_Nv2d");
	NcVar *cell_Nq2d_v = dataFile.get_var("cell_Nq2d");
	NcVar *icell_Nq_v = dataFile.get_var("InnerEdge_cell_Nq2d");

	K_v->get(&K, 1);
	Nv_v->get(&Nv, 1);
	Nlayer_v->get(&Nlayer, 1);
	N_v->get(&N, 1);
	Nz_v->get(&Nz, 1);
	Ne_inner_v->get(&Ne_inner, 1);
	Ne_boundary_v->get(&Ne_boundary, 1);
	Ne_bottom_v->get(&Ne_bottom, 1);
	Ne_bottomboundary_v->get(&Ne_bottomboundary, 1);
	Ne_surface_v->get(&Ne_surface, 1);
	Nfp_v->get(&Nfp, 1);
	NfpBottom_v->get(&NfpBottom, 1);
	Np_v->get(&Np, 1);
	cell_Nv_v->get(&cell_Nv, 1);
	cell_Nq_v->get(&cell_Nq, 1);
	cell_Nface_v->get(&cell_Nface, 1);
	Nzadd1_v->get(&Nzadd1, 1);
	//Nadd1_v->get(&Nadd1, 1);

	K2d_v->get(&K2d, 1);
	Nv2d_v->get(&Nv2d, 1);
	Ne_inner2d_v->get(&Ne_inner2d, 1);
	Ne_boundary2d_v->get(&Ne_boundary2d, 1);
	Nfp2d_v->get(&Nfp2d, 1);
	Np2d_v->get(&Np2d, 1);
	cell_Nv2d_v->get(&cell_Nv2d, 1);
	cell_Nq2d_v->get(&cell_Nq2d, 1);
	icell_Nq_v->get(&icell_Nq, 1);
}



void MeshUnion_dim::ncvar_read(double *&meshunion_data, const char* ncvarname, int &dim1, int &dim2)
{
	meshunion_data = (double*)malloc(dim1*dim2 * sizeof(double));
	//std::cout << dim1 << "*" << dim2 << std::endl;
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1, dim2);
}

void MeshUnion_dim::ncvar_read(double *&meshunion_data, const char* ncvarname, int &dim1)
{
	meshunion_data = (double*)malloc(dim1 * sizeof(double));
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1);
}

void MeshUnion_dim::ncvar_read(int *&meshunion_data, const char* ncvarname, int &dim1, int &dim2)
{
	meshunion_data = (int*)malloc(dim1*dim2 * sizeof(int));
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1, dim2);
}

void MeshUnion_dim::ncvar_read(int *&meshunion_data, const char* ncvarname, int &dim1)
{
	meshunion_data = (int*)malloc(dim1 * sizeof(int));
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1);
}

void MeshUnion_dim::ncvar_read(signed char *&meshunion_data, const char* ncvarname, int &dim1, int &dim2)
{
	meshunion_data = (signed char*)malloc(dim1*dim2 * sizeof(signed char));
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1, dim2);
}

void MeshUnion_dim::ncvar_read(signed char *&meshunion_data, const char* ncvarname, int &dim1)
{
	meshunion_data = (signed char*)malloc(dim1 * sizeof(signed char));
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1);
}