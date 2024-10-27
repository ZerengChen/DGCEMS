#pragma once
#include "netcdfcpp.h"
#include <iostream>

class MeshUnion_dim
{
public:
	MeshUnion_dim();
	~MeshUnion_dim();

	static void ncdim_read();
	static void ncvar_read(double *&meshunion_data, const char* ncvarname, int &dim1, int &dim2);
	static void ncvar_read(double *&meshunion_data, const char* ncvarname, int &dim1);
	static void ncvar_read(int *&meshunion_data, const char* ncvarname, int &dim1, int &dim2);
	static void ncvar_read(int *&meshunion_data, const char* ncvarname, int &dim1);
	static void ncvar_read(signed char *&meshunion_data, const char* ncvarname, int &dim1, int &dim2);
	static void ncvar_read(signed char *&meshunion_data, const char* ncvarname, int &dim1);

	static int K; //total elements in a layer
	static int Nv;//total vertices in a layer
	static int Nlayer;//the layer 
	static int N;//stage of horizontal
	static int Nz;//stage of vertical
	static int Ne_inner;
	static int Ne_boundary;
	static int Ne_bottom;
	static int Ne_bottomboundary;
	static int Ne_surface;
	static int Nfp;
	static int NfpBottom;
	static int Np;//in a cell
	static int cell_Nv;
	static int cell_Nq;
	static int cell_Nface;
	static int one;
	static int two;
	static int three;
	static int four;
	static int Nzadd1;
	
	//static int Nadd1;
	//2Dmesh
	static int K2d; 
	static int Nv2d;
	static int Ne_inner2d;
	static int Ne_boundary2d;
	static int Nfp2d;
	static int Np2d;
	static int cell_Nv2d;
	static int cell_Nq2d;
	static int icell_Nq;

};

