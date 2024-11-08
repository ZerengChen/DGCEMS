#include "NdgPhysMat.h"
#include "CalculateWD.h"
#include "VertLimit3d.h"

#ifdef _BAROCLINIC
#include "CalculateDensityField.h"
#include "CalculateBaroclinicTerm.h"
#endif

#ifdef COUPLING_SWAN
#include"c_coupler_cpp_interface_wrj.h"
#endif

#include <chrono>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <unistd.h>

using namespace std;
double ADVtime = 0.0;
double DIF_Htime = 0.0;
double PCEtime = 0.0;
double SOURCEtime = 0.0;
double OMGtime = 0.0;
double DIF_Vtime = 0.0;

#ifdef _BAROCLINIC
double Baroclinictime = 0.0;
#endif

char buff[255];
void newTecfile();

extern int* dg_demo_comp_id;//for coupler
extern int* decomp_id;
extern int* grid_h2d_id;//for coupler
extern int* time_step;//for coupler
extern int* coupling_freq;//for coupler
extern int First_CouplingStep;//for coupler
extern int Calculating_CouplingSteps;//for coupler
extern int  Switch_Limiter3D;
extern double gra;
extern double Hcrit;
extern int Nvar;
extern int Nfield3d;
extern int Nfield2d;
extern double dt;
extern double tidalinterval;
extern int Nout;
extern double ftime;
extern int StepNumber;
extern signed char *Status3d;
// float *HS = new float[*meshunion->Nv];
// float *DIR = new float[*meshunion->Nv];
// float *RTP = new float[*meshunion->Nv];
// float *QB = new float[*meshunion->Nv];
// float *Hm = new float[*meshunion->Nv];
// float *Um= new float[*meshunion->Nv];
// float *Vm = new float[*meshunion->Nv];

NdgPhysMat::NdgPhysMat(): 
startTime(0),
finalTime(ftime),
outputIntervalNum(10),
//tidalinterval(15.0),
//abstractoutputfile("result0001.nc", ftime / NOut, NOut),
ndgswehorizsmagrinskydiffsolver(0.25)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	if (MyID == 0) {
		getcwd(buff, 255);//将当前工作目录的绝对路径复制到参数buffer所指的内存空间中
		newTecfile();
	}
	Np = meshunion->cell_p->Np;
	K = meshunion->K;
	K2d = meshunion->mesh2d_p->K2d;
	Np2d = meshunion->mesh2d_p->mesh2dcell_p->Np2d;

	BENfp = meshunion->boundaryedge_p->Nfp;
	BENe = meshunion->boundaryedge_p->Ne;
	BENfp2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->Nfp2d;
	BENe2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->Ne2d;
	IENfp = meshunion->inneredge_p->Nfp;
	IENe = meshunion->inneredge_p->Ne;
	IENfp2d = meshunion->mesh2d_p->mesh2dinneredge_p->Nfp2d;
	IENe2d = meshunion->mesh2d_p->mesh2dinneredge_p->Ne2d;
	IEFToE3d = meshunion->inneredge_p->FToE;
	IEFToE2d = meshunion->mesh2d_p->mesh2dinneredge_p->FToE2d;
	BotENe = meshunion->bottomedge_p->Ne;
	BotENfp = meshunion->bottomedge_p->Nfp;
	BotBENe = meshunion->bottomboundaryedge_p->Ne;
	BotBENfp = meshunion->bottomboundaryedge_p->Nfp;
	SurfBENe = meshunion->surfaceboundaryedge_p->Ne;
	SurfBENfp = meshunion->surfaceboundaryedge_p->Nfp;
	Nface = meshunion->cell_p->Nface;
	Nface2d = meshunion->mesh2d_p->mesh2dcell_p->Nface2d;
	Nlayer3d = meshunion->Nlayer;

	fphys = (double *)malloc((*Np) * (*K) * Nfield3d * sizeof(double));
	fphys2d = (double *)malloc((*Np2d) * (*K2d) * Nfield2d * sizeof(double));
	fext = (double *)malloc((*BENfp) * (*BENe) * Nfield3d * sizeof(double));
	fext2d = (double *)malloc((*BENfp2d) * (*BENe2d) * Nfield2d * sizeof(double));

	memset(fphys, 0.0, (*Np) * (*K) * Nfield3d * sizeof(double));
	memset(fphys2d, 0.0, (*Np2d) * (*K2d) * Nfield2d * sizeof(double));
	memset(fext, 0.0, (*BENfp) * (*BENe) * Nfield3d * sizeof(double));
	memset(fext2d, 0.0, (*BENfp2d) * (*BENe2d) * Nfield2d * sizeof(double));

	varFieldIndex = (int *)malloc(Nvar * sizeof(int));//泥沙为2+组分，波流为2
	for (int i = 0; i < Nvar; i++)
	{
		varFieldIndex[i] = i + 1;
	}
//With T,S only
#ifdef _BAROCLINIC
	varFieldIndex[0] = 1;
	varFieldIndex[1] = 2;
	varFieldIndex[2] = 15;
	varFieldIndex[3] = 16;
//With T,S and sediment
#ifdef SEDIMENT_MODEL
	for (int i = 4; i < Nvar; i++)
	{
		varFieldIndex[i] = i + 13;
	}
#endif
#endif
//With sediment only
#ifndef _BAROCLINIC
#ifdef SEDIMENT_MODEL
	varFieldIndex[0] = 1;
	varFieldIndex[1] = 2;
	for (int i = 2; i < Nvar; i++)
	{
		varFieldIndex[i] = i + 15;
	}
#endif
#endif

/*read the fphys and fphys2d from ncfile*/

	//std::cout << "Read init_fphys.nc." << endl;

	NcFile dataFile("init_fphys.nc", NcFile::ReadOnly);
	
	NcVar *fphys_v = dataFile.get_var("fphys");
	fphys_v->get(fphys, (*Np)*(*K)*Nfield3d);
	NcVar *fphys2d_v = dataFile.get_var("fphys2d");
	fphys2d_v->get(fphys2d, (*Np2d)*(*K2d)*Nfield2d);

	//std::cout << "End reading init_fphys.nc." << endl;
/*************end reading***************/
	/////////////////////////////////////////////////找到BE所对应的节点编号，并将z值赋给fext的第四维
	double *ind = (double *)malloc((*BENfp2d) * (*BENe2d) * sizeof(double));
	double *BEFToE2d = meshunion->mesh2d_p->mesh2dboundaryedge_p->FToE2d;
	double *BEFToN12d = meshunion->mesh2d_p->mesh2dboundaryedge_p->FToN12d;
	double *fext_4 = fext2d + 3 * (*BENe2d)*(*BENfp2d);
	double *fphys_4 = fphys2d + 3 * (*Np2d)*(*K2d);
	for (int i = 0; i < *BENe2d; i++) {
		for (int j = 0; j < *BENfp2d; j++) {
			ind[i * (*BENfp2d) + j] = (BEFToE2d[*BENfp2d * i + j] - 1) * (*Np2d) + (BEFToN12d[*BENfp2d * i + j] - 1);
			fext_4[i * (*BENfp2d) + j] = fphys_4[(int)ind[i * (*BENfp2d) + j]];
		}
	}
	free(ind); ind = NULL;
	////////////////////////////////////////////////

	double *ftype = meshunion->mesh2d_p->mesh2dboundaryedge_p->ftype2d;
	for (int i = 0; i < *BENe2d; i++)
	{
		NdgEdgeType type = (NdgEdgeType)ftype[i];
		if (ftype[i] == NdgEdgeClampedDepth)
		{
			obeindex.push_back(i);
		}
		if (ftype[i] == NdgEdgeClampedVel)
		{
#ifdef _BAROCLINIC
			for (int j = 0; j < *BENfp2d; j++) {
				//For plume2023
				fext2d[i * (*BENfp2d) + j] = 0.5;//定值hu
				fext2d[(*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = 0.0;//定值hv
				//fext2d[4 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = 75.0;//定值hT
				fext2d[4 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = 40.0;
				fext2d[5 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = 0.0;//定值hS
		}
#endif
	}
		if (ftype[i] == NdgEdgeNonLinearFlatherFlow)
		{
#ifndef _BAROCLINIC
			for (int j = 0; j < *BENfp2d; j++) {
				fext2d[i * (*BENfp2d) + j] = 0.0;//定值hu
				fext2d[(*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = 0.0;//定值hv
				fext2d[2 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = -fext2d[3 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j];//定值h
			}
#endif
#ifdef _BAROCLINIC
			for (int j = 0; j < *BENfp2d; j++) {
				//For plume2023
				fext2d[i * (*BENfp2d) + j] = 0.0;//定值hu
				fext2d[(*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = 0.0;//定值hv
				fext2d[2 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = -fext2d[3 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j];//定值h
				//fext2d[4 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = fext2d[2 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] * 5.0;//定值hT
				fext2d[4 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = fext2d[2 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] * 20.0;//定值hT
				fext2d[5 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] = fext2d[2 * (*BENfp2d) * (*BENe2d) + i * (*BENfp2d) + j] * 32.0;//定值hS
			}
#endif
		}
	}

	std::ifstream data("TideElevation.txt");//read tidal data
	if (!data.is_open())
	{
		std::cout << "Error File Path !!!" << std::endl;
		// system("pause");
	}
	double point_tidal;
	while (data >> point_tidal)
		tidal.push_back(point_tidal);

	data.close();
	

	//std::cout<<"num_of_BE_OBC :"<<obeindex.size()<<endl;

	ifstream fort14("fort.14");
	double d;
	int d1;
	int size;

	fort14 >> d;
	
	int fort_Ne = (int)d; //网格
	fort14 >> d;
	
	int fort_Nv = (int)d;  //fort14点

	for (size_t i = 0; i < fort_Nv * 4; i++)
	{
		fort14 >> d;
		
	}

	for (size_t i = 0; i < fort_Ne; i++)
	{
		for (size_t j = 0; j < 2; j++)
		{
			fort14 >> d1;
		}

		for (size_t k = 0; k < 3; k++)
		{
			//while (data >> d)
			fort14 >> d1;
			d1--;
                                                //cout<<d1<<endl;
			DG_Swan_Node.push_back(d1);//������ѹ���ջ��//

		}

	}

	fort14.close();
	


	//std::cout << fort_Ne << endl << fort_Nv << endl;
	//int i = 0;
	//cout<<DG_Swan_Node.size()<<endl;
	for (size_t i = 0; i < fort_Nv; i++)
	{
		size = 0;
		for (int j = 0; j < fort_Ne * 3; j++)
		{
			if (DG_Swan_Node[j] == i) {
				Swan_DG_Node.push_back(j);
				size++;
			}
		}
		//cout << "size of size: " << size << endl;
		sizeof_PerNode.push_back(size);
	}



	vector<int>::iterator it;
	int sum = 0;
	int i = 0;
	for (it = sizeof_PerNode.begin(); it != sizeof_PerNode.end(); it++)
	{
		// cout << "Swan_DG_Node[" << i << "]=" /*<< setprecision(16)*/ << sizeof_PerNode[i] << endl;
		sum = sum + sizeof_PerNode[i];

		i++;
	}
	//std::cout << "sum: " << sum << endl;

}


NdgPhysMat::~NdgPhysMat()
{
	//free(fphys); fphys = NULL;
	//free(fphys2d); fphys2d = NULL;
	//free(fext); fext = NULL;
	//free(fext2d); fext2d = NULL;
	//free(frhs); frhs = NULL;

	//free(varFieldIndex); varFieldIndex = NULL;

}


void NdgPhysMat::matSolver(int*pE2d, int*pE3d, int*pV)
{
	matEvaluateIMEXRK222(pE2d, pE3d, pV);
}


void NdgPhysMat::matEvaluateIMEXRK222(int*pE2d, int*pE3d, int*pV)
{
	
	bool interface_status;
	double*fphys_1 = fphys2d;//h2d
	double*fphys_2 = fphys2d + 1 * (*Np2d)*(*K2d);//hu2d
	double*fphys_3 = fphys2d + 2 * (*Np2d)*(*K2d);//hv2d
	double*fphys_4 = fphys2d + 3 * (*Np2d)*(*K2d);//z2d
	signed char *status = meshunion->mesh2d_p->status;
	int Nz = *(meshunion->cell_p->Nz);
	int _Nz = Nz + 1;
	int NLayer = *(meshunion->Nlayer);
	typedef enum {
		NdgRegionNormal = 1,
		NdgRegionRefine = 2,
		NdgRegionSponge = 3,
		NdgRegionWet = 4,
		NdgRegionDry = 5,
		NdgRegionPartialWet = 6,
		NdgRegionPartialWetFlood = 7,
		NdgRegionPartialWetDamBreak = 8
	} NdgRegionType;

	NdgRegionType type = NdgRegionPartialWet;
	//auto begintime = steady_clock::now();

	const int num = (*K)*(*Np)*Nvar;
	const int num2d = (*K2d)*(*Np2d);
	const int num_1 = (*K)*(*Np);
	int Interface = *Nlayer3d + 1;

	double time = startTime;
	double ftime = finalTime;

	double *Tempfphys2d = (double *)malloc((*Np2d)*(*K2d)*sizeof(double));
	double *Tempfphys = (double *)malloc((*Np)*(*K)* Nvar * sizeof(double));
	double *Limited_huhv2D = (double *)malloc((*Np2d)*(*K2d) * 2 * sizeof(double));//For WD reconstruction,存修正前的二维动量
	EXfrhs2d = (double *)malloc((*Np2d) * (*K2d) * 2 * sizeof(double));//save h2d in two steps
	EXfrhs = (double *)malloc((*Np) * (*K) * 2 * Nvar * sizeof(double));//save hu3d and hv3d in two steps
	IMfrhs = (double *)malloc((*Np) * (*K) * Nvar * sizeof(double));

	memset(Limited_huhv2D, 0, (*Np2d) * (*K2d) * 2 * sizeof(double));
	memset(EXfrhs2d, 0, (*Np2d) * (*K2d) * 2 * sizeof(double));
	memset(EXfrhs, 0, (*Np) * (*K) * 2 * Nvar * sizeof(double));
	memset(IMfrhs, 0, (*Np) * (*K) * Nvar * sizeof(double));
	if (MyID == 0) {
		std::cout << "Allocate Memory." << endl;
	}
	/*Allocate Memory for Advection, H_Diffusion, V_Difffusion, PCE, Source term, Verticalvelocity */
	AllocateMemory.AdvMemoryAllocation(*Np, *K, Nvar, *IENfp, *IENe, *Nface, *BENfp, *BENe, *BotENfp, *BotENe, *BotBENfp, *BotBENe, *SurfBENfp, *SurfBENe);
	AllocateMemory.HorizDiffMemoryAllocation(*Np, *K, Nvar, *Nface, *BENfp, *BENe, *IENfp, *IENe);
	AllocateMemory.PCEUpdatedMemoryAllocation(*IENfp2d, *IENe2d, *Np2d, *K2d, *Nface2d, *BENe2d, *BENfp2d, *IENfp, *IENe, *BENe, *BENfp);
	AllocateMemory.VertDiffMemoryAllocation(*Np2d, *K2d, *Nlayer3d, Nvar);
	AllocateMemory.UpdatedVertVelocitySolverMemoryAllocation(*Np2d, *K2d, *IENfp2d, *IENe2d, *Nface2d, *BENe2d, *BENfp2d, *Np, *K, *IENfp, *IENe, *BENe, *BENfp);
	AllocateMemory.GotmSolverMemoryAllocation(num2d, Interface, *Np2d, *K);

#ifdef _BAROCLINIC
	AllocateMemory.BaroclinicPartMemoryAllocation(*Np, *K, *K2d, *BotENe, *BotENfp, *IENe, *IENfp, *Nface, *BENe, *BENfp);
	AllocateMemory.BaroDensityMemoryAllocation(*Np, *K);
#endif

#ifdef SEDIMENT_MODEL
	AllocateMemory.SedimentModelMemoryAllocation(*Np, *K);
#endif

#ifdef COUPLING_SWAN
	AllocateMemory.RollerWaveRadiationMemoryAllocation(*Np, *K, *IENfp, *IENe, *BENfp, *BENe, *BotENfp, *BotENe, *Nface, *BotBENfp, *BotBENe, *SurfBENfp, *SurfBENe, *Np2d, *K2d);
	float *HS_from_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *T_from_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *DIR_from_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *QB_from_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *WLEN_from_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *UBOT_from_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *TMBOT_from_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));

	float *H_to_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *U_to_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *V_to_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));
	float *test_to_swan = (float *)malloc((*meshunion->mesh2d_p->Nv2d) * sizeof(double));

	double *HS = (double *)malloc((*Np2d)*(*K2d) * sizeof(double));
	double *T = (double *)malloc((*Np2d)*(*K2d) * sizeof(double));
	double *DIR = (double *)malloc((*Np2d)*(*K2d) * sizeof(double));
	double *QB = (double *)malloc((*Np2d)*(*K2d) * sizeof(double));
	double *WLEN = (double *)malloc((*Np2d)*(*K2d) * sizeof(double));
	double *UBOT = (double *)malloc((*Np2d)*(*K2d) * sizeof(double));
	double *TMBOT = (double *)malloc((*Np2d)*(*K2d) * sizeof(double));
#endif
	/*************  end Allocation  ***************/
	if (MyID == 0) {
		std::cout << "End allocate Memory." << endl;
	}
	//abstractoutputfile.ncFile_create(Np2d, K2d, Nvar+1);
	//if (MyID == 0) {
	//	std::cout << "Output create." << endl;
	//}
#ifdef COUPLING_SWAN
	if (MyID == 0) {
		//cout << "before register_component_coupling_configuration_wrj\n";
		register_component_coupling_configuration_wrj(HS_from_swan, T_from_swan, DIR_from_swan, QB_from_swan, WLEN_from_swan, UBOT_from_swan, TMBOT_from_swan, \
			H_to_swan, U_to_swan, V_to_swan, test_to_swan);
	}
	double *HS3d = (double *)malloc((*Np) * (*K) * sizeof(double));
	double *T3d = (double *)malloc((*Np) * (*K) * sizeof(double));
	double *DIR3d = (double *)malloc((*Np) * (*K) * sizeof(double));
	double *QB3d = (double *)malloc((*Np) * (*K) * sizeof(double));
	double *WLEN3d = (double *)malloc((*Np) * (*K) * sizeof(double));
	double *UBOT3d = (double *)malloc((*Np) * (*K) * sizeof(double));
	double *TMBOT3d = (double *)malloc((*Np) * (*K) * sizeof(double));
	if (MyID == 0) {
		cout << "after register_component_coupling_configuration_wrj\n";
	}
#endif
	int control=0;
	int TEC_out_i = 0;
	//clock_t start,end;
	double dg_time=0;
	double timeRatio;
	double tloc;
	    /*******************************************************************************************************************/
	/*Get IMEXRK222 Parameter from the paper: 《Solving Unsteady Convection - Diffusion Problems in One and
		More Dimensions with Local Discontinuous Galerkin Methods
		and Implicit - Explicit Runge - Kutta Time Stepping》(2016)*/
	//Change to matEvaluateSSPRK22
	double GAMA = (2.0 - sqrt(2)) / 2.0;
	double rkb[4] = { 1, 0.5, 0.0, 0.5 };
	double rkt[2] = { 0, 1.0 };

#ifdef _BAROCLINIC
	CalculateDensityField(fphys, pE3d, MyID);//Init Density for baroclinic model
#endif

		/*******************************************************************************************************************/
	if (MyID == 0) {
		std::cout << "Enter the time discretization." << endl;
	}
	cout << "MPIrank = " << MyID << "  " << endl;

	while (time < ftime)
	{
		if (MyID == 0) {
			std::cout << "Time is " << time << " s" << endl;
		}

		if (time + dt > ftime) {
			dt = ftime - time;
		}

#ifdef COUPLING_SWAN                               
		if (control % ((*coupling_freq) / (*time_step)) == 0)  //耦合的时候的步数1/0.05
		{
			// Gather processes to 0
			GatherToZero(fphys_1, pE2d, *Np2d, *K2d, 1);//h2d
			GatherToZero(fphys_2, pE2d, *Np2d, *K2d, 1);//hu2d_Euler
			GatherToZero(fphys_3, pE2d, *Np2d, *K2d, 1);//hv2d_Euler
			GatherToZero(fphys_4, pE2d, *Np2d, *K2d, 1);//z2d

			// For process 0, calculate the vertex-averaged h2d, hu2d, hv2d and z2d first and then coupler to SWAN
			if (MyID == 0) {
				int location = 0;
				int flag = 0;
				//cout<<"DG_H ***********************************"<<endl;
				for (size_t i = 0; i < (*meshunion->mesh2d_p->Nv2d); i++)
				{
					H_to_swan[i] = 0;
					U_to_swan[i] = 0;
					V_to_swan[i] = 0;
					test_to_swan[i] = 0;
					flag = 0;
					for (size_t j = 0; j < sizeof_PerNode[i]; j++)    //这里稳定性有问题
					{
						flag = sizeof_PerNode[i];
						if (fphys_1[Swan_DG_Node[location + j]] > 0.0) {
							H_to_swan[i] = H_to_swan[i] + (fphys_1[Swan_DG_Node[location + j]] + fphys_4[Swan_DG_Node[location + j]]);
							//H_to_swan[i] = H_to_swan[i] + fphys_6[Swan_DG_Node[location + j]];
							//cout<<"Swan_DG_Node  "<<i<<":"<<(fphys_1[Swan_DG_Node[location+j]]+fphys_4[Swan_DG_Node[location+j]])<<endl;
							U_to_swan[i] = U_to_swan[i] + (fphys_2[Swan_DG_Node[location + j]] / fphys_1[Swan_DG_Node[location + j]]);
							//cout<<"U_dg  "<<U_to_swan[i]<<endl;
							V_to_swan[i] = V_to_swan[i] + (fphys_3[Swan_DG_Node[location + j]] / fphys_1[Swan_DG_Node[location + j]]);
							//if(control==0) test_to_swan[i]=test_to_swan[i]+fphys_1[Swan_DG_Node[location+j]]; //初始水深
							test_to_swan[i] = test_to_swan[i] + fphys_1[Swan_DG_Node[location + j]];
						}
						else {
							flag = flag - 1;
						}

					}

					if (flag > 0) {
						//H_to_swan[i] = H_to_swan[i] / (float)sizeof_PerNode[i];
						H_to_swan[i] = H_to_swan[i] / (float)flag;
						U_to_swan[i] = U_to_swan[i] / (float)flag;
						V_to_swan[i] = V_to_swan[i] / (float)flag;
						test_to_swan[i] = test_to_swan[i] / (float)flag;
					}
					else {
						U_to_swan[i] = 0.0;
						V_to_swan[i] = 0.0;
						test_to_swan[i] = 0.0;
					}

					if (test_to_swan[i] <= 0.0) {
						U_to_swan[i] = 0.0;
						V_to_swan[i] = 0.0;
					}

					//U_to_swan[i]=100;
					//V_to_swan[i]=100;
					//cout<<H_to_swan[i]<<endl;
					location = location + sizeof_PerNode[i];
				}

				//cout<<"DG_H ***********************************"<<endl;
				interface_status = execute_interface_using_name_wrj(dg_demo_comp_id, "send_data_to_swan", false, "execute interface for sending data to swan");
				//cout<<interface_status<<endl;
				interface_status = execute_interface_using_name_wrj(dg_demo_comp_id, "receive_data_from_swan", false, "execute interface for receiving data from swan");

				if (control < (First_CouplingStep + Calculating_CouplingSteps) * ((*coupling_freq) / (*time_step))) { //Couple when the steps are in need

					//保证SWAN传过来的波浪要素是有效的
					for (size_t i = 0; i < *K2d; i++)
					{
						for (size_t j = 0; j < *Np2d; j++)
						{
							if (HS_from_swan[DG_Swan_Node[i*(*Np2d) + j]] > 0.0 && HS_from_swan[DG_Swan_Node[i*(*Np2d) + j]] <= 99.0) {
								HS[i*(*Np2d) + j] = HS_from_swan[DG_Swan_Node[i*(*Np2d) + j]];
							}
							else {
								HS[i*(*Np2d) + j] = 0.0;
							}

							if (T_from_swan[DG_Swan_Node[i*(*Np2d) + j]] > 0.0 && T_from_swan[DG_Swan_Node[i*(*Np2d) + j]] <= 99.0) {
								T[i*(*Np2d) + j] = T_from_swan[DG_Swan_Node[i*(*Np2d) + j]];
							}
							else {
								T[i*(*Np2d) + j] = 0.0;
							}

							if (DIR_from_swan[DG_Swan_Node[i*(*Np2d) + j]] > 0.0 && DIR_from_swan[DG_Swan_Node[i*(*Np2d) + j]] <= 360.0) {
								DIR[i*(*Np2d) + j] = DIR_from_swan[DG_Swan_Node[i*(*Np2d) + j]];
							}
							else {
								DIR[i*(*Np2d) + j] = 0.0;
							}

							if (QB_from_swan[DG_Swan_Node[i*(*Np2d) + j]] > 0.0 && QB_from_swan[DG_Swan_Node[i*(*Np2d) + j]] <= 1.0) {
								QB[i*(*Np2d) + j] = QB_from_swan[DG_Swan_Node[i*(*Np2d) + j]];
							}
							else {
								QB[i*(*Np2d) + j] = 0.0;
							}

							if (WLEN_from_swan[DG_Swan_Node[i * (*Np2d) + j]] > 0.0 && WLEN_from_swan[DG_Swan_Node[i*(*Np2d) + j]] <= 999.0) {
								WLEN[i * (*Np2d) + j] = WLEN_from_swan[DG_Swan_Node[i * (*Np2d) + j]];
							}
							else {
								WLEN[i * (*Np2d) + j] = 0.0;
							}

							if (UBOT_from_swan[DG_Swan_Node[i * (*Np2d) + j]] > 0.0 && UBOT_from_swan[DG_Swan_Node[i*(*Np2d) + j]] <= 99.0) {
								UBOT[i * (*Np2d) + j] = UBOT_from_swan[DG_Swan_Node[i * (*Np2d) + j]];
							}
							else {
								UBOT[i * (*Np2d) + j] = 0.0;
							}

							if (TMBOT_from_swan[DG_Swan_Node[i*(*Np2d) + j]] > 0.0 && TMBOT_from_swan[DG_Swan_Node[i*(*Np2d) + j]] <= 99.0) {
								TMBOT[i*(*Np2d) + j] = TMBOT_from_swan[DG_Swan_Node[i*(*Np2d) + j]];
							}
							else {
								TMBOT[i*(*Np2d) + j] = 0.0;
							}
						}
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);

			// Boastcast 0 processe to all DG processes
			BroadcastToAll(HS, pE2d, *Np2d, *K2d, 1);
			BroadcastToAll(T, pE2d, *Np2d, *K2d, 1);
			BroadcastToAll(DIR, pE2d, *Np2d, *K2d, 1);
			BroadcastToAll(QB, pE2d, *Np2d, *K2d, 1);
			BroadcastToAll(WLEN, pE2d, *Np2d, *K2d, 1);
			BroadcastToAll(UBOT, pE2d, *Np2d, *K2d, 1);
			BroadcastToAll(TMBOT, pE2d, *Np2d, *K2d, 1);
			MPI_Barrier(MPI_COMM_WORLD);

			if (control < (First_CouplingStep + Calculating_CouplingSteps) * ((*coupling_freq) / (*time_step))) { //Couple when the steps are in need
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
				for (int i = 0; i < *K2d; i++) {
					if (MyID == pE2d[i]) {
						NdgExtend2dField((double*)HS3d, (double*)HS, *Np2d, i, *Np, NLayer, Nz);
						NdgExtend2dField((double*)T3d, (double*)T, *Np2d, i, *Np, NLayer, Nz);
						NdgExtend2dField((double*)DIR3d, (double*)DIR, *Np2d, i, *Np, NLayer, Nz);
						NdgExtend2dField((double*)QB3d, (double*)QB, *Np2d, i, *Np, NLayer, Nz);
						NdgExtend2dField((double*)WLEN3d, (double*)WLEN, *Np2d, i, *Np, NLayer, Nz);
						NdgExtend2dField((double*)UBOT3d, (double*)UBOT, *Np2d, i, *Np, NLayer, Nz);
						NdgExtend2dField((double*)TMBOT3d, (double*)TMBOT, *Np2d, i, *Np, NLayer, Nz);
					}
				}
				// 再限制一下保证湿的地方才传递数据,很关键!
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
				for (int k = 0; k < *K; k++) {
					if (MyID == pE3d[k]) {
						for (int n = 0; n < *Np; n++) {
							if (fphys[3 * (*K) * (*Np) + k * (*Np) + n] <= Hcrit) {
								HS3d[k * (*Np) + n] = 0.0;
							}
						}
					}
				}

				if (control < First_CouplingStep * ((*coupling_freq) / (*time_step))) { //Couple when the wave field is steady
					memset(HS3d, 0, (*Np) * (*K) * sizeof(double));
					memset(T3d, 0, (*Np) * (*K) * sizeof(double));
					memset(DIR3d, 0, (*Np) * (*K) * sizeof(double));
					memset(QB3d, 0, (*Np) * (*K) * sizeof(double));
					memset(WLEN3d, 0, (*Np) * (*K) * sizeof(double));
					memset(UBOT3d, 0, (*Np) * (*K) * sizeof(double));
					memset(TMBOT3d, 0, (*Np) * (*K) * sizeof(double));
				}
			}
		}
#endif 

		MPI_Barrier(MPI_COMM_WORLD);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < *K2d; k++) {
			if (MyID == pE2d[k]) {
				for (int n = 0; n < *Np2d; n++) {
					Tempfphys2d[k * (*Np2d) + n] = fphys2d[k * (*Np2d) + n];
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < *K; k++) {
			if (MyID == pE3d[k]) {
				for (int i = 0; i < Nvar; i++) {
					for (int n = 0; n < *Np; n++) {
						Tempfphys[i * (*K) * (*Np) + k * (*Np) + n] = fphys[(varFieldIndex[i] - 1) * (*K) * (*Np) + k * (*Np) + n];
					}
				}
			}
		}


/*****************************************************************  Start RK time step  *****************************************************************/
		//OUTPUT TECPLOT
		if (TEC_out_i % StepNumber == 0) {
			addTecdata(fphys, fphys2d, sizeof_PerNode, Swan_DG_Node, (int)(TEC_out_i*dt), pE2d, pE3d, pV);//update in every StepNumber*dt 
		}
		TEC_out_i++;

/*****  Step  1 *****/
		int intRK = 0;
		tloc = time + rkt[intRK] * dt;

		UpdateExternalField(tloc, fphys2d, fphys);//有潮位记得打开

		/********** Here the wet and dry status are updated. ***********/
		if (time == 0.0) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
			for (int k = 0; k < *K2d; k++) {
				if (MyID == pE2d[k]) {
					for (int n = 0; n < *Np2d; n++) {
						Limited_huhv2D[k * (*Np2d) + n] = fphys2d[*K2d*(*Np2d) + k * (*Np2d) + n];//hu2d
						Limited_huhv2D[*K2d * (*Np2d) + k * (*Np2d) + n] = fphys2d[*K2d*(*Np2d) * 2 + k * (*Np2d) + n];//hv2d
					}
				}
			}
			UpdateWetDryState(fphys, fphys2d, Limited_huhv2D,Nlayer3d, status, *Np, *K, *Np2d, *K2d, pE3d, pE2d, MyID);//第一步一般huhv2D=0
		}
		/****** ********************** END ********************** ******/
#ifdef _BAROCLINIC
		CalculateDensityField(fphys, pE3d, MyID);//Init/Update Density for baroclinic model
#endif
		Exchange(fphys, IEFToE3d, pE3d, *IENe, *Np, *K, Nfield3d);// use MPI to exchange 3d message
		Exchange(fphys2d, IEFToE2d, pE2d, *IENe2d, *Np2d, *K2d, 3);// use MPI to exchange 2d message

#ifndef COUPLING_SWAN
		EvaluateRHS_Nowave(fphys, EXfrhs, time, fext, varFieldIndex, fphys2d, fext2d, EXfrhs2d, pE2d, pE3d, MyID);
#endif

#ifdef COUPLING_SWAN
		EvaluateRHS(fphys, EXfrhs, time, fext, varFieldIndex, fphys2d, fext2d, EXfrhs2d, (double*)HS3d, (double*)T3d,\
			(double*)DIR3d, (double*)QB3d, (double*)WLEN3d, (double*)UBOT3d, (double*)TMBOT3d, pE2d, pE3d, MyID);
#endif

		//Update h2d,hu3d,hv3d field
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < (*K2d); k++) {
			if (MyID == pE2d[k]) {
				for (int n = 0; n < (*Np2d); n++) {
					fphys2d[k*(*Np2d) + n] = Tempfphys2d[k*(*Np2d) + n] + dt * EXfrhs2d[k*(*Np2d) + n];
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < (*K); k++) {
			if (MyID == pE3d[k]) {
				for (int i = 0; i < Nvar; i++) {
					for (int n = 0; n < (*Np); n++) {
						fphys[(varFieldIndex[i] - 1)*(*K)*(*Np) + k * (*Np) + n] = Tempfphys[i*(*K)*(*Np) + k * (*Np) + n] + dt * EXfrhs[i*(*K)*(*Np) + k * (*Np) + n];
					}
				}
			}
		}

		Exchange(fphys, IEFToE3d, pE3d, *IENe, *Np, *K, 2);// use MPI to exchange 3d message

		if ((int)Switch_Limiter3D == 1) {
			Exchange(fphys2d, IEFToE2d, pE2d, *IENe2d, *Np2d, *K2d, 1);// use MPI to exchange 2d message
			Limiter2d(fphys2d, pE2d, MyID);//h2d
			Limiter3d(fphys, pE3d, MyID);//hu
			Limiter3d(fphys + (*K)*(*Np), pE3d, MyID);//hv
		}

		//Update h3d field
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int i = 0; i < *K2d; i++) {
			if (MyID == pE2d[i]) {
				NdgExtend2dField(fphys + (*Np)*(*K) * 3, fphys2d, *Np2d, i, *Np, NLayer, Nz);
			}
		}

		//Update hu2d,hv2d field
		verticalColumnIntegralField.EvaluateVerticalIntegral(fphys2d + (*Np2d) * (*K2d), fphys, pE2d,MyID);
		verticalColumnIntegralField.EvaluateVerticalIntegral(fphys2d + (*Np2d) * (*K2d) * 2, fphys + (*Np) * (*K), pE2d, MyID);



		/********** Here the wet and dry status are updated. ***********/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < *K2d; k++) {
			if (MyID == pE2d[k]) {
				for (int n = 0; n < *Np2d; n++) {
					Limited_huhv2D[k * (*Np2d) + n] = fphys2d[*K2d*(*Np2d) + k * (*Np2d) + n];//hu2d
					Limited_huhv2D[*K2d * (*Np2d) + k * (*Np2d) + n] = fphys2d[*K2d*(*Np2d) * 2 + k * (*Np2d) + n];//hv2d
				}
			}
		}

		//UpdateWetDryState(fphys, fphys2d, Limited_huhv2D, Nlayer3d, status, *Np, *K, *Np2d, *K2d, pE2d, MyID);
		/****** ********************** END ********************** ******/


		//Update zeta field (7 = 4 + 6) if wet, and 7 = 0.0 if dry
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < (*K); k++) {
			if (MyID == pE3d[k]) {
				if ((NdgRegionType)Status3d[k] != NdgRegionDry) {
					for (int n = 0; n < (*Np); n++) {
						fphys[(*K)*(*Np) * 6 + k * (*Np) + n] = fphys[(*K)*(*Np) * 3 + k * (*Np) + n] + fphys[(*K)*(*Np) * 5 + k * (*Np) + n];

#ifdef COUPLING_SWAN
						//Update Euler vilocity hu_e and hv_e
						fphys[(*K)*(*Np) * 9 + k * (*Np) + n] = fphys[k * (*Np) + n] - fphys[(*K)*(*Np) * 11 + k * (*Np) + n];
						fphys[(*K)*(*Np) * 10 + k * (*Np) + n] = fphys[(*K)*(*Np) + k * (*Np) + n] - fphys[(*K)*(*Np) * 12 + k * (*Np) + n];
#endif
					}
				}
				else {
					continue;
				}
			}
		}

		Exchange(fphys, IEFToE3d, pE3d, *IENe, *Np, *K, Nfield3d);// use MPI to exchange 3d message
		Exchange(fphys2d, IEFToE2d, pE2d, *IENe2d, *Np2d, *K2d, 3);// use MPI to exchange 2d message

		//Updata Vertical Velocity
		std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
		calculateVerticalVelocity.EvaluateVerticalVelocity(fphys2d, fphys, fext2d, fext, pE2d, pE3d, MyID);
		std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
		std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
		OMGtime = OMGtime + time_used.count();

		/*****  End Step  1  *****/


		/*****  Step 2 *****/
		intRK = 1;
		tloc = time + rkt[intRK] * dt;

		UpdateExternalField(tloc, fphys2d, fphys);//有潮位记得打开

#ifdef _BAROCLINIC
		CalculateDensityField(fphys, pE3d, MyID);
#endif

		Exchange(fphys, IEFToE3d, pE3d, *IENe, *Np, *K, Nfield3d);// use MPI to exchange 3d message
		Exchange(fphys2d, IEFToE2d, pE2d, *IENe2d, *Np2d, *K2d, 3);// use MPI to exchange 2d message

#ifndef COUPLING_SWAN
		EvaluateRHS_Nowave(fphys, EXfrhs + Nvar * (*Np) * (*K), time, fext, varFieldIndex, fphys2d, fext2d, EXfrhs2d + (*Np2d) * (*K2d), pE2d, pE3d, MyID);
#endif

#ifdef COUPLING_SWAN
		EvaluateRHS(fphys, EXfrhs + Nvar * (*Np) * (*K), time, fext, varFieldIndex, fphys2d, fext2d, EXfrhs2d + (*Np2d) * (*K2d),\
			(double*)HS3d, (double*)T3d, (double*)DIR3d, (double*)QB3d, (double*)WLEN3d, (double*)UBOT3d, (double*)TMBOT3d, pE2d, pE3d, MyID);
#endif

		//Update h2d,hu3d,hv3d field
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < (*K2d); k++) {
			if (MyID == pE2d[k]) {
				for (int n = 0; n < (*Np2d); n++) {
					fphys2d[k*(*Np2d) + n] = Tempfphys2d[k*(*Np2d) + n] + rkb[intRK] * dt * EXfrhs2d[k*(*Np2d) + n] + \
						rkb[2 * intRK + 1] * dt * EXfrhs2d[(*Np2d) * (*K2d) + k * (*Np2d) + n];
				}
			}
		}
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < (*K); k++) {
			if (MyID == pE3d[k]) {
				for (int i = 0; i < Nvar; i++) {
					for (int n = 0; n < (*Np); n++) {
						fphys[(varFieldIndex[i] - 1)*(*K)*(*Np) + k * (*Np) + n] = Tempfphys[i*(*K)*(*Np) + k * (*Np) + n] + \
							rkb[intRK] * dt * EXfrhs[i * (*K) * (*Np) + k * (*Np) + n] + \
							rkb[2 * intRK + 1] * dt * EXfrhs[Nvar * (*K)*(*Np) + i * (*K)*(*Np) + k * (*Np) + n];
					}
				}
			}
		}

		Exchange(fphys, IEFToE3d, pE3d, *IENe, *Np, *K, 2);// use MPI to exchange 3d message

		if ((int)Switch_Limiter3D == 1) {
			Exchange(fphys2d, IEFToE2d, pE2d, *IENe2d, *Np2d, *K2d, 1);// use MPI to exchange 2d message
			Limiter2d(fphys2d, pE2d, MyID);//h2d
			Limiter3d(fphys, pE3d, MyID);//hu
			Limiter3d(fphys + (*K)*(*Np), pE3d, MyID);//hv
		}

		//Update h3d field
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int i = 0; i < *K2d; i++) {
			if (MyID == pE2d[i]) {
				NdgExtend2dField(fphys + (*Np)*(*K) * 3, fphys2d, *Np2d, i, *Np, NLayer, Nz);
			}
		}

		//Update zeta field (7 = 4 + 6) if wet, and 7 = 0.0 if dry
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < (*K); k++) {
			if (MyID == pE3d[k]) {
				if ((NdgRegionType)Status3d[k] != NdgRegionDry) {
					for (int n = 0; n < (*Np); n++) {
						fphys[(*K)*(*Np) * 6 + k * (*Np) + n] = fphys[(*K)*(*Np) * 3 + k * (*Np) + n] + fphys[(*K)*(*Np) * 5 + k * (*Np) + n];

#ifdef COUPLING_SWAN
						//Update Euler vilocity hu_e and hv_e
						fphys[(*K)*(*Np) * 9 + k * (*Np) + n] = fphys[k * (*Np) + n] - fphys[(*K)*(*Np) * 11 + k * (*Np) + n];
						fphys[(*K)*(*Np) * 10 + k * (*Np) + n] = fphys[(*K)*(*Np) + k * (*Np) + n] - fphys[(*K)*(*Np) * 12 + k * (*Np) + n];
#endif
					}
				}
				else {
					continue;
				}
			}
		}

		/*****  End Step 2  *****/
#ifdef _BAROCLINIC
		CalculateDensityField(fphys, pE3d, MyID);
#endif
		//matUpdateImplicitVerticalDiffusion垂向扩散求解
		t0 = std::chrono::steady_clock::now();
#ifndef COUPLING_SWAN
		ndgswevertgotmdiffsolver.EvaluateVertDiffRHS(fphys, IMfrhs, &time, fphys2d, 1, varFieldIndex, pE2d, pE3d, MyID);
#endif
#ifdef COUPLING_SWAN
		ndgswevertgotmdiffsolver.EvaluateVertDiffRHS_CW(fphys, IMfrhs, &time, fphys2d, 1, UBOT, TMBOT, varFieldIndex, pE2d, pE3d, MyID);
#endif
		t1 = std::chrono::steady_clock::now();
		time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
		DIF_Vtime = DIF_Vtime + time_used.count();

		//Update hu2d,hv2d field
		verticalColumnIntegralField.EvaluateVerticalIntegral(fphys2d + (*Np2d) * (*K2d), fphys, pE2d, MyID);
		verticalColumnIntegralField.EvaluateVerticalIntegral(fphys2d + (*Np2d) * (*K2d) * 2, fphys + (*Np) * (*K), pE2d, MyID);

#ifdef COUPLING_SWAN
		//Update Euler hu2d,hv2d field
		verticalColumnIntegralField.EvaluateVerticalIntegral(fphys2d + (*Np2d) * (*K2d) * 6, fphys + (*Np) * (*K) * 9);
		verticalColumnIntegralField.EvaluateVerticalIntegral(fphys2d + (*Np2d) * (*K2d) * 7, fphys + (*Np) * (*K) * 10);
#endif

		/********** Here the wet and dry status are updated. ***********/
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
		for (int k = 0; k < *K2d; k++) {
			if (MyID == pE2d[k]) {
				for (int n = 0; n < *Np2d; n++) {
					Limited_huhv2D[k * (*Np2d) + n] = fphys2d[*K2d*(*Np2d) + k * (*Np2d) + n];//hu2d
					Limited_huhv2D[*K2d * (*Np2d) + k * (*Np2d) + n] = fphys2d[*K2d*(*Np2d) * 2 + k * (*Np2d) + n];//hv2d
				}
			}
		}

		UpdateWetDryState(fphys, fphys2d, Limited_huhv2D, Nlayer3d, status, *Np, *K, *Np2d, *K2d, pE3d, pE2d, MyID);
		/****** ********************** END ********************** ******/

		Exchange(fphys, IEFToE3d, pE3d, *IENe, *Np, *K, Nfield3d);// use MPI to exchange 3d message
		Exchange(fphys2d, IEFToE2d, pE2d, *IENe2d, *Np2d, *K2d, 3);// use MPI to exchange 2d message

		//Updata Vertical Velocity
		t0 = std::chrono::steady_clock::now();
		calculateVerticalVelocity.EvaluateVerticalVelocity(fphys2d, fphys, fext2d, fext, pE2d, pE3d, MyID);
		t1 = std::chrono::steady_clock::now();
		time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
		OMGtime = OMGtime + time_used.count();

		//暂存的右端项清零
		memset(EXfrhs, 0, (*Np) * (*K) * Nvar * 2 * sizeof(double));
		memset(EXfrhs2d, 0, (*Np2d) * (*K2d) * 2 * sizeof(double));
		memset(IMfrhs, 0, (*Np) * (*K) * Nvar * sizeof(double));
/*****************************************************************  END RK time step  *****************************************************************/
		time = time + dt;

		//UpdateOutputResult( time, fphys2d, fphys );后期会改，按照需要输出的结果自行调节
		//UpdateOutputResult(time, fphys2d, fphys);

		timeRatio = time / ftime;

#ifdef COUPLING_SWAN
		if (MyID == 0) {
			advance_time_wrj(dg_demo_comp_id, "dg_demo advances time for one step");
		}
		control++;
#endif

	}

	addTecdata(fphys, fphys2d, sizeof_PerNode, Swan_DG_Node, (int)(TEC_out_i*dt), pE2d, pE3d, pV);//update in every StepNumber*dt 

	//free(Tempfphys2d); Tempfphys2d = NULL;
	//free(Tempfphys); Tempfphys = NULL;
	//free(Limited_huhv2D); Limited_huhv2D = NULL;
	//free(SystemRHS); SystemRHS = NULL;
	//free(EXfrhs); EXfrhs = NULL;
	//free(EXfrhs2d); EXfrhs2d = NULL;
	//free(IMfrhs); IMfrhs = NULL;
	//UpdateFinalResult( time, fphys2d, fphys );
	//abstractoutputfile.closencfile();
	if (MyID == 0) {
		std::cout << "Time is " << time << " s" << endl;
		std::cout << "Total advection time is :" << ADVtime << " s" << endl;
		std::cout << "Total horizontal diffusion time is :" << DIF_Htime << " s" << endl;
		std::cout << "Total PCE solver time is :" << PCEtime << " s" << endl;
		std::cout << "Total source solver time is :" << SOURCEtime << " s" << endl;
		std::cout << "Total vertical diffusion time is :" << DIF_Vtime << " s" << endl;
		std::cout << "Total omega solver time is :" << OMGtime << " s" << endl;
#ifdef _BAROCLINIC
		std::cout << "Total BaroclinicTerm solver time is :" << Baroclinictime << " s" << endl;
#endif
	}

	/*Clear Allocate Memory for Advection, H_Diffusion, V_Difffusion, PCE, Source term, Verticalvelocity */
	AllocateMemory.AdvMemoryDeAllocation();
	AllocateMemory.HorizDiffMemoryDeAllocation();
	AllocateMemory.GotmSolverMemoryDeAllocation();
	AllocateMemory.VertDiffMemoryDeAllocation();
	AllocateMemory.PCEUpdatedMemoryDeAllocation();
	AllocateMemory.UpdatedVertVelocitySolverMemoryDeAllocation();

#ifdef _BAROCLINIC
	AllocateMemory.BaroDensityMemoryDeAllocation();
	AllocateMemory.BaroclinicPartMemoryDeAllocation();
#endif

#ifdef SEDIMENT_MODEL
	AllocateMemory.SedimentModelMemoryDeAllocation();
#endif

#ifdef COUPLING_SWAN
	AllocateMemory.RollerWaveRadiationMemoryDeAllocation();
	if (MyID == 0) {
		finalize_wrj(true, "dg_demo finalizes C-Coupler2");
	}
#endif
	/*************  end DeAllocation  ***************/

#ifdef COUPLING_SWAN
	free(HS_from_swan); HS_from_swan = NULL;
	free(T_from_swan); T_from_swan = NULL;
	free(DIR_from_swan); DIR_from_swan = NULL;
	free(QB_from_swan); QB_from_swan = NULL;
	free(WLEN_from_swan); WLEN_from_swan = NULL;
	free(UBOT_from_swan); UBOT_from_swan = NULL;
	free(TMBOT_from_swan); TMBOT_from_swan = NULL;

	free(H_to_swan); H_to_swan = NULL;
	free(U_to_swan); U_to_swan = NULL;
	free(V_to_swan); V_to_swan = NULL;
	free(test_to_swan); test_to_swan = NULL;

	free(HS); HS = NULL;
	free(T); T = NULL;
	free(DIR); DIR = NULL;
	free(QB); QB = NULL;
	free(WLEN); WLEN = NULL;
	free(UBOT); UBOT = NULL;
	free(TMBOT); TMBOT = NULL;
#endif

#ifdef COUPLING_SWAN
	free(HS3d); HS3d = NULL;
	free(T3d); T3d = NULL;
	free(DIR3d); DIR3d = NULL;
	free(QB3d); QB3d = NULL;
	free(WLEN3d); WLEN3d = NULL;
	free(UBOT3d); UBOT3d = NULL;
	free(TMBOT3d); TMBOT3d = NULL;

#endif
} 

#ifndef COUPLING_SWAN
void NdgPhysMat::EvaluateRHS_Nowave(double *fphys, double *frhs, double time, double *fext, int *varFieldIndex, double *fphys2d, double *fext2d, double *frhs2d, int* pE2d, int*pE3d, int MyID)
{
	//Advection
	//startclock1 = clock();
	std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
	ndgquadfreestrongformadvsolver3d.evaluateAdvectionRHS(fphys, frhs, fext, varFieldIndex, pE3d, MyID);
	std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	ADVtime = ADVtime + time_used.count();
	//endclock1 = clock();
	//ADVtime = ADVtime + (double)(endclock1 - startclock1) / CLOCKS_PER_SEC;

	//Diffusion
	t0 = std::chrono::steady_clock::now();
	ndgswehorizsmagrinskydiffsolver.EvaluateDiffRHS_Nowave(fphys, frhs, fext, varFieldIndex, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	DIF_Htime = DIF_Htime + time_used.count();

	//PCE
	t0 = std::chrono::steady_clock::now();
	ndgquadfreestrongformPECsolver2d.evaluatePCERHSUpdated(fphys, frhs2d, fext, varFieldIndex,fphys2d, fext2d, pE2d, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	PCEtime = PCEtime + time_used.count();

	//Source Term
	t0 = std::chrono::steady_clock::now();
	ndgsourcetermsolver3d.EvaluateSourceTerm(fphys, frhs, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	SOURCEtime = SOURCEtime + time_used.count();

	//BaroclinicTerm
#ifdef _BAROCLINIC
	t0 = std::chrono::steady_clock::now();
	EvaluateBaroclinicTerm(fphys, frhs, fext, pE2d, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	Baroclinictime = Baroclinictime + time_used.count();
#endif

};
#endif

#ifdef COUPLING_SWAN
void NdgPhysMat::EvaluateRHS(double *fphys, double *frhs, double time, double *fext, int *varFieldIndex, double *fphys2d, double *fext2d, double *frhs2d,\
	double *HS3d, double *T3d, double *DIR3d, double *QB3d, double *WLEN3d, double *UBOT3d, double *TMBOT3d, int* pE2d, int*pE3d, int MyID)
{
	//Advection
	std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
	ndgquadfreestrongformadvsolver3d.evaluateAdvectionRHS(fphys, frhs, fext, varFieldIndex, pE3d, MyID);
	std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	ADVtime = ADVtime + time_used.count();
	//Diffusion
	t0 = std::chrono::steady_clock::now();
	ndgswehorizsmagrinskydiffsolver.EvaluateDiffRHS(fphys, frhs, fext, varFieldIndex, &time, (double*)HS3d, (double*)WLEN3d, (double*)UBOT3d, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	DIF_Htime = DIF_Htime + time_used.count();
	//PCE
	t0 = std::chrono::steady_clock::now();
	ndgquadfreestrongformPECsolver2d.evaluatePCERHSUpdated(fphys, frhs2d, fext, varFieldIndex, fphys2d, fext2d, pE2d, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	PCEtime = PCEtime + time_used.count();
	//Source Term
	t0 = std::chrono::steady_clock::now();
	ndgsourcetermsolver3d.EvaluateSourceTerm(fphys, frhs, &time, (double*)HS3d, (double*)T3d, (double*)DIR3d, (double*)QB3d, (double*)WLEN3d, pE2d, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	SOURCEtime = SOURCEtime + time_used.count();

	//BaroclinicTerm
#ifdef _BAROCLINIC
	t0 = std::chrono::steady_clock::now();
	EvaluateBaroclinicTerm(fphys, frhs, fext, pE2d, pE3d, MyID);
	t1 = std::chrono::steady_clock::now();
	time_used = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
	Baroclinictime = Baroclinictime + time_used.count();
#endif

};
#endif


//void NdgPhysMat::UpdateOutputResult(double time, double *fphys) {};
void NdgPhysMat::UpdateExternalField(double tloc, double *fphys2d, double *fphys)
{
	//const int benfp2d = *meshunion->mesh2d_p->mesh2dboundaryedge_p->Nfp2d;
	//const int bene2d = *meshunion->mesh2d_p->mesh2dboundaryedge_p->Ne2d;
	//const int benfp3d = *meshunion->boundaryedge_p->Nfp;
	//const int bene3d = *meshunion->boundaryedge_p->Ne;

	int benfp2d = *meshunion->mesh2d_p->mesh2dboundaryedge_p->Nfp2d;
	int bene2d = *meshunion->mesh2d_p->mesh2dboundaryedge_p->Ne2d;
	int benfp3d = *meshunion->boundaryedge_p->Nfp;
	int bene3d = *meshunion->boundaryedge_p->Ne;

	const int obnum = obeindex.size()+1;

	const double delta = tidalinterval;

	const int s1 = (int)ceil(tloc / delta);
	const double alpha1 = (delta*s1 - tloc) / delta;
	double alpha2 = (tloc - delta * (s1 - 1)) / delta;

	//const int s1 = floor(tloc / delta) + 1;
	//const int s2 = s1 + 1;
	//double alpha1 = (delta * (s2 - 1) - tloc) / delta;
	//double alpha2 = (tloc - delta * (s1 - 1)) / delta;

	std::vector<double> fnT;

	for (int i = 0; i < obnum; i++) {
		double temp = tidal[(s1 - 1)*obnum + i] * alpha1 + tidal[s1*obnum + i] * alpha2;
		fnT.push_back(temp);
	}

	double *fext_4 = fext2d + 3 * benfp2d * bene2d;
	for (int i = 0; i < obeindex.size(); i++) {
		for (int j = 0; j < benfp2d; j++)
		{
			//将开边界水位转换成边界的水深
			fext2d[2 * benfp2d * bene2d + obeindex[i] * benfp2d + j] = fmax(fnT[i + j] - fext_4[obeindex[i] * benfp2d + j], 0);
#ifdef _BAROCLINIC
			fext2d[4 * benfp2d * bene2d + obeindex[i] * benfp2d + j] = fext2d[2 * benfp2d * bene2d + obeindex[i] * benfp2d + j] * 20.0;
			fext2d[5 * benfp2d * bene2d + obeindex[i] * benfp2d + j] = fext2d[2 * benfp2d * bene2d + obeindex[i] * benfp2d + j] * 35.0;
#endif
		}
	}

	verticalrepmatfacialvalue.EvaluateRepmatFacialValue(fext2d + 2 * benfp2d * bene2d,fext + 2 * benfp3d * bene3d);
#ifdef _BAROCLINIC
	verticalrepmatfacialvalue.EvaluateRepmatFacialValue(fext2d, fext);//first for hu
	verticalrepmatfacialvalue.EvaluateRepmatFacialValue(fext2d + benfp2d * bene2d, fext + benfp3d * bene3d);//second for hv
	verticalrepmatfacialvalue.EvaluateRepmatFacialValue(fext2d + 4 * benfp2d * bene2d, fext + 3 * benfp3d * bene3d);//forth for hT
	verticalrepmatfacialvalue.EvaluateRepmatFacialValue(fext2d + 5 * benfp2d * bene2d, fext + 4 * benfp3d * bene3d);//fifth for hS
#endif
}


void NdgPhysMat::UpdateOutputResult(double &time, double *fphys2d, double *fphys)
{
	//abstractoutputfile.outputIntervalResult(time, fphys2d, Nvar+1, Np2d, K2d);
	//abstractoutputfile.outputIntervalResult(time, fphys, Nvar, Np, K);
};


void newTecfile()
{
	const int size = 255;
	char path1[size];
	strcpy(path1, buff);
	strcat(path1, "/Tidal_Current.tec");
//#ifdef _BAROCLINIC
//	char path3[size];
//	strcpy(path3, buff);
//	strcat(path3, "/VerticalField_RHO_T_S.tec");
//	ofstream outfile3;
//#endif
	ofstream outfile1;

	outfile1.open(path1, ios::trunc);

	if (outfile1.is_open())
	{
#ifndef _BAROCLINIC
		outfile1 << "Title = \"ProjectPostProcess \"\nVARIABLES =\"X\",\"Y\",\"Z\",\"h\",\"u\",\"v\",\"usurf\",\"u1\",\"u2\",\"u3\",\"u4\",\"u5\",\"u6\",\"u7\",\"u8\",\"u9\",\"ubot\",\"vsurf\",\"v1\",\"v2\",\"v3\",\"v4\",\"v5\",\"v6\",\"v7\",\"v8\",\"v9\",\"vbot\",\"zeta\"";
#endif
#ifdef _BAROCLINIC
		outfile1 << "Title = \"ProjectPostProcess \"\nVARIABLES =\"X\",\"Y\",\"Z\",\"h\",\"u\",\"v\",\"usurf\",\"u1\",\"u2\",\"u3\",\"u4\",\"u5\",\"u6\",\"u7\",\"u8\",\"u9\",\"ubot\",\"vsurf\",\"v1\",\"v2\",\"v3\",\"v4\",\"v5\",\"v6\",\"v7\",\"v8\",\"v9\",\"vbot\",\"Ssurf\",\"Sbot\",\"Tsurf\",\"Tbot\",\"zeta\"";
		//outfile1 << "Title = \"ProjectPostProcess \"\nVARIABLES =\"X\",\"Y\",\"Z\",\"Tsurf\",\"T1\",\"T2\",\"T3\",\"T4\",\"T5\",\"T6\",\"T7\",\"T8\",\"T9\",\"T10\",\"T11\",\"T12\",\"T13\",\"T14\",\"T15\",\"T16\",\"T17\",\"T18\",\"T19\",\"Tbot\"";
		//outfile3 << "Title = \"ProjectPostProcess \"\nVARIABLES =\"X\",\"Y\",\"Z\",\"Ssurf\",\"S1\",\"S2\",\"S3\",\"S4\",\"S5\",\"S6\",\"S7\",\"S8\",\"S9\",\"S10\",\"S11\",\"S12\",\"S13\",\"S14\",\"S15\",\"S16\",\"S17\",\"S18\",\"S19\",\"Sbot\"";
		//outfile3.close();
#endif
		outfile1.close();
	}
}

void NdgPhysMat::addTecdata(double *fphys, double*fphys2d, vector<int> sizeof_PerNode, vector<int> Swan_DG_Node, int n_,int*pE2d,int*pE3d,int*pV)
{
	const int size = 255;
	int location = 0;
	int *Np2d = meshunion->mesh2d_p->mesh2dcell_p->Np2d;
	int *K2d = meshunion->mesh2d_p->K2d;
	int *Np3d = meshunion->cell_p->Np;
	int *K3d = meshunion->K;
	int NLayer = *meshunion->Nlayer;
	int Nz = *meshunion->cell_p->Nz;
	double *x_out = meshunion->x;
	double *y_out = meshunion->y;
	MPI_Status status_MPI;
	MPI_Request request_MPI[99];

	ofstream outfile1;
//#ifdef _BAROCLINIC
//	ofstream outfile3;
//#endif
	double d_1;
	int d1_1;
	std::ifstream fort14_1("fort.14");
	int fort_Ne1;
	int fort_Nv1;

	if (MyID == 0) {
		//----------------------------------------------------------File 1
		char path1[size];
		strcpy(path1, buff);
		strcat(path1, "/Tidal_Current.tec");

		outfile1.open(path1, ios::app);

		outfile1 << "\nZONE T =\"P_";
		outfile1 << n_;
		outfile1 << "\",F=FEPOINT,ET=TRIANGLE,";
		outfile1 << "N=";

		fort14_1 >> d_1;
		fort_Ne1 = (int)d_1; //网格
		fort14_1 >> d_1;
		fort_Nv1 = (int)d_1;  //fort14点
		outfile1 << fort_Nv1;
		outfile1 << ",E=";
		outfile1 << fort_Ne1;
		outfile1 << "\n";
//#ifdef _BAROCLINIC
//	    char path3[size];
//	    strcpy(path3, buff);
//	    strcat(path3, "/VerticalField_RHO_T_S.tec");
//		
//		ofstream outfile3;
//		outfile3.open(path3, ios::app);
//		
//		outfile3 << "\nZONE T =\"P_";
//		outfile3 << n_;
//		outfile3 << "\",F=FEPOINT,ET=TRIANGLE,";
//		outfile3 << "N=";
//		
//		std::ifstream fort14_3("fort.14");
//		double d_3;
//		int d1_3;
//		fort14_3 >> d_3;
//		int fort_Ne3 = (int)d_3; //网格
//		fort14_3 >> d_3;
//		int fort_Nv3 = (int)d_3;  //fort14点
//		outfile3 << fort_Nv3;
//		outfile3 << ",E=";
//		outfile3 << fort_Ne3;
//	outfile3 << "\n";
//#endif
		//----------------------------------------------------------End
	}//MPI 0 run

	double *fphys_h = fphys2d + (*Np2d)*(*K2d) * 0;
	double *fphys_hu = fphys2d + (*Np2d)*(*K2d) * 1;
	double *fphys_hv = fphys2d + (*Np2d)*(*K2d) * 2;
	double *fphys_z = fphys2d + (*Np2d)*(*K2d) * 3;
	double *fphys_hc = fphys2d + (*Np2d)*(*K2d) * 4;
	double *fphys3d_hu = fphys + (*Np3d)*(*K3d) * 0;
	double *fphys3d_hv = fphys + (*Np3d)*(*K3d) * 1;

	double *ub_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *vb_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *us_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *vs_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u1_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u2_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u3_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u4_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u5_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u6_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u7_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u8_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u9_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v1_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v2_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v3_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v4_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v5_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v6_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v7_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v8_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v9_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *Surface_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *Surface_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *Bottom_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *Bottom_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U1_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U2_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U3_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U4_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U5_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U6_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U7_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U8_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *U9_hu = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V1_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V2_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V3_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V4_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V5_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V6_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V7_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V8_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *V9_hv = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));

#ifdef COUPLING_SWAN
	double *fphys_hu_Euler = fphys2d + (*Np2d)*(*K2d) * 6;
	double *fphys_hv_Euler = fphys2d + (*Np2d)*(*K2d) * 7;
	double *fphys3d_hu_Euler = fphys + (*Np3d)*(*K3d) * 9;
	double *fphys3d_hv_Euler = fphys + (*Np3d)*(*K3d) * 10;
#endif

#ifdef _BAROCLINIC
	//double *fphys3d_rho = fphys + (*Np3d)*(*K3d) * 13;
	double *fphys3d_hT = fphys + (*Np3d)*(*K3d) * 14;
	double *fphys3d_hS = fphys + (*Np3d)*(*K3d) * 15;

	double *Sb_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *Tb_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *Ss_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *Ts_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *Surface_hS = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *Surface_hT = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *Bottom_hS = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
	double *Bottom_hT = (double *)malloc(sizeof(double)*(*Np2d)*(*K2d));
#endif


//-----------------------------Output vertical fields-----------------------------//
//	double *sigma = meshunion->z;
//	double *sigma_To_z = (double *)malloc(sizeof(double)*(*K3d)*(*Np3d));
//	double *U3d_out = (double *)malloc(sizeof(double)*(*K3d)*(*Np3d));
//	double *V3d_out = (double *)malloc(sizeof(double)*(*K3d)*(*Np3d));
//#ifdef _BAROCLINIC
//	double *fphys3d_rho = fphys + (*Np3d)*(*K3d) * 13;
//	double *T_out = (double *)malloc(sizeof(double)*(*K3d)*(*Np3d));
//	double *S_out = (double *)malloc(sizeof(double)*(*K3d)*(*Np3d));
//#endif
//
//#ifndef COUPLING_SWAN
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//	for (int k = 0; k < *K3d; k++) {
//		if (MyID == pE3d[k]) {
//			for (int n = 0; n < *Np3d; n++) {
//				if (fphys[(*Np3d)*(*K3d) * 3 + k * (*Np3d) + n] > Hcrit) {
//					sigma_To_z[k * (*Np3d) + n] = sigma[k * (*Np3d) + n] * fphys[(*Np3d)*(*K3d) * 3 + k * (*Np3d) + n] + fphys[(*Np3d)*(*K3d) * 6 + k * (*Np3d) + n];
//					U3d_out[k * (*Np3d) + n] = fphys[(*Np3d)*(*K3d) * 0 + k * (*Np3d) + n] / fphys[(*Np3d)*(*K3d) * 3 + k * (*Np3d) + n];
//					V3d_out[k * (*Np3d) + n] = fphys[(*Np3d)*(*K3d) * 1 + k * (*Np3d) + n] / fphys[(*Np3d)*(*K3d) * 3 + k * (*Np3d) + n];
//				}
//				else {
//					U3d_out[k * (*Np3d) + n] = 0.0;
//					V3d_out[k * (*Np3d) + n] = 0.0;
//				}
//			}
//		}
//		if (pE3d[k] > 0 && pE3d[k] <MPIsize) {//如果这个单元不在0分区
//			if (MyID == pE3d[k]) {//对于当前进程，传递给0进程
//				MPI_Isend(sigma_To_z + k * (*Np3d), *Np3d, MPI_DOUBLE, 0, pE3d[k], MPI_COMM_WORLD, &request_MPI[0]);
//				MPI_Isend(U3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, 0, pE3d[k] + MPIsize, MPI_COMM_WORLD, &request_MPI[1]);
//				MPI_Isend(V3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, 0, pE3d[k] + 2 * MPIsize, MPI_COMM_WORLD, &request_MPI[2]);
//			}
//			else if (MyID == 0) {//0进程接收
//				MPI_Irecv(sigma_To_z + k * (*Np3d), *Np3d, MPI_DOUBLE, pE3d[k], pE3d[k], MPI_COMM_WORLD, &request_MPI[0]);
//				MPI_Wait(&request_MPI[0], &status_MPI);
//				MPI_Irecv(U3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, pE3d[k], pE3d[k] + MPIsize, MPI_COMM_WORLD, &request_MPI[1]);
//				MPI_Wait(&request_MPI[1], &status_MPI);
//				MPI_Irecv(V3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, pE3d[k], pE3d[k] + 2 * MPIsize, MPI_COMM_WORLD, &request_MPI[2]);
//				MPI_Wait(&request_MPI[2], &status_MPI);
//			}
//			else {
//			}
//		}
//
//	}
//#endif
//
//#ifdef COUPLING_SWAN
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(DG_THREADS)
//#endif
//	for (int k = 0; k < *K3d; k++) {
//		if (MyID == pE3d[k]) {
//			for (int n = 0; n < *Np3d; n++) {
//				sigma_To_z[k * (*Np3d) + n] = sigma[k * (*Np3d) + n] * fphys[(*Np3d)*(*K3d) * 3 + k * (*Np3d) + n] + fphys[(*Np3d)*(*K3d) * 6 + k * (*Np3d) + n];
//				U3d_out[k * (*Np3d) + n] = fphys3d_hu_Euler[k * (*Np3d) + n] / fphys[(*Np3d)*(*K3d) * 3 + k * (*Np3d) + n];
//				V3d_out[k * (*Np3d) + n] = fphys3d_hv_Euler[k * (*Np3d) + n] / fphys[(*Np3d)*(*K3d) * 3 + k * (*Np3d) + n];
//			}
//		}
//		if (pE3d[k] > 0 && pE3d[k] < MPIsize) {//如果这个单元不在0分区
//			if (MyID == pE3d[k]) {//对于当前进程，传递给0进程
//				MPI_Isend(sigma_To_z + k * (*Np3d), *Np3d, MPI_DOUBLE, 0, pE3d[k], MPI_COMM_WORLD, &request_MPI[0]);
//				MPI_Isend(U3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, 0, pE3d[k] + MPIsize, MPI_COMM_WORLD, &request_MPI[1]);
//				MPI_Isend(V3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, 0, pE3d[k] + 2 * MPIsize, MPI_COMM_WORLD, &request_MPI[2]);
//			}
//			else if (MyID == 0) {//0进程接收
//				MPI_Irecv(sigma_To_z + k * (*Np3d), *Np3d, MPI_DOUBLE, pE3d[k], pE3d[k], MPI_COMM_WORLD, &request_MPI[0]);
//				MPI_Wait(&request_MPI[0], &status_MPI);
//				MPI_Irecv(U3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, pE3d[k], pE3d[k] + MPIsize, MPI_COMM_WORLD, &request_MPI[1]);
//				MPI_Wait(&request_MPI[1], &status_MPI);
//				MPI_Irecv(V3d_out + k * (*Np3d), *Np3d, MPI_DOUBLE, pE3d[k], pE3d[k] + 2 * MPIsize, MPI_COMM_WORLD, &request_MPI[2]);
//				MPI_Wait(&request_MPI[2], &status_MPI);
//			}
//			else {
//			}
//		}
//	}
//#endif
//-----------------------------End Output vertical fields-----------------------------//

	//terrain change
	//double *fphys_zzz = fphys_zz;
	//double *h_out = NULL, *u_out = NULL, *v_out = NULL, *z_out = NULL, *hc_out = NULL;
	double *h_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *u_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *v_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *zeta_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));
	double *z_out = (double *)malloc(sizeof(double)*(*meshunion->mesh2d_p->Nv2d));

	//Here, we first calculate hu3d, hv3d, hT and hS under all threads.
	for (int k = 0; k < (*K2d); k++) {
		if (MyID == pE2d[k]) {
			for (int n = 0; n < (*Np2d); n++) {
#ifndef COUPLING_SWAN
				Surface_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + (*Np2d) * Nz + n];
				Surface_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + (*Np2d) * Nz + n];
				Bottom_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + (NLayer - 1) * (*Np3d) + n];
				Bottom_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + (NLayer - 1) * (*Np3d) + n];
				U1_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + (*Np3d) + (*Np2d) * Nz + n];
				U2_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 2 * (*Np3d) + (*Np2d) * Nz + n];
				U3_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 3 * (*Np3d) + (*Np2d) * Nz + n];
				U4_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 4 * (*Np3d) + (*Np2d) * Nz + n];
				U5_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 5 * (*Np3d) + (*Np2d) * Nz + n];
				U6_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 6 * (*Np3d) + (*Np2d) * Nz + n];
				U7_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 7 * (*Np3d) + (*Np2d) * Nz + n];
				U8_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 8 * (*Np3d) + (*Np2d) * Nz + n];
				U9_hu[k * (*Np2d) + n] = fphys3d_hu[k * (*Np3d) * NLayer + 9 * (*Np3d) + (*Np2d) * Nz + n];
				V1_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + (*Np3d) + (*Np2d) * Nz + n];
				V2_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 2 * (*Np3d) + (*Np2d) * Nz + n];
				V3_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 3 * (*Np3d) + (*Np2d) * Nz + n];
				V4_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 4 * (*Np3d) + (*Np2d) * Nz + n];
				V5_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 5 * (*Np3d) + (*Np2d) * Nz + n];
				V6_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 6 * (*Np3d) + (*Np2d) * Nz + n];
				V7_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 7 * (*Np3d) + (*Np2d) * Nz + n];
				V8_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 8 * (*Np3d) + (*Np2d) * Nz + n];
				V9_hv[k * (*Np2d) + n] = fphys3d_hv[k * (*Np3d) * NLayer + 9 * (*Np3d) + (*Np2d) * Nz + n];
#endif

#ifdef COUPLING_SWAN
				Surface_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + (*Np2d) * Nz + n];
				Surface_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + (*Np2d) * Nz + n];
				Bottom_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + (NLayer - 1) * (*Np3d) + n];
				Bottom_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + (NLayer - 1) * (*Np3d) + n];
				U1_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + (*Np3d) + (*Np2d) * Nz + n];
				U2_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 2 * (*Np3d) + (*Np2d) * Nz + n];
				U3_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 3 * (*Np3d) + (*Np2d) * Nz + n];
				U4_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 4 * (*Np3d) + (*Np2d) * Nz + n];
				U5_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 5 * (*Np3d) + (*Np2d) * Nz + n];
				U6_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 6 * (*Np3d) + (*Np2d) * Nz + n];
				U7_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 7 * (*Np3d) + (*Np2d) * Nz + n];
				U8_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 8 * (*Np3d) + (*Np2d) * Nz + n];
				U9_hu[k * (*Np2d) + n] = fphys3d_hu_Euler[k * (*Np3d) * NLayer + 9 * (*Np3d) + (*Np2d) * Nz + n];
				V1_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + (*Np3d) + (*Np2d) * Nz + n];
				V2_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 2 * (*Np3d) + (*Np2d) * Nz + n];
				V3_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 3 * (*Np3d) + (*Np2d) * Nz + n];
				V4_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 4 * (*Np3d) + (*Np2d) * Nz + n];
				V5_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 5 * (*Np3d) + (*Np2d) * Nz + n];
				V6_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 6 * (*Np3d) + (*Np2d) * Nz + n];
				V7_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 7 * (*Np3d) + (*Np2d) * Nz + n];
				V8_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 8 * (*Np3d) + (*Np2d) * Nz + n];
				V9_hv[k * (*Np2d) + n] = fphys3d_hv_Euler[k * (*Np3d) * NLayer + 9 * (*Np3d) + (*Np2d) * Nz + n];
#endif

#ifdef _BAROCLINIC
				Surface_hS[k * (*Np2d) + n] = fphys3d_hS[k * (*Np3d) * NLayer + (*Np2d) * Nz + n];
				Surface_hT[k * (*Np2d) + n] = fphys3d_hT[k * (*Np3d) * NLayer + (*Np2d) * Nz + n];
				Bottom_hS[k * (*Np2d) + n] = fphys3d_hS[k * (*Np3d) * NLayer + (NLayer - 1) * (*Np3d) + n];
				Bottom_hT[k * (*Np2d) + n] = fphys3d_hT[k * (*Np3d) * NLayer + (NLayer - 1) * (*Np3d) + n];
#endif
			}
		}
	}

	//Gather h2d,hu2d,hv2d,z to zero process	//Gather hu3d,hv3d,hT,hS at 10 layers to zero process
	for (int k = 0; k < *K2d; k++) {
		if (pE2d[k] > 0 && pE2d[k] < MPIsize) {
			if (MyID == pE2d[k]) {
				MPI_Isend(fphys_h + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, k, MPI_COMM_WORLD, &request_MPI[5]);
#ifndef COUPLING_SWAN
				MPI_Isend(fphys_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, *K2d + k, MPI_COMM_WORLD, &request_MPI[6]);
				MPI_Isend(fphys_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, *K2d * 2 + k, MPI_COMM_WORLD, &request_MPI[7]);
#endif
#ifdef COUPLING_SWAN
				MPI_Isend(fphys_hu_Euler + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, *K2d * 4 + k, MPI_COMM_WORLD, &request_MPI[6]);
				MPI_Isend(fphys_hv_Euler + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, *K2d * 5 + k, MPI_COMM_WORLD, &request_MPI[7]);
#endif
				MPI_Isend(fphys_z + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, *K2d * 3 + k, MPI_COMM_WORLD, &request_MPI[8]);
//---------------------------------------------------------------------------------------------------------------------
				MPI_Isend(Surface_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 1 + k, MPI_COMM_WORLD, &request_MPI[11]);
				MPI_Isend(Surface_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 2 + k, MPI_COMM_WORLD, &request_MPI[12]);
				MPI_Isend(Bottom_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 3 + k, MPI_COMM_WORLD, &request_MPI[13]);
				MPI_Isend(Bottom_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 4 + k, MPI_COMM_WORLD, &request_MPI[14]);
				MPI_Isend(U1_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 5 + k, MPI_COMM_WORLD, &request_MPI[15]);
				MPI_Isend(U2_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 6 + k, MPI_COMM_WORLD, &request_MPI[16]);
				MPI_Isend(U3_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 7 + k, MPI_COMM_WORLD, &request_MPI[17]);
				MPI_Isend(U4_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 8 + k, MPI_COMM_WORLD, &request_MPI[18]);
				MPI_Isend(U5_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 9 + k, MPI_COMM_WORLD, &request_MPI[19]);
				MPI_Isend(U6_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 10 + k, MPI_COMM_WORLD, &request_MPI[20]);
				MPI_Isend(U7_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 11 + k, MPI_COMM_WORLD, &request_MPI[21]);
				MPI_Isend(U8_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 12 + k, MPI_COMM_WORLD, &request_MPI[22]);
				MPI_Isend(U9_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 13 + k, MPI_COMM_WORLD, &request_MPI[23]);
				MPI_Isend(V1_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 14 + k, MPI_COMM_WORLD, &request_MPI[24]);
				MPI_Isend(V2_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 15 + k, MPI_COMM_WORLD, &request_MPI[25]);
				MPI_Isend(V3_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 16 + k, MPI_COMM_WORLD, &request_MPI[26]);
				MPI_Isend(V4_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 17 + k, MPI_COMM_WORLD, &request_MPI[27]);
				MPI_Isend(V5_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 18 + k, MPI_COMM_WORLD, &request_MPI[28]);
				MPI_Isend(V6_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 19 + k, MPI_COMM_WORLD, &request_MPI[29]);
				MPI_Isend(V7_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 20 + k, MPI_COMM_WORLD, &request_MPI[30]);
				MPI_Isend(V8_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 21 + k, MPI_COMM_WORLD, &request_MPI[31]);
				MPI_Isend(V9_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 22 + k, MPI_COMM_WORLD, &request_MPI[32]);
#ifdef _BAROCLINIC
				MPI_Isend(Surface_hS + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 23 + k, MPI_COMM_WORLD, &request_MPI[33]);
				MPI_Isend(Surface_hT + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 24 + k, MPI_COMM_WORLD, &request_MPI[34]);
				MPI_Isend(Bottom_hS + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 25 + k, MPI_COMM_WORLD, &request_MPI[35]);
				MPI_Isend(Bottom_hT + k * (*Np2d), *Np2d, MPI_DOUBLE, 0, 26 + k, MPI_COMM_WORLD, &request_MPI[36]);
#endif
			}
			else if (MyID == 0) {
				MPI_Irecv(fphys_h + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], k, MPI_COMM_WORLD, &request_MPI[5]);
				MPI_Wait(&request_MPI[5], &status_MPI);
				MPI_Irecv(fphys_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], *K2d + k, MPI_COMM_WORLD, &request_MPI[6]);
				MPI_Wait(&request_MPI[6], &status_MPI);
				MPI_Irecv(fphys_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], *K2d * 2 + k, MPI_COMM_WORLD, &request_MPI[7]);
				MPI_Wait(&request_MPI[7], &status_MPI);
				MPI_Irecv(fphys_z + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], *K2d * 3 + k, MPI_COMM_WORLD, &request_MPI[8]);
				MPI_Wait(&request_MPI[8], &status_MPI);
#ifdef COUPLING_SWAN
				MPI_Irecv(fphys_hu_Euler + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], *K2d * 4 + k, MPI_COMM_WORLD, &request_MPI[6]);
				MPI_Wait(&request_MPI[6], &status_MPI);
				MPI_Irecv(fphys_hv_Euler + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], *K2d * 5 + k, MPI_COMM_WORLD, &request_MPI[7]);
				MPI_Wait(&request_MPI[7], &status_MPI);
#endif
//---------------------------------------------------------------------------------------------------------------------
				MPI_Irecv(Surface_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 1 + k, MPI_COMM_WORLD, &request_MPI[11]);
				MPI_Wait(&request_MPI[11], &status_MPI);
				MPI_Irecv(Surface_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 2 + k, MPI_COMM_WORLD, &request_MPI[12]);
				MPI_Wait(&request_MPI[12], &status_MPI);
				MPI_Irecv(Bottom_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 3 + k, MPI_COMM_WORLD, &request_MPI[13]);
				MPI_Wait(&request_MPI[13], &status_MPI);
				MPI_Irecv(Bottom_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 4 + k, MPI_COMM_WORLD, &request_MPI[14]);
				MPI_Wait(&request_MPI[14], &status_MPI);
				MPI_Irecv(U1_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 5 + k, MPI_COMM_WORLD, &request_MPI[15]);
				MPI_Wait(&request_MPI[15], &status_MPI);
				MPI_Irecv(U2_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 6 + k, MPI_COMM_WORLD, &request_MPI[16]);
				MPI_Wait(&request_MPI[16], &status_MPI);
				MPI_Irecv(U3_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 7 + k, MPI_COMM_WORLD, &request_MPI[17]);
				MPI_Wait(&request_MPI[17], &status_MPI);
				MPI_Irecv(U4_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 8 + k, MPI_COMM_WORLD, &request_MPI[18]);
				MPI_Wait(&request_MPI[18], &status_MPI);
				MPI_Irecv(U5_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 9 + k, MPI_COMM_WORLD, &request_MPI[19]);
				MPI_Wait(&request_MPI[19], &status_MPI);
				MPI_Irecv(U6_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 10 + k, MPI_COMM_WORLD, &request_MPI[20]);
				MPI_Wait(&request_MPI[20], &status_MPI);
				MPI_Irecv(U7_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 11 + k, MPI_COMM_WORLD, &request_MPI[21]);
				MPI_Wait(&request_MPI[21], &status_MPI);
				MPI_Irecv(U8_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 12 + k, MPI_COMM_WORLD, &request_MPI[22]);
				MPI_Wait(&request_MPI[22], &status_MPI);
				MPI_Irecv(U9_hu + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 13 + k, MPI_COMM_WORLD, &request_MPI[23]);
				MPI_Wait(&request_MPI[23], &status_MPI);
				MPI_Irecv(V1_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 14 + k, MPI_COMM_WORLD, &request_MPI[24]);
				MPI_Wait(&request_MPI[24], &status_MPI);
				MPI_Irecv(V2_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 15 + k, MPI_COMM_WORLD, &request_MPI[25]);
				MPI_Wait(&request_MPI[25], &status_MPI);
				MPI_Irecv(V3_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 16 + k, MPI_COMM_WORLD, &request_MPI[26]);
				MPI_Wait(&request_MPI[26], &status_MPI);
				MPI_Irecv(V4_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 17 + k, MPI_COMM_WORLD, &request_MPI[27]);
				MPI_Wait(&request_MPI[27], &status_MPI);
				MPI_Irecv(V5_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 18 + k, MPI_COMM_WORLD, &request_MPI[28]);
				MPI_Wait(&request_MPI[28], &status_MPI);
				MPI_Irecv(V6_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 19 + k, MPI_COMM_WORLD, &request_MPI[29]);
				MPI_Wait(&request_MPI[29], &status_MPI);
				MPI_Irecv(V7_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 20 + k, MPI_COMM_WORLD, &request_MPI[30]);
				MPI_Wait(&request_MPI[30], &status_MPI);
				MPI_Irecv(V8_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 21 + k, MPI_COMM_WORLD, &request_MPI[31]);
				MPI_Wait(&request_MPI[31], &status_MPI);
				MPI_Irecv(V9_hv + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 22 + k, MPI_COMM_WORLD, &request_MPI[32]);
				MPI_Wait(&request_MPI[32], &status_MPI);
#ifdef _BAROCLINIC
				MPI_Irecv(Surface_hS + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 23 + k, MPI_COMM_WORLD, &request_MPI[33]);
				MPI_Wait(&request_MPI[33], &status_MPI);
				MPI_Irecv(Surface_hT + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 24 + k, MPI_COMM_WORLD, &request_MPI[34]);
				MPI_Wait(&request_MPI[34], &status_MPI);
				MPI_Irecv(Bottom_hS + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 25 + k, MPI_COMM_WORLD, &request_MPI[35]);
				MPI_Wait(&request_MPI[35], &status_MPI);
				MPI_Irecv(Bottom_hT + k * (*Np2d), *Np2d, MPI_DOUBLE, pE2d[k], 26 + k, MPI_COMM_WORLD, &request_MPI[36]);
				MPI_Wait(&request_MPI[36], &status_MPI);
#endif
			}
			else {
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (MyID == 0) {
		for (size_t i = 0; i < (*meshunion->mesh2d_p->Nv2d); i++)
		{
			h_out[i] = 0.0;
			u_out[i] = 0.0;
			v_out[i] = 0.0;
			zeta_out[i] = 0.0;
			z_out[i] = 0.0;
			ub_out[i] = 0.0;
			vb_out[i] = 0.0;
			us_out[i] = 0.0;
			vs_out[i] = 0.0;
			u1_out[i] = 0.0;
			u2_out[i] = 0.0;
			u3_out[i] = 0.0;
			u4_out[i] = 0.0;
			u5_out[i] = 0.0;
			u6_out[i] = 0.0;
			u7_out[i] = 0.0;
			u8_out[i] = 0.0;
			u9_out[i] = 0.0;
			v1_out[i] = 0.0;
			v2_out[i] = 0.0;
			v3_out[i] = 0.0;
			v4_out[i] = 0.0;
			v5_out[i] = 0.0;
			v6_out[i] = 0.0;
			v7_out[i] = 0.0;
			v8_out[i] = 0.0;
			v9_out[i] = 0.0;

#ifdef _BAROCLINIC
			Sb_out[i] = 0.0;
			Tb_out[i] = 0.0;
			Ss_out[i] = 0.0;
			Ts_out[i] = 0.0;
#endif
			for (size_t j = 0; j < sizeof_PerNode[i]; j++)
			{
				h_out[i] = h_out[i] + fphys_h[Swan_DG_Node[location + j]];
				z_out[i] = z_out[i] + fphys_z[Swan_DG_Node[location + j]];
				if (fphys_h[Swan_DG_Node[location + j]] > Hcrit) {
#ifndef COUPLING_SWAN
					u_out[i] = u_out[i] + (fphys_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v_out[i] = v_out[i] + (fphys_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
#endif

#ifdef COUPLING_SWAN
					u_out[i] = u_out[i] + (fphys_hu_Euler[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v_out[i] = v_out[i] + (fphys_hv_Euler[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
#endif
					zeta_out[i] = zeta_out[i] + fphys_h[Swan_DG_Node[location + j]] + fphys_z[Swan_DG_Node[location + j]];
					us_out[i] = us_out[i] + (Surface_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					vs_out[i] = vs_out[i] + (Surface_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					ub_out[i] = ub_out[i] + (Bottom_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					vb_out[i] = vb_out[i] + (Bottom_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u1_out[i] = u1_out[i] + (U1_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u2_out[i] = u2_out[i] + (U2_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u3_out[i] = u3_out[i] + (U3_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u4_out[i] = u4_out[i] + (U4_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u5_out[i] = u5_out[i] + (U5_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u6_out[i] = u6_out[i] + (U6_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u7_out[i] = u7_out[i] + (U7_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u8_out[i] = u8_out[i] + (U8_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					u9_out[i] = u9_out[i] + (U9_hu[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v1_out[i] = v1_out[i] + (V1_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v2_out[i] = v2_out[i] + (V2_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v3_out[i] = v3_out[i] + (V3_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v4_out[i] = v4_out[i] + (V4_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v5_out[i] = v5_out[i] + (V5_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v6_out[i] = v6_out[i] + (V6_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v7_out[i] = v7_out[i] + (V7_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v8_out[i] = v8_out[i] + (V8_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					v9_out[i] = v9_out[i] + (V9_hv[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));

#ifdef _BAROCLINIC
					Ss_out[i] = Ss_out[i] + (Surface_hS[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					Ts_out[i] = Ts_out[i] + (Surface_hT[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					Sb_out[i] = Sb_out[i] + (Bottom_hS[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
					Tb_out[i] = Tb_out[i] + (Bottom_hT[Swan_DG_Node[location + j]] / (fphys_h[Swan_DG_Node[location + j]]));
#endif
				}

			}

			h_out[i] = h_out[i] / (float)sizeof_PerNode[i];
			z_out[i] = z_out[i] / (float)sizeof_PerNode[i];

			if (h_out[i] > Hcrit) {

				u_out[i] = u_out[i] / (float)sizeof_PerNode[i];
				v_out[i] = v_out[i] / (float)sizeof_PerNode[i];
				zeta_out[i] = zeta_out[i] / (float)sizeof_PerNode[i];
				us_out[i] = us_out[i] / (float)sizeof_PerNode[i];
				vs_out[i] = vs_out[i] / (float)sizeof_PerNode[i];
				ub_out[i] = ub_out[i] / (float)sizeof_PerNode[i];
				vb_out[i] = vb_out[i] / (float)sizeof_PerNode[i];
				u1_out[i] = u1_out[i] / (float)sizeof_PerNode[i];
				u2_out[i] = u2_out[i] / (float)sizeof_PerNode[i];
				u3_out[i] = u3_out[i] / (float)sizeof_PerNode[i];
				u4_out[i] = u4_out[i] / (float)sizeof_PerNode[i];
				u5_out[i] = u5_out[i] / (float)sizeof_PerNode[i];
				u6_out[i] = u6_out[i] / (float)sizeof_PerNode[i];
				u7_out[i] = u7_out[i] / (float)sizeof_PerNode[i];
				u8_out[i] = u8_out[i] / (float)sizeof_PerNode[i];
				u9_out[i] = u9_out[i] / (float)sizeof_PerNode[i];
				v1_out[i] = v1_out[i] / (float)sizeof_PerNode[i];
				v2_out[i] = v2_out[i] / (float)sizeof_PerNode[i];
				v3_out[i] = v3_out[i] / (float)sizeof_PerNode[i];
				v4_out[i] = v4_out[i] / (float)sizeof_PerNode[i];
				v5_out[i] = v5_out[i] / (float)sizeof_PerNode[i];
				v6_out[i] = v6_out[i] / (float)sizeof_PerNode[i];
				v7_out[i] = v7_out[i] / (float)sizeof_PerNode[i];
				v8_out[i] = v8_out[i] / (float)sizeof_PerNode[i];
				v9_out[i] = v9_out[i] / (float)sizeof_PerNode[i];

#ifdef _BAROCLINIC
				Ss_out[i] = Ss_out[i] / (float)sizeof_PerNode[i];
				Ts_out[i] = Ts_out[i] / (float)sizeof_PerNode[i];
				Sb_out[i] = Sb_out[i] / (float)sizeof_PerNode[i];
				Tb_out[i] = Tb_out[i] / (float)sizeof_PerNode[i];
#endif
			}
			else {
				u_out[i] = 0.0;
				v_out[i] = 0.0;
				us_out[i] = 0.0;
				vs_out[i] = 0.0;
				ub_out[i] = 0.0;
				vb_out[i] = 0.0;
				zeta_out[i] = 0.0;
				u1_out[i] = 0.0;
				u2_out[i] = 0.0;
				u3_out[i] = 0.0;
				u4_out[i] = 0.0;
				u5_out[i] = 0.0;
				u6_out[i] = 0.0;
				u7_out[i] = 0.0;
				u8_out[i] = 0.0;
				u9_out[i] = 0.0;
				v1_out[i] = 0.0;
				v2_out[i] = 0.0;
				v3_out[i] = 0.0;
				v4_out[i] = 0.0;
				v5_out[i] = 0.0;
				v6_out[i] = 0.0;
				v7_out[i] = 0.0;
				v8_out[i] = 0.0;
				v9_out[i] = 0.0;

#ifdef _BAROCLINIC
				Ss_out[i] = 0.0;
				Ts_out[i] = 0.0;
				Sb_out[i] = 0.0;
				Tb_out[i] = 0.0;
#endif
			}

			location = location + sizeof_PerNode[i];

		}
//--------------------------------------------------------------------//
//------------------- file 1 -------------------//
		for (size_t i = 0; i < fort_Nv1 * 4; i++)
		{
			fort14_1 >> d_1;
			if (i % 4 == 0) continue;
			outfile1 << d_1;
			outfile1 << " ";
			if (i % 4 == 3)
			{

				outfile1 << h_out[i / 4];
				outfile1 << " ";
				outfile1 << u_out[i / 4];
				outfile1 << " ";
				outfile1 << v_out[i / 4];
				outfile1 << " ";
				outfile1 << us_out[i / 4];
				outfile1 << " ";
				outfile1 << u1_out[i / 4];
				outfile1 << " ";
				outfile1 << u2_out[i / 4];
				outfile1 << " ";
				outfile1 << u3_out[i / 4];
				outfile1 << " ";
				outfile1 << u4_out[i / 4];
				outfile1 << " ";
				outfile1 << u5_out[i / 4];
				outfile1 << " ";
				outfile1 << u6_out[i / 4];
				outfile1 << " ";
				outfile1 << u7_out[i / 4];
				outfile1 << " ";
				outfile1 << u8_out[i / 4];
				outfile1 << " ";
				outfile1 << u9_out[i / 4];
				outfile1 << " ";
				outfile1 << ub_out[i / 4];
				outfile1 << " ";
				outfile1 << vs_out[i / 4];
				outfile1 << " ";
				outfile1 << v1_out[i / 4];
				outfile1 << " ";
				outfile1 << v2_out[i / 4];
				outfile1 << " ";
				outfile1 << v3_out[i / 4];
				outfile1 << " ";
				outfile1 << v4_out[i / 4];
				outfile1 << " ";
				outfile1 << v5_out[i / 4];
				outfile1 << " ";
				outfile1 << v6_out[i / 4];
				outfile1 << " ";
				outfile1 << v7_out[i / 4];
				outfile1 << " ";
				outfile1 << v8_out[i / 4];
				outfile1 << " ";
				outfile1 << v9_out[i / 4];
				outfile1 << " ";
				outfile1 << vb_out[i / 4];
#ifdef _BAROCLINIC
				outfile1 << " ";
				outfile1 << Ss_out[i / 4];
				outfile1 << " ";
				outfile1 << Sb_out[i / 4];
				outfile1 << " ";
				outfile1 << Ts_out[i / 4];
				outfile1 << " ";
				outfile1 << Tb_out[i / 4];
#endif
				outfile1 << " ";
				outfile1 << zeta_out[i / 4];
				//outfile << " ";
				//outfile << z_out[i / 4];
				//outfile << " ";
				//outfile << hc_out[i / 4];
				outfile1 << "\n";
			}
		}

		//--------------------------------------------------------------------//
		//------------------- file 1 -------------------//
		for (size_t i = 0; i < fort_Ne1; i++)
		{
			for (size_t j = 0; j < 2; j++)
			{
				fort14_1 >> d1_1;
			}

			for (size_t k = 0; k < 3; k++)
			{
				//while (data >> d)
				fort14_1 >> d1_1;
				outfile1 << d1_1;
				outfile1 << " ";

			}
			outfile1 << "\n";
		}
		fort14_1.close();
		outfile1.close();


	}//MPI 0 run
	MPI_Barrier(MPI_COMM_WORLD);

	free(h_out), h_out = NULL;
	free(u_out), u_out = NULL;
	free(v_out), v_out = NULL;
	free(zeta_out), zeta_out = NULL;
	free(ub_out), ub_out = NULL;
	free(vb_out), vb_out = NULL;
	free(us_out), us_out = NULL;
	free(vs_out), vs_out = NULL;
	free(u1_out), u1_out = NULL;
	free(u2_out), u2_out = NULL;
	free(u3_out), u3_out = NULL;
	free(u4_out), u4_out = NULL;
	free(u5_out), u5_out = NULL;
	free(u6_out), u6_out = NULL;
	free(u7_out), u7_out = NULL;
	free(u8_out), u8_out = NULL;
	free(u9_out), u9_out = NULL;
	free(v1_out), v1_out = NULL;
	free(v2_out), v2_out = NULL;
	free(v3_out), v3_out = NULL;
	free(v4_out), v4_out = NULL;
	free(v5_out), v5_out = NULL;
	free(v6_out), v6_out = NULL;
	free(v7_out), v7_out = NULL;
	free(v8_out), v8_out = NULL;
	free(v9_out), v9_out = NULL;
	free(Surface_hu), Surface_hu = NULL;
	free(Surface_hv), Surface_hv = NULL;
	free(Bottom_hu), Bottom_hu = NULL;
	free(Bottom_hv), Bottom_hv = NULL;
	free(U1_hu), U1_hu = NULL;
	free(U2_hu), U2_hu = NULL;
	free(U3_hu), U3_hu = NULL;
	free(U4_hu), U4_hu = NULL;
	free(U5_hu), U5_hu = NULL;
	free(U6_hu), U6_hu = NULL;
	free(U7_hu), U7_hu = NULL;
	free(U8_hu), U8_hu = NULL;
	free(U9_hu), U9_hu = NULL;
	free(V1_hv), V1_hv = NULL;
	free(V2_hv), V2_hv = NULL;
	free(V3_hv), V3_hv = NULL;
	free(V4_hv), V4_hv = NULL;
	free(V5_hv), V5_hv = NULL;
	free(V6_hv), V6_hv = NULL;
	free(V7_hv), V7_hv = NULL;
	free(V8_hv), V8_hv = NULL;
	free(V9_hv), V9_hv = NULL;

	//-----------------------------Output vertical fields-----------------------------//
		//free(sigma_To_z), sigma_To_z = NULL;
		//free(U3d_out), U3d_out = NULL;
		//free(V3d_out), V3d_out = NULL;
#ifdef _BAROCLINIC
	free(Sb_out), Sb_out = NULL;
	free(Tb_out), Tb_out = NULL;
	free(Ss_out), Ss_out = NULL;
	free(Ts_out), Ts_out = NULL;
	free(Surface_hS), Surface_hS = NULL;
	free(Surface_hT), Surface_hT = NULL;
	free(Bottom_hS), Bottom_hS = NULL;
	free(Bottom_hT), Bottom_hT = NULL;

#endif
}