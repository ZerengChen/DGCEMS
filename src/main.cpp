//#include "MeshUnion.h"
#include"NdgPhysMat.h"
#include"GlobalVar.h"
#include <chrono>
#ifdef COUPLING_SWAN
#include"c_coupler_cpp_interface_wrj.h"
#endif
#include<mpi.h>
#include<string>
#include<iostream>

MeshUnion_dim Read_NC_dim;
MeshUnion mesh;
const MeshUnion *meshunion = &mesh;

int main(int argc, char* argv)
{
	//-------- For MPI --------//
	//MPI_Init(&argc, &argv);
	//extern int MPIrank;
	//extern int MPIsize;
	//MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	//MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	*dg_demo_comp_id = 0;
	//-------- For MPI --------//

	//Read parameter
	ifstream infile;

	infile.open("param.txt", ios::in);
	if (!infile.is_open())
	{
		cout << "Fail to read parameter!" << endl;
		return 0;
	}

	string line;
	getline(infile, line);
	infile >> ftime;//final time (s)

	getline(infile, line);
	getline(infile, line);
	infile >> *coupling_freq;// coupling freqency (steps). If dt in SWAN is 1(s),dt in DG is 0.01(s),and we meed 2(s) couling,the value is 200.

	getline(infile, line);
	getline(infile, line);
	infile >> *time_step; //Must be a int, for example:1 .We must pay attention to the dt in DG model and the dt in SWAN model.

	getline(infile, line);
	getline(infile, line);
	infile >> First_CouplingStep; //Skip x step in order to wait the wave field steady.

	getline(infile, line);
	getline(infile, line);
	infile >> Calculating_CouplingSteps; //In the model, the total coupling steps are x.

	getline(infile, line);
	getline(infile, line);
	infile >> dt;// One calculating time step of DG model,usually 0.01 (s)

	getline(infile, line);
	getline(infile, line);
	infile >> tidalinterval;// The tidal interval (s), connected to the tidal file.

	getline(infile, line);
	getline(infile, line);
	infile >> NOut; // Output times (times).If 20 means that the results will out 20 times.

	getline(infile, line);
	getline(infile, line);
	infile >> gra; //9.81

	getline(infile, line);
	getline(infile, line);
	infile >> Hcrit; //depmin (m)

	getline(infile, line);
	getline(infile, line);
	infile >> Nvar; //If only current-wave coupling:2,Hu,Hv.If current-wave coupling with sediment:3,Hu,Hv,Hc.

	getline(infile, line);
	getline(infile, line);
	infile >> Nfield3d; //16=9+4+3,{hu hv omega h nv z eta zx zy /U_Euler V_Euler U_Stokes V_Stokes / rho hT hS / hc1......}

	getline(infile, line);
	getline(infile, line);
	infile >> Nfield2d; //8, {'h'  'hu'  'hv'  'z'  'zx'  'zy' / 'U_Euler2d' 'V_Euler2d'}

	getline(infile, line);
	getline(infile, line);
	infile >> rho0;  //density of water,1000kg/m3

	//紊流模型NoGOTM
	getline(infile, line);
	getline(infile, line);
	infile >> z0b;  //Value of the physical bottom roughness length

	getline(infile, line);
	getline(infile, line);
	infile >> z0s;  //Value of the physical surface roughness length:0.02

	//风应力tau/rho
	getline(infile, line);
	getline(infile, line);
	infile >> WindTauxC;  //wind stress at x for constant

	getline(infile, line);
	getline(infile, line);
	infile >> WindTauyC;  //wind stress at y for constant

	getline(infile, line);
	getline(infile, line);
	infile >> Switch_Limiter3D;//1(T) to use Limiter3D, (0)F to close

	getline(infile, line);
	getline(infile, line);
	infile >> Switch_GOTM;//1(T) to use GOTM model, (0)F to use constant cf and nv_con

	getline(infile, line);
	getline(infile, line);
	infile >> cf;//Drag force coefficient, e.g. 0.005

	getline(infile, line);
	getline(infile, line);
	infile >> nv_con;//constant vertical diffusion

	//For tecplot output;
	getline(infile, line);
	getline(infile, line);
	infile >> StepNumber;//OutputTimeStep

	//For Baroclinic
	getline(infile, line);
	getline(infile, line);
	infile >> T0;//init temperature

	getline(infile, line);
	getline(infile, line);
	infile >> S0;//init salimity


	getline(infile, line);
	getline(infile, line);
	infile >> alphaT;//Parameter of accelerate temperature

	getline(infile, line);
	getline(infile, line);
	infile >> betaS;//Parameter of accelerate salimity

	getline(infile, line);
	getline(infile, line);
	infile >> Latitude;//The Latitude

	infile.close();

	std::cout << " Welcome to the 3D DG-FEM Hydrodynamic model by CZR in 2022." << std::endl;
	cout << "ftime = " << ftime << " s" << endl;
	cout << "coupling_freq = " << *coupling_freq << " step/s" << endl;
	cout << "time_step = " << *time_step << " s" << endl;
	cout << "First_CouplingStep = " << First_CouplingStep << " step" << endl;
	cout << "Calculating_CouplingSteps = " << Calculating_CouplingSteps << " step" << endl;
	cout << "dt = " << dt << " s" << endl;
	cout << "tidalinterval = " << tidalinterval << " s" << endl;
	cout << "NOut = " << NOut << " times" << endl;
	cout << "gra = " << gra << " m/s2" << endl;
	cout << "Hcrit = " << Hcrit << " m" << endl;
	cout << "Nvar = " << Nvar << endl;
	cout << "Nfield3d = " << Nfield3d << endl;
	cout << "Nfield2d = " << Nfield2d << endl;
	cout << "rho0 = " << rho0 << " kg/m3" << endl;
	cout << "z0b = " << z0b << " m" << endl;
	cout << "z0s = " << z0s << " m" << endl;
	cout << "WindTauxC = " << WindTauxC << " m2/s2" << endl;
	cout << "WindTauyC = " << WindTauyC << " m2/s2" << endl;
	cout << "Switch_Limiter3D = " << Switch_Limiter3D << endl;
	cout << "Switch_GOTM = " << Switch_GOTM << endl;
	cout << "cf = " << cf << " m" << endl;
	cout << "nv_con = " << nv_con << endl;
	cout << "OutputTime = " << StepNumber * dt << " s" << endl;
	cout << "Init temperature = " << T0 << " C degree" << endl;
	cout << "Init salimity = " << S0 << " g/kg" << endl;
	cout << "Accelerate temperature parameter = " << alphaT << " " << endl;
	cout << "Accelerate salimity parameter = " << betaS << " " << endl;
	cout << "Latitude = " << Latitude << " degree" << endl;
	cout << "End read parameter." << endl;

	//End read parameter.
#ifdef COUPLING_SWAN
	int mpicom = MPI_COMM_NULL ;
#endif

	NdgPhysMat Solver;
	cout << "Solver has been built!" << endl;
#ifdef COUPLING_SWAN	
	register_dg_demo_component_wrj(&mpicom);
#endif

	chrono::steady_clock::time_point t0 = chrono::steady_clock::now();
	Solver.matSolver();
	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
	chrono::duration<double> time_used = chrono::duration_cast < chrono::duration<double>>(t1 - t0);
	cout << "Total calculating time is :" << time_used.count() << " s" << endl;

	cout << "TADA!" << endl;

	//-------- For MPI --------//
	//MPI_Finalize();
	//-------- For MPI --------//
	delete[] dg_demo_comp_id;
	delete[] decomp_id;
	delete[] grid_h2d_id;
	delete[] time_step;
	delete[] coupling_freq;

	return 0;
}

