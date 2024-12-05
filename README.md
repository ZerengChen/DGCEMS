# DGCEMS

## Introduction

This software uses Nodal Discontinuous Galerkin Finite Element Methods to solve various hydraulic problems. 
The software is written in  C++ languages and uses OpenMP and MPI for parallelization.
The simulation results can be exported in the nc format.

## Code download

To download the latest version of DGCEMS:<br>
https://github.com/ZerengChen/DGCEMS

## Environment Deployment & Installation
Before running DGCEMS, make sure you have the following dependency libraries installed.
We recommend users to use.

* Intel oneAPI  v2021.6.0
* gcc v6.1.0
* NetCDF v4.1.2
* NetCDF-cxx4.2
* Cmake v3.6.3
* GOTM v4.1.0
* METIS v5.1.0
* Lapack (need OpenBlas 0.3.8)

## Input file creation
* Step 1: Get pre-processing code<br>

* Step 2: Add folder of example called '@Test' <br>
The '@Test' folder should contain 'Test.m' and grid file 'fort.14'<br>
* Step 3: Edit 'Test.m' <br>
    * Modify: Test < SWEBarotropic3d <br> 
    * Modify: SMSFile = [pwd,'\@Test\fort.14'];<br>
    * Modify: function obj = Test( N, Nz, Mz );<br>

* Step 4: Generating the solver<br>
    * Open the path of 'NDGOM-master'<br>
    * <code>NdgSetup</code><br>
    * <code>Solver=Test(1,1,Layers)</code><br>
* Step 5: <br>
    * Run 'init_fphys.m' and 'meshUnion_nc_Output.m'<br>
    * If the user uses the MPI parallel way,Run 'Geteptreind.m'<br>
* Step 6: Give 'TideElevation.txt'<br>
    * Time series of tide elevation given by the number of open boundary nodes.<br>
* Step 7: Set up 'param.txt' <br>

## DGCEMS Installation
* Step 1: Edit 'CMakeLists.txt' and select the required modules,including:<br>
    * -DDG_THREADS=<br>
    * -D_OPENMP<br>
    * -DCOUPLING_SWAN<br>
    * -D_BAROCLINIC<br>
* Step 2: Configure the project in the build/ directory. <br>    
<code>cd build</code><br>
<code>cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++</code><br>
* Step 3: Compiling projects<br>
<code>make</icode><br>
* Step 4: clean up the build file<br>
If the user wants to clean up the build file, they can use the<br>
<code>make clean</code><br>
Or delete the entire build directory and rebuild<br>
<code>rm -rf build</code><br>
<code>mkdir build</code> <br>
<code>cd build</code> <br>
<code>cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++</code><br>
<code>make</code><br>


## Set up and run
* Step 1: make a folder and copy/link the executable file 'DGCEMS' to this folder.<br>
* Step 2: prepare all required input files.<br>
* Step 3: run the model.<br>