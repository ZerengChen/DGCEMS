#pragma once
#include<netcdfcpp.h>
#include<vector>
#include<string.h>
class AbstractOutputFile
{
public:
	AbstractOutputFile(const char*, double timeInerval, int StepPerFile);
	~AbstractOutputFile();
	void ncFile_create(int *Np2d, int *K2d, int Nvar);
	void outputIntervalResult(double& time, double *field2d, int Nvar, int *Np2d, int *K2d);
	void outputResult(double time, double *field2d, int Nvar, int *Np2d, int *K2d);
	void closencfile();

	NcFile resultFile;
	NcVar *output_time;
	NcVar *output_fphys;
	char* NCfile_name;

protected:
	double timePrevious;
	double timeInerval;
	int outputStep;
	int StepPerFile;


};


