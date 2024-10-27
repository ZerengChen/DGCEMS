#include "AbstractOutputFile.h"
#include <iostream>
using namespace std;

AbstractOutputFile::AbstractOutputFile(const char* NCfile_name, double timeInerval, int StepPerFile) :timePrevious(0), outputStep(0), resultFile(NCfile_name, NcFile::Replace)
{

	this->timeInerval = timeInerval;
	this->StepPerFile = StepPerFile;
}



AbstractOutputFile::~AbstractOutputFile()
{

}

void AbstractOutputFile::closencfile()
{
	resultFile.close();
}
void AbstractOutputFile::outputIntervalResult(double &time, double *field2d, int Nvar, int *Np2d, int *K2d)
{
	if ((time - timePrevious) >= timeInerval)
	{
		outputResult(time, field2d, Nvar, Np2d, K2d);
		timePrevious = time;
	}
};

void AbstractOutputFile::outputResult(double time, double *field2d, int Nvar, int *Np2d, int *K2d)
{
	output_time->set_cur(outputStep);
	output_time->put(&time, 1);

	output_fphys->set_cur(outputStep, 0, 0, 0);
	output_fphys->put(field2d, 1, Nvar, *K2d, *Np2d);

	if (outputStep == StepPerFile + 1)
	{
		//resultFile.close();
	}
	else
	{
		outputStep++;
	}

};

void AbstractOutputFile::ncFile_create(int *Np2d, int *K2d, int Nvar)
{


	NcDim *dimNp = resultFile.add_dim("Np2d", *Np2d);
	NcDim *dimK = resultFile.add_dim("K2d", *K2d);
	NcDim *dimNfield = resultFile.add_dim("Nvar2d", Nvar);
	NcDim *dimtime = resultFile.add_dim("Nt");

	//std::vector<NcDim> dims_time;
	//dims_time.push_back(dimtime);

	//std::cout << "CZRtest1\n";

	//std::vector<NcDim> dims_fphys;
	//dims_fphys.push_back(dimtime);
	//dims_fphys.push_back(dimNfield);
	//dims_fphys.push_back(dimK);
	//dims_fphys.push_back(dimNp);

	output_time = resultFile.add_var("time", ncDouble, dimtime);
	output_time->add_att("units", "s");
	output_fphys = resultFile.add_var("fphys2d", ncDouble, dimtime, dimNfield, dimK, dimNp);

	//resultFile.enddef();
	//std::cout << "CZRtest2\n";
};