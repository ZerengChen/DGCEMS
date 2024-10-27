#include "NdgMemory.h"
#include <iostream>
using namespace std;

NdgMemory::NdgMemory()
{
}

NdgMemory::~NdgMemory()
{
}

void MemoryAllocationCheck(double *dest, int size){
    while(dest == NULL){
        dest = (double *)malloc(size); 
    }
}
void MemoryAllocationCheck(signed char *dest, int size) {
	while (dest == NULL) {
		dest = (signed char *)malloc(size);
	}
}

/*The following part is for calculation of density, and is called from mxCalculateDensityField.c*/
double *BaroclinicT = NULL, *BaroclinicS = NULL, *BaroclinicDTS = NULL;

//char *BaroDensityInitialized = "False";

void NdgMemory::BaroDensityMemoryAllocation(int Np, int K) {
	BaroclinicT = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicT, Np*K * sizeof(double));
	BaroclinicS = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicS, Np*K * sizeof(double));
	BaroclinicDTS = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicDTS, Np*K * sizeof(double));
	//BaroDensityInitialized = "True";
}

void NdgMemory::BaroDensityMemoryDeAllocation() {
	free(BaroclinicT), BaroclinicT = NULL;
	free(BaroclinicS), BaroclinicS = NULL;
	free(BaroclinicDTS), BaroclinicDTS = NULL;
	//BaroDensityInitialized = "False";
}

/*The following part is for calculation of baroclinic term, and is called from mxCalculateBaroclinicTerm.c*/

double *BaroclinicPDPX = NULL, *BaroclinicPDPY = NULL, *BaroclinicPRHOPX = NULL, *BaroclinicPRHOPY = NULL, \
*BaroclinicPRHOPS = NULL, *BaroclinicInXPartOne = NULL, *BaroclinicInXPartTwo = NULL, *BaroclinicInYPartOne = NULL, \
*BaroclinicInYPartTwo = NULL, *BaroclinicInXTempRHS = NULL, *BaroclinicInYTempRHS = NULL, *Baroclinicfmod = NULL, \
*BaroclinicBotEfm = NULL, *BaroclinicBotEfp = NULL, *BaroclinicBotEFluxM = NULL, *BaroclinicBotEFluxP = NULL, \
*BaroclinicBotEFluxS = NULL, *BaroclinicIEfm = NULL, *BaroclinicIEfp = NULL, *BaroclinicIEfluxM = NULL, \
*BaroclinicIEfluxP = NULL, *BaroclinicIEfluxS = NULL, *BaroclinicERHS = NULL, *BaroclinicTempFacialIntegral = NULL, \
*BaroclinicTempVolumeIntegral = NULL, *BaroclinicBEfm = NULL, *BaroclinicBEfp = NULL, *BaroclinicBEfluxM = NULL, \
*BaroclinicBEfluxS = NULL;

//char *BaroclinicPartInitialized = "False";

void NdgMemory::BaroclinicPartMemoryAllocation(int Np, int K, int K2d, int BotENe, int BotENfp, int IENe, \
	int IENfp, int Nface, int BENe, int BENfp) {
	BaroclinicPDPX = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPDPX, Np*K * sizeof(double));
	BaroclinicPDPY = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPDPY, Np*K * sizeof(double));
	BaroclinicPRHOPX = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPRHOPX, Np*K * sizeof(double));
	BaroclinicPRHOPY = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPRHOPY, Np*K * sizeof(double));
	BaroclinicPRHOPS = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicPRHOPS, Np*K * sizeof(double));
	BaroclinicInXPartOne = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInXPartOne, Np*K * sizeof(double));
	BaroclinicInXPartTwo = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInXPartTwo, Np*K * sizeof(double));
	BaroclinicInYPartOne = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInYPartOne, Np*K * sizeof(double));
	BaroclinicInYPartTwo = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInYPartTwo, Np*K * sizeof(double));
	BaroclinicInXTempRHS = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInXTempRHS, Np*K * sizeof(double));
	BaroclinicInYTempRHS = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicInYTempRHS, Np*K * sizeof(double));
	Baroclinicfmod = (double*)malloc(Np*K2d * sizeof(double));
	MemoryAllocationCheck(Baroclinicfmod, Np*K2d * sizeof(double));
	BaroclinicBotEfm = (double*)malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEfm, Np*K * sizeof(double));
	BaroclinicBotEfp = (double*)malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEfp, BotENe*BotENfp * sizeof(double));
	BaroclinicBotEFluxM = (double*)malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEFluxM, BotENe*BotENfp * sizeof(double));
	BaroclinicBotEFluxP = (double*)malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEFluxP, BotENe*BotENfp * sizeof(double));
	BaroclinicBotEFluxS = (double*)malloc(BotENe*BotENfp * sizeof(double));
	MemoryAllocationCheck(BaroclinicBotEFluxS, BotENe*BotENfp * sizeof(double));
	BaroclinicTempFacialIntegral = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicTempFacialIntegral, Np*K * sizeof(double));
	BaroclinicIEfm = (double*)malloc(IENe*IENfp * 2 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfm, IENe*IENfp * 2 * sizeof(double));
	BaroclinicIEfp = (double*)malloc(IENe*IENfp * 2 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfp, IENe*IENfp * 2 * sizeof(double));
	BaroclinicIEfluxM = (double*)malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfluxM, IENe*IENfp * 4 * sizeof(double));
	BaroclinicIEfluxP = (double*)malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfluxP, IENe*IENfp * 4 * sizeof(double));
	BaroclinicIEfluxS = (double*)malloc(IENe*IENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicIEfluxS, IENe*IENfp * 4 * sizeof(double));
	BaroclinicERHS = (double*)malloc(4 * Np*K*(Nface - 2) * sizeof(double));
	MemoryAllocationCheck(BaroclinicERHS, 4 * Np*K*(Nface - 2) * sizeof(double));
	BaroclinicTempVolumeIntegral = (double*)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(BaroclinicTempVolumeIntegral, Np*K * sizeof(double));
	BaroclinicBEfm = (double*)malloc(BENe*BENfp * 2 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfm, BENe*BENfp * 2 * sizeof(double));
	BaroclinicBEfp = (double*)malloc(BENe*BENfp * 3 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfp, BENe*BENfp * 3 * sizeof(double));
	BaroclinicBEfluxM = (double*)malloc(BENe*BENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfluxM, BENe*BENfp * 4 * sizeof(double));
	BaroclinicBEfluxS = (double*)malloc(BENe*BENfp * 4 * sizeof(double));
	MemoryAllocationCheck(BaroclinicBEfluxS, BENe*BENfp * 4 * sizeof(double));
	//BaroclinicPartInitialized = "True";
}

void NdgMemory::BaroclinicPartMemoryDeAllocation() {
	free(BaroclinicPDPX), BaroclinicPDPX = NULL;
	free(BaroclinicPDPY), BaroclinicPDPY = NULL;
	free(BaroclinicPRHOPX), BaroclinicPRHOPX = NULL;
	free(BaroclinicPRHOPY), BaroclinicPRHOPY = NULL;
	free(BaroclinicPRHOPS), BaroclinicPRHOPS = NULL;
	free(BaroclinicInXPartOne), BaroclinicInXPartOne = NULL;
	free(BaroclinicInXPartTwo), BaroclinicInXPartTwo = NULL;
	free(BaroclinicInYPartOne), BaroclinicInYPartOne = NULL;
	free(BaroclinicInYPartTwo), BaroclinicInYPartTwo = NULL;
	free(BaroclinicInXTempRHS), BaroclinicInXTempRHS = NULL;
	free(BaroclinicInYTempRHS), BaroclinicInYTempRHS = NULL;
	free(Baroclinicfmod), Baroclinicfmod = NULL;
	free(BaroclinicBotEfm), BaroclinicBotEfm = NULL;
	free(BaroclinicBotEfp), BaroclinicBotEfp = NULL;
	free(BaroclinicBotEFluxM), BaroclinicBotEFluxM = NULL;
	free(BaroclinicBotEFluxP), BaroclinicBotEFluxP = NULL;
	free(BaroclinicBotEFluxS), BaroclinicBotEFluxS = NULL;
	free(BaroclinicTempFacialIntegral), BaroclinicTempFacialIntegral = NULL;
	free(BaroclinicIEfm), BaroclinicIEfm = NULL;
	free(BaroclinicIEfp), BaroclinicIEfp = NULL;
	free(BaroclinicIEfluxM), BaroclinicIEfluxM = NULL;
	free(BaroclinicIEfluxP), BaroclinicIEfluxP = NULL;
	free(BaroclinicIEfluxS), BaroclinicIEfluxS = NULL;
	free(BaroclinicERHS), BaroclinicERHS = NULL;
	free(BaroclinicTempFacialIntegral), BaroclinicTempFacialIntegral = NULL;
	free(BaroclinicTempVolumeIntegral), BaroclinicTempVolumeIntegral = NULL;
	free(BaroclinicBEfm), BaroclinicBEfm = NULL;
	free(BaroclinicBEfp), BaroclinicBEfp = NULL;
	free(BaroclinicBEfluxM), BaroclinicBEfluxM = NULL;
	free(BaroclinicBEfluxS), BaroclinicBEfluxS = NULL;
	//BaroclinicPartInitialized = "False";
}

/*This is for vertical diffusion part*/
double *Tau = NULL, *BBE = NULL, *SBE = NULL, *Hhuv2d = NULL, \
*Tau_w = NULL, *Tau_CW = NULL, *fw = NULL, *Ab = NULL, *Eta2d = NULL;
//char *VertDiffInitialized = "False";

void NdgMemory::VertDiffMemoryAllocation(const int Np2d, int K2d, const int Nz){
	Tau = (double *)malloc(sizeof(double)*(Np2d*K2d*(Nz+1)));
    MemoryAllocationCheck(Tau, sizeof(double)*(Np2d*K2d*(Nz+1)));
	//memset(Tau, 0, sizeof(double)*(Np2d*K2d*(Nz + 1)));
	BBE = (double *)malloc(sizeof(double)*(Np2d*K2d*4));//Include T S
	MemoryAllocationCheck(BBE, sizeof(double)*(Np2d*K2d * 4));
	memset(BBE, 0, sizeof(double)*(Np2d*K2d * 4));
	SBE = (double *)malloc(sizeof(double)*(Np2d*K2d * 4));
	MemoryAllocationCheck(SBE, sizeof(double)*(Np2d*K2d * 4));
	memset(SBE, 0, sizeof(double)*(Np2d*K2d * 4));
	Hhuv2d = (double *)malloc(sizeof(double)*Np2d*K2d*3);
	MemoryAllocationCheck(Hhuv2d, sizeof(double)*(Np2d*K2d*3));
	//memset(Hhuv2d, 0, sizeof(double)*(Np2d*K2d * 3));
	Eta2d = (double *)malloc(sizeof(double)*(Np2d*K2d));
	MemoryAllocationCheck(Eta2d, sizeof(double)*(Np2d*K2d));
	memset(Eta2d, 0, sizeof(double)*(Np2d*K2d));
	//VertDiffInitialized = "True";
}

void NdgMemory::VertDiffMemoryDeAllocation(){
	free(Tau);
	Tau = NULL;
	//free(BBE); BBE = NULL;
	//free(SBE); SBE = NULL;
	//free(Hhuv2d); Hhuv2d = NULL;
	//VertDiffInitialized = "False";
}

/*This is for GOTM part*/
double *tkeGOTM = NULL, *epsGOTM = NULL, *LGOTM = NULL, *nuhGOTM = NULL, \
*numGOTM = NULL, *layerHeight = NULL, *huCentralDate = NULL, *hvCentralDate = NULL, \
*huVerticalLine = NULL, *hvVerticalLine = NULL, *shearFrequencyDate = NULL, *buoyanceFrequencyDate = NULL, \
*BottomFrictionLength = NULL, *BottomFrictionVelocity = NULL, *SurfaceFrictionLength = NULL, \
*SurfaceFrictionVelocity = NULL, *eddyViscosityDate = NULL, *rhoCentralDate = NULL, *rhoVerticalLine = NULL, \
*eddyDiffusionDate = NULL, *eddyTKEDate = NULL, *eddyLengthDate = NULL, *eddyEPSDate = NULL;

//char *GOTMInitialized = "False";

void NdgMemory::GotmSolverMemoryAllocation(int Num2d, int Interface, int Np2d, int K3d){
	tkeGOTM = (double *)malloc(sizeof(double)*(Num2d*Interface)); 
    MemoryAllocationCheck(tkeGOTM, sizeof(double)*(Num2d*Interface));
	memset(tkeGOTM, 0, sizeof(double)*(Num2d*Interface));
	epsGOTM = (double *)malloc(sizeof(double)*(Num2d*Interface)); 
    MemoryAllocationCheck(epsGOTM, sizeof(double)*(Num2d*Interface));
	memset(epsGOTM, 0, sizeof(double)*(Num2d*Interface));
	LGOTM = (double *)malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(LGOTM, sizeof(double)*(Num2d*Interface));
	memset(LGOTM, 0, sizeof(double)*(Num2d*Interface));
	nuhGOTM = (double *)malloc(sizeof(double)*(Num2d*Interface)); 
    MemoryAllocationCheck(nuhGOTM, sizeof(double)*(Num2d*Interface));
	memset(nuhGOTM, 0, sizeof(double)*(Num2d*Interface));
	numGOTM = (double *)malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(numGOTM, sizeof(double)*(Num2d*Interface));
	memset(numGOTM, 0, sizeof(double)*(Num2d*Interface));
	layerHeight = (double *)malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(layerHeight, sizeof(double)*(Num2d*Interface));
	memset(layerHeight, 0, sizeof(double)*(Num2d*Interface));
	huCentralDate = (double *)malloc(sizeof(double)*(Np2d*K3d));
    MemoryAllocationCheck(huCentralDate, sizeof(double)*(Np2d*K3d));
	memset(huCentralDate, 0, sizeof(double)*(Np2d*K3d));
	hvCentralDate = (double *)malloc(sizeof(double)*(Np2d*K3d));
    MemoryAllocationCheck(hvCentralDate, sizeof(double)*(Np2d*K3d));
	memset(hvCentralDate, 0, sizeof(double)*(Np2d*K3d));
	huVerticalLine = (double *)malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(huVerticalLine, sizeof(double)*(Num2d*Interface));
	memset(huVerticalLine, 0, sizeof(double)*(Num2d*Interface));
	hvVerticalLine = (double *)malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(hvVerticalLine, sizeof(double)*(Num2d*Interface));
	memset(hvVerticalLine, 0, sizeof(double)*(Num2d*Interface));
	shearFrequencyDate = (double *)malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(shearFrequencyDate, sizeof(double)*(Num2d*Interface));
	memset(shearFrequencyDate, 0, sizeof(double)*(Num2d*Interface));
	buoyanceFrequencyDate = (double *)malloc(sizeof(double)*(Num2d*Interface));
    MemoryAllocationCheck(buoyanceFrequencyDate, sizeof(double)*(Num2d*Interface));
	memset(buoyanceFrequencyDate, 0, sizeof(double)*(Num2d*Interface));
	BottomFrictionLength = (double *)malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(BottomFrictionLength, sizeof(double)*Num2d);
	memset(BottomFrictionLength, 0, sizeof(double)*Num2d);
	BottomFrictionVelocity = (double *)malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(BottomFrictionVelocity, sizeof(double)*Num2d);
	memset(BottomFrictionVelocity, 0, sizeof(double)*Num2d);
	SurfaceFrictionLength = (double *)malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(SurfaceFrictionLength, sizeof(double)*Num2d);
	memset(SurfaceFrictionLength, 0, sizeof(double)*Num2d);
	SurfaceFrictionVelocity = (double *)malloc(sizeof(double)*Num2d);
    MemoryAllocationCheck(SurfaceFrictionVelocity, sizeof(double)*Num2d);
	memset(SurfaceFrictionVelocity, 0, sizeof(double)*Num2d);
	eddyViscosityDate = (double *)malloc(sizeof(double)*(Num2d * Interface));
    MemoryAllocationCheck(eddyViscosityDate, sizeof(double)*(Num2d * Interface));
	memset(eddyViscosityDate, 0, sizeof(double)*(Num2d * Interface));
	rhoCentralDate = (double *)malloc(sizeof(double)*(Np2d*K3d));
	MemoryAllocationCheck(rhoCentralDate, sizeof(double)*(Np2d*K3d));
	rhoVerticalLine = (double *)malloc(sizeof(double)*(Num2d*Interface));
	MemoryAllocationCheck(rhoVerticalLine, sizeof(double)*(Num2d*Interface));
	eddyTKEDate = (double *)malloc(sizeof(double)*(Num2d * Interface));
	MemoryAllocationCheck(eddyTKEDate, sizeof(double)*(Num2d * Interface));
	eddyLengthDate = (double *)malloc(sizeof(double)*(Num2d * Interface));
	MemoryAllocationCheck(eddyLengthDate, sizeof(double)*(Num2d * Interface));
	eddyEPSDate = (double *)malloc(sizeof(double)*(Num2d * Interface));
	MemoryAllocationCheck(eddyEPSDate, sizeof(double)*(Num2d * Interface));
	eddyDiffusionDate = (double *)malloc(sizeof(double)*(Num2d * Interface));
	MemoryAllocationCheck(eddyDiffusionDate, sizeof(double)*(Num2d * Interface));
	//GOTMInitialized = "True";
}

void NdgMemory::GotmSolverMemoryDeAllocation(){
	free(nuhGOTM); nuhGOTM = NULL;
	free(numGOTM); numGOTM = NULL;
	free(tkeGOTM); tkeGOTM = NULL;
	free(epsGOTM); epsGOTM = NULL;
	free(LGOTM); LGOTM = NULL;
	free(layerHeight); layerHeight = NULL;
	free(huCentralDate); huCentralDate = NULL;
	free(hvCentralDate); hvCentralDate = NULL;
	free(huVerticalLine); huVerticalLine = NULL;
	free(hvVerticalLine); hvVerticalLine = NULL;
	free(shearFrequencyDate); shearFrequencyDate = NULL;
	free(buoyanceFrequencyDate); buoyanceFrequencyDate = NULL;
	free(BottomFrictionLength); BottomFrictionLength = NULL;
	free(BottomFrictionVelocity); BottomFrictionVelocity = NULL;
	free(SurfaceFrictionLength); SurfaceFrictionLength = NULL;
	free(SurfaceFrictionVelocity); SurfaceFrictionVelocity = NULL;
	free(eddyViscosityDate); eddyViscosityDate = NULL;
	//GOTMInitialized = "False";
}

/*This is for updated vertical velocity solver part*/
double *UpdatedVSrhs2d = NULL, *UpdatedVSIEfm2d = NULL, *UpdatedVSIEfp2d = NULL, *UpdatedVSIEFluxM2d = NULL, \
*UpdatedVSIEFluxP2d = NULL, *UpdatedVSIEFluxS2d = NULL, *UpdatedVSERHS2d = NULL, *UpdatedVSVolumeIntegralX = NULL, \
*UpdatedVSTempVolumeIntegralX = NULL, *UpdatedVSVolumeIntegralY = NULL, *UpdatedVSTempVolumeIntegralY = NULL, \
*UpdatedVSBEfm2d = NULL, *UpdatedVSBEzM2d = NULL, *UpdatedVSBEfp2d = NULL, *UpdatedVSBEzP2d = NULL, *UpdatedVSBEFluxS2d = NULL, \
*UpdatedVSBEFluxM2d = NULL, *UpdatedVSTempFacialIntegral = NULL, *UpdatedVSfield2d = NULL, *UpdatedVSrhs3d = NULL, \
*UpdatedVSIEfm3d = NULL, *UpdatedVSIEfp3d = NULL, *UpdatedVSIEFluxM3d = NULL, *UpdatedVSIEFluxP3d = NULL, *UpdatedVSIEFluxS3d = NULL, \
*UpdatedVSERHS3d = NULL, *UpdatedVSVolumeIntegralX3d = NULL, *UpdatedVSTempVolumeIntegralX3d = NULL, \
*UpdatedVSVolumeIntegralY3d = NULL, *UpdatedVSTempVolumeIntegralY3d = NULL, *UpdatedVSBEfm3d = NULL, \
*UpdatedVSBEzM3d = NULL, *UpdatedVSBEfp3d = NULL, *UpdatedVSBEzP3d = NULL, *UpdatedVSBEFluxS3d = NULL, *UpdatedVSBEFluxM3d = NULL, \
*UpdatedVSTempFacialIntegral3d = NULL, *UpdatedVSIEfmod = NULL, *UpdatedVSBEfmod = NULL, *Updatedfmod = NULL;

//char *UpdatedVertVelocityInitialized = "False";

void NdgMemory::UpdatedVertVelocitySolverMemoryAllocation(int Np2d, int K2d, int IENfp2d, int IENe2d, int Nface, int BENe2d, int BENfp2d, int Np3d, \
	int K3d, int IENfp3d, int IENe3d, int BENe3d, int BENfp3d){
	UpdatedVSrhs2d = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSrhs2d, Np2d*K2d*sizeof(double));
	memset(UpdatedVSrhs2d, 0, Np2d*K2d * sizeof(double));
	UpdatedVSIEfm2d = (double *)malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfm2d, IENfp2d*IENe2d * 3 * sizeof(double));
	memset(UpdatedVSIEfm2d, 0, IENfp2d*IENe2d * 3 * sizeof(double));
	UpdatedVSIEfp2d = (double *)malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfp2d, IENfp2d*IENe2d * 3 * sizeof(double));
	memset(UpdatedVSIEfp2d, 0, IENfp2d*IENe2d * 3 * sizeof(double));
	UpdatedVSIEFluxM2d = (double *)malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxM2d, IENfp2d*IENe2d*sizeof(double));
	memset(UpdatedVSIEFluxM2d, 0, IENfp2d*IENe2d * sizeof(double));
	UpdatedVSIEFluxP2d = (double *)malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxP2d, IENfp2d*IENe2d*sizeof(double));
	memset(UpdatedVSIEFluxP2d, 0, IENfp2d*IENe2d * sizeof(double));
	UpdatedVSIEFluxS2d = (double *)malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxS2d, IENfp2d*IENe2d*sizeof(double));
	memset(UpdatedVSIEFluxS2d, 0, IENfp2d*IENe2d * sizeof(double));
	UpdatedVSERHS2d = (double *)malloc(Np2d*K2d*Nface*sizeof(double));
	MemoryAllocationCheck(UpdatedVSERHS2d, Np2d*K2d*Nface*sizeof(double));
	memset(UpdatedVSERHS2d, 0, Np2d*K2d*Nface * sizeof(double));
	UpdatedVSVolumeIntegralX = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralX, Np2d*K2d*sizeof(double));
	memset(UpdatedVSVolumeIntegralX, 0, Np2d*K2d * sizeof(double));
	UpdatedVSTempVolumeIntegralX = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralX, Np2d*K2d*sizeof(double));
	memset(UpdatedVSTempVolumeIntegralX, 0, Np2d*K2d * sizeof(double));
	UpdatedVSVolumeIntegralY = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralY, Np2d*K2d*sizeof(double));
	memset(UpdatedVSVolumeIntegralY, 0, Np2d*K2d * sizeof(double));
	UpdatedVSTempVolumeIntegralY = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralY, Np2d*K2d*sizeof(double));
	memset(UpdatedVSTempVolumeIntegralY, 0, Np2d*K2d * sizeof(double));
	UpdatedVSBEfm2d = (double *)malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfm2d, BENe2d * BENfp2d * 3 * sizeof(double));
	memset(UpdatedVSBEfm2d, 0, BENe2d * BENfp2d * 3 * sizeof(double));
	UpdatedVSBEzM2d = (double *)malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzM2d, BENe2d * BENfp2d*sizeof(double));
	memset(UpdatedVSBEzM2d, 0, BENe2d * BENfp2d * sizeof(double));
	UpdatedVSBEfp2d = (double *)malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfp2d, BENe2d * BENfp2d * 3 * sizeof(double));
	memset(UpdatedVSBEfp2d, 0, BENe2d * BENfp2d * 3 * sizeof(double));
	UpdatedVSBEzP2d = (double *)malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzP2d, BENe2d * BENfp2d*sizeof(double));
	memset(UpdatedVSBEzP2d, 0, BENe2d * BENfp2d * sizeof(double));
	UpdatedVSBEFluxS2d = (double *)malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxS2d, BENe2d*BENfp2d*sizeof(double));
	memset(UpdatedVSBEFluxS2d, 0, BENe2d*BENfp2d * sizeof(double));
	UpdatedVSBEFluxM2d = (double *)malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxM2d, BENe2d*BENfp2d*sizeof(double));
	memset(UpdatedVSBEFluxM2d, 0, BENe2d*BENfp2d * sizeof(double));
	UpdatedVSTempFacialIntegral = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempFacialIntegral, Np2d*K2d*sizeof(double));
	memset(UpdatedVSTempFacialIntegral, 0, Np2d*K2d * sizeof(double));
	UpdatedVSfield2d = (double *)malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSfield2d, Np3d*K3d*sizeof(double));
	memset(UpdatedVSfield2d, 0, Np3d*K3d * sizeof(double));
	UpdatedVSrhs3d = (double *)malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSrhs3d, Np3d*K3d*sizeof(double));
	memset(UpdatedVSrhs3d, 0, Np3d*K3d * sizeof(double));
	UpdatedVSIEfm3d = (double *)malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfm3d, IENfp3d*IENe3d * 3 * sizeof(double));
	memset(UpdatedVSIEfm3d, 0, IENfp3d*IENe3d * 3 * sizeof(double));
	UpdatedVSIEfp3d = (double *)malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfp3d, IENfp3d*IENe3d * 3 * sizeof(double));
	memset(UpdatedVSIEfp3d, 0, IENfp3d*IENe3d * 3 * sizeof(double));
	UpdatedVSIEFluxM3d = (double *)malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxM3d, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSIEFluxM3d, 0, IENfp3d*IENe3d * sizeof(double));
	UpdatedVSIEFluxP3d = (double *)malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxP3d, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSIEFluxP3d, 0, IENfp3d*IENe3d * sizeof(double));
	UpdatedVSIEFluxS3d = (double *)malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEFluxS3d, IENfp3d*IENe3d*sizeof(double));
	memset(UpdatedVSIEFluxS3d, 0, IENfp3d*IENe3d * sizeof(double));
	UpdatedVSERHS3d = (double *)malloc(Np3d*K3d*Nface*sizeof(double));
	MemoryAllocationCheck(UpdatedVSERHS3d, Np3d*K3d*Nface*sizeof(double));
	memset(UpdatedVSERHS3d, 0, Np3d*K3d*Nface * sizeof(double));
	UpdatedVSVolumeIntegralX3d = (double *)malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralX3d, Np3d*K3d*sizeof(double));
	memset(UpdatedVSVolumeIntegralX3d, 0, Np3d*K3d * sizeof(double));
	UpdatedVSTempVolumeIntegralX3d = (double *)malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralX3d, Np3d*K3d*sizeof(double));
	memset(UpdatedVSTempVolumeIntegralX3d, 0, Np3d*K3d * sizeof(double));
	UpdatedVSVolumeIntegralY3d = (double *)malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSVolumeIntegralY3d, Np3d*K3d*sizeof(double));
	memset(UpdatedVSVolumeIntegralY3d, 0, Np3d*K3d * sizeof(double));
	UpdatedVSTempVolumeIntegralY3d = (double *)malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempVolumeIntegralY3d, Np3d*K3d*sizeof(double));
	memset(UpdatedVSTempVolumeIntegralY3d, 0, Np3d*K3d * sizeof(double));
	UpdatedVSBEfm3d = (double *)malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfm3d, BENe3d * BENfp3d * 3 * sizeof(double));
	memset(UpdatedVSBEfm3d, 0, BENe3d * BENfp3d * 3 * sizeof(double));
	UpdatedVSBEzM3d = (double *)malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzM3d, BENe3d * BENfp3d * sizeof(double));
	memset(UpdatedVSBEzM3d, 0, BENe3d * BENfp3d * sizeof(double));
	UpdatedVSBEfp3d = (double *)malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfp3d, BENe3d * BENfp3d * 3 * sizeof(double));
	memset(UpdatedVSBEfp3d, 0, BENe3d * BENfp3d * 3 * sizeof(double));
	UpdatedVSBEzP3d = (double *)malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEzP3d, BENe3d * BENfp3d * sizeof(double));
	memset(UpdatedVSBEzP3d, 0, BENe3d * BENfp3d * sizeof(double));
	UpdatedVSBEFluxS3d = (double *)malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxS3d, BENe3d*BENfp3d*sizeof(double));
	memset(UpdatedVSBEFluxS3d, 0, BENe3d*BENfp3d * sizeof(double));
	UpdatedVSBEFluxM3d = (double *)malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxM3d, BENe3d*BENfp3d*sizeof(double));
	memset(UpdatedVSBEFluxM3d, 0, BENe3d*BENfp3d * sizeof(double));
	UpdatedVSTempFacialIntegral3d = (double *)malloc(Np3d*K3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSTempFacialIntegral3d, Np3d*K3d*sizeof(double));
	memset(UpdatedVSTempFacialIntegral3d, 0, Np3d*K3d * sizeof(double));
	UpdatedVSBEFluxM3d = (double *)malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEFluxM3d, BENe3d*BENfp3d*sizeof(double));
	memset(UpdatedVSBEFluxM3d, 0, BENe3d*BENfp3d * sizeof(double));
	UpdatedVSIEfmod = (double *)malloc(IENe2d*IENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSIEfmod, IENe2d*IENfp3d*sizeof(double));
	memset(UpdatedVSIEfmod, 0, IENe2d*IENfp3d * sizeof(double));
	UpdatedVSBEfmod = (double *)malloc(BENe2d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(UpdatedVSBEfmod, BENe2d*BENfp3d*sizeof(double));
	memset(UpdatedVSBEfmod, 0, BENe2d*BENfp3d * sizeof(double));
	Updatedfmod = (double *)malloc(K2d*Np3d*sizeof(double));
	MemoryAllocationCheck(Updatedfmod, K2d*Np3d*sizeof(double));
	memset(Updatedfmod, 0, K2d*Np3d * sizeof(double));
	//UpdatedVertVelocityInitialized = "True";
}

void NdgMemory::UpdatedVertVelocitySolverMemoryDeAllocation(){
	//free(UpdatedVSrhs2d), UpdatedVSrhs2d = NULL;//有问题，放最后
	free(UpdatedVSIEfm2d); UpdatedVSIEfm2d = NULL;
	free(UpdatedVSIEfp2d); UpdatedVSIEfp2d = NULL;
	free(UpdatedVSIEFluxM2d); UpdatedVSIEFluxM2d = NULL;
	free(UpdatedVSIEFluxP2d); UpdatedVSIEFluxP2d = NULL;
	free(UpdatedVSIEFluxS2d); UpdatedVSIEFluxS2d = NULL;
	free(UpdatedVSERHS2d);  UpdatedVSERHS2d = NULL;
	free(UpdatedVSVolumeIntegralX); UpdatedVSVolumeIntegralX = NULL;
	free(UpdatedVSTempVolumeIntegralX); UpdatedVSTempVolumeIntegralX = NULL;
	free(UpdatedVSVolumeIntegralY); UpdatedVSVolumeIntegralY = NULL;
	free(UpdatedVSTempVolumeIntegralY); UpdatedVSTempVolumeIntegralY = NULL;
	free(UpdatedVSBEfm2d); UpdatedVSBEfm2d = NULL;
	free(UpdatedVSBEzM2d); UpdatedVSBEzM2d = NULL;
	free(UpdatedVSBEfp2d); UpdatedVSBEfp2d = NULL;
	free(UpdatedVSBEzP2d); UpdatedVSBEzP2d = NULL;
	free(UpdatedVSBEFluxS2d); UpdatedVSBEFluxS2d = NULL;
	free(UpdatedVSBEFluxM2d); UpdatedVSBEFluxM2d = NULL;
	free(UpdatedVSTempFacialIntegral); UpdatedVSTempFacialIntegral = NULL;
	free(UpdatedVSfield2d); UpdatedVSfield2d = NULL;
	free(UpdatedVSrhs3d); UpdatedVSrhs3d = NULL;
	free(UpdatedVSIEfm3d); UpdatedVSIEfm3d = NULL;
	free(UpdatedVSIEfp3d); UpdatedVSIEfp3d = NULL;
	free(UpdatedVSIEFluxM3d); UpdatedVSIEFluxM3d = NULL;
	free(UpdatedVSIEFluxP3d); UpdatedVSIEFluxP3d = NULL;
	free(UpdatedVSIEFluxS3d); UpdatedVSIEFluxS3d = NULL;
	free(UpdatedVSERHS3d); UpdatedVSERHS3d = NULL;
	free(UpdatedVSVolumeIntegralX3d); UpdatedVSVolumeIntegralX3d = NULL;
	free(UpdatedVSTempVolumeIntegralX3d); UpdatedVSTempVolumeIntegralX3d = NULL;
	free(UpdatedVSVolumeIntegralY3d); UpdatedVSVolumeIntegralY3d = NULL;
	free(UpdatedVSTempVolumeIntegralY3d); UpdatedVSTempVolumeIntegralY3d = NULL;
	free(UpdatedVSBEfm3d); UpdatedVSBEfm3d = NULL;
	free(UpdatedVSBEzM3d); UpdatedVSBEzM3d = NULL;
	free(UpdatedVSBEfp3d); UpdatedVSBEfp3d = NULL;
	free(UpdatedVSBEzP3d); UpdatedVSBEzP3d = NULL;
	free(UpdatedVSBEFluxS3d); UpdatedVSBEFluxS3d = NULL;
	free(UpdatedVSBEFluxM3d); UpdatedVSBEFluxM3d = NULL;
	free(UpdatedVSTempFacialIntegral3d); UpdatedVSTempFacialIntegral3d = NULL;
	free(UpdatedVSIEfmod); UpdatedVSIEfmod = NULL;
	free(UpdatedVSBEfmod); UpdatedVSBEfmod = NULL;
	free(Updatedfmod); Updatedfmod = NULL;
	//free(UpdatedVSrhs2d); UpdatedVSrhs2d = NULL;
	//UpdatedVertVelocityInitialized = "False";
}

/*This is for Updated PCE Solver part*/
double *PCEUpdatedIEfm2d = NULL, *PCEUpdatedIEfp2d = NULL, *PCEUpdatedIEFluxM2d = NULL, *PCEUpdatedIEFluxP2d = NULL, *PCEUpdatedIEFluxS2d = NULL, \
*PCEUpdatedERHS2d = NULL, *PCEUpdatedVolumeIntegralX = NULL, *PCEUpdatedTempVolumeIntegralX = NULL, *PCEUpdatedVolumeIntegralY = NULL, \
*PCEUpdatedTempVolumeIntegralY = NULL, *PCEUpdatedBEfm2d = NULL, *PCEUpdatedBEzM2d = NULL, *PCEUpdatedBEfp2d = NULL, *PCEUpdatedBEzP2d = NULL, \
*PCEUpdatedBEFluxS2d = NULL, *PCEUpdatedBEFluxM2d = NULL, *PCEUpdatedPCETempFacialIntegral = NULL, *PCEUpdatedIEfmod = NULL, *PCEUpdatedBEfmod = NULL, \
*PCEUpdatedIEfm3d = NULL, *PCEUpdatedIEfp3d = NULL, *PCEUpdatedIEFluxS3d = NULL, *PCEUpdatedBEfm3d = NULL, *PCEUpdatedBEfp3d = NULL, *PCEUpdatedBEFluxS3d = NULL, \
*PCEUpdatedBEzM3d = NULL, *PCEUpdatedBEzP3d = NULL;

//char *PCEUpdatedInitialized = "False";

void NdgMemory::PCEUpdatedMemoryAllocation(int IENfp2d, int IENe2d, int Np2d, int K2d, int Nface, int BENe2d, int BENfp2d, \
	int IENfp3d, int IENe3d, int BENe3d, int BENfp3d){
	PCEUpdatedIEfm2d = (double *)malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfm2d, IENfp2d*IENe2d * 3 * sizeof(double));
	memset(PCEUpdatedIEfm2d, 0, IENfp2d*IENe2d * 3 * sizeof(double));
	PCEUpdatedIEfp2d = (double *)malloc(IENfp2d*IENe2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfp2d, IENfp2d*IENe2d * 3 * sizeof(double));
	memset(PCEUpdatedIEfp2d, 0, IENfp2d*IENe2d * 3 * sizeof(double));
	PCEUpdatedIEFluxM2d = (double *)malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxM2d, IENfp2d*IENe2d*sizeof(double));
	memset(PCEUpdatedIEFluxM2d, 0, IENfp2d*IENe2d * sizeof(double));
	PCEUpdatedIEFluxP2d = (double *)malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxP2d, IENfp2d*IENe2d*sizeof(double));
	memset(PCEUpdatedIEFluxP2d, 0, IENfp2d*IENe2d * sizeof(double));
	PCEUpdatedIEFluxS2d = (double *)malloc(IENfp2d*IENe2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxS2d, IENfp2d*IENe2d*sizeof(double));
	memset(PCEUpdatedIEFluxS2d, 0, IENfp2d*IENe2d * sizeof(double));
	PCEUpdatedERHS2d = (double *)malloc(Np2d*K2d*Nface*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedERHS2d, Np2d*K2d*Nface*sizeof(double));
	memset(PCEUpdatedERHS2d, 0, Np2d*K2d*Nface * sizeof(double));
	PCEUpdatedVolumeIntegralX = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedVolumeIntegralX, Np2d*K2d*sizeof(double));
	memset(PCEUpdatedVolumeIntegralX, 0, Np2d*K2d * sizeof(double));
	PCEUpdatedTempVolumeIntegralX = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedTempVolumeIntegralX, Np2d*K2d*sizeof(double));
	memset(PCEUpdatedTempVolumeIntegralX, 0, Np2d*K2d * sizeof(double));
	PCEUpdatedVolumeIntegralY = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedVolumeIntegralY, Np2d*K2d*sizeof(double));
	memset(PCEUpdatedVolumeIntegralY, 0, Np2d*K2d * sizeof(double));
	PCEUpdatedTempVolumeIntegralY = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedTempVolumeIntegralY, Np2d*K2d*sizeof(double));
	memset(PCEUpdatedTempVolumeIntegralY, 0, Np2d*K2d * sizeof(double));
	PCEUpdatedBEfm2d = (double *)malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfm2d, BENe2d * BENfp2d * 3 * sizeof(double));
	memset(PCEUpdatedBEfm2d, 0, BENe2d * BENfp2d * 3 * sizeof(double));
	PCEUpdatedBEzM2d = (double *)malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzM2d, BENe2d * BENfp2d*sizeof(double));
	memset(PCEUpdatedBEzM2d, 0, BENe2d * BENfp2d * sizeof(double));
	PCEUpdatedBEfp2d = (double *)malloc(BENe2d * BENfp2d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfp2d, BENe2d * BENfp2d * 3 * sizeof(double));
	memset(PCEUpdatedBEfp2d, 0, BENe2d * BENfp2d * 3 * sizeof(double));
	PCEUpdatedBEzP2d = (double *)malloc(BENe2d * BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzP2d, BENe2d * BENfp2d*sizeof(double));
	memset(PCEUpdatedBEzP2d, 0, BENe2d * BENfp2d * sizeof(double));
	PCEUpdatedBEFluxS2d = (double *)malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEFluxS2d, BENe2d*BENfp2d*sizeof(double));
	memset(PCEUpdatedBEFluxS2d, 0, BENe2d * BENfp2d * sizeof(double));
	PCEUpdatedBEFluxM2d = (double *)malloc(BENe2d*BENfp2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEFluxM2d, BENe2d*BENfp2d*sizeof(double));
	memset(PCEUpdatedBEFluxM2d, 0, BENe2d * BENfp2d * sizeof(double));
	PCEUpdatedPCETempFacialIntegral = (double *)malloc(Np2d*K2d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedPCETempFacialIntegral, Np2d*K2d*sizeof(double));
	memset(PCEUpdatedPCETempFacialIntegral, 0, Np2d*K2d * sizeof(double));

	PCEUpdatedIEfmod = (double *)malloc(IENe2d*IENfp3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfmod, IENe2d*IENfp3d*sizeof(double));
	memset(PCEUpdatedIEfmod, 0, IENe2d*IENfp3d * sizeof(double));
	PCEUpdatedBEfmod = (double *)malloc(BENe2d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfmod, BENe2d*BENfp3d*sizeof(double));
	memset(PCEUpdatedBEfmod, 0, BENe2d*BENfp3d * sizeof(double));

	PCEUpdatedIEfm3d = (double *)malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfm3d, IENfp3d*IENe3d * 3 * sizeof(double));
	memset(PCEUpdatedIEfm3d, 0, IENfp3d*IENe3d * 3 * sizeof(double));

	PCEUpdatedIEfp3d = (double *)malloc(IENfp3d*IENe3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEfp3d, IENfp3d*IENe3d * 3 * sizeof(double));
	memset(PCEUpdatedIEfp3d, 0, IENfp3d*IENe3d * 3 * sizeof(double));

	PCEUpdatedIEFluxS3d = (double *)malloc(IENfp3d*IENe3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedIEFluxS3d, IENfp3d*IENe3d*sizeof(double));
	memset(PCEUpdatedIEFluxS3d, 0, IENfp3d*IENe3d * sizeof(double));

	PCEUpdatedBEfm3d = (double *)malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfm3d, BENe3d * BENfp3d * 3 * sizeof(double));
	memset(PCEUpdatedBEfm3d, 0, BENe3d * BENfp3d * 3 * sizeof(double));

	PCEUpdatedBEfp3d = (double *)malloc(BENe3d * BENfp3d * 3 * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEfp3d, BENe3d * BENfp3d * 3 * sizeof(double));
	memset(PCEUpdatedBEfp3d, 0, BENe3d * BENfp3d * 3 * sizeof(double));

	PCEUpdatedBEzM3d = (double *)malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzM3d, BENe3d * BENfp3d * sizeof(double));
	memset(PCEUpdatedBEzM3d, 0, BENe3d * BENfp3d * sizeof(double));

	PCEUpdatedBEzP3d = (double *)malloc(BENe3d * BENfp3d * sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEzP3d, BENe3d * BENfp3d * sizeof(double));
	memset(PCEUpdatedBEzP3d, 0, BENe3d * BENfp3d * sizeof(double));

	PCEUpdatedBEFluxS3d = (double *)malloc(BENe3d*BENfp3d*sizeof(double));
	MemoryAllocationCheck(PCEUpdatedBEFluxS3d, BENe3d*BENfp3d*sizeof(double));
	memset(PCEUpdatedBEFluxS3d, 0, BENe3d * BENfp3d * sizeof(double));

	//PCEUpdatedInitialized = "True";
}

void NdgMemory::PCEUpdatedMemoryDeAllocation(){
	free(PCEUpdatedIEfm2d); PCEUpdatedIEfm2d = NULL;
	free(PCEUpdatedIEfp2d); PCEUpdatedIEfp2d = NULL;
	free(PCEUpdatedIEFluxM2d); PCEUpdatedIEFluxM2d = NULL;
	free(PCEUpdatedIEFluxP2d); PCEUpdatedIEFluxP2d = NULL;
	free(PCEUpdatedIEFluxS2d); PCEUpdatedIEFluxS2d = NULL;
	free(PCEUpdatedERHS2d); PCEUpdatedERHS2d = NULL;
	free(PCEUpdatedVolumeIntegralX); PCEUpdatedVolumeIntegralX = NULL;
	free(PCEUpdatedTempVolumeIntegralX); PCEUpdatedTempVolumeIntegralX = NULL;
	free(PCEUpdatedVolumeIntegralY); PCEUpdatedVolumeIntegralY = NULL;
	free(PCEUpdatedTempVolumeIntegralY); PCEUpdatedTempVolumeIntegralY = NULL;
	free(PCEUpdatedBEfm2d); PCEUpdatedBEfm2d = NULL;
	free(PCEUpdatedBEzM2d); PCEUpdatedBEzM2d = NULL;
	free(PCEUpdatedBEfp2d); PCEUpdatedBEfp2d = NULL;
	free(PCEUpdatedBEzP2d); PCEUpdatedBEzP2d = NULL;
	free(PCEUpdatedBEFluxS2d); PCEUpdatedBEFluxS2d = NULL;
	free(PCEUpdatedBEFluxM2d); PCEUpdatedBEFluxM2d = NULL;
	free(PCEUpdatedPCETempFacialIntegral); PCEUpdatedPCETempFacialIntegral = NULL;
	free(PCEUpdatedIEfmod); PCEUpdatedIEfmod = NULL;
	free(PCEUpdatedBEfmod); PCEUpdatedBEfmod = NULL;
	free(PCEUpdatedIEfm3d); PCEUpdatedIEfm3d = NULL;
	free(PCEUpdatedIEfp3d); PCEUpdatedIEfp3d = NULL;
	free(PCEUpdatedBEfm3d); PCEUpdatedBEfm3d = NULL;
	free(PCEUpdatedBEfp3d); PCEUpdatedBEfp3d = NULL;
	free(PCEUpdatedBEzM3d); PCEUpdatedBEzM3d = NULL;
	free(PCEUpdatedBEzP3d); PCEUpdatedBEzP3d = NULL;
	free(PCEUpdatedBEFluxS3d); PCEUpdatedBEFluxS3d = NULL;
	//PCEUpdatedInitialized = "False";
}

/*This is for advection memory part*/
double *TempFacialIntegral = NULL, *IEfm = NULL, *IEfp = NULL, *IEFluxM = NULL, *IEFluxP = NULL, \
*IEFluxS = NULL, *ERHS = NULL, *BEfm = NULL, *BEfp = NULL, *AdvzM = NULL, *AdvzP = NULL, *BEFluxM = NULL, \
*BEFluxS = NULL, *BotEfm = NULL, *BotEfp = NULL, *BotEFluxM = NULL, *BotEFluxP = NULL, *BotEFluxS = NULL, \
*BotBEfm = NULL, *BotBEFluxM = NULL, *BotBEFluxS = NULL, *SurfBEfm = NULL, *SurfBEFluxM = NULL, \
*SurfBEFluxS = NULL, *E = NULL, *G = NULL, *H = NULL, *TempVolumeIntegral = NULL;

signed char *Status3d = NULL;
//char *AdvInitialized = "False";

void NdgMemory::AdvMemoryAllocation(int Np, int K, int Nvar, int IENfp, int IENe, int Nface, int BENfp, int BENe, int BotENfp,\
                        int BotENe, int BotBENfp, int BotBENe, int SurfBENfp, int SurfBENe){
	TempFacialIntegral = (double *)malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(TempFacialIntegral,Np*K*Nvar*sizeof(double));
	memset(TempFacialIntegral, 0, Np*K*Nvar * sizeof(double));
	IEfm = (double *)malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(IEfm,IENfp*IENe*(Nvar + 1)*sizeof(double));
	memset(IEfm, 0, IENfp*IENe*(Nvar + 1) * sizeof(double));
	IEfp = (double *)malloc(IENfp*IENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(IEfp,IENfp*IENe*(Nvar + 1)*sizeof(double));
	memset(IEfp, 0, IENfp*IENe*(Nvar + 1) * sizeof(double));
	IEFluxM = (double *)malloc(IENfp*IENe*Nvar*sizeof(double));
    MemoryAllocationCheck(IEFluxM,IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxM, 0, IENfp*IENe*Nvar * sizeof(double));
	IEFluxP = (double *)malloc(IENfp*IENe*Nvar*sizeof(double));
    MemoryAllocationCheck(IEFluxP,IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxP, 0, IENfp*IENe*Nvar * sizeof(double));
	IEFluxS = (double *)malloc(IENfp*IENe*Nvar*sizeof(double));
    MemoryAllocationCheck(IEFluxS,IENfp*IENe*Nvar*sizeof(double));
	memset(IEFluxS, 0, IENfp*IENe*Nvar * sizeof(double));
	ERHS = (double *)malloc(Np*K*Nvar*Nface*sizeof(double));
    MemoryAllocationCheck(ERHS,Np*K*Nvar*Nface*sizeof(double));
	memset(ERHS, 0, Np*K*Nvar*Nface * sizeof(double));
	BEfm = (double *)malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(BEfm,BENfp*BENe*(Nvar + 1)*sizeof(double));
	memset(BEfm, 0, BENfp*BENe*(Nvar + 1) * sizeof(double));
	BEfp = (double *)malloc(BENfp*BENe*(Nvar + 1)*sizeof(double));
    MemoryAllocationCheck(BEfp,BENfp*BENe*(Nvar + 1)*sizeof(double));
	memset(BEfp, 0, BENfp*BENe*(Nvar + 1) * sizeof(double));
	AdvzM = (double *)malloc(BENfp*BENe*sizeof(double));
    MemoryAllocationCheck(AdvzM,BENfp*BENe*sizeof(double));
	memset(AdvzM, 0, BENfp*BENe * sizeof(double));
	AdvzP = (double *)malloc(BENfp*BENe*sizeof(double));
    MemoryAllocationCheck(AdvzP,BENfp*BENe*sizeof(double));
	memset(AdvzP, 0, BENfp*BENe * sizeof(double));
	BEFluxM = (double *)malloc(BENfp*BENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BEFluxM,BENfp*BENe*Nvar*sizeof(double));
	memset(BEFluxM, 0, BENfp*BENe*Nvar * sizeof(double));
	BEFluxS = (double *)malloc(BENfp*BENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BEFluxS,BENfp*BENe*Nvar*sizeof(double));
	memset(BEFluxS, 0, BENfp*BENe*Nvar * sizeof(double));
	BotEfm = (double *)malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(BotEfm,BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	memset(BotEfm, 0, BotENfp*BotENe*(Nvar + 2) * sizeof(double));
	BotEfp = (double *)malloc(BotENfp*BotENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(BotEfp,BotENfp*BotENe*(Nvar + 2)*sizeof(double));
	memset(BotEfp, 0, BotENfp*BotENe*(Nvar + 2) * sizeof(double));
	BotEFluxM = (double *)malloc(BotENfp*BotENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotEFluxM,BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxM, 0, BotENfp*BotENe*Nvar * sizeof(double));
	BotEFluxP = (double *)malloc(BotENfp*BotENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotEFluxP,BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxP, 0, BotENfp*BotENe*Nvar * sizeof(double));
	BotEFluxS = (double *)malloc(BotENfp*BotENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotEFluxS,BotENfp*BotENe*Nvar*sizeof(double));
	memset(BotEFluxS, 0, BotENfp*BotENe*Nvar * sizeof(double));
	BotBEfm = (double *)malloc(BotBENfp*BotBENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(BotBEfm,BotBENfp*BotBENe*(Nvar + 2)*sizeof(double));
	memset(BotBEfm, 0, BotBENfp*BotBENe*(Nvar + 2) * sizeof(double));
	BotBEFluxM = (double *)malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotBEFluxM,BotBENfp*BotBENe*Nvar*sizeof(double));
	memset(BotBEFluxM, 0, BotBENfp*BotBENe*Nvar * sizeof(double));
	BotBEFluxS = (double *)malloc(BotBENfp*BotBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(BotBEFluxS,BotBENfp*BotBENe*Nvar*sizeof(double));
	memset(BotBEFluxS, 0, BotBENfp*BotBENe*Nvar * sizeof(double));
	SurfBEfm = (double *)malloc(SurfBENfp*SurfBENe*(Nvar + 2)*sizeof(double));
    MemoryAllocationCheck(SurfBEfm,SurfBENfp*SurfBENe*(Nvar + 2)*sizeof(double));
	memset(SurfBEfm, 0, SurfBENfp*SurfBENe*(Nvar + 2) * sizeof(double));
	SurfBEFluxM = (double *)malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(SurfBEFluxM,SurfBENfp*SurfBENe*Nvar*sizeof(double));
	memset(SurfBEFluxM, 0, SurfBENfp*SurfBENe*Nvar * sizeof(double));
	SurfBEFluxS = (double *)malloc(SurfBENfp*SurfBENe*Nvar*sizeof(double));
    MemoryAllocationCheck(SurfBEFluxS,SurfBENfp*SurfBENe*Nvar*sizeof(double));
	memset(SurfBEFluxS, 0, SurfBENfp*SurfBENe*Nvar * sizeof(double));
	/*Allocate memory for E, G and H, and calculate these volume flux term*/
	E = (double *)malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(E,Np*K*Nvar*sizeof(double));
	memset(E, 0, Np*K*Nvar * sizeof(double));
	G = (double *)malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(G,Np*K*Nvar*sizeof(double));
	memset(G, 0, Np*K*Nvar * sizeof(double));
	H = (double *)malloc(Np*K*Nvar*sizeof(double));
    MemoryAllocationCheck(H,Np*K*Nvar*sizeof(double));
	memset(H, 0, Np*K*Nvar * sizeof(double));
	TempVolumeIntegral = (double *)malloc(Np*K*sizeof(double));
    MemoryAllocationCheck(TempVolumeIntegral,Np*K*sizeof(double));
	memset(TempVolumeIntegral, 0, Np*K * sizeof(double));

	Status3d = (signed char *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Status3d, Np*K * sizeof(double));
	memset(Status3d, 0, Np*K * sizeof(double));


	//AdvInitialized = "True";
	//cout<<"Adv is OK"<<endl;
}

void NdgMemory::AdvMemoryDeAllocation(){
	//free(TempFacialIntegral); TempFacialIntegral = NULL;
	free(IEfm); IEfm = NULL;
	free(IEfp); IEfp = NULL;
	free(IEFluxM); IEFluxM = NULL;
	free(IEFluxP); IEFluxP = NULL;
	free(IEFluxS); IEFluxS = NULL;
	free(ERHS); ERHS = NULL;
	free(BEfm); BEfm = NULL;
	free(BEfp); BEfp = NULL;
	free(AdvzM); AdvzM = NULL;
	free(AdvzP); AdvzP = NULL;
	free(BEFluxM); BEFluxM = NULL;
	free(BEFluxS); BEFluxS = NULL;
	free(BotEfm); BotEfm = NULL;
	free(BotEfp); BotEfp = NULL;
	free(BotEFluxM); BotEFluxM = NULL;
	free(BotEFluxP); BotEFluxP = NULL;
	free(BotEFluxS); BotEFluxS = NULL;
	free(BotBEfm); BotBEfm = NULL;
	free(BotBEFluxM); BotBEFluxM = NULL;
	free(BotBEFluxS); BotBEFluxS = NULL;
	free(SurfBEfm); BotBEFluxS = NULL;
	free(SurfBEFluxM); SurfBEFluxM = NULL;
	free(SurfBEFluxS); SurfBEFluxS = NULL;
	free(E); E = NULL;
	free(G); G = NULL;
	free(H); H = NULL;
	free(TempVolumeIntegral); TempVolumeIntegral = NULL;
	free(Status3d); Status3d = NULL;
	//AdvInitialized = "False";
}


/*This is for horizontal diffusion memory part*/
double *HorDiffnv = NULL, *HorDiffvariable = NULL, *HorDiffBEfp = NULL, *HorDiffzM = NULL, \
*HorDiffzP = NULL, *HorDiffTempBEfp = NULL, *HorDiffTempBEfm = NULL, *HorDiffAVx = NULL, \
*HorDiffAVy = NULL, *HorDiffVx = NULL, *HorDiffTempVx = NULL, *HorDiffVy = NULL, *HorDiffTempVy = NULL, \
*HorDiffIEfm = NULL, *HorDiffAVIEfm = NULL, *HorDiffIEfp = NULL, *HorDiffAVIEfp = NULL, *HorDiffIEFluxM = NULL, \
*HorDiffIEFluxP = NULL, *HorDiffBEfm = NULL, *HorDiffIEFluxS = NULL, *HorDiffAVBEfm = NULL, *HorDiffBEFluxM = NULL, \
*HorDiffBEFluxS = NULL, *HorDiffERHSX = NULL, *HorDiffERHSY = NULL, *HorDiffLocalPrimitiveDiffTermX = NULL, \
*HorDiffLocalPrimitiveDiffTermY = NULL, *HorDiffLPDTIEfm = NULL, *HorDiffLPDTIEfp = NULL, *HorDiffLPDTBEfm = NULL, \
*HorDiffTempFacialIntegralX = NULL, *HorDiffTempFacialIntegralY = NULL, *HorDiffInnerEdgeTau = NULL, \
*HorDiffBoundaryEdgeTau = NULL, *HorDiffIEnvfm = NULL, *HorDiffIEnvfp = NULL, *HorDiffBEnvfm = NULL, \
*u_u = NULL, *v_v = NULL, *rx_Dr_u = NULL, *rx_Dr_v = NULL, *sx_Ds_u = NULL, *sx_Ds_v = NULL, \
*ry_Dr_u = NULL, *ry_Dr_v = NULL, *sy_Ds_u = NULL, *sy_Ds_v = NULL, \
*Hrms = NULL, *WaveNumber_ = NULL, *nv_c = NULL, *nv_w = NULL, *nv_cw = NULL, \
*BBE3d = NULL;

//char *HorDiffInitialized = "False";

void NdgMemory::HorizDiffMemoryAllocation(/*NdgMeshType type,*/ int Np, int K, int Nvar, int tempNface, int BENfp, int BENe, int IENfp, int IENe){
	int Nfield;
	int Nface;
	//if (type == Two){
	//	/*For 2d shallow water problem, no horizontal diffusion terms are included in the governing equation for water depth $H$*/
	//	Nfield = Nvar - 1;
	//	/*For 2d shallow water problem, the face number is equal to TempNface, since there is no surface edge and bottom edge*/
	//	Nface = tempNface;
	//}
	//else if (type == Three){
		Nfield = Nvar;
		/*For 3d shallow water problem, the face number is equal to TempNface - 2, since there
		* is surface edge and bottom edge is not considered for horizontal diffusion term*/
		Nface = tempNface - 2;
	//}
	HorDiffzM = (double *)malloc(BENfp*BENe*sizeof(double));
    MemoryAllocationCheck(HorDiffzM,BENfp*BENe*sizeof(double));
	memset(HorDiffzM, 0, BENfp*BENe * sizeof(double));
	HorDiffzP = (double *)malloc(BENfp*BENe*sizeof(double));
    MemoryAllocationCheck(HorDiffzP,BENfp*BENe*sizeof(double));
	memset(HorDiffzP, 0, BENfp*BENe * sizeof(double));
	HorDiffnv = (double *)malloc(Np*K*sizeof(double));
    MemoryAllocationCheck(HorDiffnv,Np*K*sizeof(double));
	memset(HorDiffnv, 0, Np*K * sizeof(double));
	HorDiffTempBEfp = (double *)malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
    MemoryAllocationCheck(HorDiffTempBEfp,BENe*BENfp*(Nfield + 1)*sizeof(double));
	memset(HorDiffTempBEfp, 0, BENe*BENfp*(Nfield + 1) * sizeof(double));
	HorDiffTempBEfm = (double *)malloc(BENe*BENfp*(Nfield + 1)*sizeof(double));
    MemoryAllocationCheck(HorDiffTempBEfm,BENe*BENfp*(Nfield + 1)*sizeof(double));
	memset(HorDiffTempBEfm, 0, BENe*BENfp*(Nfield + 1) * sizeof(double));
	/*Allocate memory for the original variable over boundary edge*/
	HorDiffBEfp = (double *)malloc(BENfp*BENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEfp,BENfp*BENe*Nfield*sizeof(double));
	memset(HorDiffBEfp, 0, BENfp*BENe*Nfield * sizeof(double));
	/*Allocate memory for the local face value at boundary edge*/
	HorDiffBEfm = (double *)malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEfm,BENe*BENfp*Nfield*sizeof(double));
	memset(HorDiffBEfm, 0, BENfp*BENe*Nfield * sizeof(double));
	/*Allocate memory for the original variable $u,v$ and $\theta$*/
	HorDiffvariable = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffvariable,Np*K*Nfield*sizeof(double));
	memset(HorDiffvariable, 0, Np*K*Nfield * sizeof(double));
	/*Allocate memory for auxiallary variable $q_x=\frac{\partial u(v,\theta)}{\partial x}$*/
	HorDiffAVx = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVx,Np*K*Nfield*sizeof(double));
	memset(HorDiffAVx, 0, Np*K*Nfield * sizeof(double));
	/*Allocate memory for auxiallary variable $q_y=\frac{\partial u(v,\theta)}{\partial y}$*/
	HorDiffAVy = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVy,Np*K*Nfield*sizeof(double));
	memset(HorDiffAVy, 0, Np*K*Nfield * sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{x1}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffVx = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffVx,Np*K*Nfield*sizeof(double));
	memset(HorDiffVx, 0, Np*K*Nfield * sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{x2}=\bold{r_x}\cdot (D_r*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffTempVx = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempVx,Np*K*Nfield*sizeof(double));
	memset(HorDiffTempVx, 0, Np*K*Nfield * sizeof(double));
	/*Allocate memory for part of the volumn integral of auxiallary variable $q_{y1}=\bold{r_y}\cdot (D_r*u(v,\theta))$,
	and for part of the volumn integral part of the second order operator $\frac{\partial q}{\partial x}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$ */
	HorDiffVy = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffVy,Np*K*Nfield*sizeof(double));
	memset(HorDiffVy, 0, Np*K*Nfield * sizeof(double));
	/*Allocate memory for the rest part of the volumn integral of auxiallary variable $q_{y2}=\bold{s_y}\cdot (D_s*u(v,\theta))$,
	and for the rest part of the volumn integral part of the second order operator $\frac{\partial q}{\partial y}$. Here, $q$ can
	be both $\frac{\partial u(v,\theta)}{\partial x}$ and $\frac{\partial u(v,\theta)}{\partial y}$*/
	HorDiffTempVy = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempVy,Np*K*Nfield*sizeof(double));
	memset(HorDiffTempVy, 0, Np*K*Nfield * sizeof(double));
	/*Allocate memory for the local face value at inner edge*/
	HorDiffIEfm = (double *)malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEfm,IENe*IENfp*Nfield*sizeof(double));
	memset(HorDiffIEfm, 0, Np*K * sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at inner edge*/
	HorDiffAVIEfm = (double *)malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVIEfm,IENe*IENfp*Nfield*sizeof(double));
	memset(HorDiffAVIEfm, 0, IENe*IENfp*Nfield * sizeof(double));
	/*Allocate memory for the adjacent face value at inner edge*/
	HorDiffIEfp = (double *)malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEfp,IENe*IENfp*Nfield*sizeof(double));
	memset(HorDiffIEfp, 0, IENe*IENfp*Nfield * sizeof(double));
	/*Allocate memory for the adjacent face value of the auxialary variable at inner edge*/
	HorDiffAVIEfp = (double *)malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVIEfp,IENe*IENfp*Nfield*sizeof(double));
	memset(HorDiffAVIEfp, 0, IENe*IENfp*Nfield * sizeof(double));
	/*Allocate memory for the local flux term at inner edge*/
	HorDiffIEFluxM = (double *)malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEFluxM,IENe*IENfp*Nfield*sizeof(double));
	memset(HorDiffIEFluxM, 0, IENe*IENfp*Nfield * sizeof(double));
	/*Allocate memory for the adjacent flux term at inner edge*/
	HorDiffIEFluxP = (double *)malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEFluxP,IENe*IENfp*Nfield*sizeof(double));
	memset(HorDiffIEFluxP, 0, IENe*IENfp*Nfield * sizeof(double));
	/*Allocate memory for the numerical flux term at inner edge*/
	HorDiffIEFluxS = (double *)malloc(IENe*IENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffIEFluxS,IENe*IENfp*Nfield*sizeof(double));
	memset(HorDiffIEFluxS, 0, IENe*IENfp*Nfield * sizeof(double));
	/*Allocate memory for the local face value of the auxialary variable at boundary edge*/
	HorDiffAVBEfm = (double *)malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffAVBEfm,BENe*BENfp*Nfield*sizeof(double));
	memset(HorDiffAVBEfm, 0, BENe*BENfp*Nfield * sizeof(double));
	/*Allocate memory for the local flux term at boundary edge*/
	HorDiffBEFluxM = (double *)malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEFluxM,BENe*BENfp*Nfield*sizeof(double));
	memset(HorDiffBEFluxM, 0, BENe*BENfp*Nfield * sizeof(double));
	/*Allocate memory for the numerical flux term at boundary edge*/
	HorDiffBEFluxS = (double *)malloc(BENe*BENfp*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffBEFluxS,BENe*BENfp*Nfield*sizeof(double));
	memset(HorDiffBEFluxS, 0, BENe*BENfp*Nfield * sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in x direction*/
	HorDiffERHSX = (double *)malloc(Np*K*Nfield*Nface*sizeof(double));
    MemoryAllocationCheck(HorDiffERHSX,Np*K*Nfield*Nface*sizeof(double));
	memset(HorDiffERHSX, 0, Np*K*Nfield*Nface * sizeof(double));
	/*Allocate memory for interior edge contribution to right hand side in y direction*/
	HorDiffERHSY = (double *)malloc(Np*K*Nfield*Nface*sizeof(double));
    MemoryAllocationCheck(HorDiffERHSY,Np*K*Nfield*Nface*sizeof(double));
	memset(HorDiffERHSY, 0, Np*K*Nfield*Nface * sizeof(double));

	/*Allocate memory for local primitive diffusion term in both x and y direction, $\nu\nabla u(v,\theta)$*/
	/*This part is used because we need to extract the trace of the local derivative operator to compute the numerical flux*/
	HorDiffLocalPrimitiveDiffTermX = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLocalPrimitiveDiffTermX,Np*K*Nfield*sizeof(double));
	memset(HorDiffLocalPrimitiveDiffTermX, 0, Np*K*Nfield * sizeof(double));
	HorDiffLocalPrimitiveDiffTermY = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLocalPrimitiveDiffTermY,Np*K*Nfield*sizeof(double));
	memset(HorDiffLocalPrimitiveDiffTermY, 0, Np*K*Nfield * sizeof(double));
	HorDiffLPDTIEfm = (double *)malloc(IENfp*IENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLPDTIEfm,IENfp*IENe*Nfield*sizeof(double));
	memset(HorDiffLPDTIEfm, 0, IENfp*IENe*Nfield * sizeof(double));
	HorDiffLPDTIEfp = (double *)malloc(IENfp*IENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLPDTIEfp,IENfp*IENe*Nfield*sizeof(double));
	memset(HorDiffLPDTIEfp, 0, IENfp*IENe*Nfield * sizeof(double));
	HorDiffLPDTBEfm = (double *)malloc(BENfp*BENe*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffLPDTBEfm,BENfp*BENe*Nfield*sizeof(double));
	memset(HorDiffLPDTBEfm, 0, BENfp*BENe*Nfield * sizeof(double));
	HorDiffTempFacialIntegralX = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempFacialIntegralX,Np*K*Nfield*sizeof(double));
	memset(HorDiffTempFacialIntegralX, 0, Np*K*Nfield * sizeof(double));
	HorDiffTempFacialIntegralY = (double *)malloc(Np*K*Nfield*sizeof(double));
    MemoryAllocationCheck(HorDiffTempFacialIntegralY,Np*K*Nfield*sizeof(double));
	memset(HorDiffTempFacialIntegralY, 0, Np*K*Nfield * sizeof(double));
	HorDiffInnerEdgeTau = (double *)malloc(IENe*IENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffInnerEdgeTau,IENe*IENfp*sizeof(double));
	memset(HorDiffInnerEdgeTau, 0, IENe*IENfp * sizeof(double));
	HorDiffBoundaryEdgeTau = (double *)malloc(BENe*BENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffBoundaryEdgeTau,BENe*BENfp*sizeof(double));
	memset(HorDiffBoundaryEdgeTau, 0, BENe*BENfp * sizeof(double));
	HorDiffIEnvfm = (double *)malloc(IENe*IENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffIEnvfm,IENe*IENfp*sizeof(double));
	memset(HorDiffIEnvfm, 0, IENe*IENfp * sizeof(double));
	HorDiffIEnvfp = (double *)malloc(IENe*IENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffIEnvfp,IENe*IENfp*sizeof(double));
	memset(HorDiffIEnvfp, 0, IENe*IENfp * sizeof(double));
	HorDiffBEnvfm = (double *)malloc(BENe*BENfp*sizeof(double));
    MemoryAllocationCheck(HorDiffBEnvfm,BENe*BENfp*sizeof(double));
	memset(HorDiffBEnvfm, 0, BENe*BENfp * sizeof(double));
	//HorDiffInitialized = "True";
	u_u = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(u_u, Np * K * sizeof(double));
	memset(u_u, 0, Np*K * sizeof(double));
	v_v = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(v_v, Np * K * sizeof(double));
	memset(v_v, 0, Np*K * sizeof(double));
	rx_Dr_u = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(rx_Dr_u, Np * K * sizeof(double));
	memset(rx_Dr_u, 0, Np*K * sizeof(double));
	rx_Dr_v = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(rx_Dr_v, Np * K * sizeof(double));
	memset(rx_Dr_v, 0, Np*K * sizeof(double));
	sx_Ds_u = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(sx_Ds_u, Np * K * sizeof(double));
	memset(sx_Ds_u, 0, Np*K * sizeof(double));
	sx_Ds_v = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(sx_Ds_v, Np * K * sizeof(double));
	memset(sx_Ds_v, 0, Np*K * sizeof(double));
	ry_Dr_u = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(ry_Dr_u, Np * K * sizeof(double));
	memset(ry_Dr_u, 0, Np*K * sizeof(double));
	ry_Dr_v = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(ry_Dr_v, Np * K * sizeof(double));
	memset(ry_Dr_v, 0, Np*K * sizeof(double));
	sy_Ds_u = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(sy_Ds_u, Np * K * sizeof(double));
	memset(sy_Ds_u, 0, Np*K * sizeof(double));
	sy_Ds_v = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(sy_Ds_v, Np * K * sizeof(double));
	memset(sy_Ds_v, 0, Np*K * sizeof(double));
	Hrms = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Hrms, Np * K * sizeof(double));
	memset(Hrms, 0, Np*K * sizeof(double));
	WaveNumber_ = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(WaveNumber_, Np * K * sizeof(double));
	memset(WaveNumber_, 0, Np*K * sizeof(double));
	nv_c = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(nv_c, Np * K * sizeof(double));
	memset(nv_c, 0, Np*K * sizeof(double));
	nv_w = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(nv_w, Np * K * sizeof(double));
	memset(nv_w, 0, Np*K * sizeof(double));
	nv_cw = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(nv_cw, Np * K * sizeof(double));
	memset(nv_cw, 0, Np*K * sizeof(double));
	BBE3d = (double *)malloc(Np * K * 2 * sizeof(double));
	MemoryAllocationCheck(BBE3d, Np * K * 2 * sizeof(double));
	memset(BBE3d, 0, Np*K*2 * sizeof(double));
}

void NdgMemory::HorizDiffMemoryDeAllocation()
{
	free(HorDiffnv); HorDiffnv = NULL;
	free(HorDiffvariable); HorDiffvariable = NULL;
	free(HorDiffBEfp); HorDiffBEfp = NULL;
	free(HorDiffzM); HorDiffzM = NULL;
	free(HorDiffzP); HorDiffzP = NULL;
	free(HorDiffTempBEfp); HorDiffTempBEfp = NULL;
	free(HorDiffTempBEfm); HorDiffTempBEfm = NULL;
	free(HorDiffAVx); HorDiffAVx = NULL;
	free(HorDiffAVy); HorDiffAVy = NULL;
	free(HorDiffVx); HorDiffVx = NULL;
	free(HorDiffTempVx); HorDiffVx = NULL;
	free(HorDiffVy); HorDiffVy = NULL;
	free(HorDiffTempVy); HorDiffTempVy = NULL;
	free(HorDiffIEfm); HorDiffIEfm = NULL;
	free(HorDiffAVIEfm); HorDiffAVIEfm = NULL;
	free(HorDiffIEfp); HorDiffIEfp = NULL;
	free(HorDiffAVIEfp); HorDiffAVIEfp = NULL;
	free(HorDiffIEFluxM); HorDiffIEFluxM = NULL;
	free(HorDiffIEFluxP); HorDiffIEFluxP = NULL;
	free(HorDiffBEfm); HorDiffBEfm = NULL;
	free(HorDiffIEFluxS); HorDiffIEFluxS = NULL;
	free(HorDiffAVBEfm); HorDiffAVBEfm = NULL;
	free(HorDiffBEFluxM); HorDiffBEFluxM = NULL;
	free(HorDiffBEFluxS); HorDiffBEFluxS = NULL;
	free(HorDiffERHSX); HorDiffERHSX = NULL;
	free(HorDiffERHSY); HorDiffERHSY = NULL;
	free(HorDiffLocalPrimitiveDiffTermX); HorDiffLocalPrimitiveDiffTermX = NULL;
	free(HorDiffLocalPrimitiveDiffTermY); HorDiffLocalPrimitiveDiffTermY = NULL;
	free(HorDiffLPDTIEfm); HorDiffLPDTIEfm = NULL;
	free(HorDiffLPDTIEfp); HorDiffLPDTIEfp = NULL;
	free(HorDiffLPDTBEfm); HorDiffLPDTBEfm = NULL;
	free(HorDiffTempFacialIntegralX); HorDiffTempFacialIntegralX = NULL;
	free(HorDiffTempFacialIntegralY); HorDiffTempFacialIntegralY = NULL;
	free(HorDiffInnerEdgeTau); HorDiffInnerEdgeTau = NULL;
	free(HorDiffBoundaryEdgeTau); HorDiffBoundaryEdgeTau = NULL;
	free(HorDiffIEnvfm); HorDiffIEnvfm = NULL;
	free(HorDiffIEnvfp); HorDiffIEnvfp = NULL;
	free(HorDiffBEnvfm); HorDiffBEnvfm = NULL;

	free(u_u); u_u = NULL;
	free(v_v); v_v = NULL;
	free(rx_Dr_u); rx_Dr_u = NULL;
	free(rx_Dr_v); rx_Dr_v = NULL;
	free(sx_Ds_u); sx_Ds_u = NULL;
	free(sx_Ds_v); sx_Ds_v = NULL;
	free(ry_Dr_u); ry_Dr_u = NULL;
	free(ry_Dr_v); ry_Dr_v = NULL;
	free(sy_Ds_u); sy_Ds_u = NULL;
	free(sy_Ds_v); sy_Ds_v = NULL;

	free(Hrms); Hrms = NULL;
	free(WaveNumber_); WaveNumber_ = NULL;
	free(nv_c); nv_c = NULL;
	free(nv_w); nv_w = NULL;
	free(nv_cw); nv_cw = NULL;
	free(BBE3d); BBE3d = NULL;
	//HorDiffInitialized = "False";
}

	/*This is for wave radiation and wave surface roller*/
double *SIN_DIR = NULL, *COS_DIR = NULL, *WaveNumber = NULL, *RSIEfm = NULL, *RSIEfp = NULL, \
*RSIEFluxM = NULL, *RSIEFluxP = NULL, *RSIEFluxS = NULL, *RSERHS = NULL, *RSBotEfm = NULL, \
*RSBotEfp = NULL, *RSBotEFluxM = NULL, *RSBotEFluxP = NULL, *RSBotEFluxS = NULL, \
*RSBEfm = NULL, *RSBEfp = NULL, *RSBEFluxM = NULL, *RSBEFluxS = NULL, \
*WaveEnergy = NULL, *KD = NULL, *FSS = NULL, *FCS = NULL, *FSC = NULL, *FCC = NULL, *M2 = NULL, *N1 = NULL, \
*N2 = NULL, *Ar = NULL, *Rz = NULL, *R = NULL, *CFF1 = NULL, *CFF2 = NULL, *CFF3 = NULL, *CFF4 = NULL, \
*Roller1 = NULL, *Roller2 = NULL, *Roller3 = NULL, *H_Radiation1 = NULL, *H_Radiation2 = NULL, \
*H_Radiation3 = NULL, *V_Radiation = NULL, *RS_E = NULL, *RS_G = NULL, *Hx = NULL, *Hy = NULL, *tempRHSx = NULL, \
*tempRHSy = NULL, *midtempRHSx = NULL, *midtempRHSy = NULL, *VtempRHSx = NULL, *VtempRHSy = NULL, \
*midERHS = NULL, *RSBotBEfm = NULL, *RSBotBEFluxM = NULL, *RSBotBEFluxS = NULL, *RSSurfBEfm = NULL, \
*RSSurfBEFluxM = NULL, *RSSurfBEFluxS = NULL, *CFF5 = NULL, *CFF6 = NULL, *CFF7 = NULL, *CFF8 = NULL, *wave_C = NULL,\
*Rzn = NULL, *Rz1 = NULL;

void NdgMemory::RollerWaveRadiationMemoryAllocation(int Np, int K, int IENfp, int IENe, int BENfp, int BENe, int BotENfp, int BotENe, int Nface, int BotBENfp, int BotBENe, int SurfBENfp, int SurfBENe,int Np2d,int K2d) {
	SIN_DIR = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(SIN_DIR, Np*K * sizeof(double));
	memset(SIN_DIR, 0, Np*K * sizeof(double));
	COS_DIR = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(COS_DIR, Np*K * sizeof(double));
	memset(COS_DIR, 0, Np*K * sizeof(double));
	WaveNumber = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(WaveNumber, Np*K * sizeof(double));
	memset(WaveNumber, 0, Np*K * sizeof(double));
	RSIEfm = (double *)malloc(IENfp*IENe * 3 * sizeof(double));
	MemoryAllocationCheck(RSIEfm, IENfp*IENe * 3 * sizeof(double));
	memset(RSIEfm, 0, IENfp*IENe * 3 * sizeof(double));
	RSIEfp = (double *)malloc(IENfp*IENe * 3 * sizeof(double));
	MemoryAllocationCheck(RSIEfp, IENfp*IENe * 3 * sizeof(double));
	memset(RSIEfp, 0, IENfp*IENe * 3 * sizeof(double));
	RSBEfm = (double *)malloc(BENfp*BENe * 3 * sizeof(double));
	MemoryAllocationCheck(RSBEfm, BENfp*BENe * 3 * sizeof(double));
	memset(RSBEfm, 0, BENfp*BENe * 3 * sizeof(double));
	RSBEfp = (double *)malloc(BENfp*BENe * 3 * sizeof(double));
	MemoryAllocationCheck(RSBEfp, BENfp*BENe * 3 * sizeof(double));
	memset(RSBEfp, 0, BENfp*BENe * 3 * sizeof(double));

	RSIEFluxM = (double *)malloc(IENfp*IENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSIEFluxM, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxM, 0, IENfp*IENe * 2 * sizeof(double));
	RSIEFluxP = (double *)malloc(IENfp*IENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSIEFluxP, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxP, 0, IENfp*IENe * 2 * sizeof(double));
	RSIEFluxS = (double *)malloc(IENfp*IENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSIEFluxS, IENfp*IENe * 2 * sizeof(double));
	memset(RSIEFluxS, 0, IENfp*IENe * 2 * sizeof(double));
	RSBEFluxM = (double *)malloc(BENfp*BENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBEFluxM, BENfp*BENe * 2 * sizeof(double));
	memset(RSBEFluxM, 0, BENfp*BENe * 2 * sizeof(double));
	RSBEFluxS = (double *)malloc(BENfp*BENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBEFluxS, BENfp*BENe * 2 * sizeof(double));
	memset(RSBEFluxS, 0, BENfp*BENe * 2 * sizeof(double));

	RSBotEfm = (double *)malloc(BotENfp*BotENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotEfm, BotENfp*BotENe * 2 * sizeof(double));
	memset(RSBotEfm, 0, BotENfp*BotENe * 2 * sizeof(double));
	RSBotEfp = (double *)malloc(BotENfp*BotENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotEfp, BotENfp*BotENe * 2 * sizeof(double));
	memset(RSBotEfp, 0, BotENfp*BotENe * 2 * sizeof(double));
	RSBotEFluxM = (double *)malloc(BotENfp*BotENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotEFluxM, BotENfp*BotENe * 2 * sizeof(double));
	memset(RSBotEFluxM, 0, BotENfp*BotENe * 2 * sizeof(double));
	RSBotEFluxP = (double *)malloc(BotENfp*BotENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotEFluxP, BotENfp*BotENe * 2 * sizeof(double));
	memset(RSBotEFluxP, 0, BotENfp*BotENe * 2 * sizeof(double));
	RSBotEFluxS = (double *)malloc(BotENfp*BotENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotEFluxS, BotENfp*BotENe * 2 * sizeof(double));
	memset(RSBotEFluxS, 0, BotENfp*BotENe * 2 * sizeof(double));

	WaveEnergy = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(WaveEnergy, Np*K * sizeof(double));
	memset(WaveEnergy, 0, Np*K * sizeof(double));
	KD = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(KD, Np*K * sizeof(double));
	memset(KD, 0, Np*K * sizeof(double));
	FSS = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(FSS, Np*K * sizeof(double));
	memset(FSS, 0, Np*K * sizeof(double));
	FCS = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(FCS, Np*K * sizeof(double));
	memset(FCS, 0, Np*K * sizeof(double));
	FSC = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(FSC, Np*K * sizeof(double));
	memset(FSC, 0, Np*K * sizeof(double));
	FCC = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(FCC, Np*K * sizeof(double));
	memset(FCC, 0, Np*K * sizeof(double));
	M2 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(M2, Np*K * sizeof(double));
	memset(M2, 0, Np*K * sizeof(double));
	N1 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(N1, Np*K * sizeof(double));
	memset(N1, 0, Np*K * sizeof(double));
	N2 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(N2, Np*K * sizeof(double));
	memset(N2, 0, Np*K * sizeof(double));
	Ar = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Ar, Np*K * sizeof(double));
	memset(Ar, 0, Np*K * sizeof(double));
	Rz = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Rz, Np*K * sizeof(double));
	memset(Rz, 0, Np*K * sizeof(double));
	R = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(R, Np*K * sizeof(double));
	memset(R, 0, Np*K * sizeof(double));
	wave_C = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(wave_C, Np*K * sizeof(double));
	memset(wave_C, 0, Np*K * sizeof(double));
	Rzn = (double *)malloc(Np2d * K2d * sizeof(double));
	MemoryAllocationCheck(Rzn, Np2d*K2d * sizeof(double));
	memset(Rzn, 0, Np2d*K2d * sizeof(double));
	Rz1 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Rz1, Np*K * sizeof(double));
	memset(Rz1, 0, Np*K * sizeof(double));
	CFF1 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF1, Np*K * sizeof(double));
	memset(CFF1, 0, Np*K * sizeof(double));
	CFF2 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF2, Np*K * sizeof(double));
	memset(CFF2, 0, Np*K * sizeof(double));
	CFF3 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF3, Np*K * sizeof(double));
	memset(CFF3, 0, Np*K * sizeof(double));
	CFF4 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF4, Np*K * sizeof(double));
	memset(CFF4, 0, Np*K * sizeof(double));
	CFF5 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF5, Np*K * sizeof(double));
	memset(CFF5, 0, Np*K * sizeof(double));
	CFF6 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF6, Np*K * sizeof(double));
	memset(CFF6, 0, Np*K * sizeof(double));
	CFF7 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF7, Np*K * sizeof(double));
	memset(CFF7, 0, Np*K * sizeof(double));
	CFF8 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(CFF8, Np*K * sizeof(double));
	memset(CFF8, 0, Np*K * sizeof(double));
	Roller1 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Roller1, Np*K * sizeof(double));
	memset(Roller1, 0, Np*K * sizeof(double));
	Roller2 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Roller2, Np*K * sizeof(double));
	memset(Roller2, 0, Np*K * sizeof(double));
	Roller3 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(Roller3, Np*K * sizeof(double));
	memset(Roller3, 0, Np*K * sizeof(double));
	H_Radiation1 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(H_Radiation1, Np*K * sizeof(double));
	memset(H_Radiation1, 0, Np*K * sizeof(double));
	H_Radiation2 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(H_Radiation2, Np*K * sizeof(double));
	memset(H_Radiation2, 0, Np*K * sizeof(double));
	H_Radiation3 = (double *)malloc(Np * K * sizeof(double));
	MemoryAllocationCheck(H_Radiation3, Np*K * sizeof(double));
	memset(H_Radiation3, 0, Np*K * sizeof(double));
	V_Radiation = (double *)malloc(2 * Np * K * sizeof(double));
	MemoryAllocationCheck(V_Radiation, 2 * Np*K * sizeof(double));
	memset(V_Radiation, 0, 2 * Np*K * sizeof(double));
	RSERHS = (double *)malloc(Np*K * 2 * Nface * sizeof(double));
	MemoryAllocationCheck(RSERHS, Np*K * 2 * Nface * sizeof(double));
	memset(RSERHS, 0, Np*K * 2 * Nface * sizeof(double));
	RS_E = (double *)malloc(Np*K * 2 * sizeof(double));
	MemoryAllocationCheck(RS_E, Np*K * 2 * sizeof(double));
	memset(RS_E, 0, Np*K * 2 * sizeof(double));
	RS_G = (double *)malloc(Np*K * 4 * sizeof(double));
	MemoryAllocationCheck(RS_G, Np*K * 4 * sizeof(double));
	memset(RS_G, 0, Np*K * 4 * sizeof(double));
	Hx = (double *)malloc(Np*K * 2 * sizeof(double));
	MemoryAllocationCheck(Hx, Np*K * 2 * sizeof(double));
	memset(Hx, 0, Np*K * 2 * sizeof(double));
	Hy = (double *)malloc(Np*K * 2 * sizeof(double));
	MemoryAllocationCheck(Hy, Np*K * 2 * sizeof(double));
	memset(Hy, 0, Np*K * 2 * sizeof(double));
	tempRHSx = (double *)malloc(Np*K * 2 * sizeof(double));
	MemoryAllocationCheck(tempRHSx, Np*K * 2 * sizeof(double));
	memset(tempRHSx, 0, Np*K * 2 * sizeof(double));
	tempRHSy = (double *)malloc(Np*K * 2 * sizeof(double));
	MemoryAllocationCheck(tempRHSy, Np*K * 2 * sizeof(double));
	memset(tempRHSy, 0, Np*K * 2 * sizeof(double));
	midtempRHSx = (double *)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(midtempRHSx, Np*K * sizeof(double));
	memset(midtempRHSx, 0, Np*K * sizeof(double));
	midtempRHSy = (double *)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(midtempRHSy, Np*K * sizeof(double));
	memset(midtempRHSy, 0, Np*K * sizeof(double));
	VtempRHSx = (double *)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(VtempRHSx, Np*K * sizeof(double));
	memset(VtempRHSx, 0, Np*K * sizeof(double));
	VtempRHSy = (double *)malloc(Np*K * sizeof(double));
	MemoryAllocationCheck(VtempRHSy, Np*K * sizeof(double));
	memset(VtempRHSy, 0, Np*K * sizeof(double));
	midERHS = (double *)malloc(2 * Np * K * Nface * sizeof(double));
	MemoryAllocationCheck(midERHS, 2 * Np*K * Nface * sizeof(double));
	memset(midERHS, 0, 2 * Np*K * Nface * sizeof(double));
	RSBotBEfm = (double *)malloc(BotBENfp*BotBENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotBEfm, BotBENfp*BotBENe * 2 * sizeof(double));
	memset(RSBotBEfm, 0, BotBENfp*BotBENe * 2 * sizeof(double));
	RSBotBEFluxM = (double *)malloc(BotBENfp*BotBENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotBEFluxM, BotBENfp*BotBENe * 2 * sizeof(double));
	memset(RSBotBEFluxM, 0, BotBENfp*BotBENe * 2 * sizeof(double));
	RSBotBEFluxS = (double *)malloc(BotBENfp*BotBENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSBotBEFluxS, BotBENfp*BotBENe * 2 * sizeof(double));
	memset(RSBotBEFluxS, 0, BotBENfp*BotBENe * 2 * sizeof(double));
	RSSurfBEfm = (double *)malloc(SurfBENfp*SurfBENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSSurfBEfm, SurfBENfp*SurfBENe * 2 * sizeof(double));
	memset(RSSurfBEfm, 0, SurfBENfp*SurfBENe * 2 * sizeof(double));
	RSSurfBEFluxM = (double *)malloc(SurfBENfp*SurfBENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSSurfBEFluxM, SurfBENfp*SurfBENe * 2 * sizeof(double));
	memset(RSSurfBEFluxM, 0, SurfBENfp*SurfBENe * 2 * sizeof(double));
	RSSurfBEFluxS = (double *)malloc(SurfBENfp*SurfBENe * 2 * sizeof(double));
	MemoryAllocationCheck(RSSurfBEFluxS, SurfBENfp*SurfBENe * 2 * sizeof(double));
	memset(RSSurfBEFluxS, 0, SurfBENfp*SurfBENe * 2 * sizeof(double));
}

void NdgMemory::RollerWaveRadiationMemoryDeAllocation()
{
	free(SIN_DIR); SIN_DIR = NULL;
	free(COS_DIR); COS_DIR = NULL;
	free(WaveNumber); WaveNumber = NULL;
	free(RSIEfm); RSIEfm = NULL;
	free(RSIEfp); RSIEfp = NULL;
	free(RSBEfm); RSBEfm = NULL;
	free(RSBEfp); RSBEfp = NULL;
	free(RSIEFluxM); RSIEFluxM = NULL;
	free(RSIEFluxP); RSIEFluxP = NULL;
	free(RSIEFluxS); RSIEFluxS = NULL;
	free(RSBEFluxM); RSBEFluxM = NULL;
	free(RSBEFluxS); RSBEFluxS = NULL;
	free(RSBotEfm); RSBotEfm = NULL;
	free(RSBotEfp); RSBotEfp = NULL;
	free(RSBotEFluxM); RSBotEFluxM = NULL;
	free(RSBotEFluxP); RSBotEFluxP = NULL;
	free(RSBotEFluxS); RSBotEFluxS = NULL;
	free(WaveEnergy); WaveEnergy = NULL;
	free(KD); KD = NULL;
	free(FSS); FSS = NULL;
	free(FCS); FCS = NULL;
	free(FSC); FSC = NULL;
	free(FCC); FCC = NULL;
	free(M2); M2 = NULL;
	free(N1); N1 = NULL;
	free(N2); N2 = NULL;
	free(Ar); Ar = NULL;
	free(Rz); Rz = NULL;
	free(R); R = NULL;
	free(CFF1); CFF1 = NULL;
	free(CFF2); CFF2 = NULL;
	free(CFF3); CFF3 = NULL;
	free(CFF4); CFF4 = NULL;
	free(Roller1); Roller1 = NULL;
	free(Roller2); Roller2 = NULL;
	free(Roller3); Roller3 = NULL;
	free(H_Radiation1); H_Radiation1 = NULL;
	free(H_Radiation2); H_Radiation2 = NULL;
	free(H_Radiation3); H_Radiation3 = NULL;
	free(V_Radiation); V_Radiation = NULL;
	free(RSERHS); RSERHS = NULL;
	free(RS_E); RS_E = NULL;
	free(RS_G); RS_G = NULL;
	free(Hx); Hx = NULL;
	free(Hy); Hy = NULL;
	free(tempRHSx); tempRHSx = NULL;
	free(tempRHSy); tempRHSy = NULL;
	free(midtempRHSx); midtempRHSx = NULL;
	free(midtempRHSy); midtempRHSy = NULL;
	free(VtempRHSx); VtempRHSx = NULL;
	free(VtempRHSy); VtempRHSy = NULL;
	free(midERHS); midERHS = NULL;
	free(RSBotBEfm); RSBotBEfm = NULL;
	free(RSBotBEFluxM); RSBotBEFluxM = NULL;
	free(RSBotBEFluxS); RSBotBEFluxS = NULL;
	free(RSSurfBEfm); RSBotBEFluxS = NULL;
	free(RSSurfBEFluxM); RSSurfBEFluxM = NULL;
	free(RSSurfBEFluxS); RSSurfBEFluxS = NULL;
	free(CFF5); CFF5 = NULL;
	free(CFF6); CFF6 = NULL;
	free(CFF7); CFF7 = NULL;
	free(CFF8); CFF8 = NULL;
	free(wave_C); wave_C = NULL;
	free(Rzn); Rzn = NULL;
	free(Rz1); Rz1 = NULL;
}

void NdgMemory::SedimentModelMemoryAllocation(int Np,int K)
{


}

void NdgMemory::SedimentModelMemoryDeAllocation()
{

}