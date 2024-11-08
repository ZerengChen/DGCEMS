//#if !defined(_WIN32)
//#define dgesv dgesv_
//#endif

#include "NdgMath.h"
#include "NdgMemory.h"
#include "UpdateImplicitRHS.h"
#include <stdio.h>

#include <stdlib.h>

#include <string.h>

#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102

UpdateImplicitRHS::UpdateImplicitRHS()
{
}

UpdateImplicitRHS::~UpdateImplicitRHS()
{
}

extern double *Tau;
extern double dt;
extern int Nvar;
extern signed char *Status3d;

void FetchBoundaryData(double *dest, double *source, const int Np2d, double *Eid)
{
	for (int i = 0; i < Np2d; i++) {
		dest[i] = source[(int)(Eid[i]) - 1];
	}
}

void CalculatePenaltyParameter(double *dest, const int Np2d, const int Np3d, double *UpEidM, double *BotEidM, double *nv, \
	int Nz, int P, int Nface)
{
	//double nvM[Np2d], nvP[Np2d];
	double *nvM = (double *)malloc(Np2d*sizeof(double));
	double *nvP = (double *)malloc(Np2d * sizeof(double));
	//Note: the penalty parameter for the topmost face of each column is not needed, so we leave them undefined
	for (int Layer = 1; Layer < Nz; Layer++)
	{
		FetchBoundaryData(nvM, nv + (Layer - 1)*Np3d, Np2d, BotEidM);
		FetchBoundaryData(nvP, nv + Layer * Np3d, Np2d, UpEidM);
		for (int p = 0; p < Np2d; p++)
		{
			dest[Layer*Np2d + p] = ((P + 1)*(P + 3) / 3.0)*(Nface / 2.0)*Nz*fmax(nvM[p], nvP[p]);
		}
	}
	FetchBoundaryData(nvM, nv + (Nz - 1)*Np3d, Np2d, BotEidM);
	for (int p = 0; p < Np2d; p++)
	{
		dest[Nz*Np2d + p] = ((P + 1)*(P + 3) / 3.0)*(Nface / 2.0)*Nz*nvM[p];
	}

	free(nvM);
	free(nvP);
}

/*Note: This function is used to assemble the local columns index and local rows index of the studied cell, and has been verified*/
void GetLocalRowsAndColumns(double *LocalRows, double *LocalColumns, int Index, int Np)
{
	for (int i = 0; i < Np; i++) {
		LocalRows[i] = Index * Np + i + 1;
		LocalColumns[i] = Index * Np + i + 1;
	}
}
/*Note: This function is used to calculate the volumn integral contained in primal form*/
void VolumnIntegral(double *dest, double *Dz, double *Mass3d, double *diff, int Np)
{
	double zero = 0.0;
	double Alpha1 = -1.0;
	double Alpha2 = 1.0;
	double *Tempdest = (double *)malloc(Np*Np * sizeof(double));
	MatrixMultiplyTN(Dz, Mass3d, Tempdest, Np, Np, Np, Alpha1);
	MatrixMultiply(Tempdest, diff, dest, Np, Np, Np, Alpha2);
	/*dgemm("t", "n", &Np, &Np, &Np, &Alpha1, Dz, &Np, Mass3d, &Np, &zero, Tempdest, &Np);
	dgemm("n", "n", &Np, &Np, &Np, &Alpha2, Tempdest, &Np, diff, &Np, &zero, dest, &Np);*/
	free(Tempdest);
}

/*This function is used to impose the Newmann boundary at the surface and bottom boundary*/
void ImposeNewmannBoundary(double *Eid, double *mass2d, double* InvMassMatrix3d, double dt, double Imparam, \
	 double* BoundNewmannData, int Np2d, int K2d, double *SystemRHS, int Np3d, int K3d, double *StiffMatrix, int Nvar, char *ch, int *varIndex_)
{
	double *tempRHS = (double *)malloc(Np3d * sizeof(double));
	double *RHS2d = (double *)malloc(Np2d * sizeof(double));
	for (int i = 0; i < Nvar; i++) {// For TOP's hu,hv
		memset(tempRHS, 0, Np3d * sizeof(double));
		memset(RHS2d, 0, Np2d * sizeof(double));
		int colB = 1;
		double Alpha = 1.0;
		double Beta = 0.0;
		MatrixMultiply(mass2d, BoundNewmannData + i * Np2d*K2d, RHS2d, Np2d, colB, Np2d, Alpha);
		//dgemm("n", "n", &Np2d, &colB, &Np2d, &Alpha, mass2d, &Np2d, BoundNewmannData + i * Np2d*K2d, &Np2d, &Beta, RHS2d, &Np2d);
		for (int j = 0; j < Np2d; j++)
			tempRHS[(int)Eid[j] - 1] = RHS2d[j];
		MatrixMultiply(InvMassMatrix3d, tempRHS, StiffMatrix + i * Np3d, Np3d, colB, Np3d, Alpha);
		//dgemm("n", "n", &Np3d, &colB, &Np3d, &Alpha, InvMassMatrix3d, &Np3d, tempRHS, &Np3d, &Beta, StiffMatrix + i * Np3d, &Np3d);
		for (int j = 0; j < Np3d; j++)
			SystemRHS[(varIndex_[i] - 1)*K3d*Np3d + j] += dt * Imparam*(*(StiffMatrix + i * Np3d + j));		

	}

	free(tempRHS);
	free(RHS2d);
}
/*Note: This function is used to get the row index of cell bellow the studied one*/
void GetBottomRows(double *BottomRows, double *LocalRows, int Np3d)
{
	for (int i = 0; i < Np3d; i++)
		BottomRows[i] = LocalRows[i] + Np3d;
}
/*Note: This function is used to get the row index of cell on top of the studied one*/
void GetUpRows(double *UpRows, double *LocalRows, int Np3d)
{
	for (int i = 0; i < Np3d; i++)
		UpRows[i] = LocalRows[i] - Np3d;
}
/*Note: This function is used to get the facial diff matrix, and has been checked*/
void AssembleFacialDiffMatrix(double *dest, double *source, double *Eid, int Np2d, int Np3d)
{
	int Index = 0;
	for (int colI = 0; colI < Np3d; colI++) {
		for (int RowI = 0; RowI < Np2d; RowI++)
		{
			dest[Index] = source[colI*Np3d + (int)Eid[RowI] - 1];
			Index++;
		}
	}
}

/*Note: This function is used to assemble the global stiff matrix*/
void AssembleGlobalStiffMatrix(double *dest, double *invMass, double *OP, double *rows, double *columns, double alpha, int StiffRow, int Np)
{
	double *TempDest = (double *)malloc(Np * Np * sizeof(double));
	double Beta = 0;
	MatrixMultiply(invMass, OP, TempDest, Np, Np, Np, alpha);
	//dgemm("n", "n", &Np, &Np, &Np, &alpha, invMass, &Np, OP, &Np, &Beta, TempDest, &Np);
	//here Np stands Np3d, StiffRow is the row of the stiffmatrix
	AssembleContributionIntoRowAndColumn(dest, TempDest, rows, columns, StiffRow, (int)Np, 1);
	free(TempDest);
}
/*Note: This function is used to solve the equation AX=b to get term X, and has been checked*/
void EquationSolve(double *dest, int m, int n, double *CoeMatrix, double dt, double ImplicitParam, int Nvar, int K3d, int Np, int *varIndex)
{
	double *TempCoeMatrix = (double *)malloc(m*m * sizeof(double));
	int *iPivot = (int *)malloc(m * sizeof(int));
	int info;
	int p = m;
	for (int var = 0; var < Nvar; var++) {
		memcpy(TempCoeMatrix, CoeMatrix + var * m*m, m*m * sizeof(double));
		for (int col = 0; col < m; col++) {
			for (int row = 0; row < m; row++) {
				TempCoeMatrix[col*m + row] = -1 * dt*ImplicitParam*TempCoeMatrix[col*m + row];
			}
			TempCoeMatrix[col*m + col] += 1;
		}
		//dgesv(&m, &n, TempCoeMatrix, &m, iPivot, dest + var * K3d*Np, &p, &info);//求逆
		LAPACKE_dgesv(LAPACK_COL_MAJOR, m, n, TempCoeMatrix, m, iPivot, dest + (varIndex[var] - 1) * K3d*Np, m);
	}
	free(iPivot);
	free(TempCoeMatrix);
}
/*Note: This function is used to calculate the local boundary contribution to the local stiff operator OP11*/
void LocalBoundaryIntegral(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, int Np3d, int Np2d, int Flag, double epsilon)
{
	double *FDiffMatrix = (double *)malloc(Np2d*Np3d * sizeof(double));
	AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, (int)Np2d, (int)Np3d);
	double Alpha = -1 * epsilon*Flag*0.5, Beta = 0.0;
	double *EdgeContribution = (double *)malloc(Np2d*Np3d * sizeof(double));
	MatrixMultiplyTN(FDiffMatrix, mass2d, EdgeContribution, Np3d, Np2d, Np2d, Alpha);
	//dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
	AssembleContributionIntoColumn(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
	Alpha = 0.5*Flag;
	MatrixMultiply(mass2d, FDiffMatrix, EdgeContribution, Np2d, Np3d, Np2d, Alpha);
	//dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
	AssembleContributionIntoRow(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
	double *DoubleJump = (double *)malloc(Np2d*Np2d * sizeof(double));
	DiagMultiply(DoubleJump, mass2d, TempTau, (int)Np2d);
	AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, (int)Np3d, (int)Np2d, -1);
	free(FDiffMatrix);
	free(EdgeContribution);
	free(DoubleJump);
}
/*Note: This function is used to impose the homogeneous Dirichlet boundary condition at the bottom boundary*/
void ImposeDirichletBoundary(double *eid, double *DiffMatrix, double *mass2d, double *TempTau, double *OP11, int Np3d, int Np2d, int Flag, const double epsilon)
{
	/*	double *Tau = (double *)malloc(Np2d*sizeof(double));
	for (int i = 0; i < Np2d; i++)
		Tau[i] = 2 * TempTau[i];
	 */
	 //LocalBoundaryIntegral(eid, DiffMatrix, mass2d, Tau, OP11, Np3d, Np2d, Flag, epsilon);
	double *FDiffMatrix = (double *)malloc(Np2d*Np3d * sizeof(double));
	AssembleFacialDiffMatrix(FDiffMatrix, DiffMatrix, eid, (int)Np2d, (int)Np3d);
	double Alpha = -1 * epsilon*Flag*0.5 * 2, Beta = 0.0;
	double *EdgeContribution = (double *)malloc(Np2d*Np3d * sizeof(double));
	MatrixMultiplyTN(FDiffMatrix, mass2d, EdgeContribution, Np3d, Np2d, Np2d, Alpha);
	//dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
	AssembleContributionIntoColumn(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
	Alpha = 0.5*Flag * 2;
	MatrixMultiply(mass2d, FDiffMatrix, EdgeContribution, Np2d, Np3d, Np2d, Alpha);
	//dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
	AssembleContributionIntoRow(OP11, EdgeContribution, eid, (int)Np3d, (int)Np2d);
	double *DoubleJump = (double *)malloc(Np2d*Np2d * sizeof(double));
	//DiagMultiply(DoubleJump, mass2d, Tau, (int)Np2d);
	DiagMultiply(DoubleJump, mass2d, TempTau, (int)Np2d);
	AssembleContributionIntoRowAndColumn(OP11, DoubleJump, eid, eid, (int)Np3d, (int)Np2d, -1);
	free(FDiffMatrix);
	free(EdgeContribution);
	free(DoubleJump);
	//	free(Tau);
}
/*Note: This function is used to calculate the adjacent boundary contribution to the adjacent stiff operator OP12, here adjacent means the test function is defined over the adjacent cell and not the local one*/
void AdjacentBoundaryIntegral(double *eidM, double *eidP, double *LocalDiff, double *AdjacentDiff, double *mass2d, double *TempTau, double *OP12, int Np3d, \
	int Np2d, int Flag, const double epsilon)
{
	double *FDiffMatrix = (double *)malloc(Np2d*Np3d * sizeof(double));
	AssembleFacialDiffMatrix(FDiffMatrix, AdjacentDiff, eidP, (int)Np2d, (int)Np3d);
	double Alpha = -1 * epsilon*Flag*0.5, Beta = 0;
	double *EdgeContribution = (double *)malloc(Np2d*Np3d * sizeof(double));
	MatrixMultiplyTN(FDiffMatrix, mass2d, EdgeContribution, Np3d, Np2d, Np2d, Alpha);
	//dgemm("t", "n", &Np3d, &Np2d, &Np2d, &Alpha, FDiffMatrix, &Np2d, mass2d, &Np2d, &Beta, EdgeContribution, &Np3d);
	AssembleContributionIntoColumn(OP12, EdgeContribution, eidM, (int)Np3d, (int)Np2d);
	//Alpha = 0.5*Flag;
	Alpha = -1 * 0.5*Flag;
	AssembleFacialDiffMatrix(FDiffMatrix, LocalDiff, eidM, (int)Np2d, (int)Np3d);
	MatrixMultiply(mass2d, FDiffMatrix, EdgeContribution, Np2d, Np3d, Np2d, Alpha);
	//dgemm("n", "n", &Np2d, &Np3d, &Np2d, &Alpha, mass2d, &Np2d, FDiffMatrix, &Np2d, &Beta, EdgeContribution, &Np2d);
	AssembleContributionIntoRow(OP12, EdgeContribution, eidP, (int)Np3d, (int)Np2d);
	double *DoubleJump = (double *)malloc(Np2d*Np2d * sizeof(double));
	DiagMultiply(DoubleJump, mass2d, TempTau, (int)Np2d);
	AssembleContributionIntoRowAndColumn(OP12, DoubleJump, eidP, eidM, (int)Np3d, (int)Np2d, 1);
	free(FDiffMatrix);
	free(EdgeContribution);
	free(DoubleJump);
}
/*Note: This function is used to add the Newmman boundary part back to the implicit right hand side*/
void AssembleBoundaryContribution(double *dest, double *source, int Np, int K3d, int Nvar)
{
	for (int var = 0; var < Nvar; var++) {
		for (int i = 0; i < Np; i++)
			//dest[2 * var*K3d*Np + i] += source[var*Np + i];
		    dest[var*K3d*Np + i] += source[var*Np + i];
	}
}

void UpdateImplicitRHS::EvaluateupdateimplicitRHS(double *fphys_, double *nv_v,double *frhs_, double ImplicitA, double *BBE, double *SBE,int *varIndex,int*pE2d,int*pE3d,int MyID)
{
	//mexAtExit(&MyExit);
	double *J2d = meshunion->mesh2d_p->J2d;
	double *J3d = meshunion->J;
	double *M2d = meshunion->mesh2d_p->mesh2dcell_p->M2d;
	double *M3d = meshunion->cell_p->M;
	double *tz = meshunion->tz;
	double *Dt = meshunion->cell_p->Dt;
	//double *Diff = fphys_ + 4 * Np*K3d;
	int Nface = *(meshunion->cell_p->Nface);
	int BENfp = (*meshunion->boundaryedge_p->Nfp);
	int BotBENfp = *(meshunion->bottomboundaryedge_p->Nfp);
	int Node = fmax(BENfp, BotBENfp);
	double *bot = BBE;
	double *surf = SBE;
	double ImplicitParam = ImplicitA;
	//double *RHS = SystemRHS_;
	double Prantl = 1.0;
	int K2d = *(meshunion->mesh2d_p->K2d);
	int K3d = *(meshunion->K);
	int Nz = K3d / K2d;
	int Np = *(meshunion->cell_p->Np);
	int N = *(meshunion->cell_p->N);
	int Np2d = *(meshunion->mesh2d_p->mesh2dcell_p->Np2d);
	double *UpEidM = meshunion->cell_p->Fmask + (Nface - 1)*Node;
	double *BotEidM = meshunion->cell_p->Fmask + (Nface - 2)*Node;
	int P = *(meshunion->cell_p->Nz);
	
	double *Diff = nv_v;
	char* BoundaryType;
	BoundaryType = (char*)meshunion->bottomboundaryedge_p->ftype;
	signed char *status = meshunion->mesh2d_p->status;
	//    printf("%s\n", BoundaryType);
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

	double *fphys = fphys_;


#ifdef COUPLING_SWAN
	//Change to Euler vilocity
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		if (MyID == pE3d[k]) {
			for (int n = 0; n < Np; n++) {
				fphys[k * Np + n] = fphys[K3d * Np * 9 + k * Np + n];
				fphys[K3d * Np + k * Np + n] = fphys[K3d * Np * 10 + k * Np + n];
			}
		}
	}
#endif

	double *ImplicitRHS = frhs_;
	/*Here, epsilon is set to be -1.0, and this correspoinding to the SIPG. This value can be set to 1.0(IIPG or NIPG) or 0.0(IIPG or NIPG)*/
	const double epsilon = -1.0;


#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {
		if (MyID == pE2d[i]) {
			CalculatePenaltyParameter(Tau + i * Np2d*(Nz + 1), Np2d, Np, UpEidM, BotEidM, Diff + i * Np*Nz, Nz, P, Nface);
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int i = 0; i < K2d; i++) {
		if (MyID == pE2d[i]) {
			NdgRegionType type = (NdgRegionType)status[i];
			if (type == NdgRegionWet) {
				double *StiffMatrix = (double *)malloc(Np*Nz*Np*Nz*Nvar * sizeof(double));
				memset(StiffMatrix, 0, Np*Nz*Np*Nz*Nvar * sizeof(double));
				double *EleMass3d = (double *)malloc(Np*Np * sizeof(double));
				double *InvEleMass3d = (double *)malloc(Np*Np * sizeof(double));
				double *EleMass2d = (double *)malloc(Np2d*Np2d * sizeof(double));
				double *LocalPhysicalDiffMatrix = (double *)malloc(Np*Np * sizeof(double));
				double *Dz = (double *)malloc(Np*Np * sizeof(double));
				double *OP11 = (double *)malloc(Np*Np * sizeof(double));
				double *LocalRows = (double *)malloc(Np * sizeof(double));
				double *LocalColumns = (double *)malloc(Np * sizeof(double));
				DiagMultiply(EleMass3d, M3d, J3d + i * Nz*Np, Np);
				memcpy(InvEleMass3d, EleMass3d, Np*Np * sizeof(double));
				//MatrixInverseTest();
				MatrixInverse(InvEleMass3d, (int)Np);

				DiagMultiply(EleMass2d, M2d, J2d + i * Np2d, Np2d);
				DiagMultiply(Dz, Dt, tz + i * Nz*Np, Np);
				DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i * Nz*Np, Np);
				GetLocalRowsAndColumns(LocalRows, LocalColumns, 0, Np);
				/*Calculate the volumn integral for the first cell, and impose the surface Newmann boundary condition*/
				VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
				double *SurfBoundStiffTerm = (double *)malloc(Np*Nvar * sizeof(double));
				double *BotBoundStiffTerm = (double *)malloc(Np*Nvar * sizeof(double));
				char ch = 'T';
				ImposeNewmannBoundary(UpEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, surf + i * Np2d, (int)Np2d, K2d, fphys + i * Np*Nz, (int)Np, K3d, SurfBoundStiffTerm, Nvar, &ch, varIndex);
				ch = 'B';
				//有温盐泥沙时Nvar不等于2，边界条件得修改

				if (Nz != 1) {
					/*When vertical layers greater than one, we calculate the bottom boundary integral part */
					double *BottomAdjacentRows = (double *)malloc(Np * sizeof(double));
					double *UpAdjacentRows = (double *)malloc(Np * sizeof(double));
					GetBottomRows(BottomAdjacentRows, LocalRows, Np);
					double *BottomPhysicalDiffMatrix = (double *)malloc(Np*Np * sizeof(double));
					double *UpPhysicalDiffMatrix = (double *)malloc(Np*Np * sizeof(double));
					DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i * Nz*Np + Np, Np);
					LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + 2 - 1), OP11, (int)Np, (int)Np2d, -1, epsilon);
					double *OP12 = (double *)malloc(Np*Np * sizeof(double));
					memset(OP12, 0, Np*Np * sizeof(double));
					AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + 2 - 1), OP12, Np, Np2d, -1, epsilon);
					for (int var = 0; var < 2; var++) {
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1.0, Np*Nz, Np);
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, 1.0, Np*Nz, Np);
					}
					for (int var = 2; var < Nvar; var++) {
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
					}
					for (int j = 1; j < Nz - 1; j++) {
						/*Calculate both the volume and surface integral for the second to the last but one cell*/
						GetLocalRowsAndColumns(LocalRows, LocalColumns, j, Np);
						GetBottomRows(BottomAdjacentRows, LocalRows, Np);
						GetUpRows(UpAdjacentRows, LocalRows, Np);
						DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i * Nz*Np + (j - 1)*Np, Np);
						DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i * Nz*Np + j * Np, Np);
						DiagMultiply(BottomPhysicalDiffMatrix, Dz, Diff + i * Nz*Np + (j + 1)*Np, Np);
						VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
						LocalBoundaryIntegral(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + j + 1), OP11, (int)Np, (int)Np2d, -1, epsilon);
						LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + j), OP11, (int)Np, (int)Np2d, 1, epsilon);
						memset(OP12, 0, Np*Np * sizeof(double));
						AdjacentBoundaryIntegral(UpEidM, BotEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + j), OP12, Np, Np2d, 1, epsilon);

						for (int var = 0; var < 2; var++) {
							AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1, Np*Nz, Np);
							AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, 1, Np*Nz, Np);
						}
						for (int var = 2; var < Nvar; var++) {
							AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
							AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
						}
						memset(OP12, 0, Np*Np * sizeof(double));
						AdjacentBoundaryIntegral(BotEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + j + 1), OP12, Np, Np2d, -1, epsilon);
						for (int var = 0; var < 2 && var < Nvar; var++) {
							AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, 1, Np*Nz, Np);
						}
						for (int var = 2; var < Nvar; var++) {
							AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, BottomAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
						}
					}
					/*For the bottom most cell, calculate the volumn integral and impose Neumann boundary condition or the Dirichlet boundary condition*/
					GetLocalRowsAndColumns(LocalRows, LocalColumns, Nz - 1, Np);
					GetUpRows(UpAdjacentRows, LocalRows, Np);
					DiagMultiply(LocalPhysicalDiffMatrix, Dz, Diff + i * Nz*Np + (Nz - 1)*Np, Np);
					DiagMultiply(UpPhysicalDiffMatrix, Dz, Diff + i * Nz*Np + (Nz - 1 - 1)*Np, Np);
					VolumnIntegral(OP11, Dz, EleMass3d, LocalPhysicalDiffMatrix, Np);
					LocalBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + Nz - 1), OP11, (int)Np, (int)Np2d, 1, epsilon);
					/*For passive transport substances, we only impose Neumann boundary condition, so this has no effect on the stiff matrix*/
					for (int var = 2; var < Nvar; var++) {
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
					}
					/*We note that, the Neumann data is zero by default, so the impositon of Neumann BC for hu and hv will not affect the Dirichlet boundary condition for hu and hv*/
					ImposeNewmannBoundary(BotEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, bot + i * Np2d, (int)Np2d, K2d, fphys + i * Np*Nz + (Nz - 1)*Np, (int)Np, K3d, BotBoundStiffTerm, Nvar, &ch, varIndex);
					/*The following is used to add homogeneous dirichlet boundary for hu and hv*/
					//if (!strcmp(BoundaryType, "Dirichlet")) {
					//	ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + Nz), OP11, (int)Np, (int)Np2d, -1, epsilon);
					//}
					memset(OP12, 0, Np*Np * sizeof(double));
					AdjacentBoundaryIntegral(UpEidM, BotEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + Nz - 1), OP12, (int)Np, (int)Np2d, 1, epsilon);
					for (int var = 0; var < 2; var++) {
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, 1, Np*Nz, Np);
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1, Np*Nz, Np);
					}
					for (int var = 2; var < Nvar; var++) {
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP12, UpAdjacentRows, LocalColumns, Prantl, Np*Nz, Np);
					}

					free(BottomAdjacentRows);
					free(UpAdjacentRows);
					free(BottomPhysicalDiffMatrix);
					free(UpPhysicalDiffMatrix);
					free(OP12);
				}
				else {/*If only one layer included in vertical direction, we need to consider the bottom face, and the Upper face has been considered in the first part, from 291-312*/
					/*For passive transport substances, we only impose Neumann boundary condition*/
					for (int var = 2; var < Nvar; var++) {
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, Prantl, Np*Nz, Np);
					}
					/*We note that, the Neumann data is zero by default, so the impositon of Neumann BC for hu and hv will not affect the Dirichlet boundary condition for hu and hv*/
					ImposeNewmannBoundary(BotEidM, EleMass2d, InvEleMass3d, dt, ImplicitParam, bot + i * Np2d, (int)Np2d, K2d, fphys + i * Np*Nz + (Nz - 1)*Np, (int)Np, K3d, BotBoundStiffTerm, Nvar, &ch, varIndex);
					/* The following is used to add homogeneous dirichlet boundary for hu and hv*/
					//if (!strcmp(BoundaryType, "Dirichlet")) {
					//	ImposeDirichletBoundary(BotEidM, LocalPhysicalDiffMatrix, EleMass2d, Tau + Np2d * (i*(Nz + 1) + Nz), OP11, (int)Np, (int)Np2d, -1, epsilon);
					//}
					for (int var = 0; var < 2; var++) {
						AssembleGlobalStiffMatrix(StiffMatrix + var * Np*Nz*Np*Nz, InvEleMass3d, OP11, LocalRows, LocalColumns, 1, Np*Nz, Np);
					}
				}

				EquationSolve(fphys + i * Nz*Np, Nz*Np, 1, StiffMatrix, dt, ImplicitParam, Nvar, K3d, Np, varIndex);
				//EquationSolve(fphys + i * Nz*Np, Nz*Np, 1, StiffMatrix, dt, ImplicitParam, 2, K3d, Np, varIndex);

				const int dimension = Np * Nz;
				const int colB = 1;
				double Alpha = 1.0, Beta = 0.0;
				for (int var = 0; var < Nvar; var++) {
					MatrixMultiply(StiffMatrix + var * Np*Nz*Np*Nz, fphys + (varIndex[var] - 1) * K3d*Np + i * Nz*Np, ImplicitRHS + var * K3d*Np + i * Nz*Np, dimension, colB, dimension, Alpha);
					//MatrixMultiply(StiffMatrix + var * Np*Nz*Np*Nz, fphys + var * K3d*Np + i * Nz*Np, ImplicitRHS + var * K3d*Np + i * Nz*Np, dimension, colB, dimension, Alpha);
					//dgemm("n", "n", &dimension, &colB, &dimension, &Alpha, StiffMatrix + var * Np*Nz*Np*Nz, &dimension, fphys + var * K3d*Np + i * Nz*Np, &dimension, &Beta, ImplicitRHS + var * K3d*Np + i * Nz*Np, &dimension);
				}
				AssembleBoundaryContribution(ImplicitRHS + i * Nz*Np, SurfBoundStiffTerm, Np, K3d, Nvar);
				AssembleBoundaryContribution(ImplicitRHS + i * Nz*Np + (Nz - 1)*Np, BotBoundStiffTerm, Np, K3d, Nvar);

				free(EleMass3d);
				free(InvEleMass3d);
				free(EleMass2d);
				free(LocalPhysicalDiffMatrix);
				free(Dz);
				free(OP11);
				free(LocalRows);
				free(LocalColumns);
				free(SurfBoundStiffTerm);
				free(BotBoundStiffTerm);
				free(StiffMatrix);

			}
			else {
				continue;
			}
		}
	}

#ifdef COUPLING_SWAN
	//Update Lagrange hu hv
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
	for (int k = 0; k < K3d; k++) {
		if (MyID == pE3d[k]) {
			for (int n = 0; n < Np; n++) {
				fphys[K3d * Np * 9 + k * Np + n] = fphys[k * Np + n];
				fphys[K3d * Np * 10 + k * Np + n] = fphys[K3d * Np + k * Np + n];
				fphys[k * Np + n] = fphys[k * Np + n] + fphys[K3d * Np * 11 + k * Np + n];
				fphys[K3d * Np + k * Np + n] = fphys[K3d * Np + k * Np + n] + fphys[K3d * Np * 12 + k * Np + n];
			}
		}
	}

#endif

}
