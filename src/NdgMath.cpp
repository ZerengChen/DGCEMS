#include "NdgMath.h"
#include "cblas.h"
#include <iostream>
#include<mpi.h>
using namespace std;
//#include <Accelerate.h>
//#include <omp.h>
//extern "C" {
//	// LU decomoposition of a general matrix 
//	void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
//	lapack_int LAPACKE_dgetrf(int matrix_layout, lapack_int m, lapack_int n,
//		double* a, lapack_int lda, lapack_int* ipiv);
//	// generate inverse of a matrix given its LU decomposition 
//	void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
//	lapack_int LAPACKE_dgetri(int matrix_layout, lapack_int n, double* a,
//		lapack_int lda, const lapack_int* ipiv);
//}

#ifndef lapack_int
#define lapack_int     int
#endif
#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102

void Add(double *dest, double *sourcea, double *sourceb, int size) {
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] + sourceb[i];
}

void AddByConstant(double *dest, double *sourcea, double ConstData, int size) {
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] + ConstData;
}

/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the column index, and has been checked*/
void AssembleContributionIntoColumn(double *dest, double *source, double *column, int Np3d, int Np2d)
{
	for (int colI = 0; colI < Np2d; colI++) {
		for (int RowI = 0; RowI < Np3d; RowI++)
			dest[((int)column[colI] - 1)*Np3d + RowI] += source[colI*Np3d + RowI];
	}
}
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index, and has been checked*/
void AssembleContributionIntoRow(double *dest, double *source, double *Row, int Np3d, int Np2d)
{
	for (int colI = 0; colI < Np3d; colI++) {
		for (int RowI = 0; RowI < Np2d; RowI++)
			dest[colI*Np3d + (int)Row[RowI] - 1] += source[colI*Np2d + RowI];
	}
}
/*Note: This function is used to assemble the facial integral term into the local stiff operator according to the row index and the column index, and has been checked*/
void AssembleContributionIntoRowAndColumn(double *dest, double *source, double *Row, double *column, int Np3d, int Np2d, int Flag)
{
	for (int colI = 0; colI < Np2d; colI++)
	{
		for (int RowI = 0; RowI < Np2d; RowI++)
			dest[((int)column[colI] - 1)*Np3d + (int)Row[RowI] - 1] += Flag * source[colI*Np2d + RowI];
	}
}

void AssembleDataIntoPoint(double *dest, double *source, double *PIndex, int Size) {
	for (int p = 0; p < Size; p++) {
		dest[(int)PIndex[p] - 1] += source[p];
	}
}

/*Note: this function is used to assemble the element mass matrix and the physical diff matrix, and has been verified.
Left multiply the matrix source with a diagonal matrix composed of element contained in coe.
*/
void DiagMultiply(double *dest, const double *source, const double *coe, int Np)
{
	for (int colI = 0; colI < Np; colI++) {
		for (int RowI = 0; RowI < Np; RowI++)
			dest[colI*Np + RowI] = coe[RowI] * source[colI*Np + RowI];
	}
}

void DiagRightMultiply(double *dest, const double *source, const double *coe, int Np)
{
	for (int colI = 0; colI < Np; colI++) {
		for (int RowI = 0; RowI < Np; RowI++)
			dest[colI*Np + RowI] = coe[colI] * source[colI*Np + RowI];
	}
}

void DotCriticalDivide(double *dest, double *source, double *criticalValue, double *Depth, int size) {
	for (int i = 0; i < size; i++) {
		if (Depth[i] >= *criticalValue)
			dest[i] = source[i] / Depth[i];
		else
			dest[i] = 0;
	}
}

void DotDivide(double *dest, double *source, double *Coefficient, int size) {
	for (int i = 0; i < size; i++)
		dest[i] = source[i] / Coefficient[i];
}

void DotDivideByConstant(double *dest, double *Source, double Coefficient, int Np) {
	for (int i = 0; i < Np; i++)
		dest[i] = Source[i] / Coefficient;
}

void DotProduct(double *dest, double *sourcea, double *sourceb, int size) {
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] * sourceb[i];
}

void Exchange(double *fphys, double *IEFToE, int *pE, int IENe, int Np, int K, int Nfield)
{
	//对2D,3D均适用
	int MyID, MPIsize;
	int ThisElement, AdjacentElement;
	MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	MPI_Status status_MPI; // variable to contain status information
	//MPI_Request request_MPI[MPIsize];
	for (int i = 0; i < IENe; i++) {
		ThisElement = (int)IEFToE[2 * i];//本侧的单元编号
		AdjacentElement = (int)IEFToE[2 * i + 1];//相邻的单元编号
		if ((pE[AdjacentElement - 1] != pE[ThisElement - 1])) {//如果两者的分区不一致，调出这两个分区对应的进程交换数据
			//if (MyID == pE[ThisElement - 1]){
			//	for (int field = 0; field < Nfield; field++) {
			//		int tag1 = i + field;
			//		int tag2 = i + Nfield + field;
			//		MPI_Sendrecv(fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag1, \
			//			fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag2, \
			//			MPI_COMM_WORLD, &status_MPI);
			//	}
   //         }
			//else if (MyID == pE[AdjacentElement - 1]) {
			//	for (int field = 0; field < Nfield; field++) {
			//		int tag3 = i + Nfield + field;
			//		int tag4 = i + field;
			//		MPI_Sendrecv(fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag3, \
			//			fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag4, \
			//			MPI_COMM_WORLD, &status_MPI);
			//	}
			//}
			//if (MyID == pE[ThisElement - 1]){
			//	for (int field = 0; field < Nfield; field++) {
			//		int tag1 = field * IENe * MPIsize + i * MPIsize + pE[ThisElement - 1];
			//		int tag2 = field * IENe * MPIsize + i * MPIsize + pE[AdjacentElement - 1];
			//		MPI_Isend(fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag1, MPI_COMM_WORLD, &request_MPI[0]);
			//		MPI_Irecv(fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag2, MPI_COMM_WORLD, &request_MPI[1]);
			//		MPI_Wait(&request_MPI[1], &status_MPI);
			//	}
			//}
			//else if (MyID == pE[AdjacentElement - 1]) {
			//	for (int field = 0; field < Nfield; field++) {
			//		int tag1 = field * IENe * MPIsize + i * MPIsize + pE[AdjacentElement - 1];
			//		int tag2 = field * IENe * MPIsize + i * MPIsize + pE[ThisElement - 1];
			//		MPI_Irecv(fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag2, MPI_COMM_WORLD, &request_MPI[0]);
			//		MPI_Isend(fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag1, MPI_COMM_WORLD, &request_MPI[1]);
			//		MPI_Wait(&request_MPI[0], &status_MPI);
			//	}
			//}
			if (MyID == pE[ThisElement - 1]) {
				for (int field = 0; field < Nfield; field++) {
					int tag1 = pE[ThisElement - 1] + field;
					int tag2 = pE[AdjacentElement - 1] + field;
					MPI_Send(fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag1, MPI_COMM_WORLD);
					//std::cout << "ThisID "<< pE[ThisElement - 1] <<" send to ID "<< pE[AdjacentElement - 1] << std::endl;
					MPI_Recv(fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag2, MPI_COMM_WORLD, &status_MPI);
					//std::cout << "ThisID " << pE[ThisElement - 1] << " recv from ID " << pE[AdjacentElement - 1] << std::endl;
				}
			}
			else if (MyID == pE[AdjacentElement - 1]) {
				for (int field = 0; field < Nfield; field++) {
					int tag1 = pE[AdjacentElement - 1] + field;
					int tag2 = pE[ThisElement - 1] + field;
					MPI_Recv(fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag2, MPI_COMM_WORLD, &status_MPI);
					//std::cout << "AdjacentID " << pE[ThisElement - 1] << " send to ID " << pE[AdjacentElement - 1] << std::endl;
					MPI_Send(fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag1, MPI_COMM_WORLD);
					//std::cout << "AdjacentID " << pE[ThisElement - 1] << " recv from ID " << pE[AdjacentElement - 1] << std::endl;
				}
			}
			//else{}
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Exchange(signed char *fphys, double *IEFToE, int *pE, int IENe, int Np, int K, int Nfield)
{
	//对2D,3D均适用
	int MyID, MPIsize;
	int ThisElement, AdjacentElement;
	MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	MPI_Status status_MPI; // variable to contain status information
	//MPI_Request request_MPI[MPIsize];
	for (int i = 0; i < IENe; i++) {
		ThisElement = (int)IEFToE[2 * i];//本侧的单元编号
		AdjacentElement = (int)IEFToE[2 * i + 1];//相邻的单元编号
		if ((pE[AdjacentElement - 1] != pE[ThisElement - 1])) {	
			if (MyID == pE[ThisElement - 1]) {
				for (int field = 0; field < Nfield; field++) {
					int tag1 = pE[ThisElement - 1] + field;
					int tag2 = pE[AdjacentElement - 1] + field;
					MPI_Send(fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag1, MPI_COMM_WORLD);
					//std::cout << "ThisID "<< pE[ThisElement - 1] <<" send to ID "<< pE[AdjacentElement - 1] << std::endl;
					MPI_Recv(fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[AdjacentElement - 1], tag2, MPI_COMM_WORLD, &status_MPI);
					//std::cout << "ThisID " << pE[ThisElement - 1] << " recv from ID " << pE[AdjacentElement - 1] << std::endl;
				}
			}
			else if (MyID == pE[AdjacentElement - 1]) {
				for (int field = 0; field < Nfield; field++) {
					int tag1 = pE[AdjacentElement - 1] + field;
					int tag2 = pE[ThisElement - 1] + field;
					MPI_Recv(fphys + Np * (ThisElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag2, MPI_COMM_WORLD, &status_MPI);
					//std::cout << "AdjacentID " << pE[ThisElement - 1] << " send to ID " << pE[AdjacentElement - 1] << std::endl;
					MPI_Send(fphys + Np * (AdjacentElement - 1) + field * Np * K, Np, MPI_DOUBLE, pE[ThisElement - 1], tag1, MPI_COMM_WORLD);
					//std::cout << "AdjacentID " << pE[ThisElement - 1] << " recv from ID " << pE[AdjacentElement - 1] << std::endl;
				}
			}
			//else{}
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void GatherToZero(double *fphys, int *pE, int Np, int K, int Nfield) {
	//对2D,3D均适用
	int MyID, MPIsize;
	MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	MPI_Status status_MPI; // variable to contain status information
	MPI_Request request_MPI;
	for (int k = 0; k < K; k++) {
		if (pE[k] > 0) {
			if (MyID == pE[k]) {
				//MPI_Send(fphys + Np * k, Np, MPI_DOUBLE, 0, pE[k], MPI_COMM_WORLD);
				MPI_Isend(fphys + Np * k, Np, MPI_DOUBLE, 0, pE[k], MPI_COMM_WORLD, &request_MPI);
			}
			else{}
		}
		else {
			//MPI_Recv(fphys + Np * k, Np, MPI_DOUBLE, pE[k], pE[k], MPI_COMM_WORLD, &status_MPI);
			MPI_Irecv(fphys + Np * k, Np, MPI_DOUBLE, pE[k], pE[k], MPI_COMM_WORLD, &request_MPI);
			MPI_Wait(&request_MPI, &status_MPI);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void BroadcastToAll(double *fphys, int *pE, int Np, int K, int Nfield) {
	//对2D,3D均适用
	int MyID, MPIsize;
	MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
	MPI_Status status_MPI; // variable to contain status information
	MPI_Request request_MPI;
	if (MyID == 0) {
		for (int i = 1; i < MPIsize - 1; i++) {
			MPI_Isend(fphys, Np * K * Nfield, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request_MPI);
		}
	}
	else {
		MPI_Irecv(fphys, Np * K * Nfield, MPI_DOUBLE, 0, MyID, MPI_COMM_WORLD, &request_MPI);
		MPI_Wait(&request_MPI, &status_MPI);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}


void FetchInnerEdgeFacialValue(double *fm, double *fp, double *source, \
	double *FToE, double *FToN1, double *FToN2, int Np, int Nfp){
	int ind1 = ((int)FToE[0] - 1)*Np;
	int ind2 = ((int)FToE[1] - 1)*Np;
	for (int i = 0; i < Nfp; i++){
		fm[i] = source[ind1 + (int)FToN1[i] - 1];
		fp[i] = source[ind2 + (int)FToN2[i] - 1];
	}
}

void FetchBoundaryEdgeFacialValue(double *fm, double *source, \
	double *FToE, double *FToN1, int Np, int Nfp){
	int ind1 = ((int)FToE[0] - 1)*Np;
	for (int i = 0; i < Nfp; i++){
		fm[i] = source[ind1 + (int)FToN1[i] - 1];
	}
}

void Flip(double *dest, int size){
	double Temp;
	int i, j;
	i = 0;
	j = size - 1;
	while (i < j){
		Temp = dest[i];
		dest[i] = dest[j];
		dest[j] = Temp;
		++i;
		--j;
	}
}


void GetFacialFluxTerm2d(double *dest, double *hu, double *hv, double *nx, double *ny, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = hu[i] * nx[i] + hv[i] * ny[i];
}

void GetScalarFluxTerm2d(double *dest, double *C, double *Vector, int Nfp){
	for (int i = 0; i < Nfp; i++)
		dest[i] = C[i] * Vector[i];
}

void GetMeshAverageValue(double *dest, double *LAV, int *ROPA, int *COPB, int *COPA, double *A, double *fphys, double *Jacobian, int *LDC, double *wq){

	GetMeshIntegralValue(dest, ROPA, COPB, COPA, A, fphys, Jacobian, LDC, wq);

	(*dest) = (*dest) / (*LAV);
}

void GetMeshIntegralValue(double *dest, int *ROPA, int *COPB, int *COPA, double *A, double *fphys, double *Jacobian, int *LDC, double *wq){

	double *Jq = (double *)malloc((*LDC)*sizeof(double));
	double *fq = (double *)malloc((*LDC) * sizeof(double));
	//double Jq[(*LDC)], fq[(*LDC)];

	// map the node values fvar to quadrature nodes by
	// \f$ fq = Vq * fvar \f$
	MatrixMultiply(A, fphys, fq, *ROPA, *COPB, *COPA, 1.0);
	MatrixMultiply(A, Jacobian, Jq, *ROPA, *COPB, *COPA, 1.0);
	/*dgemm(transA, transB, ROPA, COPB, COPA, Alpha, A,
		LDA, fphys, LDB, Beta, fq, LDC);
	dgemm(transA, transB, ROPA, COPB, COPA, Alpha, A,
		LDA, Jacobian, LDB, Beta, Jq, LDC);*/

	for (int n = 0; n < (*LDC); n++) {
		(*dest) += wq[n] * Jq[n] * fq[n];
	}

	free(Jq);
	free(fq);
}

void GetVolumnIntegral1d(double *dest, int *RowOPA, int *ColOPB, int *ColOPA, double *alpha, \
	double *Dt, int *LDA, double *B, int *LDB, double *Beta, int *LDC, double *tz){

	//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dt, LDA, B, LDB, Beta, dest, LDC);
	MatrixMultiply(Dt, B, dest, *RowOPA, *ColOPB, *ColOPA, *alpha);

	DotProduct(dest, dest, tz, (int)(*LDC));


}

void GetVolumnIntegral2d(double *dest, double *tempdest, int *RowOPA, int *ColOPB, int *ColOPA, double *alpha, \
	double *Dr, double *Ds, int *LDA, double *B, int *LDB, double *Beta, int *LDC, double *rx, double *sx){
	/*$Dr*u(v,\theta)$*/
	MatrixMultiply(Dr, B, dest, *RowOPA, *ColOPB, *ColOPA, *alpha);
	//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dr, LDA, B, LDB, Beta, dest, LDC);
	/*$Ds*u(v,\theta)$*/
	MatrixMultiply(Ds, B, tempdest, *RowOPA, *ColOPB, *ColOPA, *alpha);
	//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Ds, LDA, B, LDB, Beta, tempdest, LDC);
	/*$rx\cdot Dr*u(v,\theta)$*/
	DotProduct(dest, dest, rx, (int)(*LDC));
	/*$sx\cdot Ds*u(v,\theta)$*/
	DotProduct(tempdest, tempdest, sx, (int)(*LDC));
	/*$rx\cdot Dr*u(v,\theta) + sx\cdot Ds*u(v,\theta)$*/
	Add(dest, dest, tempdest, (int)(*LDC));
}

void GetVolumnIntegralOnlyVertical(double *dest, int *RowOPA, int *ColOPB, int *ColOPA, double *alpha, \
	double *Dt, int *LDA, double *B, int *LDB, double *Beta, int *LDC, double *tz) {
	/*$Dt*Spx$*/
	MatrixMultiply(Dt, B, dest, *RowOPA, *ColOPB, *ColOPA, *alpha);
	//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dt, LDA, H + n*Np*K, LDB, Beta, tempdest, LDC);
	/*$tz\cdot Dt*H$*/
	DotProduct(dest, dest, tz, (int)(*LDC));
}

void GetVolumnIntegral3d(double *dest, double *tempdest, int *RowOPA, int *ColOPB, int *ColOPA, double *alpha, \
	double *Dr, double *Ds, double *Dt, double *E, double *G, double *H, int *LDA, int *LDB, double *Beta, int *LDC, \
	double *rx, double *sx, double *ry, double *sy, double *tz, int Nvar, int Np, int K)
{
	for (int n = 0; n < Nvar; n++){
		/*$Dr*E$*/
		MatrixMultiply(Dr, E + n * Np*K, dest + n * Np*K, *RowOPA, *ColOPB, *ColOPA, *alpha);
		//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dr, LDA, E + n*Np*K, LDB, Beta, dest + n*Np*K, LDC);
		/*$rx\cdot Dr*E$*/
		DotProduct(dest + n*Np*K, dest + n*Np*K, rx, (int)(*LDC));
		/*$Ds*E$*/
		MatrixMultiply(Ds, E + n * Np*K, tempdest, *RowOPA, *ColOPB, *ColOPA, *alpha);
		//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Ds, LDA, E + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$sx\cdot Ds*E$*/
		DotProduct(tempdest, tempdest, sx, (int)(*LDC));
		/*rx\cdot Dr*E + sx\cdot Ds*E*/
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
		/*$Dr*G$*/
		MatrixMultiply(Dr, G + n * Np*K, tempdest, *RowOPA, *ColOPB, *ColOPA, *alpha);
		//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dr, LDA, G + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$ry\cdot Dr*G$*/
		DotProduct(tempdest, tempdest, ry, (int)(*LDC));
		/*rx\cdot Dr*E + sx\cdot Ds*E + ry\cdot Dr*G*/
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
		/*$Ds*G$*/
		MatrixMultiply(Ds, G + n * Np*K, tempdest, *RowOPA, *ColOPB, *ColOPA, *alpha);
		//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Ds, LDA, G + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$sy\cdot Ds*G$*/
		DotProduct(tempdest, tempdest, sy, (int)(*LDC));
		/*rx\cdot Dr*E + sx\cdot Ds*E + ry\cdot Dr*G + sy\cdot Ds*G */
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
		/*$Dt*H$*/
		MatrixMultiply(Dt, H + n * Np*K, tempdest, *RowOPA, *ColOPB, *ColOPA, *alpha);
		//dgemm("N", "N", RowOPA, ColOPB, ColOPA, alpha, Dt, LDA, H + n*Np*K, LDB, Beta, tempdest, LDC);
		/*$tz\cdot Dt*H$*/
		DotProduct(tempdest, tempdest, tz, (int)(*LDC));
		/*$rx\cdot Dr*E + sx\cdot Ds*E + ry\cdot Dr*G + sy\cdot Ds*G + tz\cdot Dt*H$*/
		Add(dest + n*Np*K, dest + n*Np*K, tempdest, (int)(*LDC));
	}
}

/*This function is used to get the inverse matrix of the given matrix, and has been verified*/

//void MatrixInverse(double *dest, ptrdiff_t Np)  //chenzereng int or long int or longlongint
void MatrixInverse(double *dest, int Np)
{
	int info;
	int *IPIV = new int[Np]; //(int *)malloc(Np*sizeof(int));
	double *tempdest = new double[Np];
	//int LWORK = Np * Np;
	//int IPIV[50];
	//int WORK[2500];
	//char trans = 'N';
	//int nrhs = 1;
	//double *WORK = (double *)malloc(LWORK*sizeof(double));
	//dgetrf_(&Np, &Np, dest, &Np, IPIV, &info);
	//dgetri_(&Np, dest, &Np, IPIV, WORK, &LWORK, &info);

	//cout << "In the main function, before calling LAPACK_dgetri" << endl;
	//LAPACKE_dgetrf(LAPACK_COL_MAJOR, Np, Np, dest, Np, IPIV);
	LAPACK_dgetrf(&Np, &Np, dest, &Np, IPIV, &info);
	/*for (int i = 0; i < Np*Np; i++) {
		cout << dest[i] << endl;
	}
	cout << "In the main function, after calling LAPACK_dgetri" << endl;*/
	//info = LAPACKE_dgetri(LAPACK_COL_MAJOR, Np, dest, Np, IPIV);

	LAPACK_dgetri(&Np, dest, &Np, IPIV, tempdest,&Np, &info);
	/*for (int i = 0; i < Np*Np; i++) {
		cout << dest[i] << endl;
	}*/
	delete[] IPIV;
	delete[] tempdest;
	//free(IPIV);
	//free(WORK);
}

//void MatrixInverseTest()
//{
//	setlocale(LC_ALL, "");
//	double a[] =
//	{
//		3,-1,-1,
//		4,-2,-1,
//		-3,2,1
//	};
//	int m = 3;
//	int n = 3;
//	int lda = 3;
//	int ipiv[3];
//	int info;
//	cout << "Before calling LAPACK_dgetri" << endl;
//	//info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, a, lda, ipiv);
//	LAPACK_dgetrf(&m, &n, a, &lda, ipiv,&info);
//	for (int i = 0; i < 9; i++) {
//		cout << a[i] << endl;
//	}
//	cout << "After calling LAPACK_dgetri" << endl;
//	//info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, m, a, lda, ipiv);
//	double *b = new double[m];
//	LAPACK_dgetri(&m, a, &lda, ipiv, b, &n, &info);
//	for (int i = 0; i < 9; i++) {
//		cout << a[i] << endl;
//	}
//	delete[] b;
//}

/* This function invokes dgemm implemented by blas to calculate C = alpha*A*B + beta*C */
//void MatrixMultiply(char* TRANSA, char* TRANSB, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, double ALPHA, double *A,
//	ptrdiff_t LDA, double *B, ptrdiff_t LDB, double BETA, double *C, ptrdiff_t LDC)	{
//
//	dgemm(TRANSA, TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
//
//}

void MatrixMultiply(double *matrix1, double *matrix2, double *matrix3, const int M_, const int N_, const int K_, const double Alpha)
{
	//C=AB,A=M*K,B=K*N,C=M*N
	//const enum CBLAS_ORDER Order = CblasColMajor;
	//const enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	//const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
	//const int MM = M_;
	//const int NN = N_;
	//const int KK = K_;
	//const double alpha = Alpha;
	//const double beta = 0.0;
	//const int lda = MM;//max(MM,KK)
	//const int ldb = KK;//max(NN,KK)
	//const int ldc = MM;//max(MM,NN)
	//cblas_dgemm(Order, TransA, TransB, MM, NN, KK, alpha, matrix1, lda, matrix2, ldb, beta, matrix3, ldc);
//-------------------------------------------------------------------------------------------------------------------------
	int i, j, k;
	double *A1 = (double *)malloc((int)M_ * (int)K_ * sizeof(double));
	double *B1 = (double *)malloc((int)K_ * (int)N_ * sizeof(double));
	double *C1 = (double *)malloc((int)M_ * (int)N_ * sizeof(double));
	memset(A1, 0, (int)M_ * (int)K_ * sizeof(double));
	memset(B1, 0, (int)K_ * (int)N_ * sizeof(double));
	memset(C1, 0, (int)M_ * (int)N_ * sizeof(double));

	for (i = 0; i < (int)M_; i++) {
		for (k = 0; k < (int)K_; k++) {
			A1[i * K_ + k] = matrix1[k * M_ + i];
		}
	}

	for (k = 0; k < (int)K_; k++) {
		for (j = 0; j < (int)N_; j++) {
			B1[k * N_ + j] = Alpha * matrix2[j * K_ + k];
		}
	}

	for (i = 0; i < (int)M_; i++) {
		for (j = 0; j < (int)N_; j++) {
			for (k = 0; k < (int)K_; k++) {
				C1[j * M_ + i] = C1[j * M_ + i] + A1[i * K_ + k] * B1[k * N_ + j];
			}
		}
	}

	for (i = 0; i < (int)M_; i++) {
		for (j = 0; j < (int)N_; j++) {
			matrix3[i * N_ + j] = C1[i * N_ + j];
		}
	}

	free(A1); A1 = NULL;
	free(B1); B1 = NULL;
	free(C1); C1 = NULL;
}

void MatrixMultiplyTN(double *matrix1, double *matrix2, double *dest, const int M_, const int N_, const int K_, const double Alpha)
{
	//C=AB,A=M*K,B=K*N,C=M*N
	//const enum CBLAS_ORDER Order = CblasColMajor;
	//const enum CBLAS_TRANSPOSE TransA = CblasTrans;
	//const enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
	//const int MM = M_;
	//const int NN = N_;
	//const int KK = K_;
	//const double alpha = Alpha;
	//const double beta = 0.0;
	//const int lda = KK;//max(MM,KK)
	//const int ldb = KK;//max(NN,KK)
	//const int ldc = MM;//max(MM,NN)
	//cblas_dgemm(Order, TransA, TransB, MM, NN, KK, alpha, matrix1, lda, matrix2, ldb, beta, dest, ldc);
//-------------------------------------------------------------------------------------------------------------------------
	int i, j, k;
	double *A2 = (double *)malloc((int)M_ * (int)K_ * sizeof(double));
	double *B2 = (double *)malloc((int)K_ * (int)N_ * sizeof(double));
	double *C2 = (double *)malloc((int)M_ * (int)N_ * sizeof(double));
	memset(A2, 0, (int)M_ * (int)K_ * sizeof(double));
	memset(B2, 0, (int)K_ * (int)N_ * sizeof(double));
	memset(C2, 0, (int)M_ * (int)N_ * sizeof(double));

	for (i = 0; i < (int)M_; i++) {
		for (k = 0; k < (int)K_; k++) {
			A2[k * M_ + i] = matrix1[k * M_ + i];
		}
	}

	for (k = 0; k < (int)K_; k++) {
		for (j = 0; j < (int)N_; j++) {
			B2[k * N_ + j] = Alpha * matrix2[j * K_ + k];
		}
	}

	for (i = 0; i < (int)M_; i++) {
		for (j = 0; j < (int)N_; j++) {
			for (k = 0; k < (int)K_; k++) {
				C2[j * M_ + i] = C2[j * M_ + i] + A2[i * K_ + k] * B2[k * N_ + j];
			}
		}
	}

	for (i = 0; i < (int)M_; i++) {
		for (j = 0; j < (int)N_; j++) {
			dest[i * N_ + j] = C2[i * N_ + j];
		}
	}

	free(A2); A2 = NULL;
	free(B2); B2 = NULL;
	free(C2); C2 = NULL;
}

void MatrixMultiplyTT(double *matrix1, double *matrix2, double *dest, const int M_, const int N_, const int K_, const double Alpha)
{
	//C=AB,A=M*K,B=K*N,C=M*N
	//const enum CBLAS_ORDER Order = CblasColMajor;
	//const enum CBLAS_TRANSPOSE TransA = CblasTrans;
	//const enum CBLAS_TRANSPOSE TransB = CblasTrans;
	//const int MM = M_;
	//const int NN = N_;
	//const int KK = K_;
	//const double alpha = Alpha;
	//const double beta = 0.0;
	//const int lda = KK;//max(MM,KK)
	//const int ldb = NN;//max(NN,KK)
	//const int ldc = MM;//max(MM,NN)
	//cblas_dgemm(Order, TransA, TransB, MM, NN, KK, alpha, matrix1, lda, matrix2, ldb, beta, dest, ldc);
//-------------------------------------------------------------------------------------------------------------------------
	int i, j, k;
	double *A3 = (double *)malloc(M_ * K_ * sizeof(double));
	double *B3 = (double *)malloc(K_ * N_ * sizeof(double));
	memset(dest, 0, M_ * N_ * sizeof(double));

	for (i = 0; i < M_; i++) {
		for (k = 0; k < K_; k++) {
			A3[k * M_ + i] = Alpha * matrix1[k * M_ + i];
		}
	}

	for (k = 0; k < K_; k++) {
		for (j = 0; j < N_; j++) {
			B3[j * K_ + k] = matrix2[j * K_ + k];
		}
	}

	for (i = 0; i < M_; i++) {
		for (j = 0; j < N_; j++) {
			for (k = 0; k < K_; k++) {
				dest[j * M_ + i] = dest[j * M_ + i] + A3[i * K_ + k] * B3[k * N_ + j];
			}
		}
	}
	free(A3); A3 = NULL;
	free(B3); B3 = NULL;
}

void Minus(double *dest, double *sourcea, double *sourceb, int size){
	for (int i = 0; i < size; i++)
		dest[i] = sourcea[i] - sourceb[i];
}

void MultiEdgeContributionByLiftOperator(double *SrcAndDest, double *TempSource, int *RowOPA, int *ColOPB, int *ColOPA, double *Alpha, \
	double *A, int *LDA, int *LDB, double *Beta, int *LDC, double *J, int Np){

	MatrixMultiply(A, SrcAndDest, TempSource, *RowOPA, *ColOPB, *ColOPA, *Alpha);
	/*dgemm("N", "N", RowOPA, ColOPB, ColOPA, Alpha, A, LDA, SrcAndDest, LDB, Beta, TempSource, \
		LDC);*/
	DotDivide(SrcAndDest, TempSource, J, Np);
}

void MultiplyByConstant(double *dest, double *Source, double Coefficient, int Np){
	for (int i = 0; i < Np; i++)
		dest[i] = Source[i] * Coefficient;
}

void NdgExtend2dField(double *dest, double *source, int Np2d, int Index, int Np3d, int NLayer, int Nz){
	//直接把第一层循环中的k提出来，放置到调用层，利用openmp并行，传入的K2d替换成序号k
	//for (int k = 0; k < K2d; k++){
	for (int Layer = 0; Layer < NLayer; Layer++){
		for (int N = 0; N < Nz + 1; N++){
			for (int i = 0; i < Np2d; i++){
				dest[Index*NLayer*Np3d + Layer*Np3d + N*Np2d + i] = \
					source[Index*Np2d + i];
			}
		}
	}
}

void RepmatValue(double *dest, double *source, int size, int Npz){
	for (int i = 0; i < Npz; i++){
		for (int j = 0; j < size; j++){
			dest[i*size + j] = source[j];
		}
	}
}
/*Note: source[i] equals to zero is not excluded*/
void ReverseValue(double *dest, double *source, int size){
	for (int i = 0; i < size; i++)
		dest[i] = 1.0 / source[i];
}

void StrongFormBoundaryEdgeRHS(int edgeIndex, double *FToE, double *FToF, int Np, int K, \
	int Nfp, double *FToN1, double *fluxM_, double *fluxS_, double *Js, double *Mb, double *rhs){
	const int e1 = (int)FToE[2 * edgeIndex] - 1;
	const int f1 = (int)FToF[2 * edgeIndex] - 1;
	const int ind1 = e1 * Np - 1 + f1 * Np * K;
	const int ind = edgeIndex * Nfp;
	double *rhsM = (double *)malloc(Nfp*sizeof(double));
	memset(rhsM, 0, Nfp*sizeof(double));
	for (int n = 0; n < Nfp; n++) {
		const int sk = n + ind;
		double dfM = fluxM_[sk] - fluxS_[sk];
		double j = Js[sk];
		double *mb = Mb + n * Nfp;
		for (int m = 0; m < Nfp; m++) {
			rhsM[m] += mb[m] * j * dfM;
		}
	}
	for (int n = 0; n < Nfp; n++) {
		const int sk = n + ind;
		const int m1 = (int)FToN1[sk] + ind1;
		rhs[m1] += rhsM[n];
	}
	free(rhsM);
}

void StrongFormInnerEdgeRHS(int edgeIndex, double *FToE, double *FToF, int Np, int K,\
	int Nfp, double *FToN1, double *FToN2, double *fluxM_, double *fluxP_,\
	                  double *fluxS_, double *Js, double *Mb, double *rhs){
	const int e1 = (int)FToE[2 * edgeIndex] - 1;
	const int e2 = (int)FToE[2 * edgeIndex + 1] - 1;
    const int f1 = (int)FToF[2 * edgeIndex] - 1;
    const int f2 = (int)FToF[2 * edgeIndex + 1] - 1;
	const int ind1 = e1 * Np - 1 + f1 * Np * K;
	const int ind2 = e2 * Np - 1 + f2 * Np * K;
	const int ind = edgeIndex * Nfp;
	double *rhsM = (double *)malloc(Nfp*sizeof(double));
	double *rhsP = (double *)malloc(Nfp*sizeof(double));
	memset(rhsM, 0, Nfp*sizeof(double));
	memset(rhsP, 0, Nfp*sizeof(double));
	for (int n = 0; n < Nfp; n++) {
		const int sk = n + ind;
		double dfM = fluxM_[sk] - fluxS_[sk];
		double dfP = fluxP_[sk] - fluxS_[sk];
		double j = Js[sk];
		double *mb = Mb + n * Nfp;
		for (int m = 0; m < Nfp; m++) {
			rhsM[m] += mb[m] * j * dfM;
			rhsP[m] -= mb[m] * j * dfP;
		}
	}
	for (int n = 0; n < Nfp; n++) {
	const int sk = n + ind;
	const int m1 = (int)FToN1[sk] + ind1;
	const int m2 = (int)FToN2[sk] + ind2;
	rhs[m1] += rhsM[n];
	rhs[m2] += rhsP[n];
	}
	free(rhsM);
	free(rhsP);
}