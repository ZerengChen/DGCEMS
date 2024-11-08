#ifndef _NdgMath_H
#define _NdgMath_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lapacke.h"
#include <string.h>
#include <cblas.h>
#include <stdlib.h>
#include <cmath>

/*Note: this function is used to assemble the element mass matrix and the physical diff matrix, and has been verified*/
void DiagMultiply(double *, const double *, const double *, int );

void DiagRightMultiply(double *, const double *, const double *, int );

void Exchange(double *fphys, double *IEFToE, int *pE, int IENe, int Np, int K, int Nfield);//For MPI exchange

void Exchange(signed char *Status, double *IEFToE, int *pE, int IENe, int Np, int K, int Nfield);//For MPI exchange

void GatherToZero(double *fphys, int *pE, int Np, int K, int Nfield);//For MPI exchange

void BroadcastToAll(double *fphys, int *pE, int Np, int K, int Nfield);//For MPI exchange

void Add(double *, double *, double *, int );

void AddByConstant(double *, double *, double, int);

void AssembleDataIntoPoint(double *, double *, double *, int);

void Minus(double *, double *, double *, int);

void MatrixMultiply(double *, double *, double *, const int , const int, const int, const double);

void MatrixMultiplyTN(double *, double *, double *, const int, const int, const int, const double );

void MatrixMultiplyTT(double *, double *, double *, const int, const int, const int, const double);

void MatrixInverse(double *,int );

void DotProduct(double *, double *, double *, int );

void DotDivide(double *, double *, double *, int );

void DotDivideByConstant(double *, double *, double , int);

void MultiplyByConstant(double *, double *, double, int);

void DotCriticalDivide(double *, double *, double *, double *, int);

void StrongFormInnerEdgeRHS(int , double *, double *, int , int ,int , double *, double *, double *, \
	double *,double *, double *, double *, double *);

void StrongFormBoundaryEdgeRHS(int , double *, double *, int , int ,  int , double *, double *, double *, \
	double *, double *, double *);

void FetchInnerEdgeFacialValue(double *, double *, double *, double *, double *, double *, int , int );

void FetchBoundaryEdgeFacialValue(double *, double *, double *, double *, int , int );

void Flip(double *, int);

void ReverseValue(double *, double *, int );

void RepmatValue(double *, double *, int , int );

void AssembleContributionIntoColumn(double *, double *, double *, int, int);

void AssembleContributionIntoRow(double *, double *, double *, int, int);

void AssembleContributionIntoRowAndColumn(double *, double *, double *, double *, int, int, int);

void NdgExtend2dField(double *, double *, int, int, int, int, int);

void GetVolumnIntegral1d(double *, int *, int *, int *, double *, \
	double *, int *, double *, int *, double *, int *, double *);

void GetVolumnIntegral2d(double *, double *, int *, int *, int *, double *, double *, \
	double *, int *, double *, int *, double *, int *, double *, double *);

void GetVolumnIntegralOnlyVertical(double *, int *, int *, int *, double *, \
	double *, int *, double *, int *, double *, int *, double *);

void GetVolumnIntegral3d(double *, double *, int *, int *, int *, double *, \
	double *, double *, double *, double *, double *, double *, int *, int *, double *, int *, \
	double *, double *, double *, double *, double *, int , int , int );

void GetFacialFluxTerm2d(double *, double *, double *, double *, double *, int);

void GetScalarFluxTerm2d(double *, double *, double *, int);

void MultiEdgeContributionByLiftOperator(double *, double *, int *, int *, int *, double *, \
	double *, int *, int *, double *, int *, double *, int);

void GetMeshIntegralValue(double *, int *, int *, int *, double *, double *, double *, int *, double *);

void GetMeshAverageValue(double *, double *, int *, int *, int *, double *, double *, double *, int *, double *);

//void MatrixInverseTest();

#endif