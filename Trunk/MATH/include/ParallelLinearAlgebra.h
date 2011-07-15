/*******************************************************************************
 * File:        ParallelLinearAlgebra.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _PARALLELLINEARALGEBRA_H
#define _PARALLELLINEARALGEBRA_H

#ifdef HAVE_MPI

struct double_int{
  double d;
  int    i;
};

/* Parallel Space Matrix Vector Multiply */
void ParallelSparseMatrixVectorMultiple(double *A, double *x, double *y, int *IA, int *JA, int *IAU, int NRows, int NCols, int Rank, int NProc);
/* Parallel Gaussain PA = LU with Block-Cyclic Partition */
void ParallelGaussainLUDecompose(double **A, int *P, int *Pinv, int Dim, int ROOT, int Rank, int NProc, int RecursionLevel);
/* Parallel Sparse Arnodi process for Eigenvalue Computation */
void ParallelSparseArnoldi(double *CRSMatrix, int *IA, int *JA, int *IAU, int NRows, int NCols, int KEnd, int Iterations, int Rank, int ROOT, int NProc);
/* Parallel Sparse GMRES for Linear System Equation Ax = b */
void ParallelSparseGMRES(double *CRSMatrix, double *X, double *RHS, int *IA, int *JA, int *IAU, int Dim, int NRows, int NCols, int KEnd, int NRestarts, int Precondition, int Rank, int ROOT, int NProc);

#endif /* HAVE_MPI */

#endif /* _PARALLELLINEARALGEBRA_H */
