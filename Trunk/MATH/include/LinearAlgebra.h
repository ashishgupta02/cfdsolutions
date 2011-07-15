/*******************************************************************************
 * File:        LinearAlgebra.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _LINEARALGEBRA_H
#define _LINEARALGEBRA_H

/*******************************************************************************
 * Keep C++ compilers from getting confused
 *******************************************************************************/
#if defined __cplusplus
extern "C" {
#endif

    /* Ly = Pb */
    void Solve_PivotForwardSubstitution(double **L, double *y, double *b, int *P, int size);
    /* Pivot Ux = y */
    void Solve_PivotBackSubstitution(double **U, double *x, double *y, int *P, int size);
    /* Ly = b */
    void Solve_ForwardSubstitution(double **L, double *y, double *b, int size);
    /* Ux = y */
    void Solve_BackSubstitution(double **U, double *x, double *y, int size);
    /* Cmxn = Amxp*Bpxn */
    void MatrixMultiplyMatrix(double **A, double **B, double **C, int m, int p, int n);
    /* C = A*B */
    void MatrixMultiplyMatrixSquare(double **A, double **B, double **C, int size);
    /* Transpose of A Matrix */
    void MatrixTransposeSquare(double **A, int size);
    /* Gaussain PA = LU */
    void GaussainLUDecompose(double **A, int *P, int size, int pivot_flag);
    /* Householder A = QR */
    void HouseholderQRDecompose(double **A, double *h, int m, int n, double *b);
    /* GramSchmidt A = QR */
    void GramSchmidtQRDecompose(double **A, int m, int n, double *b, double *c);
    void GramSchmidtQRDecomposeFull(double **A, double **Q, double **R, int m, int n);
    /* Eigenvalue QR Iteration */
    void QR_Iteration(double **A, int size, int MaxIter);
    /* Hankel Matrix */
    void CreateHankelMatrix(double **A, int size);
    /* Sparse Matrix Vector Multiple y = Ax */
    void SparseMatrixVectorMultiply(double *A, double *x, double *y, int *IA, int *JA, int *IAU, int NRows, int NCols) ;

/*******************************************************************************
 * Keep C++ compilers from getting confused
 *******************************************************************************/
#if defined __cplusplus
}
#endif

#endif /* _LINEARALGEBRA_H */
