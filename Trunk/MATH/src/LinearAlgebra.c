/*******************************************************************************
 * File:        LinearAlgebra.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Application Header Files */
#include "Trim_Utils.h"
#include "LinearAlgebra.h"

//------------------------------------------------------------------------------
//! Ly = Pb
//------------------------------------------------------------------------------
void Solve_PivotForwardSubstitution(double **L, double *y, double *b, int *P, int size) {
    int i, j;
    y[0] = b[P[0]];
    for (i = 1; i < size; i++) {
        y[i] = b[P[i]];
        for (j = 0; j < i; j++)
            y[i] -= y[j] * ((i == j) ? 1.0 : L[P[i]][j]);
    }
}

//------------------------------------------------------------------------------
//! Pivot Ux = y
//------------------------------------------------------------------------------
void Solve_PivotBackSubstitution(double **U, double *x, double *y, int *P, int size) {
    int i, j;
    x[size - 1] = y[size - 1] / U[P[size - 1]][size - 1];
    for (i = size - 2; i >= 0; i--) {
        x[i] = y[i];
        for (j = size - 1; j > i; j--)
            x[i] -= x[j] * U[P[i]][j];
        x[i] /= U[P[i]][i];
    }
}

//------------------------------------------------------------------------------
//! Ly = b
//------------------------------------------------------------------------------
void Solve_ForwardSubstitution(double **L, double *y, double *b, int size) {
  int i, j;
  y[0] = b[0];
  for (i = 1; i < size; i++) {
    y[i] = b[i];
    for (j = 0; j < i; j++)
      y[i] -= y[j] * ((i == j) ? 1.0 : L[i][j]);
  }
}

//------------------------------------------------------------------------------
//! Ux = y
//------------------------------------------------------------------------------
void Solve_BackSubstitution(double **U, double *x, double *y, int size) {
  int i, j;
  x[size - 1] = y[size - 1] / U[size - 1][size - 1];
  for (i = size - 2; i >= 0; i--) {
    x[i] = y[i];
    for (j = size - 1; j > i; j--)
      x[i] -= x[j] * U[i][j];
    x[i] /= U[i][i];
  }
}

//------------------------------------------------------------------------------
//! Cmxn = Amxp*Bpxn
//------------------------------------------------------------------------------
void MatrixMultiplyMatrix(double **A, double **B, double **C, int m, int p, int n) {
    int i, j, k;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (k = 0; k < p; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
}

//------------------------------------------------------------------------------
//! C = A*B
//------------------------------------------------------------------------------
void MatrixMultiplyMatrixSquare(double **A, double **B, double **C, int size) {
    int i, j, k;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            C[i][j] = 0.0;
            for (k = 0; k < size; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
    }
}

//------------------------------------------------------------------------------
//! Transpose of A Matrix
//------------------------------------------------------------------------------
void MatrixTransposeSquare(double **A, int size) {
    int i, j;
    double temp;
    for (i = 0; i < size; i++) {
        for(j = i + 1; j < size; j++) {
            temp = A[j][i];
            A[j][i] = A[i][j];
            A[i][j] = temp;
        }
    }
}

//------------------------------------------------------------------------------
//! Gaussain PA = LU
//------------------------------------------------------------------------------
void GaussainLUDecompose(double **A, int *P, int size, int pivot_flag) {
    int i, j, c, tmp, pivot_index;
    double pivot, pivot_max;

    // Initialize Permutation Vector P
    for (i = 0; i < size; i++)
        P[i] = i;

    // Looping on Columns of A
    for (c = 0; c < size; c++) {
        if (pivot_flag) {
            pivot = A[P[c]][c];
            pivot_index = c;
            pivot_max = pivot;
            // Search Row for Pivot
            for (i = c + 1; i < size; i++) {
                if (fabs(A[P[i]][c]) > fabs(pivot_max)) {
                    pivot_max = A[P[i]][c];
                    pivot_index = i;
                }
            }
            // Check for Pivot Swap
            if (pivot_index != c) {
                tmp = P[c];
                P[c] = P[pivot_index];
                P[pivot_index] = tmp;
            }
        }

        // Compute multipliers
        for (i = c + 1; i < size; i++)
            A[P[i]][c] /= A[P[c]][c];

        for (i = c + 1; i < size; i++) {
            for (j = c + 1; j < size; j++)
                A[P[i]][j] -= A[P[i]][c] * A[P[c]][j];
        }
    }
}

//------------------------------------------------------------------------------
//! Householder A = QR
// A will be overwritten with R in upper triangle (+ diagonal),
// nonzero portions of householder vectors in lower triangle.
// h will hold first nonzero entry of each Householder vector
// (to avoid overwriting diagonal of R).
//------------------------------------------------------------------------------
void HouseholderQRDecompose(double **A, double *h, int m, int n, double *b) {
    int i, j, k;
    double *v = NULL;
    double alpha, beta, gamma, trans_scalar;

    v = (double *) malloc(m*sizeof(double));
    for (i = 0; i < m; i++)
        h[i] = 0.0;

    // loop over columns
    for (k = 0; k < n; k++){
        alpha = 0.0;
        for (i = 0; i < m; i++)
            v[i] = 0.0;

        for (i = k; i < m; i++)
            alpha += A[i][k]*A[i][k];
        alpha = sqrt(alpha);

        if (A[k][k] > 0)
            alpha *= -1.0;
        v[k] = A[k][k] - alpha;

        for (i = (k+1); i < m; i++)
            v[i] = A[i][k];

        beta = 0.0;
        for (i = k; i < m; i++)
            beta += v[i]*v[i];
        if (fabs(beta) < DBL_TOLERANCE)
            continue;

        // apply to rest of the columns of A
        for (j = k; j < n; j++) {
            gamma = 0.0;
            for(i = k; i < m; i++)
                gamma += v[i]*A[i][j];
            trans_scalar = (2.0*gamma/beta);
            for(i = k; i < m; i++)
                A[i][j] -= trans_scalar * v[i];
        }

        // apply to b
        if (b != NULL) {
            gamma = 0.0;
            for(i = k; i < m; i++)
                gamma += v[i]*b[i];
            trans_scalar = (2.0*gamma/beta);
            for (i = k; i < m; i++)
                b[i] -= trans_scalar * v[i];
        }

        h[k] = v[k];
        for (i = k+1; i < m; i++)
            A[i][k] = v[i];
    }

    if (v != NULL)
        free(v);
}

//------------------------------------------------------------------------------
//! GramSchmidt A = QR
//! A is replace with R
//! Return R overwritten in A, and destroy b!
//! c will contain orthogonalized b for Rx=c solve
//------------------------------------------------------------------------------
void GramSchmidtQRDecompose(double **A, int m, int n, double *b, double *c) { 
    int i, j, k;
    double **R = NULL;

    R = (double **) malloc(m*sizeof(double*));;
    for (i = 0; i < m; i++) {
        R[i] = (double *) malloc(n*sizeof(double));
        for (j = 0; j < n; j++)
            R[i][j] = 0.0;
    }

    for (k = 0; k < n; k++) {
        R[k][k] = 0.0;
        for (i = 0; i < m; i++)
            R[k][k] += pow(A[i][k],2);
        R[k][k] = sqrt(R[k][k]);

        // Check if A is not Linearly Independent
        if (R[k][k] < DBL_TOLERANCE)
            return;

        for (i = 0; i < m; i++)
            A[i][k] = A[i][k]/R[k][k];

        for (j = k+1; j < n; j++) {
            // R[k][j] is dot product of this column of Q (transposed) and A
            for (i = 0; i < m; i++)
                R[k][j] += A[i][k]*A[i][j];

            for (i = 0; i < m; i++)
                A[i][j] -= R[k][j]*A[i][k];
        }

        if (b != NULL && c != NULL) {
            c[k] = 0.0;
            for (i = 0; i < m; i++)
                c[k] += A[i][k]*b[i];

            // b_n = b_n-1 - q_k'*b_n-1*q_k
            for (i = 0; i < m; i++)
                b[i] -= c[k]*A[i][k];
        }
    }

    for (i = 0; i < m; i++) {
        for(j = 0; j < n; j++)
            A[i][j] = R[i][j];
    }

    if (R != NULL) {
        for (i = 0; i < m; i++)
            if (R[i] != NULL)
                free(R[i]);
        free(R);
    }
}

//------------------------------------------------------------------------------
//! GramSchmidt A = QR
//! A is reconstructed
//------------------------------------------------------------------------------
void GramSchmidtQRDecomposeFull(double **A, double **Q, double **R, int m, int n) {
    int i, j, k;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            R[i][j] = 0.0;

    for (k = 0; k < n; k++) {
        R[k][k] = 0.0;
        for (i = 0; i < m; i++)
            R[k][k] += pow(A[i][k], 2);
        R[k][k] = sqrt(R[k][k]);

        // Check if A is not Linearly Independent
        if (R[k][k] < DBL_TOLERANCE)
            return;

        for (i = 0; i < m; i++) {
            A[i][k] = A[i][k] / R[k][k];
            Q[i][k] = A[i][k];
        }

        for (j = k + 1; j < n; j++) {
            // R[k][j] is dot product of this column of Q (transposed) and A
            for (i = 0; i < m; i++)
                R[k][j] += A[i][k] * A[i][j];

            for (i = 0; i < m; i++)
                A[i][j] -= R[k][j] * A[i][k];
        }
    }

    // Reconstruct A
    MatrixMultiplyMatrix(Q, R, A, m, n, n);
}

//------------------------------------------------------------------------------
//! Eigenvalue Computation: QR Iteration
//------------------------------------------------------------------------------
void QR_Iteration(double **A, int size, int MaxIter) {
    if (size <= 1)
        return;

    int i, iter, wn, converged;
    double shift;
    double **Q = NULL;
    double **R = NULL;
    Q = (double **) malloc(size*sizeof(double *));
    R = (double **) malloc(size*sizeof(double *));
    for (i = 0; i < size; i++) {
        Q[i] = (double *) malloc(size*sizeof(double));
        R[i] = (double *) malloc(size*sizeof(double));
    }

    iter = 0;
    wn   = size;
    converged = 0;
    while (iter < MaxIter && !converged) {
        shift = A[wn - 1][wn - 1];
        for (i = 0; i < wn; i++)
            A[i][i] -= shift;

        GramSchmidtQRDecomposeFull(A, Q, R, wn, wn);
        MatrixMultiplyMatrixSquare(R, Q, A, wn);

        for (i = 0; i < wn; i++)
            A[i][i] += shift;

        // Check for convergence: Monitor A[size][size-1]
        if (fabs(A[wn - 1][wn - 2]) < DBL_ZERO) {
            wn--;
            if (wn <= 1)
                converged = 1;
        }

        iter++;
    }

    // Check if QR Iteration Fails to converge
    if (converged == 0)
        warn("QR_Iteration: Convegence Failed at Max Iteration at size %d", wn);

    for (i = 0; i < size; i++) {
        free(Q[i]);
        free(R[i]);
    }
    free(Q);
    free(R);
}

//------------------------------------------------------------------------------
//! Hankel Matrix
//------------------------------------------------------------------------------
void CreateHankelMatrix(double **A, int size) {
    int i, j;
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            A[i][j] = 1.0 / (2.0 * (size - (i + 1) - (j + 1) + 3.0/2.0));
}

//------------------------------------------------------------------------------
//! Sparse Matrix Vector Multiple
//! y = Ax
//------------------------------------------------------------------------------
void SparseMatrixVectorMultiply(double *A, double *x, double *y, int *IA, int *JA, int *IAU, int NRows, int NCols) {
    int i, idx, j;

    for (i = 0; i < NRows; i++) {
        y[i] = 0;
        for (idx = IA[i]; idx <= IAU[i]; idx++) {
            j = JA[idx];
            if (j < NCols)
                y[i] += A[idx] * x[j];
        }
    }
}
