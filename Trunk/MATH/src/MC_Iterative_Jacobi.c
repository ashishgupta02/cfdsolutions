/*******************************************************************************
 * File:        MC_Iterative_Jacobi.c
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "MC.h"

//------------------------------------------------------------------------------
//! Currently Inverse of 2x2, 3x3 and 4x4 is supported
//------------------------------------------------------------------------------
double MC_Determinant(int nRow, int nCol, double **A) {
    double det = 0.0;
    if ((nRow <= 0) || (nCol <= 0) || (nRow != nCol) || (A == NULL))
        return det;
    
    switch (nRow) {
        case 2:
            det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
            break;
        case 3:
            det = A[0][0]*(A[1][1] * A[2][2] - A[1][2] * A[2][1])
                    + A[0][1]*(A[1][2] * A[2][0] - A[1][0] * A[2][2])
                    + A[0][2]*(A[1][0] * A[2][1] - A[1][1] * A[2][0]);
            break;
        case 4:
            det =   (A[0][0] * (A[1][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2])
                            - A[1][2] * (A[2][1] * A[3][3] - A[2][3] * A[3][1])
                            + A[1][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1]))
                    - A[0][1] * (A[1][0] * (A[2][2] * A[3][3] - A[2][3] * A[3][2])
                            - A[1][2] * (A[2][0] * A[3][3] - A[2][3] * A[3][0])
                            + A[1][3] * (A[2][0] * A[3][2] - A[2][2] * A[3][0]))
                    + A[0][2] * (A[1][0] * (A[2][1] * A[3][3] - A[2][3] * A[3][1])
                            - A[1][1] * (A[2][0] * A[3][3] - A[2][3] * A[3][0])
                            + A[1][3] * (A[2][0] * A[3][1] - A[2][1] * A[3][0]))
                    - A[0][3] * (A[1][0] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])
                            - A[1][1] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])
                            + A[1][2] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])));
            break;
        default:
            break;
    }

    return det;
}

//------------------------------------------------------------------------------
//! Currently Inverse of 2x2, 3x3 and 4x4 is supported
//------------------------------------------------------------------------------
void MC_Inverse(int nRow, int nCol, double **A, double **InvA) {
    double det = 0.0;
    if ((nRow <= 0) || (nCol <= 0) || (nRow != nCol) ||
            (A == NULL) || (InvA == NULL))
        return;
    
    det = MC_Determinant(nRow, nCol, A);

    switch (nRow) {
        case 2:
            InvA[0][0] = +A[1][1] / det;
            InvA[0][1] = -A[0][1] / det;
            InvA[1][0] = -A[1][0] / det;
            InvA[1][1] = +A[0][0] / det;
            break;
        case 3:
            InvA[0][0] = +(A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
            InvA[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / det;
            InvA[0][2] = +(A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;
            InvA[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / det;
            InvA[1][1] = +(A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
            InvA[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / det;
            InvA[2][0] = +(A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
            InvA[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / det;
            InvA[2][2] = +(A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;
            break;
        case 4:
            InvA[0][0] = +(
                            +A[1][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2])
                            +A[1][2] * (A[2][3] * A[3][1] - A[2][1] * A[3][3])
                            +A[1][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])
                        ) / det;
            InvA[1][0] = -(
                            +A[1][0] * (A[2][2] * A[3][3] - A[2][3] * A[3][2])
                            +A[1][2] * (A[2][3] * A[3][0] - A[2][0] * A[3][3])
                            +A[1][3] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])
                        ) / det;
            InvA[2][0] = +(
                            +A[1][0] * (A[2][1] * A[3][3] - A[2][3] * A[3][1])
                            +A[1][1] * (A[2][3] * A[3][0] - A[2][0] * A[3][3])
                            +A[1][3] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])
                        ) / det;
            InvA[3][0] = -(
                            +A[1][0] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])
                            +A[1][1] * (A[2][2] * A[3][0] - A[2][0] * A[3][2])
                            +A[1][2] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])
                        ) / det;
            InvA[0][1] = -(
                            +A[0][1] * (A[2][2] * A[3][3] - A[2][3] * A[3][2])
                            +A[0][2] * (A[2][3] * A[3][1] - A[2][1] * A[3][3])
                            +A[0][3] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])
                        ) / det;
            InvA[1][1] = +(
                            +A[0][0] * (A[2][2] * A[3][3] - A[2][3] * A[3][2])
                            +A[0][2] * (A[2][3] * A[3][0] - A[2][0] * A[3][3])
                            +A[0][3] * (A[2][0] * A[3][2] - A[2][2] * A[3][0])
                        ) / det;
            InvA[2][1] = -(
                            +A[0][0] * (A[2][1] * A[3][3] - A[2][3] * A[3][1])
                            +A[0][1] * (A[2][3] * A[3][0] - A[2][0] * A[3][3])
                            +A[0][3] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])
                        ) / det;
            InvA[3][1] = +(
                            +A[0][0] * (A[2][1] * A[3][2] - A[2][2] * A[3][1])
                            +A[0][1] * (A[2][2] * A[3][0] - A[2][0] * A[3][2])
                            +A[0][2] * (A[2][0] * A[3][1] - A[2][1] * A[3][0])
                        ) / det;
            InvA[0][2] = +(
                            +A[0][1] * (A[1][2] * A[3][3] - A[1][3] * A[3][2])
                            +A[0][2] * (A[1][3] * A[3][1] - A[1][1] * A[3][3])
                            +A[0][3] * (A[1][1] * A[3][2] - A[1][2] * A[3][1])
                        ) / det;
            InvA[1][2] = -(
                            +A[0][0] * (A[1][2] * A[3][3] - A[1][3] * A[3][2])
                            +A[0][2] * (A[1][3] * A[3][0] - A[1][0] * A[3][3])
                            +A[0][3] * (A[1][0] * A[3][2] - A[1][2] * A[3][0])
                        ) / det;
            InvA[2][2] = +(
                            +A[0][0] * (A[1][1] * A[3][3] - A[1][3] * A[3][1])
                            +A[0][1] * (A[1][3] * A[3][0] - A[1][0] * A[3][3])
                            +A[0][3] * (A[1][0] * A[3][1] - A[1][1] * A[3][0])
                        ) / det;
            InvA[3][2] = -(
                            +A[0][0] * (A[1][1] * A[3][2] - A[1][2] * A[3][1])
                            +A[0][1] * (A[1][2] * A[3][0] - A[1][0] * A[3][2])
                            +A[0][2] * (A[1][0] * A[3][1] - A[1][1] * A[3][0])
                        ) / det;
            InvA[0][3] = -(
                            +A[0][1] * (A[1][2] * A[2][3] - A[1][3] * A[2][2])
                            +A[0][2] * (A[1][3] * A[2][1] - A[1][1] * A[2][3])
                            +A[0][3] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                        ) / det;
            InvA[1][3] = +(
                            +A[0][0] * (A[1][2] * A[2][3] - A[1][3] * A[2][2])
                            +A[0][2] * (A[1][3] * A[2][0] - A[1][0] * A[2][3])
                            +A[0][3] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                        ) / det;
            InvA[2][3] = -(
                            +A[0][0] * (A[1][1] * A[2][3] - A[1][3] * A[2][1])
                            +A[0][1] * (A[1][3] * A[2][0] - A[1][0] * A[2][3])
                            +A[0][3] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
                        ) / det;
            InvA[3][3] = +(
                            +A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                            +A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2])
                            +A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])
                        ) / det;
            break;
        default:
            break;
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void MC_Matrix_Mul_Vector(int nRow, int nCol, double **Mat, double *B, double *X) {
    int iRow, iCol;

    /* Basic Checking before proceeding */
    if ((Mat == NULL) || (B == NULL) || (X == NULL) ||
            (nRow <= 0) || (nCol <= 0) || (nRow != nCol))
        return;
    
    for (iRow = 0; iRow < nRow; iRow++) {
        X[iRow] = 0.0;
        for (iCol = 0; iCol < nCol; iCol++)
            X[iRow] += Mat[iRow][iCol]*B[iCol];
    }
}

//------------------------------------------------------------------------------
//! Currently Supports Square Matrix
//------------------------------------------------------------------------------
void MC_Matrix_Mul_Matrix(int nRow, int nCol, double **Mat1, double **Mat2, double **Out) {
    int iRow, iCol, k;

    /* Basic Checking before proceeding */
    if ((Mat1 == NULL) || (Mat2 == NULL) || (Out == NULL) ||
            (nRow <= 0) || (nCol <= 0) || (nRow != nCol))
        return;

    for (iRow = 0; iRow < nRow; iRow++) {
        for (iCol = 0; iCol < nCol; iCol++) {
            Out[iRow][iCol] = 0.0;
            for (k = 0; k < nCol; k++)
                Out[iRow][iCol] += Mat1[iRow][k] * Mat2[k][iCol];
        }
    }
}

//------------------------------------------------------------------------------
//! Currently Supports Upto 4x4 Block Matrix
//------------------------------------------------------------------------------
double MC_Iterative_Block_Jacobi_CRS (int Iteration, double Relax, MC_CRS Object) {

    int i, col, iRow, iCol, iRMS, iLoop, iStart, iEnd, MaxIter = 10000;
    double RMS = 0.0;
    double *Block_X = NULL;
    double **RHS = NULL;
    double ***InvBlock = NULL;
    
    /* Basic Checking before proceeding */
    if ((Object.IA == NULL) || (Object.JA == NULL) || (Object.IAU == NULL) ||
            (Object.A == NULL) || (Object.B == NULL) || (Object.X == NULL) ||
            (Object.Block_nRow <= 0) || (Object.Block_nCol <= 0) ||
            (Object.DIM <= 0) || (Object.nROW <= 0) || (Object.nCOL <= 0) ||
            (Object.nCOL != Object.nROW) || (Object.Block_nCol != Object.Block_nRow))
        return RMS;

    // Set the Initial Value for Solution Vector
    for (iRow = 0; iRow < Object.nROW; iRow++) {
        RMS = 0.0;
        for (iCol = 0; iCol < Object.Block_nCol; iCol++)
            RMS += Object.B[iRow][iCol];
        RMS /= Object.Block_nCol;
        for (iCol = 0; iCol < Object.Block_nCol; iCol++)
            Object.X[iRow][iCol] = RMS;
    }

    // Allocate Memory to Store Inverse of Each Diagonal Block Matrix
    InvBlock = (double ***) malloc(Object.nCOL*sizeof(double**));
    for (iCol = 0; iCol < Object.nCOL; iCol++) {
        InvBlock[iCol] = (double **) malloc(Object.Block_nRow*sizeof(double*));
        for (iRow = 0; iRow < Object.Block_nRow; iRow++) {
            InvBlock[iCol][iRow] = (double *) malloc(Object.Block_nCol*sizeof(double));
            // Initialize
            for (i = 0; i <Object.Block_nCol; i++)
                InvBlock[iCol][iRow][i] = 0.0;
        }
    }

    // Compute the Inverse of Each Diagonal Block Matrix
    for (iCol = 0; iCol < Object.nCOL; iCol++)
        MC_Inverse(Object.Block_nRow, Object.Block_nCol, Object.A[Object.IAU[iCol]], InvBlock[iCol]);
    
    // Allocate Memory to store a Copy of RHS
    RHS = (double **) malloc(Object.nCOL * sizeof (double*));
    for (iCol = 0; iCol < Object.nCOL; iCol++) {
        RHS[iCol] = (double *) malloc(Object.Block_nCol * sizeof (double));
    }

    // Allocate Memory to store Block Vector
    Block_X = (double *) malloc (Object.Block_nCol*sizeof(double));
    
    // Iterate Until Convergence or Max Iteration
    if (Iteration <= 0)
        Iteration = MaxIter;
    
    for (iLoop = 0; iLoop < Iteration; iLoop++) {
        RMS  = 0.0;
        iRMS = 0;
        // Create a Copy of RHS
        for (iCol = 0; iCol < Object.nCOL; iCol++) {
            for (iRow = 0; iRow < Object.Block_nCol; iRow++)
                RHS[iCol][iRow] = Object.B[iCol][iRow];
        }
        
        for (iRow = 0; iRow < Object.nROW; iRow++) {
            iStart = Object.IA[iRow];
            iEnd   = Object.IA[iRow+1];
            // Transfer Off-diagonal elements to RHS
            for (iCol = iStart; iCol < iEnd; iCol++) {
                col = Object.JA[iCol];
                if (col == iRow)
                    continue;
                // Matrix-Vector Multiply of Off Diagonal Elements
                MC_Matrix_Mul_Vector(Object.Block_nRow, Object.Block_nCol, Object.A[iCol], Object.X[col], Block_X);
                // Subtract Value form RHS current Row
                for (i = 0; i < Object.Block_nCol; i++)
                    RHS[iRow][i] -= Block_X[i];
            }

            // Multiply RHS with Inverse of Diagonal Matrix
            MC_Matrix_Mul_Vector(Object.Block_nRow, Object.Block_nCol, InvBlock[iRow], RHS[iRow], Block_X);
            // Update the RHS which is Solution for Loop +1
            for (i = 0; i < Object.Block_nCol; i++)
                RHS[iRow][i] = Block_X[i];

            // Compute the RMS
            for (i = 0; i < Object.Block_nCol; i++)
                RMS += (RHS[iRow][i] - Object.X[iRow][i])*(RHS[iRow][i] - Object.X[iRow][i]);
            iRMS++;
        }
        RMS /= iRMS;
        RMS = sqrt(RMS);

        if (RMS < DBL_EPSILON)
            iLoop = MaxIter;

        // Update the Value of Solution
        for (iRow = 0; iRow < Object.nROW; iRow++) {
            for (iCol = 0; iCol < Object.Block_nCol; iCol++)
                Object.X[iRow][iCol] += Relax*(RHS[iRow][iCol] - Object.X[iRow][iCol]);
        }
    }

    // Free Memory
    if (Block_X != NULL)
        free(Block_X);

    if (RHS != NULL) {
        for (iRow = 0; iRow < Object.nROW; iRow++)
            if (RHS[iRow] != NULL)
                free(RHS[iRow]);
        free(RHS);
    }
    
    if (InvBlock != NULL) {
        for (iRow = 0; iRow < Object.nROW; iRow++) {
            if (InvBlock[iRow] != NULL) {
                for (iCol = 0; iCol < Object.Block_nCol; iCol++) {
                    if (InvBlock[iRow][iCol] != NULL)
                        free(InvBlock[iRow][iCol]);
                }
                free(InvBlock[iRow]);
            }
        }
        free(InvBlock);
    }

    return RMS;
}

// *****************************************************************************
// Gauss-Seidel
// Direction = (0 for forwards, 1 for backwards, 2 for alternating)
// *****************************************************************************
double MC_Iterative_Block_LU_Jacobi_CRS(int Iteration, int Direction, MC_CRS Object) {
    int i, ii, j, jj, k, m, iRow, iLoop, iStart, iEnd, MaxIter = 10000;
    double *tempvect = NULL;
    double **RHS     = NULL;
    double ***LU     = NULL;
    double RMS       = 0.0;
    
    /* Basic Checking before proceeding */
    if ((Object.IA == NULL) || (Object.JA == NULL) || (Object.IAU == NULL) ||
            (Object.A == NULL) || (Object.B == NULL) || (Object.X == NULL) ||
            (Object.Block_nRow <= 0) || (Object.Block_nCol <= 0) ||
            (Object.DIM <= 0) || (Object.nROW <= 0) || (Object.nCOL <= 0) ||
            (Object.nCOL != Object.nROW) || (Object.Block_nCol != Object.Block_nRow))
        return RMS;
    
    // Allocate Memory
    tempvect = (double *)   malloc(Object.Block_nRow*sizeof(double));
    RHS      = (double **)  malloc(Object.nROW*sizeof(double*));
    LU       = (double ***) malloc(Object.nROW*sizeof(double**));
    for (iRow = 0; iRow < Object.nROW; iRow++) {
        RHS[iRow] = (double *)  malloc(Object.Block_nRow*sizeof(double));
        LU[iRow]  = (double **) malloc(Object.Block_nRow*sizeof(double*));
        for (j = 0; j < Object.Block_nRow; j++)
            LU[iRow][j]  = (double *) malloc(Object.Block_nCol*sizeof(double));
    }
    
    // Initialize the Solution Vector
    for (iRow = 0; iRow < Object.nROW; iRow++) {
        for (j = 0; j < Object.Block_nRow; j++)
            Object.X[iRow][j] = 0.0;
    }
    
    // LU factorizations with no pivot of diagonal elements
    for (iRow = 0; iRow < Object.nROW; iRow++) {
        // Copy A in LU
        for (i = 0; i < Object.Block_nRow; i++) {
            for (j = 0; j < Object.Block_nCol; j++)
                LU[iRow][i][j] = Object.A[Object.IAU[iRow]][i][j];
        }

        for (k = 0; k < Object.Block_nRow - 1; k++) {
            for (j = k + 1; j < Object.Block_nRow; j++)
                LU[iRow][j][k] = LU[iRow][j][k] / LU[iRow][k][k];
            
            for (j = k + 1; j < Object.Block_nCol; j++) {
                for (i = k + 1; i < Object.Block_nRow; i++)
                    LU[iRow][i][j] -= LU[iRow][i][k] * LU[iRow][k][j];
            }
        }
    }

    // Iterate Until Convergence or Max Iteration
    if (Iteration <= 0)
        Iteration = MaxIter;

    // Progressively iterate to approximate solutions
    for (iLoop = 0; iLoop < Iteration; iLoop++) {
        RMS = 0.0;

        // initializing right hand side
        for (iRow = 0; iRow < Object.nROW; iRow++) {
            for (j = 0; j < Object.Block_nRow; j++)
                RHS[iRow][j] = Object.B[iRow][j];
        }

        if ((Direction == 0) || (Direction == 2 && iLoop % 2 == 0)) {
            for (iRow = 0; iRow < Object.nROW; iRow++) {
                iStart = Object.IA[iRow];
                iEnd   = Object.IA[iRow + 1];
                for (j = iStart; j < iEnd; j++) {  
                    // Matrix Vector Multiply
                    for (ii = 0; ii < Object.Block_nRow; ii++)
                        tempvect[ii] = 0.0;
                    
                    for (ii = 0; ii < Object.Block_nRow; ii++) {
                        for (jj = 0; jj < Object.Block_nCol; jj++)
                            tempvect[ii] += Object.A[j][ii][jj] * Object.X[Object.JA[j]][jj];
                    }

                    for (m = 0; m < Object.Block_nRow; m++)
                        RHS[iRow][m] -= tempvect[m];
                }

                for (m = 0; m < Object.Block_nRow; m++)
                    RMS += RHS[iRow][m] * RHS[iRow][m];
                
                // Matrix Vector Multiply
                for (ii = 0; ii < Object.Block_nRow; ii++)
                    tempvect[ii] = 0.0;
                for (ii = 0; ii < Object.Block_nRow; ii++) {
                    for (jj = 0; jj < Object.Block_nCol; jj++)
                        tempvect[ii] += Object.A[Object.IAU[iRow]][ii][jj] * Object.X[iRow][jj];
                }

                for (m = 0; m < Object.Block_nRow; m++)
                    RHS[iRow][m] += tempvect[m];
                
                // Solve LU
                for (ii = 0; ii < Object.Block_nRow; ii++) {
                    tempvect[ii] = RHS[iRow][ii];
                    for (jj = 0; jj < ii; jj++)
                        tempvect[ii] -= LU[iRow][ii][jj] * tempvect[jj];
                }
                for (ii = Object.Block_nRow - 1; ii >= 0; ii--) {
                    Object.X[iRow][ii] = tempvect[ii];
                    for (jj = Object.Block_nCol - 1; jj > ii; jj--)
                        Object.X[iRow][ii] -= LU[iRow][ii][jj] * Object.X[iRow][jj];

                    Object.X[iRow][ii] = Object.X[iRow][ii] / LU[iRow][ii][ii];
                }
            }
        } else {
            for (iRow = Object.nROW - 1; iRow >= 0; iRow--) {
                iStart = Object.IA[iRow];
                iEnd   = Object.IA[iRow + 1];
                for (j = iStart; j < iEnd; j++) {
                    // Matrix Vector Multiply
                    for (ii = 0; ii < Object.Block_nRow; ii++)
                        tempvect[ii] = 0.0;
                    for (ii = 0; ii < Object.Block_nRow; ii++) {
                        for (jj = 0; jj < Object.Block_nCol; jj++)
                            tempvect[ii] += Object.A[j][ii][jj] * Object.X[Object.JA[j]][jj];
                    }
                    
                    for (m = 0; m < Object.Block_nRow; m++)
                        RHS[iRow][m] -= tempvect[m];
                }

                for (m = 0; m < Object.Block_nRow; m++)
                    RMS += RHS[iRow][m] * RHS[iRow][m];
                
                // Matrix Vector Multiply
                for (ii = 0; ii < Object.Block_nRow; ii++)
                    tempvect[ii] = 0.0;
                for (ii = 0; ii < Object.Block_nRow; ii++) {
                    for (jj = 0; jj < Object.Block_nCol; jj++)
                        tempvect[ii] += Object.A[Object.IAU[iRow]][ii][jj] * Object.X[iRow][jj];
                }
                
                for (m = 0; m < Object.Block_nRow; m++)
                    RHS[iRow][m] += tempvect[m];
                
                // Solve LU
                for (ii = 0; ii < Object.Block_nRow; ii++) {
                    tempvect[ii] = RHS[iRow][ii];
                    for (jj = 0; jj < ii; jj++)
                        tempvect[ii] -= LU[iRow][ii][jj] * tempvect[jj];
                }
                for (ii = Object.Block_nRow - 1; ii >= 0; ii--) {
                    Object.X[iRow][ii] = tempvect[ii];
                    for (jj = Object.Block_nCol - 1; jj > ii; jj--)
                        Object.X[iRow][ii] -= LU[iRow][ii][jj] * Object.X[iRow][jj];
                    
                    Object.X[iRow][ii] = Object.X[iRow][ii] / LU[iRow][ii][ii];
                }
            }
        }

        RMS = sqrt(RMS / ((double) Object.nROW));
        
        if (RMS < 100.0*DBL_EPSILON)
            iLoop = MaxIter;
    }

    // Free Memory
    if (tempvect != NULL)
        free(tempvect);
    
    if (RHS != NULL) {
        for (iRow = 0; iRow < Object.nROW; iRow++)
            if (RHS[iRow] != NULL)
                free(RHS[iRow]);
        free(RHS);
    }

    if (LU != NULL) {
        for (iRow = 0; iRow < Object.nROW; iRow++) {
            if (LU[iRow] != NULL) {
                for (jj = 0; jj < Object.Block_nRow; jj++) {
                    if (LU[iRow][jj] != NULL)
                        free(LU[iRow][jj]);
                }
                free(LU[iRow]);
            }
        }
        free(LU);
    }
    
    return RMS;
}

