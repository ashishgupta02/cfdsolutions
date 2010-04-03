/*
 * File:   MC_Iterative_Jacobi.c
 * Author: Ashish Gupta
 *
 * Created on February 10, 2010, 7:09 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "MC.h"

// *****************************************************************************
/* Currently Inverse of 2x2, 3x3 and 4x4 is supported */
// *****************************************************************************
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

// *****************************************************************************
/* Currently Inverse of 2x2, 3x3 and 4x4 is supported */
// *****************************************************************************
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

// *****************************************************************************
// *****************************************************************************
void MC_Matrix_Mul_Vector(int nRow, int nCol, double **Mat, double *B, double *X) {
    int iRow, iCol;

    /* Basic Checking before proceeding */
    if ((Mat == NULL) || (B == NULL) || (X == NULL) ||
            (nRow <= 0) || (nCol <= 0) || (nRow != nCol))
        return;

    for (iRow = 0; iRow < nRow; iRow++)
        X[iRow] = 0.0;
    
    for (iRow = 0; iRow < nRow; iRow++) {
        for (iCol = 0; iCol < nCol; iCol++)
            X[iRow] += Mat[iRow][iCol]*B[iCol];
    }
}

// *****************************************************************************
// *****************************************************************************
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

            // Mutiply RHS with Inverse of Diagonal Matrix
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

