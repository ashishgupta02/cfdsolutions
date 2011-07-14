/*******************************************************************************
 * File:        MC.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _MC_H
#define	_MC_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct _MC_CRS {
        int Block_nRow;
        int Block_nCol;
        int DIM;
        int nROW;
        int nCOL;
        int *IA;
        int *JA;
        int *IAU;
        double ***A;
        double **B;
        double **X;
    } MC_CRS;
    
    double MC_Determinant(int nRow, int nCol, double **A);
    void   MC_Inverse(int nRow, int nCol, double **A, double **InvA);
    void   MC_Matrix_Mul_Vector(int nRow, int nCol, double **Mat, double *B, double *X);
    void   MC_Matrix_Mul_Matrix(int nRow, int nCol, double **Mat1, double **Mat2, double **Out);
    double MC_Iterative_Block_Jacobi_CRS (int Iteration, double Relax, MC_CRS Object);
    double MC_Iterative_Block_LU_Jacobi_CRS(int Iteration, int Direction, MC_CRS Object);

#ifdef	__cplusplus
}
#endif

#endif	/* _MC_H */

