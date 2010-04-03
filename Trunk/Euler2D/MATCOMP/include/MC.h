/* 
 * File:   MC.h
 * Author: Ashish Gupta
 *
 * Created on February 11, 2010, 2:20 PM
 */

#ifndef _MC_H
#define	_MC_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct {
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
    void MC_Inverse(int nRow, int nCol, double **A, double **InvA);
    void MC_Matrix_Mul_Vector(int nRow, int nCol, double **Mat, double *B, double *X);
    double MC_Iterative_Block_Jacobi_CRS (int Iteration, double Relax, MC_CRS Object);

#ifdef	__cplusplus
}
#endif

#endif	/* _MC_H */

