/*******************************************************************************
 * File:        TriangleWinslowSmoother.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _TRIANGLEWINSLOWSMOOTHER_H
#define	_TRIANGLEWINSLOWSMOOTHER_H

#include "TriangleMeshDB.h"

class TriangleWinslowSmoother: virtual public TriangleMeshDB {
protected:
    // Data Structure for Boundary Rotation
    int    RState;
    double RAngle;
    double CGx;
    double CGy;
    List RBoundary;
    // Compressed Row Storage
    int CRS_DIM;
    int *CRS_IA;
    int *CRS_IAU;
    int *CRS_JA;
    double ***CRS_MATRIX;
public:
    TriangleWinslowSmoother();
    virtual ~TriangleWinslowSmoother();
    void WinslowSmooth();
    void Rotate_Boundary(int Mode);
private:
    void Create_CRS();
    void Compute_Gauss_Gradient(double *F, double *Fx, double *Fy);
    void Compute_CRS_Matrix(double *F, double *Fx, double *Fy, double *G, double *Gx, double *Gy);
    double Solve_CRS_Linear_Equation(int Iteration, double Relax, double *F, double *G);
private:
    void Init();
};


#endif	/* _TRIANGLEWINSLOWSMOOTHER_H */

