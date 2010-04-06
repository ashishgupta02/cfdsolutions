/* 
 * File:   TriangleLinearElasticSmoother.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 * Modified on April 6, 2010
 */

#ifndef _TRIANGLELINEARELASTICSMOOTHER_H
#define	_TRIANGLELINEARELASTICSMOOTHER_H

#include "TriangleMeshDB.h"

class TriangleLinearElasticSmoother: virtual public TriangleMeshDB {
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
    TriangleLinearElasticSmoother();
    virtual ~TriangleLinearElasticSmoother();
    void LESmooth();
    void Rotate_Boundary(int Mode);
protected:
    void Create_CRS();
    void Compute_CRS_Matrix();
    double Solve_CRS_Linear_Equation(int Iteration, double Relax, double *F, double *G);
    double Get_Poission_Ratio(int TriID);
    double Get_Youngs_Modulus(int TriID);
    double Compute_Tri_ConditionNumber(int n0, int n1, int n2);
    double Compute_Tri_AspectRatio(int n0, int n1, int n2);
    double Compute_Tri_Area(int n0, int n1, int n2);
private:
    void Init();
};


#endif	/* _TRIANGLELINEARELASTICSMOOTHER_H */

