/* 
 * File:   TriangleLinearElasticSmoother.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#ifndef _TRIANGLELINEARELASTICSMOOTHER_H
#define	_TRIANGLELINEARELASTICSMOOTHER_H

#include "List.h"

class TriangleLinearElasticSmoother {
public:
    TriangleLinearElasticSmoother();
    virtual ~TriangleLinearElasticSmoother();
    void LESmooth();
    void SimCenterMeshReader(const char* FileName);
    void SimCenterMeshWriter(const char* FileName);
    void GnuplotWriter(const char* FileName);
    void Rotate_Boundary(int Mode);
private:
    void Init();
    void Create_Interior_Boundary_Tag();
    void Create_Connectivity();
    void Create_Node2Cell_Connectivity();
    void Create_Cell2Cell_Connectivity();
    void Create_Node2Node_Connectivity();
    void Create_CRS();
    void Compute_CRS_Matrix();
    double Solve_CRS_Linear_Equation(int Iteration, double Relax, double *F, double *G);
    double Get_Poission_Ratio(int TriID);
    double Get_Youngs_Modulus(int TriID);
    double Compute_Tri_ConditionNumber(int n0, int n1, int n2);
    double Compute_Tri_AspectRatio(int n0, int n1, int n2);
    double Compute_Tri_Area(int n0, int n1, int n2);
private:
    int State;
    int NNode;
    int NBlock;
    int NTri;
    int NQuad;
    int NBoundary;
    int *NBoundarySegments;
    int ***BoundarySegments;
    // Coordinates
    double *X;
    double *Y;
    // Cell2Node Connectivity
    int (*Tri)[3];
    int (*Quad)[4];
    // Connectivity Data Structure
    int  **Cell2Cell;
    List **Node2Cell;
    List **Node2Node;
    // Tag Interior and Boundary Nodes
    int *IBTag;
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
};


#endif	/* _TRIANGLELINEARELASTICSMOOTHER_H */

