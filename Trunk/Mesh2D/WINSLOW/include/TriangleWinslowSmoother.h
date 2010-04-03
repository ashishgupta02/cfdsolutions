/* 
 * File:   TriangleWinslowSmoother.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#ifndef _TRIANGLEWINSLOWSMOOTHER_H
#define	_TRIANGLEWINSLOWSMOOTHER_H

#include "List.h"

class TriangleWinslowSmoother {
public:
    TriangleWinslowSmoother();
    virtual ~TriangleWinslowSmoother();
    void WinslowSmooth();
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
    void Compute_Gauss_Gradient(double *F, double *Fx, double *Fy);
    void Compute_CRS_Matrix(double *F, double *Fx, double *Fy, double *G, double *Gx, double *Gy);
    double Solve_CRS_Linear_Equation(int Iteration, double Relax, double *F, double *G);
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


#endif	/* _TRIANGLEWINSLOWSMOOTHER_H */

