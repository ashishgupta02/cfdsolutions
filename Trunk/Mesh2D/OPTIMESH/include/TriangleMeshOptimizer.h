/* 
 * File:   TriangleMeshOptimizer.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#ifndef _TRIANGLEMESHOPTIMIZER_H
#define	_TRIANGLEMESHOPTIMIZER_H

#include "List.h"

class TriangleMeshOptimizer {
public:
    TriangleMeshOptimizer();
    virtual ~TriangleMeshOptimizer();
    void Optimize();
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
    double Compute_Node_Cost(int iNode);
    double Compute_Tri_ConditionNumber(int n0, int n1, int n2);
    int Opposite_Edge_Normal_Smoothing(int iNode, double *cost, double *U, double *V);
    int Neigbour_Edge_Smoothing(int iNode, double *cost, double *U, double *V);
    int X_Y_Smoothing(int iNode, double *cost, double *U, double *V);
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
};

#endif	/* _TRIANGLEMESHOPTIMIZER_H */

