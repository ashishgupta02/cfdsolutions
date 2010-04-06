/* 
 * File:   TriangleMeshOptimizer.h
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 * Modified on April 6, 2010
 */

#ifndef _TRIANGLEMESHOPTIMIZER_H
#define	_TRIANGLEMESHOPTIMIZER_H

#include "TriangleMeshDB.h"

class TriangleMeshOptimizer: virtual public TriangleMeshDB {
protected:
    // Data Structure for Boundary Rotation
    int    RState;
    double RAngle;
    double CGx;
    double CGy;
    List   RBoundary;
public:
    TriangleMeshOptimizer();
    virtual ~TriangleMeshOptimizer();
    void Optimize();
    void Rotate_Boundary(int Mode);
protected:
    double Compute_Node_Cost(int iNode);
    double Compute_Tri_ConditionNumber(int n0, int n1, int n2);
    int Opposite_Edge_Normal_Smoothing(int iNode, double *cost, double *U, double *V);
    int Neigbour_Edge_Smoothing(int iNode, double *cost, double *U, double *V);
    int X_Y_Smoothing(int iNode, double *cost, double *U, double *V);
private:
    void Init();
};

#endif	/* _TRIANGLEMESHOPTIMIZER_H */

