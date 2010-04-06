/* 
 * File:   TriangleMeshOptimizer.cpp
 * Author: Ashish Gupta
 * 
 * Created on February 5, 2010, 7:22 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctime>
#include <iostream>

#include "Utils.h"
#include "MUtils.h"
#include "TriangleMeshOptimizer.h"

// *****************************************************************************
// *****************************************************************************
TriangleMeshOptimizer::TriangleMeshOptimizer() {

    printf("\n====================================================");
    printf("\n    Mesh Module : Optimization                      ");
    printf("\n====================================================\n");

    // Initialize the Data
    Init();
    
    return;
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Init() {
    RState  = 0;
    RAngle  = 0.0;
    CGx     = 0.0;
    CGy     = 0.0;
}

// *****************************************************************************
// *****************************************************************************
TriangleMeshOptimizer::~TriangleMeshOptimizer() {
    
}

// *****************************************************************************
// *****************************************************************************
int TriangleMeshOptimizer::Opposite_Edge_Normal_Smoothing(int iNode, double* cost, double *U, double *V) {
    int cellid, improved;
    int n0, n1, n2;
    double nx, ny, mag, newcost;
    double xorg, yorg, dx, dy;

    improved = 0;
    for (int iCell = 0; iCell < Node2Cell[iNode]->max; iCell++) {
        cellid = Node2Cell[iNode]->list[iCell];
        
        // Get the Opposite Nodes Edge
        for (int i = 0; i < 3; i++) {
            n0 = Tri[cellid][i];
            n1 = Tri[cellid][(i+1)%3];
            n2 = Tri[cellid][(i+2)%3];
            if (n0 == iNode)
                break;
        }

        // Get the unit Normal to the Opposte Edge
        nx = Y[n1] - Y[n2];
        ny = X[n2] - X[n1];
        mag = sqrt(nx*nx + ny*ny);
        nx /= mag;
        ny /= mag;
        
        // Get the purturbation magnitude
        mag = sqrt(U[iNode]*U[iNode] + V[iNode]*V[iNode]);
        dx = nx*mag;
        dy = ny*mag;

        // Save the Original Position
        xorg = X[iNode];
        yorg = Y[iNode];
        // Perturbe the Node Position
        X[iNode] += dx;
        Y[iNode] += dy;

        // Compute the New Cost
        newcost = Compute_Node_Cost(iNode);

        if (newcost < *cost) {
            U[iNode] *= 1.1;
            V[iNode] *= 1.1;
            *cost = newcost;
            improved = 1;
        } else {
            X[iNode] = xorg;
            Y[iNode] = yorg;
        }
    }

    return (improved);
}

// *****************************************************************************
// *****************************************************************************
int TriangleMeshOptimizer::Neigbour_Edge_Smoothing(int iNode, double *cost, double *U, double *V) {
    int nodeid, improved;
    double ex, ey, mag, newcost;
    double xorg, yorg, dx, dy;
    
    improved = 0;
    for (int i = 0; i < Node2Node[iNode]->max; i++) {
        nodeid = Node2Node[iNode]->list[i];
        ex = X[nodeid] - X[iNode];
        ey = Y[nodeid] - Y[iNode];
        mag = sqrt(ex*ex + ey*ey);
        ex /= mag;
        ey /= mag;

        // Get the purturbation magnitude
        mag = sqrt(U[iNode]*U[iNode] + V[iNode]*V[iNode]);
        dx = ex*mag;
        dy = ey*mag;
        
        // Save the Original Position
        xorg = X[iNode];
        yorg = Y[iNode];
        // Perturbe the Node Position
        X[iNode] += dx;
        Y[iNode] += dy;

        // Compute the New Cost
        newcost = Compute_Node_Cost(iNode);

        if (newcost < *cost) {
            U[iNode] *= 1.1;
            V[iNode] *= 1.1;
            *cost = newcost;
            improved = 1;
        } else {
            X[iNode] = xorg;
            Y[iNode] = yorg;
        }
    }
    
    return (improved);
}

// *****************************************************************************
// *****************************************************************************
int TriangleMeshOptimizer::X_Y_Smoothing(int iNode, double* cost, double* U, double* V) {
    int improved;
    double xorg, yorg, newcost;

    improved = 0;

    // Save the Original Position
    xorg = X[iNode];
    yorg = Y[iNode];

    // Perturbe the Node Position in U direction
    X[iNode] += U[iNode];
    
    // Compute the New Cost
    newcost = Compute_Node_Cost(iNode);

    if (newcost < *cost) {
        U[iNode] *= 1.1;
        *cost = newcost;
        improved = 1;
    } else {
        X[iNode] = xorg;
        U[iNode] *= -0.5;
    }

    // Perturbe the Node Position in V direction
    Y[iNode] += V[iNode];

    // Compute the New Cost
    newcost = Compute_Node_Cost(iNode);

    if (newcost < *cost) {
        V[iNode] *= 1.1;
        *cost = newcost;
        improved = 1;
    } else {
        Y[iNode] = yorg;
        V[iNode] *= -0.5;
    }
    
    return (improved);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Optimize() {
    int i, iNode, nodeid, count, improved;
    int invertNode;
    double dx, dy, len;
    double cost, maxcost, mincost;
    double eps = DBL_MAX;
    double *U = NULL;
    double *V = NULL;

    // Create Interior and Boundary Node Tagging
    Create_Interior_Boundary_Tag();

    // Rotate Required Boundaries
    Rotate_Boundary(1);
    
    // Create All Possible Connectivities
    Create_Connectivity(0);

    int MaxIter = 500;
    std::cout << "Input Optimization Iterations : ";
    std::cin  >> MaxIter;
    
    U = new double[NNode];
    V = new double[NNode];
    
    for (iNode = 0; iNode < NNode; iNode++) {
        U[iNode] = 0.0;
        V[iNode] = 0.0;
        len = 0.0;
        count = Node2Node[iNode]->max;
        for (i = 0; i < count; i++) {
           nodeid = Node2Node[i]->list[i];
           dx = fabs(X[nodeid] - X[iNode]);
           dy = fabs(Y[nodeid] - Y[iNode]);
           len += sqrt(dx*dx + dy*dy);
           U[iNode] += dx;
           V[iNode] += dy;
        }
        count = MAX(1, count);
        U[iNode] /= (double) count;
        V[iNode] /= (double) count;
        len      /= (double) count;
        eps = MIN(eps, len);
    }
    
    eps = MAX(DBL_EPSILON, eps*1.0e-08);
    info("Computed Epsilion: %lg", eps);
    
    // Optimization Loop Starts
    info("Mesh Optimization Starts");
    printf("---------------------------------------------------\n");
    info("Iteration   InvNode      MaxCost      MinCost");
    
    for (int iter = 0; iter < MaxIter; iter++) {
        invertNode = 0;
        maxcost = DBL_MIN;
        mincost = DBL_MAX;
        for (iNode = 0; iNode < NNode; iNode++) {
            // Continue if Boundary Node
            if (IBTag[iNode] > -1)
                continue;

            // Ensure U & V are Non-Zero
            U[iNode] = SIGN(MAX(fabs(U[iNode]), eps), U[iNode]);
            V[iNode] = SIGN(MAX(fabs(V[iNode]), eps), V[iNode]);

            // Compute Current Cost
            cost = Compute_Node_Cost(iNode);

            // Check if Inverted Element
            if (cost > 1.0)
                invertNode++;
            
            // Start Opposite Edge Normal Smoothing
            improved = Opposite_Edge_Normal_Smoothing(iNode, &cost, U, V);
            // Start Neighbour Edge Smoothing
            improved = Neigbour_Edge_Smoothing(iNode, &cost, U, V);
            // X-Y Smoothing
            if (improved != 1)
                improved = X_Y_Smoothing(iNode, &cost, U, V);
            maxcost = MAX(maxcost, cost);
            mincost = MIN(mincost, cost);
        }
        if (!(iter%10))
            info("%9d %9d %12.5e %12.5e", iter, invertNode, maxcost, mincost);
    }
    printf("---------------------------------------------------\n");
    delete []U;
    delete []V;
}

// *****************************************************************************
// *****************************************************************************
double TriangleMeshOptimizer::Compute_Node_Cost(int iNode) {
    int iCell, cellid;
    int n0, n1, n2;
    double J, cost, avgcost, maxcost;
    int count = Node2Cell[iNode]->max;

    maxcost = 0.0;
    avgcost = 0.0;
    for (iCell = 0; iCell < count; iCell++) {
        cellid = Node2Cell[iNode]->list[iCell];
        n0 = Tri[cellid][0];
        n1 = Tri[cellid][1];
        n2 = Tri[cellid][2];

        // Compute the Jacobian of Triangle
        J = (X[n1] - X[n0])*(Y[n2]-Y[n0]) - (X[n2] - X[n0])*(Y[n1]-Y[n0]);

        if (J < 0.0)
            cost = 1.0 - J;
        else
            cost = 1.0 - (1.0/Compute_Tri_ConditionNumber(n0, n1, n2));
        
        cost *= cost;
        avgcost += cost;
        maxcost = MAX(cost, maxcost);
    }

    count = MAX(1, count);
    avgcost /= count;
    
    if (maxcost > 1.0) {
        cost = maxcost;
    } else {
        double ratio = maxcost*maxcost*maxcost;
        cost = ratio*maxcost + MAX(0.0, (1.0 - ratio))*avgcost;
    }
    
    return (cost);
}

// *****************************************************************************
// *****************************************************************************
double TriangleMeshOptimizer::Compute_Tri_ConditionNumber(int n0, int n1, int n2) {
    double W[4], InW[4], A[4], InA[4], tmp[4];
    double ux, uy, vx, vy, cn;

    ux = X[n1] - X[n0];
    uy = Y[n1] - Y[n0];
    vx = X[n2] - X[n0];
    vy = Y[n2] - Y[n0];

    W[0] = 1.0;
    W[1] = 0.5;
    W[2] = 0.0;
    W[3] = 0.5 * sqrt(3.0);
    InvMat2x2(W, InW);

    A[0] = ux;
    A[1] = vx;
    A[2] = uy;
    A[3] = vy;
    InvMat2x2(A, InA);
    MulMat2x2(A, InW, tmp);
    cn = sqrt((tmp[0] * tmp[0])+(tmp[1] * tmp[1])+(tmp[2] * tmp[2])+(tmp[3] * tmp[3]));
    MulMat2x2(W, InA, tmp);
    cn = cn * sqrt((tmp[0] * tmp[0])+(tmp[1] * tmp[1])+(tmp[2] * tmp[2])+(tmp[3] * tmp[3])) / 2.0;

    return (cn);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Rotate_Boundary(int Mode) {
    int bndid;
    double rx, ry, rmag, theta;

    // Input Mode
    if (Mode == 0) {
        printf("---------------------------------------------------\n");
        std::cout << "Rotate Boundary (Yes = 1; No = 0): ";
        std::cin >> RState;
        if (RState == 1) {
            std::cout << "Input ID Boundary (0 ... " << (NBoundary-1) << ")" << std::endl;
            bndid = -2;
            do {
                std::cout << "Boundary To Rotate (-1 = Done) : ";
                std::cin >> bndid;
                if (bndid != -1) {
                    RBoundary.Check_List(bndid);
                }
            } while (bndid != -1);
            std::cout << "Input Angle of Rotation (0 - 90): ";
            std::cin  >> RAngle;
            RAngle *= -1.0;
            std::cout << "Input Rotation CG-X: ";
            std::cin  >> CGx;
            std::cout << "Input Rotation CG-Y: ";
            std::cin  >> CGy;
        } else
            RState = 0;
        printf("---------------------------------------------------\n");
    } else {
        if (!RState)
            return;
        printf("---------------------------------------------------\n");
        info("START: Boundary Curve Rotation");
        RAngle = RAngle * M_PI / 180.0;

        for (int ib = 0; ib < RBoundary.max; ib++) {
            bndid = RBoundary.list[ib];
            info("ROTATE: Boundary: %d", bndid);
            for (int inode = 0; inode < NNode; inode++) {
                if (IBTag[inode] != bndid)
                    continue;
                rx = X[inode] - CGx;
                ry = Y[inode] - CGy;
                rmag = sqrt(rx * rx + ry * ry);
                theta = atan2(ry, rx);
                X[inode] = CGx + rmag * cos(theta + RAngle);
                Y[inode] = CGy + rmag * sin(theta + RAngle);
            }
        }
        info("END: Boundary Cure Rotation");
        printf("---------------------------------------------------\n");
    }
    
}

