/*******************************************************************************
 * File:        TriangleWinslowSmoother.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctime>
#include <iostream>

#include "Utils.h"
#include "MUtils.h"
#include "TriangleWinslowSmoother.h"

// *****************************************************************************
// *****************************************************************************
TriangleWinslowSmoother::TriangleWinslowSmoother() {

    printf("\n====================================================");
    printf("\n    Mesh Module : Winslow Smoothing                 ");
    printf("\n====================================================\n");

    // Initialize the Data
    Init();

    return;
}

// *****************************************************************************
// *****************************************************************************
void TriangleWinslowSmoother::Init() {
    CRS_IA            = NULL;
    CRS_JA            = NULL;
    CRS_IAU           = NULL;
    CRS_MATRIX        = NULL;
    CRS_DIM           = 0;
    RState            = 0;
    RAngle            = 0.0;
    CGx               = 0.0;
    CGy               = 0.0;
}

// *****************************************************************************
// *****************************************************************************
TriangleWinslowSmoother::~TriangleWinslowSmoother() {
    int i;

    if (CRS_IA != NULL)
        delete [] CRS_IA;

    if (CRS_IAU != NULL)
        delete [] CRS_IAU;

    if (CRS_JA != NULL)
        delete [] CRS_JA;

    if (CRS_MATRIX != NULL) {
        for (i = 0; i < CRS_DIM; i++) {
            if (CRS_MATRIX[i] != NULL) {
                if (CRS_MATRIX[i][0] != NULL)
                    delete [] CRS_MATRIX[i][0];
                if (CRS_MATRIX[i][1] != NULL)
                    delete [] CRS_MATRIX[i][1];
                delete [] CRS_MATRIX[i];
            }
        }
        delete [] CRS_MATRIX;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleWinslowSmoother::WinslowSmooth() {
    int iNode, iter, nwrms;
    double wrms, lrms, relax, dx, dy, ds;
    
    double *U  = NULL;
    double *V  = NULL;
    double *Ux = NULL;
    double *Uy = NULL;
    double *Vx = NULL;
    double *Vy = NULL;
    
    // Allocate Memory for Computations
    U  = new double[NNode];
    V  = new double[NNode];
    Ux = new double[NNode];
    Vx = new double[NNode];
    Uy = new double[NNode];
    Vy = new double[NNode];

#ifdef DEBUG
    if ((U == NULL) || (V == NULL) || (Ux == NULL) || (Uy == NULL) ||
            (Vx == NULL) || (Vy == NULL))
        error("WinslowSmooth: %s\n", "Error Allocating Memory");
#endif
    
    // Get the Computational Mesh and Gradient
    // Note Computational Mesh is X, Y
    for (iNode = 0; iNode < NNode;  iNode++) {
        U[iNode]  = X[iNode];
        V[iNode]  = Y[iNode];
        Ux[iNode] = 0.0;
        Uy[iNode] = 0.0;
        Vx[iNode] = 0.0;
        Vy[iNode] = 0.0;
    }

    // Create Interior and Boundary Node Tagging
    Create_Interior_Boundary_Tag();

    // Rotate Required Boundaries
    Rotate_Boundary(1);

    // Compute Swap U, V and X, Y to Store Original Control Volume
    // This is done to Incoporate Boundary Rotation
    for (iNode = 0; iNode < NNode; iNode++) {
        Ux[iNode] = U[iNode];
        Uy[iNode] = V[iNode];
        U[iNode]  = X[iNode];
        V[iNode]  = Y[iNode];
        X[iNode]  = Ux[iNode];
        Y[iNode]  = Uy[iNode];
        Ux[iNode] = 0.0;
        Uy[iNode] = 0.0;
    }

    // Create All Possible Connectivities
    Create_Connectivity(0);

    // Create Compress Row Storage
    Create_CRS();

    int MaxIter = 1000;
    std::cout << "Input Winslow Smoothing Iterations : ";
    std::cin  >> MaxIter;
    
    // Optimization Loop Starts
    info("Winslow Mesh Smoothing Starts");
    printf("---------------------------------------------------\n");
    info("Iteration  Winslow RMS  LinearSolve RMS");

    // Winslow Smoothing Outer Loop
    for (int iws = 0; iws < MaxIter; iws++) {
        // Compute the Gradient Ux, Uy, Vx, Yy
        Compute_Gauss_Gradient(U, Ux, Uy);
        Compute_Gauss_Gradient(V, Vx, Vy);
#ifdef DEBUG
        for (iNode = 0; iNode < NNode;  iNode++) {
            info("Ux[%d] = %12.5e Uy[%d] = %12.5e Vx[%d] = %12.5e Vy[%d] = %12.5e",
                    iNode, Ux[iNode], iNode, Uy[iNode], iNode, Vx[iNode], iNode, Vy[iNode]);
        }
#endif
        // Build Global Matrix
        Compute_CRS_Matrix(U, Ux, Uy, V, Vx, Vy);

        // Copy U, V for Winslow Residual Computations
        for (iNode = 0; iNode < NNode; iNode++) {
            Ux[iNode] = U[iNode];
            Uy[iNode] = V[iNode];
        }

        // Linear Solver
        iter  = 10;
        relax = 0.5;
        lrms = Solve_CRS_Linear_Equation(iter, relax, U, V);

        // Compute Winslow RMS
        wrms  = 0.0;
        nwrms = 0;
        for (iNode = 0; iNode < NNode; iNode++) {
            dx = U[iNode] - Ux[iNode];
            dy = V[iNode] - Uy[iNode];
            ds = dx*dx + dy*dy;
            wrms +=ds;
            nwrms++;
        }
        wrms /= nwrms;
        wrms = sqrt(wrms);
        
        if (!(iws%10))
            info("%9d %12.5e %12.5e", iws, wrms, lrms);

        if (wrms < (DBL_EPSILON*10.0)) {
            iws = MaxIter;
        }
    }

    // Copy the Final Node Positions
    for (iNode = 0; iNode < NNode; iNode++) {
        X[iNode] = U[iNode];
        Y[iNode] = V[iNode];
    }
    
    printf("---------------------------------------------------\n");
    
    // Free the Used Memory
    delete []U;
    delete []V;
    delete []Ux;
    delete []Vx;
    delete []Uy;
    delete []Vy;
}

// *****************************************************************************
// *****************************************************************************
void TriangleWinslowSmoother::Rotate_Boundary(int Mode) {
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
        } else {
            RState = 0;
        }
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

// *****************************************************************************
// *****************************************************************************
void TriangleWinslowSmoother::Create_CRS() {
    int iNode;

    if (NNode <= 0)
        return;

    // Get the dimension for CRS_JA
    CRS_DIM = 0;
    for (iNode = 0; iNode < NNode; iNode++)
        CRS_DIM += Node2Node[iNode]->max;
    CRS_DIM += NNode;

    // Allocate Memory
    CRS_IA  = new int[NNode+1];
    CRS_IAU = new int[NNode];
    CRS_JA  = new int[CRS_DIM];
#ifdef DEBUG
    if ((CRS_IA == NULL) || (CRS_IAU == NULL) || (CRS_JA == NULL))
        error("Create_CRS: %s\n", "Error Allocating Memory");
#endif

    // Fill Values in the Arrays
    // Get IA
    CRS_IA[0] = 0;
    for (iNode = 1; iNode < (NNode+1); iNode++)
        CRS_IA[iNode] = CRS_IA[iNode-1] + Node2Node[iNode-1]->max+1;
    
    // Get JA and IAU
    int count = 0;
    for (iNode = 0; iNode < NNode ; iNode++) {
        // Location of Diagonal Elements is IAU
        CRS_IAU[iNode] = count;
        CRS_JA[count++] = iNode;
        for (int j = 0; j < Node2Node[iNode]->max; j++)
            CRS_JA[count++] = Node2Node[iNode]->list[j];
    }

    // Allocate Memory for Matrix
    CRS_MATRIX = new double**[CRS_DIM];
#ifdef DEBUG
    if (CRS_MATRIX == NULL)
        error("Create_CRS: %s\n", "Error Allocating Memory 1");
#endif
    for (iNode = 0; iNode < CRS_DIM; iNode++) {
        CRS_MATRIX[iNode] = new double*[2];
#ifdef DEBUG
        if (CRS_MATRIX[iNode] == NULL)
            error("Create_CRS: %s\n", "Error Allocating Memory 2");
#endif
        for (int i = 0; i < 2; i++) {
            CRS_MATRIX[iNode][i] = new double[2];
#ifdef DEBUG
            if (CRS_MATRIX[iNode][i] == NULL)
                error("Create_CRS: %s\n", "Error Allocating Memory 3");
#endif
            CRS_MATRIX[iNode][i][0] = 0.0;
            CRS_MATRIX[iNode][i][1] = 0.0;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleWinslowSmoother::Compute_Gauss_Gradient(double *F, double *Fx, double *Fy) {
    int iNode, iTri, TriID, n0, n1, n2;
    double nx0, ny0, nx1, ny1, nx2, ny2, area;

    if ((F == NULL) || (Fx == NULL) || (Fy == NULL) || (NNode <= 0))
        return;

    for (iNode = 0; iNode < NNode; iNode++) {
        // Ignore Boundary Nodes
        if (IBTag[iNode] > -1)
            continue;

        area  = 0.0;
        // Initialize Gradients
        Fx[iNode] = Fy[iNode] = 0.0;
        // Loop over all Neigbouring Triangles
        for (iTri = 0; iTri < Node2Cell[iNode]->max; iTri++) {
            TriID = Node2Cell[iNode]->list[iTri];

            n0 = Tri[TriID][0];
            n1 = Tri[TriID][1];
            n2 = Tri[TriID][2];

            // Compute the Normals
            nx0 = (Y[n2] - Y[n1]);
            ny0 = -(X[n2] - X[n1]);
            nx1 = (Y[n0] - Y[n2]);
            ny1 = -(X[n0] - X[n2]);
            nx2 = (Y[n1] - Y[n0]);
            ny2 = -(X[n1] - X[n0]);
            
            // Compute Area of the Triangle
            area += fabs(0.5*(nx1*ny2 - nx2*ny1));

            // Compute the Gradient
            Fx[iNode] -= 0.5*(F[n0]*nx0 + F[n1]*nx1 + F[n2]*nx2);
            Fy[iNode] -= 0.5*(F[n0]*ny0 + F[n1]*ny1 + F[n2]*ny2);
        }
        Fx[iNode] /= area;
        Fy[iNode] /= area;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleWinslowSmoother::Compute_CRS_Matrix(double *F, double *Fx, double *Fy, double *G, double *Gx, double *Gy) {
    int iNode, iTri, n0, n1, n2, i, istart, iend;
    double alpha, beta, gamma;
    double nx0, ny0, nx1, ny1, nx2, ny2, tx, ty, area;
    double w00, w01, w02, w10, w11, w12;
    
    if ((F == NULL) || (Fx == NULL) || (Fy == NULL) ||
            (G == NULL) || (Gx == NULL) || (Gy == NULL) || (NNode <= 0))
        return;

    // Initialize
    n0 = -1;
    n1 = -1;
    n2 = -1;
    
    // Initialize the CRS_Matrix
    for (i = 0; i < CRS_DIM; i++) {
        // First Row
        CRS_MATRIX[i][0][0] = 0.0;
        CRS_MATRIX[i][0][1] = 0.0;
        // Second Row
        CRS_MATRIX[i][1][0] = 0.0;
        CRS_MATRIX[i][1][1] = 0.0;
    }

    // For Triangles and computing all node contribution except Boundary Nodes
    for (iTri = 0; iTri < NTri; iTri++) {
        // For all Nodes in the Triangle
        for (iNode = 0; iNode < 3; iNode++) {
            switch (iNode) {
                case 0:
                    n0 = Tri[iTri][0];
                    n1 = Tri[iTri][1];
                    n2 = Tri[iTri][2];
                    break;
                case 1:
                    n0 = Tri[iTri][1];
                    n1 = Tri[iTri][2];
                    n2 = Tri[iTri][0];
                    break;
                case 2:
                    n0 = Tri[iTri][2];
                    n1 = Tri[iTri][0];
                    n2 = Tri[iTri][1];
                    break;
                default:
                    break;
            }
            
            // Check if Node0 is Boundary
            if (IBTag[n0] > -1) {
                // Diagonal Entry and Ingore offdiagonal as they will be zero
                // Since Boundary Nodes are fixed and diagonal will be Identity 'I'
                CRS_MATRIX[CRS_IAU[n0]][0][0] = 1.0;
                CRS_MATRIX[CRS_IAU[n0]][1][1] = 1.0;
            } else {
                // Compute the Normals
                nx0 = (Y[n2] - Y[n1]);
                ny0 = -(X[n2] - X[n1]);
                nx1 = (Y[n0] - Y[n2]);
                ny1 = -(X[n0] - X[n2]);
                nx2 = (Y[n1] - Y[n0]);
                ny2 = -(X[n1] - X[n0]);

                // Compute the area
                area = fabs(0.5 * (nx1 * ny2 - nx2 * ny1));

                // Compute alpha, beta, gamma, tx, ty
                alpha = Fy[n0] * Fy[n0] + Gy[n0] * Gy[n0];
                beta  = Fx[n0] * Fy[n0] + Gx[n0] * Gy[n0];
                gamma = Fx[n0] * Fx[n0] + Gx[n0] * Gx[n0];
                tx = nx0;
                ty = ny0;
                // Calulate W00, W01, W02
                // W00 - Contribution of Node0 itself
                w00 = -0.5 * (alpha * nx0 * tx - 2.0 * beta * ny0 * tx + gamma * ny0 * ty) / area;
                // W01 - Contribution of Node1 to Node 0
                w01 = -0.5 * (alpha * nx1 * tx - 2.0 * beta * ny1 * tx + gamma * ny1 * ty) / area;
                // W02 - Contribution of Node2 to Node 0
                w02 = -0.5 * (alpha * nx2 * tx - 2.0 * beta * ny2 * tx + gamma * ny2 * ty) / area;
                //-- Calulate W10, W11, W12
                // W10 - Contribution of Node0 itself
                w10 = w00;
                // W11 - Contribution of Node1 to Node 0
                w11 = w01;
                // W12 - Contribution of Node2 to Node 0
                w12 = w02;

                // Put values in CRS Matrix
                istart = CRS_IA[n0];
                iend   = CRS_IA[n0 + 1];
                // Diagonal Entry
                CRS_MATRIX[CRS_IAU[n0]][0][0] += w00;
                CRS_MATRIX[CRS_IAU[n0]][1][1] += w10;
                // Non Diagonal Entries
                for (i = istart; i < iend; i++) {
                    if (CRS_JA[i] == n1) {
                        CRS_MATRIX[i][0][0] += w01;
                        CRS_MATRIX[i][1][1] += w11;
                    }
                    if (CRS_JA[i] == n2) {
                        CRS_MATRIX[i][0][0] += w02;
                        CRS_MATRIX[i][1][1] += w12;
                    }
                }
            }
        } // Triangle Node Loop
    } // Triangle Loop
}

// *****************************************************************************
// *****************************************************************************
double TriangleWinslowSmoother::Solve_CRS_Linear_Equation(int Iteration, double Relax, double* F, double* G) {
    int ni, iNode, i, istart, iend, col, nrms;
    double RHSf, RHSg, det, df, dg, ds, rms;
    double InvMat[2][2];

    // Initialize
    nrms = 0;
    rms  = 0.0;
    
    // Compute the Inverse of Matrix
    for (iNode = 0; iNode < NNode; iNode++) {
        if (IBTag[iNode] > -1)
            continue;
        col = CRS_IAU[iNode];
        det = CRS_MATRIX[col][0][0]*CRS_MATRIX[col][1][1] -
                CRS_MATRIX[col][0][1]*CRS_MATRIX[col][1][0];
        
        InvMat[0][0] = +CRS_MATRIX[col][1][1] / det;
        InvMat[0][1] = -CRS_MATRIX[col][0][1] / det;
        InvMat[1][0] = -CRS_MATRIX[col][1][0] / det;
        InvMat[1][1] = +CRS_MATRIX[col][0][0] / det;
        
        CRS_MATRIX[col][0][0] = InvMat[0][0];
        CRS_MATRIX[col][0][1] = InvMat[0][1];
        CRS_MATRIX[col][1][0] = InvMat[1][0];
        CRS_MATRIX[col][1][1] = InvMat[1][1];
    }

    
    for (ni = 0; ni < Iteration; ni++) {
        rms = 0.0;
        nrms = 0;
        for (iNode = 0; iNode < NNode; iNode++) {
            if (IBTag[iNode] > -1)
                continue;
            
            RHSf = 0.0;
            RHSg = 0.0;
            istart = CRS_IA[iNode];
            iend   = CRS_IA[iNode+1];
            // Transfer Off-diagonal elements to RHS
            for (i = istart; i < iend; i++) {
                col = CRS_JA[i];
                if (col == iNode)
                    continue;
                RHSf -= (CRS_MATRIX[i][0][0]*F[col] + CRS_MATRIX[i][0][1]*G[col]);
                RHSg -= (CRS_MATRIX[i][1][0]*F[col] + CRS_MATRIX[i][1][1]*G[col]);
            }
            
            // Mutiply RHS with Inverse of Diagonal Matrix
            // And Update F and G
            col = CRS_IAU[iNode];
            df = (CRS_MATRIX[col][0][0]*RHSf + CRS_MATRIX[col][0][1]*RHSg) - F[iNode];
            dg = (CRS_MATRIX[col][1][0]*RHSf + CRS_MATRIX[col][1][1]*RHSg) - G[iNode];
            ds = df*df + dg*dg;
            rms +=ds;
            nrms++;
            F[iNode] += Relax*df;
            G[iNode] += Relax*dg;
        }
        rms /= ((double)nrms);
        rms = sqrt(rms);
    }
    return rms;
}

