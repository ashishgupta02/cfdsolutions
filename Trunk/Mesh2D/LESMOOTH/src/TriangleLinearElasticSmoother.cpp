/*******************************************************************************
 * File:        TriangleLinearElasticSmoother.cpp
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
#include "TriangleLinearElasticSmoother.h"

// *****************************************************************************
// *****************************************************************************
TriangleLinearElasticSmoother::TriangleLinearElasticSmoother() {

    printf("\n====================================================");
    printf("\n    Mesh Module : Linear Elastic Smoothing          ");
    printf("\n====================================================\n");

    // Initialize the Data
    Init();

    return;
}

// *****************************************************************************
// *****************************************************************************
void TriangleLinearElasticSmoother::Init() {
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
TriangleLinearElasticSmoother::~TriangleLinearElasticSmoother() {
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
void TriangleLinearElasticSmoother::LESmooth() {
    int iNode;
    double lrms, relax, dx, dy;
    double *U  = NULL;
    double *V  = NULL;
    
    // Allocate Memory for Computations of Perturbation Variables
    U  = new double[NNode];
    V  = new double[NNode];

#ifdef DEBUG
    if ((U == NULL) || (V == NULL))
        error("LESmooth: %s\n", "Error Allocating Memory");
#endif
    
    // Store Original Mesh X, Y to Compute Perturbation Distance After Boundary Rotate
    for (iNode = 0; iNode < NNode;  iNode++) {
        U[iNode]  = X[iNode];
        V[iNode]  = Y[iNode];
    }

    // Create Interior and Boundary Node Tagging
    Create_Interior_Boundary_Tag();

    // Rotate Required Boundaries
    Rotate_Boundary(1);
    
    // Compute the Boundary Rotation Purturbance
    // And Retore the Original Mesh
    for (iNode = 0; iNode < NNode; iNode++) {
        if (IBTag[iNode] > -1) {
            dx = X[iNode] - U[iNode];
            dy = Y[iNode] - V[iNode];
            X[iNode] = U[iNode];
            Y[iNode] = V[iNode];
            U[iNode] = dx;
            V[iNode] = dy;
        } else {
            U[iNode] = 0.0;
            V[iNode] = 0.0;
        }
    }
    
    // Create All Possible Connectivities
    Create_Connectivity(0);

    // Create Compress Row Storage
    Create_CRS();

    int MaxIter = 0;
    std::cout << "Input Linear Elasticity Smoothing Iterations : ";
    std::cin  >> MaxIter;
    
    // Optimization Loop Starts
    info("Linear Elasticity Mesh Smoothing Starts");
    printf("---------------------------------------------------\n");

    // Build Global Matrix
    Compute_CRS_Matrix();
    
    // Linear Solver
    relax = 0.5;
    lrms = Solve_CRS_Linear_Equation(MaxIter, relax, U, V);
    
    // Update the Final Node Positions
    for (iNode = 0; iNode < NNode; iNode++) {
        X[iNode] += U[iNode];
        Y[iNode] += V[iNode];
    }
    
    printf("---------------------------------------------------\n");
    
    // Free the Used Memory
    delete []U;
    delete []V;
}

// *****************************************************************************
// *****************************************************************************
void TriangleLinearElasticSmoother::Rotate_Boundary(int Mode) {
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
void TriangleLinearElasticSmoother::Create_CRS() {
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
void TriangleLinearElasticSmoother::Compute_CRS_Matrix() {
    int iNode, iTri, n0, n1, n2, i, istart, iend;
    double nx0, ny0, nx1, ny1, nx2, ny2, tx, ty, area;
    double M0_00, M0_01, M0_10, M0_11;
    double M1_00, M1_01, M1_10, M1_11;
    double M2_00, M2_01, M2_10, M2_11;
    double alpha11, alpha12, alpha21, alpha22;
    double theta11, theta12, theta21, theta22;
    double YoungsModulus, PoissionRatio;
    
    if (NNode <= 0)
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

        // Get Youngs Modulus and Poission Ratio of Triangle
        YoungsModulus = Get_Youngs_Modulus(iTri);
        PoissionRatio = Get_Poission_Ratio(iTri);

        // Compute alphas, thetas as the are constant for a given triangle
        alpha11 = YoungsModulus*(1.0 - PoissionRatio)/((1.0 + PoissionRatio) * (1.0 - 2.0*PoissionRatio));
        alpha12 = YoungsModulus/(2.0*(1.0 + PoissionRatio));
        alpha21 = alpha12;
        alpha22 = alpha11;
        theta11 = YoungsModulus*PoissionRatio/((1.0 + PoissionRatio) * (1.0 - 2.0*PoissionRatio));;
        theta12 = YoungsModulus/(2.0*(1.0 + PoissionRatio));
        theta21 = theta12;
        theta22 = theta11;
        
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
                // Diagonal Entry and Ignore offdiagonal as they will be zero
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
                
                // Assign tx, ty
                tx = nx0;
                ty = ny0;

                // Calculate Diagonal Terms
                // M0_00, M0_01, M0_10, M0_11 Contribution of Node0 itself
                M0_00 = -0.5*(alpha11*nx0*tx + alpha12*ny0*ty)/area;
                M0_01 = -0.5*(theta11*ny0*tx + theta12*nx0*ty)/area;
                M0_10 = -0.5*(theta21*ny0*tx + theta22*nx0*ty)/area;
                M0_11 = -0.5*(alpha21*nx0*tx + alpha22*ny0*ty)/area;

                // Calculate Off-Diagonal Terms
                // M1_00, M1_01, M1_10, M1_11 Contribution of Node1 to Node0
                M1_00 = -0.5*(alpha11*nx1*tx + alpha12*ny1*ty)/area;
                M1_01 = -0.5*(theta11*ny1*tx + theta12*nx1*ty)/area;
                M1_10 = -0.5*(theta21*ny1*tx + theta22*nx1*ty)/area;
                M1_11 = -0.5*(alpha21*nx1*tx + alpha22*ny1*ty)/area;

                // M2_00, M2_01, M2_10, M2_11 Contribution of Node2 to Node0
                M2_00 = -0.5*(alpha11*nx2*tx + alpha12*ny2*ty)/area;
                M2_01 = -0.5*(theta11*ny2*tx + theta12*nx2*ty)/area;
                M2_10 = -0.5*(theta21*ny2*tx + theta22*nx2*ty)/area;
                M2_11 = -0.5*(alpha21*nx2*tx + alpha22*ny2*ty)/area;

                // Put values in CRS Matrix
                istart = CRS_IA[n0];
                iend   = CRS_IA[n0 + 1];
                // Diagonal Entry
                CRS_MATRIX[CRS_IAU[n0]][0][0] += M0_00;
                CRS_MATRIX[CRS_IAU[n0]][0][1] += M0_01;
                CRS_MATRIX[CRS_IAU[n0]][1][0] += M0_10;
                CRS_MATRIX[CRS_IAU[n0]][1][1] += M0_11;
                // Non Diagonal Entries
                for (i = istart; i < iend; i++) {
                    if (CRS_JA[i] == n1) {
                        CRS_MATRIX[i][0][0] += M1_00;
                        CRS_MATRIX[i][0][1] += M1_01;
                        CRS_MATRIX[i][1][0] += M1_10;
                        CRS_MATRIX[i][1][1] += M1_11;
                    }
                    if (CRS_JA[i] == n2) {
                        CRS_MATRIX[i][0][0] += M2_00;
                        CRS_MATRIX[i][0][1] += M2_01;
                        CRS_MATRIX[i][1][0] += M2_10;
                        CRS_MATRIX[i][1][1] += M2_11;
                    }
                }
            }
        } // Triangle Node Loop
    } // Triangle Loop
}

// *****************************************************************************
// *****************************************************************************
double TriangleLinearElasticSmoother::Solve_CRS_Linear_Equation(int Iteration, double Relax, double* F, double* G) {
    int ni, iNode, i, istart, iend, col, nrms, MaxIter = 100000;
    double RHSf, RHSg, det, df, dg, ds, rms;
    double InvMat[2][2];

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

    if (Iteration <= 0)
        Iteration = MaxIter;

    info("Iteration  RMS");
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
        rms /= nrms;
        rms = sqrt(rms);

        if (!(ni%10))
            info("%9d %12.5e", ni, rms);

        if (rms < (DBL_EPSILON*10.0)) {
            info("%9d %12.5e", ni, rms);
            ni = MaxIter;
        }
    }
    return rms;
}

// *****************************************************************************
// *****************************************************************************
double TriangleLinearElasticSmoother::Get_Poission_Ratio(int TriID) {
    double Poission = 0.25;

    if (TriID < 0)
        Poission = 0.0;

    return Poission;
}

// *****************************************************************************
// *****************************************************************************
double TriangleLinearElasticSmoother::Get_Youngs_Modulus(int TriID) {
    double YoungsModulus = 0.0;

    if (TriID < 0) {
        YoungsModulus = 0.0;
    } else {
        double AspectRatio = 0.0;
        double ConditionNum = 0.0;
        double Area = 0.0;
        AspectRatio  = Compute_Tri_AspectRatio(Tri[TriID][0], Tri[TriID][1], Tri[TriID][2]);
        ConditionNum = Compute_Tri_ConditionNumber(Tri[TriID][0], Tri[TriID][1], Tri[TriID][2]);
        Area = Compute_Tri_Area(Tri[TriID][0], Tri[TriID][1], Tri[TriID][2]);
        YoungsModulus = MAX(AspectRatio, ConditionNum)/(Area*Area);
    }
    
    return YoungsModulus;
}

// *****************************************************************************
// *****************************************************************************
double TriangleLinearElasticSmoother::Compute_Tri_ConditionNumber(int n0, int n1, int n2) {
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
double TriangleLinearElasticSmoother::Compute_Tri_AspectRatio(int n0, int n1, int n2) {
    double a, b, c, s;
    double AspectRatio;

    a = sqrt((X[n1]-X[n0])*(X[n1]-X[n0])+(Y[n1]-Y[n0])*(Y[n1]-Y[n0]));
    b = sqrt((X[n2]-X[n1])*(X[n2]-X[n1])+(Y[n2]-Y[n1])*(Y[n2]-Y[n1]));
    c = sqrt((X[n0]-X[n2])*(X[n0]-X[n2])+(Y[n0]-Y[n2])*(Y[n0]-Y[n2]));
    s = 0.5*(a+b+c);
    AspectRatio = (a*b*c)/(8.0*(s-a)*(s-b)*(s-c));

    return AspectRatio;
}

// *****************************************************************************
// *****************************************************************************
double TriangleLinearElasticSmoother::Compute_Tri_Area(int n0, int n1, int n2) {
    double Area;

    Area = 0.5*((X[n1]-X[n0])*(Y[n2]-Y[n0])-(Y[n1]-Y[n0])*(X[n2]-X[n0]));

    return Area;
}

