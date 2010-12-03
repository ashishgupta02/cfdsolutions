/*******************************************************************************
 * File:        Euler2D_Mesh_LinearElasticSmoother.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "Utils.h"
#include "MUtils.h"
#include "Euler2D_Mesh_LinearElasticSmoother.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Mesh_LinearElasticSmoother::Euler2D_Mesh_LinearElasticSmoother() {
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Init() {
    cell                           = NULL;
    edge                           = NULL;
    node                           = NULL;
    boundaryEdge                   = NULL;
    boundaryNode                   = NULL;
    BNTag                          = NULL;
    LESmoothBlockMatrix.VectorSize = 0;
    LESmoothBlockMatrix.BlockSize  = 0;
    LESmoothBlockMatrix.CRSSize    = 0;
    LESmoothBlockMatrix.A          = NULL;
    LESmoothBlockMatrix.B          = NULL;
    LESmoothBlockMatrix.IA         = NULL;
    LESmoothBlockMatrix.IAU        = NULL;
    LESmoothBlockMatrix.JA         = NULL;
    LESmoothBlockMatrix.X          = NULL;
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Reset() {
    int i, j;

    // Un-Share the Mesh Data Structure
    cell         = NULL;
    edge         = NULL;
    node         = NULL;
    boundaryEdge = NULL;
    boundaryNode = NULL;
    BNTag        = NULL;

    // Un-Share the common data in the CRS Matrix
    LESmoothBlockMatrix.IA  = NULL;
    LESmoothBlockMatrix.JA  = NULL;
    LESmoothBlockMatrix.IAU = NULL;
    // Free the memory
    if (LESmoothBlockMatrix.A != NULL) {
        for (i = 0; i < LESmoothBlockMatrix.CRSSize; i++) {
            if (LESmoothBlockMatrix.A[i] != NULL) {
                for (j = 0; j < LESmoothBlockMatrix.BlockSize; j++) {
                    if (LESmoothBlockMatrix.A[i][j] != NULL)
                        free(LESmoothBlockMatrix.A[i][j]);
                }
                free(LESmoothBlockMatrix.A[i]);
            }
        }
        free(LESmoothBlockMatrix.A);
    }
    if (LESmoothBlockMatrix.B != NULL) {
        for (i = 0; i < LESmoothBlockMatrix.VectorSize ; i++) {
            if (LESmoothBlockMatrix.B[i] != NULL)
                free(LESmoothBlockMatrix.B[i]);
        }
        free(LESmoothBlockMatrix.B);
    }
    if (LESmoothBlockMatrix.X != NULL) {
        for (i = 0; i < LESmoothBlockMatrix.VectorSize ; i++) {
            if (LESmoothBlockMatrix.X[i] != NULL)
                free(LESmoothBlockMatrix.X[i]);
        }
        free(LESmoothBlockMatrix.X);
    }

    // Reset the values
    Init();
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Mesh_LinearElasticSmoother::~Euler2D_Mesh_LinearElasticSmoother() {
    // Free the Resource and Set default value
    Reset();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Initialize_Mesh_Smoother(MESH inputMesh,
        CELL* ptrCell, EDGE* ptrEdge, NODE* ptrNode, BOUNDARYEDGE* ptrBoundaryEdge,
        BOUNDARYNODE* ptrBoundaryNode, int* ptrBNTag, MC_CRS *Object) {
    // Share the data with the input mesh
    mesh.nnodes  = inputMesh.nnodes;
    mesh.ncells  = inputMesh.ncells;
    mesh.inside  = inputMesh.inside;
    mesh.nedges  = inputMesh.nedges;
    mesh.nbedges = inputMesh.nbedges;
    mesh.nbnodes = inputMesh.nbnodes;
    // Share the Mesh Data Structure
    cell         = ptrCell;
    edge         = ptrEdge;
    node         = ptrNode;
    boundaryEdge = ptrBoundaryEdge;
    boundaryNode = ptrBoundaryNode;
    BNTag        = ptrBNTag;
    // Set the Dimension of CRS Matrix
    VectorSize   = mesh.nnodes;
    BlockSize    = 2;

    // Create the CRS Matrix
    Create_CRS_LESmoothBlockMatrix(Object);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Mesh_Smoother_Prepare() {
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Mesh_Smoother_Finalize() {
    // Free Resource and Initialize
    Reset();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Mesh_Smoother(int MaxIter, double Relaxation, double *U, double *V) {
    double lrms;
    
    // Basic Checks
    if ((U == NULL) || (V == NULL))
        error("Mesh_Smoother: %s\n", "Invalid Input Arrays");
    
    // Build Global Matrix
    Compute_CRS_LESmoothBlockMatrix();
    
    // Linear Solver
    lrms = Solve_CRS_LESmoothBlockMatrix(MaxIter, Relaxation, U, V);
    printf("LESmooth LRMS = %10.5e\n", lrms);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Create_CRS_LESmoothBlockMatrix(MC_CRS *Object) {
    int i, j, k;

    // Get the value from Solver Block Matrix
    LESmoothBlockMatrix.VectorSize = VectorSize;
    LESmoothBlockMatrix.BlockSize  = BlockSize;
    LESmoothBlockMatrix.CRSSize    = Object->CRSSize;
    // Share the memory with Solver Block Matrix
    LESmoothBlockMatrix.IA         = Object->IA;
    LESmoothBlockMatrix.IAU        = Object->IAU;
    LESmoothBlockMatrix.JA         = Object->JA;

    // Allocate Memory for Design Specific Computations
    // Allocate Memory for CRS Matrix
    LESmoothBlockMatrix.A = NULL;
    LESmoothBlockMatrix.A = (double ***) malloc (LESmoothBlockMatrix.CRSSize*sizeof(double**));
#ifdef DEBUG
    if (LESmoothBlockMatrix.A == NULL)
        error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < LESmoothBlockMatrix.CRSSize; i++) {
        LESmoothBlockMatrix.A[i] = NULL;
        LESmoothBlockMatrix.A[i] = (double **) malloc (LESmoothBlockMatrix.BlockSize*sizeof(double*));
#ifdef DEBUG
        if (LESmoothBlockMatrix.A[i] == NULL)
            error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 2");
#endif
        for (j = 0; j < LESmoothBlockMatrix.BlockSize; j++) {
            LESmoothBlockMatrix.A[i][j] = NULL;
            LESmoothBlockMatrix.A[i][j] = (double *) malloc (LESmoothBlockMatrix.BlockSize*sizeof(double));
#ifdef DEBUG
            if (LESmoothBlockMatrix.A[i] == NULL)
                error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 3");
#endif
            for (k = 0; k < LESmoothBlockMatrix.BlockSize; k++)
                LESmoothBlockMatrix.A[i][j][k] = 0.0;
        }
    }

    // Allocate Memory of RHS
    LESmoothBlockMatrix.B = NULL;
    LESmoothBlockMatrix.B = (double **) malloc (LESmoothBlockMatrix.VectorSize*sizeof(double*));
#ifdef DEBUG
    if (LESmoothBlockMatrix.B == NULL)
        error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 4");
#endif
    for (i = 0; i < LESmoothBlockMatrix.VectorSize; i++) {
        LESmoothBlockMatrix.B[i] = NULL;
        LESmoothBlockMatrix.B[i] = (double *) malloc (LESmoothBlockMatrix.BlockSize*sizeof(double));
#ifdef DEBUG
        if (LESmoothBlockMatrix.B[i] == NULL)
            error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 5");
#endif
        for (j = 0; j < LESmoothBlockMatrix.BlockSize; j++)
            LESmoothBlockMatrix.B[i][j] = 0.0;
    }

    // Allocate Memory for X
    LESmoothBlockMatrix.X = NULL;
    LESmoothBlockMatrix.X = (double **) malloc (LESmoothBlockMatrix.VectorSize*sizeof(double*));
#ifdef DEBUG
    if (LESmoothBlockMatrix.X == NULL)
        error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 6");
#endif
    for (i = 0; i < LESmoothBlockMatrix.VectorSize; i++) {
        LESmoothBlockMatrix.X[i] = NULL;
        LESmoothBlockMatrix.X[i] = (double *) malloc (LESmoothBlockMatrix.BlockSize*sizeof(double));
#ifdef DEBUG
        if (LESmoothBlockMatrix.X[i] == NULL)
            error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 7");
#endif
        for (j = 0; j < LESmoothBlockMatrix.BlockSize; j++)
            LESmoothBlockMatrix.X[i][j] = 0.0;
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Compute_CRS_LESmoothBlockMatrix() {
    int i, j, k;
    int iNode, iCell, n0, n1, n2, istart, iend, idgn;
    double nx0, ny0, nx1, ny1, nx2, ny2, tx, ty, area;
    double M0_00, M0_01, M0_10, M0_11;
    double M1_00, M1_01, M1_10, M1_11;
    double M2_00, M2_01, M2_10, M2_11;
    double alpha11, alpha12, alpha21, alpha22;
    double theta11, theta12, theta21, theta22;
    double YoungsModulus, PoissionRatio;

    // Initialize
    n0 = -1;
    n1 = -1;
    n2 = -1;

    // Initialize the CRS Matrix
    for (i = 0; i < LESmoothBlockMatrix.CRSSize; i++) {
        for (j = 0; j < LESmoothBlockMatrix.BlockSize; j++)
            for (k = 0; k < LESmoothBlockMatrix.BlockSize; k++)
                LESmoothBlockMatrix.A[i][j][k] = 0.0;
    }
    
    // For Triangles and computing all node contribution except Boundary Nodes
    for (iCell = 0; iCell < mesh.inside; iCell++) {

        // Get Youngs Modulus and Poission Ratio of Triangle
        YoungsModulus = Get_Youngs_Modulus(iCell);
        PoissionRatio = Get_Poission_Ratio(iCell);

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
                    n0 = cell[iCell].node1;
                    n1 = cell[iCell].node2;
                    n2 = cell[iCell].node3;
                    break;
                case 1:
                    n0 = cell[iCell].node2;
                    n1 = cell[iCell].node3;
                    n2 = cell[iCell].node1;
                    break;
                case 2:
                    n0 = cell[iCell].node3;
                    n1 = cell[iCell].node1;
                    n2 = cell[iCell].node2;
                    break;
                default:
                    break;
            }

            // Check if Node0 is Boundary
            if (BNTag[n0] > -1) {
                // Diagonal Entry and Ignore offdiagonal as they will be zero
                // Since Boundary Nodes are fixed and diagonal will be Identity 'I'
                idgn = LESmoothBlockMatrix.IAU[n0];
                LESmoothBlockMatrix.A[idgn][0][0] = 1.0;
                LESmoothBlockMatrix.A[idgn][1][1] = 1.0;
            } else {
                // Compute the Normals
                nx0 =  (node[n2].y - node[n1].y);
                ny0 = -(node[n2].x - node[n1].x);
                nx1 =  (node[n0].y - node[n2].y);
                ny1 = -(node[n0].x - node[n2].x);
                nx2 =  (node[n1].y - node[n0].y);
                ny2 = -(node[n1].x - node[n0].x);

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
                istart = LESmoothBlockMatrix.IA[n0];
                iend   = LESmoothBlockMatrix.IA[n0 + 1];
                // Diagonal Entry
                idgn = LESmoothBlockMatrix.IAU[n0];
                LESmoothBlockMatrix.A[idgn][0][0] += M0_00;
                LESmoothBlockMatrix.A[idgn][0][1] += M0_01;
                LESmoothBlockMatrix.A[idgn][1][0] += M0_10;
                LESmoothBlockMatrix.A[idgn][1][1] += M0_11;
                // Non Diagonal Entries
                for (i = istart; i < iend; i++) {
                    if (LESmoothBlockMatrix.JA[i] == n1) {
                        LESmoothBlockMatrix.A[i][0][0] += M1_00;
                        LESmoothBlockMatrix.A[i][0][1] += M1_01;
                        LESmoothBlockMatrix.A[i][1][0] += M1_10;
                        LESmoothBlockMatrix.A[i][1][1] += M1_11;
                    }
                    if (LESmoothBlockMatrix.JA[i] == n2) {
                        LESmoothBlockMatrix.A[i][0][0] += M2_00;
                        LESmoothBlockMatrix.A[i][0][1] += M2_01;
                        LESmoothBlockMatrix.A[i][1][0] += M2_10;
                        LESmoothBlockMatrix.A[i][1][1] += M2_11;
                    }
                }
            }
        } // Triangle Node Loop
    } // Triangle Loop
}

// *****************************************************************************
// *****************************************************************************
double Euler2D_Mesh_LinearElasticSmoother::Solve_CRS_LESmoothBlockMatrix(int Iteration, double Relax, double* U, double* V) {
    int ni, iNode, i, istart, iend, col, nrms, MaxIter = 100000;
    double RHSf, RHSg, det, du, dv, ds, rms;
    double InvMat[2][2];

    // Compute the Inverse of Matrix
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        if (BNTag[iNode] > -1)
            continue;
        col = LESmoothBlockMatrix.IAU[iNode];
        det = LESmoothBlockMatrix.A[col][0][0]*LESmoothBlockMatrix.A[col][1][1] -
                LESmoothBlockMatrix.A[col][0][1]*LESmoothBlockMatrix.A[col][1][0];

        InvMat[0][0] = +LESmoothBlockMatrix.A[col][1][1] / det;
        InvMat[0][1] = -LESmoothBlockMatrix.A[col][0][1] / det;
        InvMat[1][0] = -LESmoothBlockMatrix.A[col][1][0] / det;
        InvMat[1][1] = +LESmoothBlockMatrix.A[col][0][0] / det;

        LESmoothBlockMatrix.A[col][0][0] = InvMat[0][0];
        LESmoothBlockMatrix.A[col][0][1] = InvMat[0][1];
        LESmoothBlockMatrix.A[col][1][0] = InvMat[1][0];
        LESmoothBlockMatrix.A[col][1][1] = InvMat[1][1];
    }

    if (Iteration <= 0)
        Iteration = MaxIter;

#ifdef DEBUG
    info("Iteration  RMS");
#endif
    for (ni = 0; ni < Iteration; ni++) {
        rms = 0.0;
        nrms = 0;
        for (iNode = 0; iNode < mesh.nnodes; iNode++) {
            if (BNTag[iNode] > -1)
                continue;

            RHSf = 0.0;
            RHSg = 0.0;
            istart = LESmoothBlockMatrix.IA[iNode];
            iend   = LESmoothBlockMatrix.IA[iNode+1];
            // Transfer Off-diagonal elements to RHS
            for (i = istart; i < iend; i++) {
                col = LESmoothBlockMatrix.JA[i];
                if (col == iNode)
                    continue;
                RHSf -= (LESmoothBlockMatrix.A[i][0][0]*U[col] + LESmoothBlockMatrix.A[i][0][1]*V[col]);
                RHSg -= (LESmoothBlockMatrix.A[i][1][0]*U[col] + LESmoothBlockMatrix.A[i][1][1]*V[col]);
            }

            // Mutiply RHS with Inverse of Diagonal Matrix
            // And Update U and V
            col = LESmoothBlockMatrix.IAU[iNode];
            du = (LESmoothBlockMatrix.A[col][0][0]*RHSf + LESmoothBlockMatrix.A[col][0][1]*RHSg) - U[iNode];
            dv = (LESmoothBlockMatrix.A[col][1][0]*RHSf + LESmoothBlockMatrix.A[col][1][1]*RHSg) - V[iNode];
            ds = du*du + dv*dv;
            rms +=ds;
            nrms++;
            U[iNode] += Relax*du;
            V[iNode] += Relax*dv;
        }
        rms /= nrms;
        rms = sqrt(rms);

#ifdef DEBUG
        if (!(ni%10))
            info("%9d %12.5e", ni, rms);
#endif

        if (rms < (DBL_EPSILON*10.0)) {
#ifdef DEBUG
            info("%9d %12.5e", ni, rms);
#endif
            ni = MaxIter;
        }
    }
    return rms;
}

// *****************************************************************************
// *****************************************************************************
double Euler2D_Mesh_LinearElasticSmoother::Get_Poission_Ratio(int CellID) {
    double Poission = 0.25;

    if (CellID < 0)
        Poission = 0.0;

    return Poission;
}

// *****************************************************************************
// *****************************************************************************
double Euler2D_Mesh_LinearElasticSmoother::Get_Youngs_Modulus(int CellID) {
    double YoungsModulus = 0.0;

    if (CellID < 0) {
        YoungsModulus = 0.0;
    } else {
        double AspectRatio = 0.0;
        double ConditionNum = 0.0;
        double Area = 0.0;
        AspectRatio  = Compute_Tri_AspectRatio(cell[CellID].node1, cell[CellID].node2, cell[CellID].node3);
        ConditionNum = Compute_Tri_ConditionNumber(cell[CellID].node1, cell[CellID].node2, cell[CellID].node3);
        Area = Compute_Tri_Area(cell[CellID].node1, cell[CellID].node2, cell[CellID].node3);
        YoungsModulus = MAX(AspectRatio, ConditionNum)/(Area*Area);
    }

    return YoungsModulus;
}

// *****************************************************************************
// *****************************************************************************
double Euler2D_Mesh_LinearElasticSmoother::Compute_Tri_ConditionNumber(int n0, int n1, int n2) {
    double W[4], InW[4], A[4], InA[4], tmp[4];
    double ux, uy, vx, vy, cn;

    ux = node[n1].x - node[n0].x;
    uy = node[n1].y - node[n0].y;
    vx = node[n2].x - node[n0].x;
    vy = node[n2].y - node[n0].y;

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
double Euler2D_Mesh_LinearElasticSmoother::Compute_Tri_AspectRatio(int n0, int n1, int n2) {
    double a, b, c, s;
    double AspectRatio;

    a = sqrt((node[n1].x-node[n0].x)*(node[n1].x-node[n0].x)+(node[n1].y-node[n0].y)*(node[n1].y-node[n0].y));
    b = sqrt((node[n2].x-node[n1].x)*(node[n2].x-node[n1].x)+(node[n2].y-node[n1].y)*(node[n2].y-node[n1].y));
    c = sqrt((node[n0].x-node[n2].x)*(node[n0].x-node[n2].x)+(node[n0].y-node[n2].y)*(node[n0].y-node[n2].y));
    s = 0.5*(a+b+c);
    AspectRatio = (a*b*c)/(8.0*(s-a)*(s-b)*(s-c));

    return AspectRatio;
}

// *****************************************************************************
// *****************************************************************************
double Euler2D_Mesh_LinearElasticSmoother::Compute_Tri_Area(int n0, int n1, int n2) {
    double Area;

    Area = 0.5*((node[n1].x-node[n0].x)*(node[n2].y-node[n0].y)
            -(node[n1].y-node[n0].y)*(node[n2].x-node[n0].x));
    
    return Area;
}

