/*******************************************************************************
 * File:        Euler2D_Mesh_LinearElasticSmoother.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "List.h"
#include "Utils.h"
#include "MUtils.h"
#include "Euler2D_Mesh_LinearElasticSmoother.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Mesh_LinearElasticSmoother::Euler2D_Mesh_LinearElasticSmoother() {
#ifdef VERBOSE
    printf("=============================================================================\n");
    printf("      Euler2D : Linear Elastic Mesh Smoother                                 \n");
    printf("=============================================================================\n");
#endif
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Init() {
    BlockSize                      = 0;
    VectorSize                     = 0;
    changeIndex                    = 0;
    Parent                         = 0;
    mesh.inside                    = 0;
    mesh.nbedges                   = 0;
    mesh.nbnodes                   = 0;
    mesh.ncells                    = 0;
    mesh.nedges                    = 0;
    mesh.nnodes                    = 0;
    cell                           = NULL;
    edge                           = NULL;
    node                           = NULL;
    boundaryEdge                   = NULL;
    boundaryNode                   = NULL;
    BNTag                          = NULL;
    LESmoothBlockMatrix.nROW       = 0;
    LESmoothBlockMatrix.nCOL       = 0;
    LESmoothBlockMatrix.Block_nRow = 0;
    LESmoothBlockMatrix.Block_nCol = 0;
    LESmoothBlockMatrix.DIM        = 0;
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
    if (Parent == 0) {
        cell         = NULL;
        edge         = NULL;
        node         = NULL;
        boundaryEdge = NULL;
        boundaryNode = NULL;
        BNTag        = NULL;
    } else {
        if (cell != NULL)
        free(cell);

        if (edge != NULL)
            free(edge);

        if (node != NULL)
            free(node);

        if (boundaryEdge != NULL)
            free(boundaryEdge);

        if (boundaryNode != NULL)
            free(boundaryNode);

        if (BNTag != NULL)
            free(BNTag);
    }

    // Un-Share the common data in the CRS Matrix
    LESmoothBlockMatrix.IA  = NULL;
    LESmoothBlockMatrix.JA  = NULL;
    LESmoothBlockMatrix.IAU = NULL;
    // Free the memory
    if (LESmoothBlockMatrix.A != NULL) {
        for (i = 0; i < LESmoothBlockMatrix.DIM; i++) {
            if (LESmoothBlockMatrix.A[i] != NULL) {
                for (j = 0; j < LESmoothBlockMatrix.Block_nRow; j++) {
                    if (LESmoothBlockMatrix.A[i][j] != NULL)
                        free(LESmoothBlockMatrix.A[i][j]);
                }
                free(LESmoothBlockMatrix.A[i]);
            }
        }
        free(LESmoothBlockMatrix.A);
    }
    if (LESmoothBlockMatrix.B != NULL) {
        for (i = 0; i < LESmoothBlockMatrix.nROW ; i++) {
            if (LESmoothBlockMatrix.B[i] != NULL)
                free(LESmoothBlockMatrix.B[i]);
        }
        free(LESmoothBlockMatrix.B);
    }
    if (LESmoothBlockMatrix.X != NULL) {
        for (i = 0; i < LESmoothBlockMatrix.nROW ; i++) {
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
void Euler2D_Mesh_LinearElasticSmoother::Initialize_Mesh_Smoother(const char* FileName, MC_CRS* Object) {
    Parent = 1;

    // Read Mesh File
    WKA_MeshReader(FileName);

    // Get Boundary Nodes and Tag Nodes
    Tag_Boundary_Nodes();

    // Extract Connectivity
    WKA_ExtractCells();

    // Compute the Expensive Geometric Properties
    Compute_Geometric_Properties();

    // Set the Dimension of CRS Matrix
    VectorSize   = mesh.nnodes;
    BlockSize    = 2;
    
    // Create the CRS Matrix
    Create_CRS_LESmoothBlockMatrix(Object);
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
    Parent       = 0;
    
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
    
#ifdef VERBOSE
    printf("LESmooth LRMS = %10.5e\n", lrms);
#endif
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Create_CRS_LESmoothBlockMatrix(MC_CRS *Object) {
    int i, j, k;

    // Get the value from Solver Block Matrix
    LESmoothBlockMatrix.nROW = VectorSize;
    LESmoothBlockMatrix.nCOL = VectorSize;
    LESmoothBlockMatrix.Block_nRow = BlockSize;
    LESmoothBlockMatrix.Block_nCol = BlockSize;
    LESmoothBlockMatrix.DIM        = Object->DIM;
    // Share the memory with Solver Block Matrix
    LESmoothBlockMatrix.IA         = Object->IA;
    LESmoothBlockMatrix.IAU        = Object->IAU;
    LESmoothBlockMatrix.JA         = Object->JA;

    // Allocate Memory for Design Specific Computations
    // Allocate Memory for CRS Matrix
    LESmoothBlockMatrix.A = NULL;
    LESmoothBlockMatrix.A = (double ***) malloc (LESmoothBlockMatrix.DIM*sizeof(double**));
#ifdef DEBUG
    if (LESmoothBlockMatrix.A == NULL)
        error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < LESmoothBlockMatrix.DIM; i++) {
        LESmoothBlockMatrix.A[i] = NULL;
        LESmoothBlockMatrix.A[i] = (double **) malloc (LESmoothBlockMatrix.Block_nRow*sizeof(double*));
#ifdef DEBUG
        if (LESmoothBlockMatrix.A[i] == NULL)
            error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 2");
#endif
        for (j = 0; j < LESmoothBlockMatrix.Block_nRow; j++) {
            LESmoothBlockMatrix.A[i][j] = NULL;
            LESmoothBlockMatrix.A[i][j] = (double *) malloc (LESmoothBlockMatrix.Block_nCol*sizeof(double));
#ifdef DEBUG
            if (LESmoothBlockMatrix.A[i] == NULL)
                error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 3");
#endif
            for (k = 0; k < LESmoothBlockMatrix.Block_nCol; k++)
                LESmoothBlockMatrix.A[i][j][k] = 0.0;
        }
    }

    // Allocate Memory of RHS
    LESmoothBlockMatrix.B = NULL;
    LESmoothBlockMatrix.B = (double **) malloc (LESmoothBlockMatrix.nROW*sizeof(double*));
#ifdef DEBUG
    if (LESmoothBlockMatrix.B == NULL)
        error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 4");
#endif
    for (i = 0; i < LESmoothBlockMatrix.nROW; i++) {
        LESmoothBlockMatrix.B[i] = NULL;
        LESmoothBlockMatrix.B[i] = (double *) malloc (LESmoothBlockMatrix.Block_nRow*sizeof(double));
#ifdef DEBUG
        if (LESmoothBlockMatrix.B[i] == NULL)
            error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 5");
#endif
        for (j = 0; j < LESmoothBlockMatrix.Block_nRow; j++)
            LESmoothBlockMatrix.B[i][j] = 0.0;
    }

    // Allocate Memory for X
    LESmoothBlockMatrix.X = NULL;
    LESmoothBlockMatrix.X = (double **) malloc (LESmoothBlockMatrix.nROW*sizeof(double*));
#ifdef DEBUG
    if (LESmoothBlockMatrix.X == NULL)
        error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 6");
#endif
    for (i = 0; i < LESmoothBlockMatrix.nROW; i++) {
        LESmoothBlockMatrix.X[i] = NULL;
        LESmoothBlockMatrix.X[i] = (double *) malloc (LESmoothBlockMatrix.Block_nRow*sizeof(double));
#ifdef DEBUG
        if (LESmoothBlockMatrix.X[i] == NULL)
            error("Create_CRS_LESmoothBlockMatrix: %s\n", "Error Allocating Memory 7");
#endif
        for (j = 0; j < LESmoothBlockMatrix.Block_nRow; j++)
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
    for (i = 0; i < LESmoothBlockMatrix.DIM; i++) {
        for (j = 0; j < LESmoothBlockMatrix.Block_nRow; j++)
            for (k = 0; k < LESmoothBlockMatrix.Block_nCol; k++)
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

// *****************************************************************************
/*                                                                         */
/* Reads mesh using edge pointers and initializes solution to linear       */
/*                                                                         */
// *****************************************************************************
/* File Discription                                                        */
/* NNode NEdges NCells NBEdges Dummy Dummy                                 */
/* Node1 Node2 Cell1 Cell2                                                 */
/* .......................                                                 */
/* sedges                                                                  */
/* BEdge BType Const1 Const2                                               */
/* .......................                                                 */
/* coordinates                                                             */
/* X Y                                                                     */
/* .......................                                                 */
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::WKA_MeshReader(const char* FileName) {
    int i, icount;
    int idum, iret;
    FILE *inputMesh;
    char dumstring[100];
    char *cdum;

#ifdef VERBOSE
    printf("=============================================================================\n");
    info("Reading Mesh File %s", FileName);
#endif

    // Open Mesh File
    if ((inputMesh = fopen(FileName, "r")) == (FILE *) NULL)
        error("Unable to Open Mesh File %s", FileName);

    // Read mesh sizes
    iret = fscanf(inputMesh, "%d %d %d %d %d %d", &mesh.nnodes, &mesh.nedges, &mesh.ncells, &mesh.nbedges, &idum, &idum);

#ifdef VERBOSE
    info("NNodes = %d NEdges = %d NCells = %d NBEdge = %d", mesh.nnodes, mesh.nedges, mesh.ncells, mesh.nbedges);
#endif
    
    mesh.inside = mesh.ncells - mesh.nbedges;

    // Allocate Memory to store Connectivity
    edge = (EDGE *) calloc(mesh.nedges, sizeof (EDGE));
    cell = (CELL *) calloc(mesh.ncells, sizeof (CELL));
    boundaryEdge = (BOUNDARYEDGE *) calloc(mesh.nbedges, sizeof (BOUNDARYEDGE));

    icount = 0;
    for (i = 0; i < mesh.nedges; i++) {
        iret = fscanf(inputMesh, "%d %d %d %d", &edge[i].node1, &edge[i].node2, &edge[i].cell1, &edge[i].cell2);
        /* Subtract 1; mesh written assuming number starts at 1 */
        /* but number really starts at 0 because this is C      */
        if (changeIndex == 1) {
            edge[i].node1 -= 1;
            edge[i].node2 -= 1;
            edge[i].cell1 -= 1;
            edge[i].cell2 -= 1;
        }
    }

    /* Read boundary edges */
    /* Types: 1000 del(q).n = 0 */
    /*        2000 q = x        */
    /*        2001 q = c1       */
    iret = fscanf(inputMesh, "\n");
    cdum = fgets(dumstring, 100, inputMesh);
    for (i = 0; i < mesh.nbedges; i++) {
        iret = fscanf(inputMesh, "%d %d %lf %lf", &idum, &(boundaryEdge[i].bcType), &(boundaryEdge[i].c1), &(boundaryEdge[i].c2));
        if (changeIndex == 1)
            boundaryEdge[i].edgeNumber = idum - 1;
        else
            boundaryEdge[i].edgeNumber = idum;
    }

    /* Now read the mesh coordinates */
    iret = fscanf(inputMesh, "\n");
    cdum = fgets(dumstring, 100, inputMesh);
    node = (NODE *) calloc(mesh.nnodes, sizeof (NODE));
    for (i = 0; i < mesh.nnodes; i++)
        iret = fscanf(inputMesh, "%le %le", &node[i].x, &node[i].y);

    // Close file
    fclose(inputMesh);

#ifdef VERBOSE
    printf("=============================================================================\n");
#endif
}

// *****************************************************************************
/*                                                                         */
/* Extract cell-to-node pointers from edge pointers                        */
/*                                                                         */
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::WKA_ExtractCells(void) {
    int i;
    int imin;
    int *work = NULL;
    int node1, node2, node3, cell1, cell2;
    int edgeNumber;
    double areamin, areamax;

    /* Allocate a work array */
    /* These should also be preset to zero by calloc */
    work = (int *) calloc(mesh.ncells, sizeof (int));

    /* Loop over cells and tag the nodes */
    for (i = 0; i < mesh.ncells; i++) {
        cell[i].node1 = -99;
        cell[i].node2 = -99;
        cell[i].node3 = -99;
        work[i] = 0;
    }

    /* Set work array to -1 for ghost cells */
    for (i = 0; i < mesh.nbedges; i++) {
        edgeNumber = boundaryEdge[i].edgeNumber;
        cell2 = edge[edgeNumber].cell2;
        work[cell2] = -1;
    }

    /* Loop over all the edges and set the first two nodes in c2n */
    for (i = 0; i < mesh.nedges; i++) {
        node1 = edge[i].node1;
        node2 = edge[i].node2;
        cell1 = edge[i].cell1;
        cell2 = edge[i].cell2;

        if (work[cell1] != -1) {
            if (work[cell1] != 1) {
                cell[cell1].node1 = node1;
                cell[cell1].node2 = node2;
                work[cell1] = 1;
            }
        }
        if (work[cell2] != -1) {
            if (work[cell2] != 1) {
                cell[cell2].node1 = node2;
                cell[cell2].node2 = node1;
                work[cell2] = 1;
            }
        }
    }

    /* Loop over the edges again to get the third node */
    for (i = 0; i < mesh.nedges; i++) {
        node1 = edge[i].node1;
        node2 = edge[i].node2;
        cell1 = edge[i].cell1;
        cell2 = edge[i].cell2;

        if (cell[cell1].node3 == -99) {
            if ((node1 != cell[cell1].node1) && (node1 != cell[cell1].node2)) {
                cell[cell1].node3 = node1;
            } else if ((node2 != cell[cell1].node1) && (node2 != cell[cell1].node2)) {
                cell[cell1].node3 = node2;
            }
        }
        if (cell2 < mesh.inside) {
            if (cell[cell2].node3 == -99) {
                if ((node1 != cell[cell2].node1) && (node1 != cell[cell2].node2)) {
                    cell[cell2].node3 = node1;
                } else if ((node2 != cell[cell2].node1) && (node2 != cell[cell2].node2)) {
                    cell[cell2].node3 = node2;
                }
            }
        }
    }

    if (work != NULL)
        free(work);

    areamin = DBL_MAX;
    areamax = DBL_MIN;
    imin = -1;
    for (i = 0; i < mesh.inside; i++) {
        double dx1, dy1;
        double dx2, dy2;
        double dx3, dy3;
        node1 = cell[i].node1;
        node2 = cell[i].node2;
        node3 = cell[i].node3;
        dx1 = -(node[node3].x - node[node2].x);
        dy1 = node[node3].y - node[node2].y;
        dx2 = -(node[node1].x - node[node3].x);
        dy2 = node[node1].y - node[node3].y;
        dx3 = -(node[node2].x - node[node1].x);
        dy3 = node[node2].y - node[node1].y;
        cell[i].area = .5 * (dx3 * dy2 - dx2 * dy3);
        if (cell[i].area >= areamax) areamax = cell[i].area;
        if (cell[i].area <= areamin) {
            areamin = cell[i].area;
            imin = i;
        }
    }

#ifdef VERBOSE
    info("Cell Area: MAX = %lf MIN = %lf", areamax, areamin);
    info("Cell ID Min Area: %d Area = %le", imin, cell[imin].area);
    printf("=============================================================================\n");
#endif
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Tag_Boundary_Nodes() {
    int i, n1, n2, bctype, ibn;
    List bnode;

    // Tag the Boundary Nodes Points
    // Method of Tagging is not robust as boundary connecting points have Tag
    // which are tagged last.
    BNTag = (int *) calloc(mesh.nnodes, sizeof (int));
    for (i = 0; i < mesh.nnodes; i++)
        BNTag[i] = -1;

    for (i = 0; i < mesh.nbedges; i++) {
        n1 = edge[boundaryEdge[i].edgeNumber].node1;
        n2 = edge[boundaryEdge[i].edgeNumber].node2;
        bctype = -1;
        // BCType 0: 0-999
        if ((boundaryEdge[i].bcType >= 0) && (boundaryEdge[i].bcType <= 999))
            bctype = 0;
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[i].bcType >= 1000) && (boundaryEdge[i].bcType <= 1999))
            bctype = 1;
        // BCType 2: 2000-2999 : Dirchlet Boundary Condition
        if ((boundaryEdge[i].bcType >= 2000) && (boundaryEdge[i].bcType <= 2999))
            bctype = 2;
        // BCType 2: 3000-3999 : Freestream
        if ((boundaryEdge[i].bcType >= 3000) && (boundaryEdge[i].bcType <= 3999))
            bctype = 3;

        if (bctype != -1) {
            BNTag[n1] = bctype;
            BNTag[n2] = bctype;
        }
    }

    mesh.nbnodes = 0;
    for (i = 0; i < mesh.nnodes; i++) {
        if (BNTag[i] != -1)
            mesh.nbnodes++;
    }

#ifdef VERBOSE
    info("NBNodes = %d", mesh.nbnodes);
#endif
    
    boundaryNode = (BOUNDARYNODE *) calloc(mesh.nbnodes, sizeof (BOUNDARYNODE));
    ibn = 0;
    for (i = 0; i < mesh.nbedges; i++) {
        n1 = edge[boundaryEdge[i].edgeNumber].node1;
        n2 = edge[boundaryEdge[i].edgeNumber].node2;

        // BCType 0: 0-999
        if ((boundaryEdge[i].bcType >= 0) && (boundaryEdge[i].bcType <= 999)) {
            if (!bnode.Is_In_List(n1)) {
                bnode.Add_To_List(n1);
                boundaryNode[ibn].bcType     = BNTag[n1];
                boundaryNode[ibn].nodeNumber = n1;
                boundaryNode[ibn].constant   = 0.0;
                ibn++;
            }
            if (!bnode.Is_In_List(n2)) {
                bnode.Add_To_List(n2);
                boundaryNode[ibn].bcType     = BNTag[n2];
                boundaryNode[ibn].nodeNumber = n2;
                boundaryNode[ibn].constant   = 0.0;
                ibn++;
            }
        }
        // BCType 1: 1000-1999 : Solid Wall
        if ((boundaryEdge[i].bcType >= 1000) && (boundaryEdge[i].bcType <= 1999)) {
            if (!bnode.Is_In_List(n1)) {
                bnode.Add_To_List(n1);
                boundaryNode[ibn].bcType     = BNTag[n1];
                boundaryNode[ibn].nodeNumber = n1;
                boundaryNode[ibn].constant   = 0.0;
                ibn++;
            }
            if (!bnode.Is_In_List(n2)) {
                bnode.Add_To_List(n2);
                boundaryNode[ibn].bcType     = BNTag[n2];
                boundaryNode[ibn].nodeNumber = n2;
                boundaryNode[ibn].constant   = 0.0;
                ibn++;
            }
        }
        // BCType 2: 2000-2999 : Dirchlet Boundary Condition
        if ((boundaryEdge[i].bcType >= 2000) && (boundaryEdge[i].bcType <= 2999)) {
            if (BNTag[n1] != 2) {
                if (!bnode.Is_In_List(n1)) {
                    bnode.Add_To_List(n1);
                    boundaryNode[ibn].bcType     = BNTag[n1];
                    boundaryNode[ibn].nodeNumber = n1;
                    if (boundaryEdge[i].bcType == 2000) {
                        boundaryNode[ibn].constant = node[n1].x;
                    } else {
                        boundaryNode[ibn].constant = boundaryEdge[i].c1;
                    }
                    ibn++;
                }
            }
            if (BNTag[n2] != 2) {
                if (!bnode.Is_In_List(n2)) {
                    bnode.Add_To_List(n2);
                    boundaryNode[ibn].bcType     = BNTag[n2];
                    boundaryNode[ibn].nodeNumber = n2;
                    if (boundaryEdge[i].bcType == 2000) {
                        boundaryNode[ibn].constant = node[n2].x;
                    } else {
                        boundaryNode[ibn].constant = boundaryEdge[i].c1;
                    }
                    ibn++;
                }
            }
        }
        // BCType 2: 3000-3999 : Freestream
        if ((boundaryEdge[i].bcType >= 3000) && (boundaryEdge[i].bcType <= 3999)) {
            if (!bnode.Is_In_List(n1)) {
                bnode.Add_To_List(n1);
                boundaryNode[ibn].bcType     = BNTag[n1];
                boundaryNode[ibn].nodeNumber = n1;
                boundaryNode[ibn].constant   = 0.0;
                ibn++;
            }
            if (!bnode.Is_In_List(n2)) {
                bnode.Add_To_List(n2);
                boundaryNode[ibn].bcType     = BNTag[n2];
                boundaryNode[ibn].nodeNumber = n2;
                boundaryNode[ibn].constant   = 0.0;
                ibn++;
            }
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh_LinearElasticSmoother::Compute_Geometric_Properties() {
    int icell, inode, iedge;
    int n1, n2, n3;
    double ux, uy, vx, vy;

    // Initialize
    for (icell = 0; icell < mesh.ncells; icell++) {
        cell[icell].xc     = 0.0;
        cell[icell].yc     = 0.0;
        cell[icell].area   = 0.0;
        cell[icell].mag12  = 0.0;
        cell[icell].mag23  = 0.0;
        cell[icell].mag31  = 0.0;
        cell[icell].unx12  = 0.0;
        cell[icell].unx23  = 0.0;
        cell[icell].unx31  = 0.0;
        cell[icell].uny12  = 0.0;
        cell[icell].uny23  = 0.0;
        cell[icell].uny31  = 0.0;
        cell[icell].magc12 = 0.0;
        cell[icell].magc23 = 0.0;
        cell[icell].magc31 = 0.0;
        cell[icell].unxc12 = 0.0;
        cell[icell].unxc23 = 0.0;
        cell[icell].unxc31 = 0.0;
        cell[icell].unyc12 = 0.0;
        cell[icell].unyc23 = 0.0;
        cell[icell].unyc31 = 0.0;
    }

    for (inode = 0; inode < mesh.nnodes; inode++)
        node[inode].area = 0.0;

    for (iedge = 0; iedge < mesh.nedges; iedge++) {
        edge[iedge].unx = 0.0;
        edge[iedge].uny = 0.0;
        edge[iedge].mag = 0.0;
    }

    // Now Compute
    for (icell = 0; icell < mesh.ncells; icell++) {
        // Check for Ghost Cells
        if ((cell[icell].node1 < 0) || (cell[icell].node2 < 0) || (cell[icell].node2 < 0))
            continue;

        // Get the Nodes
        n1 = cell[icell].node1;
        n2 = cell[icell].node2;
        n3 = cell[icell].node3;

        // Get the Median Dual Properties
        // Compute the Centriod of Cell
        cell[icell].xc = (node[n1].x + node[n2].x + node[n3].x)/3.0;
        cell[icell].yc = (node[n1].y + node[n2].y + node[n3].y)/3.0;

        // Edge 1-2
        cell[icell].unxc12 = +(cell[icell].yc - 0.5*(node[n1].y + node[n2].y));
        cell[icell].unyc12 = -(cell[icell].xc - 0.5*(node[n1].x + node[n2].x));
        cell[icell].magc12 = sqrt(cell[icell].unxc12*cell[icell].unxc12 + cell[icell].unyc12*cell[icell].unyc12);
        cell[icell].unxc12 /= cell[icell].magc12;
        cell[icell].unyc12 /= cell[icell].magc12;

        // Edge 2-3
        cell[icell].unxc23 = +(cell[icell].yc - 0.5*(node[n2].y + node[n3].y));
        cell[icell].unyc23 = -(cell[icell].xc - 0.5*(node[n2].x + node[n3].x));
        cell[icell].magc23 = sqrt(cell[icell].unxc23*cell[icell].unxc23 + cell[icell].unyc23*cell[icell].unyc23);
        cell[icell].unxc23 /= cell[icell].magc23;
        cell[icell].unyc23 /= cell[icell].magc23;

        // Edge 3-1
        cell[icell].unxc31 = +(cell[icell].yc - 0.5*(node[n3].y + node[n1].y));
        cell[icell].unyc31 = -(cell[icell].xc - 0.5*(node[n3].x + node[n1].x));
        cell[icell].magc31 = sqrt(cell[icell].unxc31*cell[icell].unxc31 + cell[icell].unyc31*cell[icell].unyc31);
        cell[icell].unxc31 /= cell[icell].magc31;
        cell[icell].unyc31 /= cell[icell].magc31;

        // Compute Median Dual Area - Cell and Node Property
        ux = node[n2].x - node[n1].x;
        uy = node[n2].y - node[n1].y;
        vx = node[n3].x - node[n1].x;
        vy = node[n3].y - node[n1].y;
        cell[icell].area = 0.5*(ux*vy - uy*vx);
        node[n1].area += (1.0/3.0)*cell[icell].area;
        node[n2].area += (1.0/3.0)*cell[icell].area;
        node[n3].area += (1.0/3.0)*cell[icell].area;

        // Cell Property - Compute Edge Normals
        // Edge 1-2
        cell[icell].unx12 = +(node[n2].y - node[n1].y);
        cell[icell].uny12 = -(node[n2].x - node[n1].x);
        cell[icell].mag12 = sqrt(cell[icell].unx12*cell[icell].unx12 + cell[icell].uny12*cell[icell].uny12);
        cell[icell].unx12 /= cell[icell].mag12;
        cell[icell].uny12 /= cell[icell].mag12;

        // Edge 2-3
        cell[icell].unx23 = +(node[n3].y - node[n2].y);
        cell[icell].uny23 = -(node[n3].x - node[n2].x);
        cell[icell].mag23 = sqrt(cell[icell].unx23*cell[icell].unx23 + cell[icell].uny23*cell[icell].uny23);
        cell[icell].unx23 /= cell[icell].mag23;
        cell[icell].uny23 /= cell[icell].mag23;

        // Edge 3-1
        cell[icell].unx31 = +(node[n1].y - node[n3].y);
        cell[icell].uny31 = -(node[n1].x - node[n3].x);
        cell[icell].mag31 = sqrt(cell[icell].unx31*cell[icell].unx31 + cell[icell].uny31*cell[icell].uny31);
        cell[icell].unx31 /= cell[icell].mag31;
        cell[icell].uny31 /= cell[icell].mag31;
    }

    for (iedge = 0; iedge < mesh.nedges; iedge++) {
        // Get the Nodes
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;

        edge[iedge].unx  = +(node[n2].y - node[n1].y);
        edge[iedge].uny  = -(node[n2].x - node[n1].x);
        edge[iedge].mag = sqrt(edge[iedge].unx*edge[iedge].unx + edge[iedge].uny*edge[iedge].uny);
        edge[iedge].unx /= edge[iedge].mag;
        edge[iedge].uny /= edge[iedge].mag;
    }
}

