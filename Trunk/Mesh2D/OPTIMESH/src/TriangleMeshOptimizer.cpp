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

#ifdef HAVE_VERDICT
#include "verdict.h"
#endif
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
    State             = 0;
    NNode             = 0;
    NBlock            = 0;
    NTri              = 0;
    NQuad             = 0;
    NBoundary         = 0;
    Tri               = NULL;
    Quad              = NULL;
    NBoundarySegments = NULL;
    BoundarySegments  = NULL;
    X                 = NULL;
    Y                 = NULL;
    Cell2Cell         = NULL;
    Node2Cell         = NULL;
    Node2Node         = NULL;
    IBTag             = NULL;
    RState            = 0;
    RAngle            = 0.0;
    CGx               = 0.0;
    CGy               = 0.0;
}

// *****************************************************************************
// *****************************************************************************
TriangleMeshOptimizer::~TriangleMeshOptimizer() {
    int i;

    if (Tri != NULL)
        delete[] Tri;
    
    if (Quad != NULL)
        delete[] Quad;

    if (BoundarySegments != NULL)
        delete[] BoundarySegments;

    if (NBoundarySegments != NULL)
        delete[] NBoundarySegments;

    if (X != NULL)
        delete[] X;

    if (Y != NULL)
        delete[] Y;

    if (Node2Cell != NULL) {
        for (i = 0; i < NNode; i++) {
            if (Node2Cell[i] != NULL)
                delete Node2Cell[i];
        }
        delete Node2Cell;
    }
    
    if (Cell2Cell != NULL) {
        for (i = 0; i < NTri; i++) {
            if (Cell2Cell[i] != NULL)
                delete Cell2Cell[i];
        }
        delete Cell2Cell;
    }

    if (Node2Node != NULL) {
        for (i = 0; i < NNode; i++) {
            if (Node2Node[i] != NULL)
                delete Node2Node[i];
        }
        delete Node2Node;
    }

    if (IBTag != NULL)
        delete []IBTag;
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
    Create_Connectivity();

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
void TriangleMeshOptimizer::SimCenterMeshReader(const char* FileName) {
    int i, b;
    const int bfsize = 132;
    char buff[bfsize];
    FILE *fp;
    
    if ((fp = fopen(FileName, "r")) == 0)
        error("SimCenterMeshReader: %s %s\n", "Couldn't Open Mesh File:", FileName);

    info("Reading Mesh File: %s", FileName);
    
    // Read number of nodes
    fgets(buff, bfsize, fp);
    fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NNode);
    info("Number of points = %d", NNode);
    
    // Allocate Memory for coordinates
    X = new double[NNode];
    Y = new double[NNode];
    
    // Read in coordinates
    for (i = 0; i < NNode; i++) {
        fgets(buff, bfsize, fp);
        sscanf(buff, "%lg %lg", &(X[i]), &(Y[i]));
    }

    // Read in Number of Blocks
    fgets(buff, bfsize, fp);
    fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NBlock);
    info("Number of blocks = %d", NBlock);

    // Read in Number of Triangles
    fgets(buff, bfsize, fp);
    fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NTri);
    info("Number of Triangles = %d", NTri);

    // Allocate Memory for Triangles
    if (NTri > 0)
        Tri = new int[NTri][3];
    
    // Read in Triangles
    for (i = 0; i < NTri; i++) {
        fgets(buff, bfsize, fp);
        sscanf(buff, "%d %d %d", &(Tri[i][0]), &(Tri[i][1]), &Tri[i][2]);
        Tri[i][0]--;
        Tri[i][1]--;
        Tri[i][2]--;
    }
    
    // Read in Number of Quadrilaterals
    fgets(buff, bfsize, fp);
    fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NQuad);
    info("Number of Quadrilaterals = %d", NQuad);
    
    // Allocate Memory for Quadrilaterals
    if (NQuad > 0)
        Quad = new int[NQuad][4];

    // Read in Quadrilaterals
    for (i = 0; i < NQuad; i++) {
        fgets(buff, bfsize, fp);
        sscanf(buff, "%d %d %d %d", &(Quad[i][0]), &(Quad[i][1]), &Quad[i][2], &Quad[i][3]);
        Quad[i][0]--;
        Quad[i][1]--;
        Quad[i][2]--;
        Quad[i][3]--;
    }
    
    // Read in Number of Boundary
    fgets(buff, bfsize, fp);
    fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NBoundary);
    info("Number of boundaries = %d", NBoundary);
    
    // Allocate for Number of Boundary Segments per Boundary
    NBoundarySegments = new int[NBoundary];
    BoundarySegments = (int***) malloc(NBoundary * sizeof (int**));
    for (b = 0; b < NBoundary; b++) {
        fgets(buff, bfsize, fp);
        fgets(buff, bfsize, fp);
        sscanf(buff, "%d", &NBoundarySegments[b]);
        info("Boundary %d has %d segments", b, NBoundarySegments[b]);

        BoundarySegments[b] = (int**) malloc(NBoundarySegments[b] * sizeof (int*));
        for (i = 0; i < NBoundarySegments[b]; i++) {
            BoundarySegments[b][i] = (int*) malloc(2 * sizeof (int));
            fgets(buff, bfsize, fp);
            sscanf(buff, "%d %d", &BoundarySegments[b][i][0], &BoundarySegments[b][i][1]);
            // Decrement for C Indexing Convention
            BoundarySegments[b][i][0]--;
            BoundarySegments[b][i][1]--;
        }
    }

    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::SimCenterMeshWriter(const char* FileName) {
    int i, j;
    FILE *fp;
    
    // Open file for write
    if ((fp = fopen(FileName, "w")) == 0)
        error("SimCenterMeshWriter: %s %s\n", "Couldn't Open Optimized Mesh File:", FileName);

    info("Writing Optimized Mesh File: %s", FileName);
    
    // Write out nodes
    fprintf(fp, "#Number of grid points\n");
    fprintf(fp, "%d\n", NNode);
    for (i = 0; i < NNode; i++)
        fprintf(fp, "%19.10e %19.10e\n", X[i], Y[i]);
    
    fprintf(fp, "#Number of blocks\n");
    fprintf(fp, "1\n");

    fprintf(fp, "#Number of triangular elements\n");
    fprintf(fp, "%d\n", NTri);
    for (i = 0; i < NTri; i++)
        fprintf(fp, "%d %d %d\n", Tri[i][0] + 1, Tri[i][1] + 1, Tri[i][2] + 1);
    
    fprintf(fp, "#Number of quadrilateral elements\n");
    fprintf(fp, "%d\n", NQuad);
    for (i = 0; i < NQuad; i++)
        fprintf(fp, "%d %d %d %d\n", Quad[i][0] + 1, Quad[i][1] + 1, Quad[i][2] + 1, Quad[i][3] + 1);
    
    fprintf(fp, "#Number of boundaries\n");
    fprintf(fp, "%d\n", NBoundary);

    for (i = 0; i < NBoundary; i++) {
        fprintf(fp, "#Number of edges for boundary %d\n", i + 1);
        fprintf(fp, "%d\n", NBoundarySegments[i]);

        for (j = 0; j < NBoundarySegments[i]; j++)
            fprintf(fp, "%d %d\n", BoundarySegments[i][j][0] + 1, BoundarySegments[i][j][1] + 1);
    }
    
    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::GnuplotWriter(const char* FileName) {
    int i, n0, n1, n2;
    FILE *fp;
    
    // Open file for write
    if ((fp = fopen(FileName, "w")) == 0)
        error("GnuplotWriter: %s %s\n", "Couldn't Open Gnuplot File:", FileName);

    info("Writing Gnuplot File: %s", FileName);

    for (i = 0; i < NTri; i++) {
        n0 = Tri[i][0];
        n1 = Tri[i][1];
        n2 = Tri[i][2];
        fprintf(fp, "%19.10e %19.10e 0.0\n", X[n0], Y[n0]);
        fprintf(fp, "%19.10e %19.10e 0.0\n", X[n1], Y[n1]);
        fprintf(fp, "%19.10e %19.10e 0.0\n", X[n2], Y[n2]);
        fprintf(fp, "%19.10e %19.10e 0.0\n\n", X[n0], Y[n0]);
    }

    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Create_Interior_Boundary_Tag() {
    int iNode, n0, n1;

    // Allocate the memory
    IBTag = new int[NNode];
#ifdef DEBUG
    if (IBTag == NULL)
        error("Create_Interior_Boundary_Tag: %s", "Error Allocating Memory");
#endif

    // Tag all Nodes to Interior = -1
    for (iNode = 0; iNode < NNode; iNode++)
        IBTag[iNode] = -1;

    // For all Boundaries
    for (int ib = 0; ib < NBoundary; ib++) {
        for (int ibs = 0; ibs < NBoundarySegments[ib]; ibs++) {
            n0 = BoundarySegments[ib][ibs][0];
            n1 = BoundarySegments[ib][ibs][1];
            IBTag[n0] = IBTag[n1] = ib;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Create_Connectivity() {
    Create_Node2Cell_Connectivity();
    Create_Cell2Cell_Connectivity();
    Create_Node2Node_Connectivity();
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Create_Node2Cell_Connectivity() {
    int i, icell;

    // Allocate the memory
    Node2Cell = new List*[NNode];
#ifdef DEBUG
    if (Node2Cell == NULL)
        error("Create_Node2Cell_Connectivity: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < NNode; i++) {
        Node2Cell[i] = new List();
#ifdef DEBUG
        if (Node2Cell[i] == NULL)
            error("Create_Node2Cell_Connectivity: %s\n", "Error Allocating Memory 2");
#endif
    }
    
    // Now Add Triangles connected to the Node
    for (icell = 0; icell < NTri; icell++) {
        for (i = 0; i < 3; i++) {
#ifdef DEBUG
            if (Tri[icell][i] >= 0)
#endif
                Node2Cell[Tri[icell][i]]->Check_List(icell);
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Create_Cell2Cell_Connectivity() {
    // Allocate memory
    int i, icell, nd, nc;
    int n0, n1, cellid;
    
    // Allocate the memory
    Cell2Cell = new int*[NTri];
#ifdef DEBUG
    if (Cell2Cell == NULL)
        error("Create_Cell2Cell_Connectivity: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < (NTri); i++) {
        Cell2Cell[i] = new int[3];
#ifdef DEBUG
        if (Cell2Cell[i] == NULL)
            error("Create_Cell2Cell_Connectivity: %s\n", "Error Allocating Memory 2");
#endif
    }

    // For Each Triangle Edge Search if they have Common Cell
    // Except the cell which whose Edge is being searched
    for (icell = 0; icell < NTri; icell++) {
        for (nd = 0; nd < 3; nd++) {
            // Initialize the Neighbour Cells to -1
            Cell2Cell[icell][nd] = -1;
            // Get the Edge Nodes of the Triangle
            n0 = Tri[icell][nd];
            n1 = Tri[icell][(nd + 1) % 3];
            // For all Node0 connected Triangles
            for (nc = 0; nc < Node2Cell[n0]->max; nc++) {
                // Get the triangle id of Node0 connected cell from list
                cellid = Node2Cell[n0]->list[nc];
                // Check if Node0 and Node1 have common triangles
                // Max of two cells can be common for an Edge and
                // Min of One Cell if Edge is on Boundary
                if ((icell != cellid) && (Node2Cell[n1]->Is_In_List(cellid))) {
                    Cell2Cell[icell][nd] = cellid;
                }
            }
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshOptimizer::Create_Node2Node_Connectivity() {
    int i, iNode, nc, cellid, nodeid;

    // Allocate the memory
    Node2Node = new List*[NNode];
#ifdef DEBUG
    if (Node2Node == NULL)
        error("Create_Node2Node_Connectivity: %s\n", "Error Allocating Memory 1");
#endif
    for (i = 0; i < NNode; i++) {
        Node2Node[i] = new List();
#ifdef DEBUG
        if (Node2Node[i] == NULL)
            error("Create_Node2Node_Connectivity: %s\n", "Error Allocating Memory 2");
#endif
    }

    // For each node use Node2Cell connectivity to get the connected triangles
    // And using cell2node connectivity get the node list connected to triangle
    for (iNode = 0; iNode < NNode; iNode++) {
        for (nc = 0; nc < Node2Cell[iNode]->max; nc++) {
            cellid = Node2Cell[iNode]->list[nc];
            for (i = 0; i < 3; i++) {
                nodeid = Tri[cellid][i];
                // Check and add the other nodes connected to the triangle
                // Add only if not added already
                if (iNode != nodeid)
                    Node2Node[iNode]->Check_List(nodeid);
            }
        }
    }
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
        else {
#ifdef HAVE_VERDICT
            double coordinates[3][3];
            coordinates[0][0] = X[n0];
            coordinates[0][1] = Y[n0];
            coordinates[0][2] = 0.0;
            coordinates[1][0] = X[n1];
            coordinates[1][1] = Y[n1];
            coordinates[1][2] = 0.0;
            coordinates[2][0] = X[n2];
            coordinates[2][1] = Y[n2];
            coordinates[2][2] = 0.0;
            cost = 1.0 - (1.0/v_tri_condition(3, coordinates));
#else
            cost = 1.0 - (1.0/Compute_Tri_ConditionNumber(n0, n1, n2));
#endif
        }

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

