/*
 * File:   TriangleWinslowSmoother.cpp
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
    Create_Connectivity();

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
void TriangleWinslowSmoother::SimCenterMeshReader(const char* FileName) {
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
#ifdef DEBUG
    if ((X == NULL) || (Y == NULL))
        error("SimCenterMeshReader: %s\n", "Error Allocating Memory");
#endif
    
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
    if (NTri > 0) {
        Tri = new int[NTri][3];
#ifdef DEBUG
    if ((Tri == NULL))
        error("SimCenterMeshReader: %s\n", "Error Allocating Memory 1");
#endif
    }
    
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
void TriangleWinslowSmoother::SimCenterMeshWriter(const char* FileName) {
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
void TriangleWinslowSmoother::GnuplotWriter(const char* FileName) {
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
void TriangleWinslowSmoother::Create_Interior_Boundary_Tag() {
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
void TriangleWinslowSmoother::Create_Connectivity() {
    Create_Node2Cell_Connectivity();
    Create_Cell2Cell_Connectivity();
    Create_Node2Node_Connectivity();
}

// *****************************************************************************
// *****************************************************************************
void TriangleWinslowSmoother::Create_Node2Cell_Connectivity() {
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
void TriangleWinslowSmoother::Create_Cell2Cell_Connectivity() {
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
void TriangleWinslowSmoother::Create_Node2Node_Connectivity() {
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
        rms /= nrms;
        rms = sqrt(rms);
    }
    return rms;
}

