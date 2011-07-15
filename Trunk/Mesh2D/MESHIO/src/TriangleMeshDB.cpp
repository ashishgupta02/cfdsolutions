/*******************************************************************************
 * File:        TriangleMeshDB.cpp
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
#include "TriangleMeshDB.h"

// *****************************************************************************
// *****************************************************************************
TriangleMeshDB::TriangleMeshDB() {
    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::Init() {
    NNode             = 0;
    NBlock            = 0;
    NTri              = 0;
    NConstant         = 0;
    NVariable         = 0;
    NBoundary         = 0;
    NBoundarySegments = NULL;
    BoundarySegments  = NULL;
    X                 = NULL;
    Y                 = NULL;
    Tri               = NULL;
    Cell2Cell         = NULL;
    Node2Cell         = NULL;
    Node2Node         = NULL;
    IBTag             = NULL;
    Variable          = NULL;
    VariableName      = NULL;
    Constant          = NULL;
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::Delete_Connectivity() {
    int i;
    
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

    Node2Cell = NULL;
    Cell2Cell = NULL;
    Node2Node = NULL;
}

// *****************************************************************************
// *****************************************************************************
TriangleMeshDB::~TriangleMeshDB() {
    int i;
    if (Tri != NULL)
        delete[] Tri;

    if (BoundarySegments != NULL)
        delete[] BoundarySegments;

    if (NBoundarySegments != NULL)
        delete[] NBoundarySegments;

    if (X != NULL)
        delete[] X;

    if (Y != NULL)
        delete[] Y;

    Delete_Connectivity();

    if (Constant != NULL)
        delete []Constant;
    
    if (Variable != NULL) {
        for (i = 0; i < NVariable; i++) {
            if (Variable[i] != NULL)
                delete Variable[i];
        }
        delete Variable;
    }
    
    if (VariableName != NULL) {
        for (i = 0; i < NVariable; i++) {
            if (VariableName[i] != NULL)
                delete VariableName[i];
        }
        delete VariableName;
    }
    
    if (IBTag != NULL)
        delete []IBTag;
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::SLK_MeshReader(const char* FileName) {
    int i, b, NQuad = 0;
    int idum;
    const int bfsize = 132;
    char buff[bfsize];
    char *cdum;
    FILE *fp;

    if ((fp = fopen(FileName, "r")) == 0)
        error("SLK_MeshReader: %s %s\n", "Couldn't Open Mesh File:", FileName);

    info("Reading Mesh File: %s", FileName);

    // Read number of nodes
    cdum = fgets(buff, bfsize, fp);
    cdum = fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NNode);
    info("Number of points = %d", NNode);

    // Allocate Memory for coordinates
    X = new double[NNode];
    Y = new double[NNode];
#ifdef DEBUG
    if ((X == NULL) || (Y == NULL))
        error("SLK_MeshReader: %s\n", "Error Allocating Memory");
#endif

    // Read in coordinates
    for (i = 0; i < NNode; i++) {
        cdum = fgets(buff, bfsize, fp);
        sscanf(buff, "%lg %lg", &(X[i]), &(Y[i]));
    }

    // Read in Number of Blocks
    cdum = fgets(buff, bfsize, fp);
    cdum = fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NBlock);
    info("Number of blocks = %d", NBlock);

    // Read in Number of Triangles
    cdum = fgets(buff, bfsize, fp);
    cdum = fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NTri);
    info("Number of Triangles = %d", NTri);

    // Allocate Memory for Triangles
    if (NTri > 0) {
        Tri = new int[NTri][3];
#ifdef DEBUG
    if ((Tri == NULL))
        error("SLK_MeshReader: %s\n", "Error Allocating Memory 1");
#endif
    }

    // Read in Triangles
    for (i = 0; i < NTri; i++) {
        cdum = fgets(buff, bfsize, fp);
        sscanf(buff, "%d %d %d", &(Tri[i][0]), &(Tri[i][1]), &Tri[i][2]);
        Tri[i][0]--;
        Tri[i][1]--;
        Tri[i][2]--;
    }

    // Read in Number of Quadrilaterals
    cdum = fgets(buff, bfsize, fp);
    cdum = fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NQuad);
    info("Number of Quadrilaterals = %d", NQuad);

    // Allocate Memory for Quadrilaterals
    if (NQuad > 0)
        error("SLK_MeshReader: %s %s\n", "Input file contains Quads - Use Hybrid Reader", FileName);

    // Read in Number of Boundary
    cdum = fgets(buff, bfsize, fp);
    cdum = fgets(buff, bfsize, fp);
    sscanf(buff, "%d", &NBoundary);
    info("Number of boundaries = %d", NBoundary);

    // Allocate for Number of Boundary Segments per Boundary
    NBoundarySegments = new int[NBoundary];
    BoundarySegments = (int***) malloc(NBoundary * sizeof (int**));
    for (b = 0; b < NBoundary; b++) {
        cdum = fgets(buff, bfsize, fp);
        cdum = fgets(buff, bfsize, fp);
        sscanf(buff, "%d", &NBoundarySegments[b]);
        info("Boundary %d has %d segments", b, NBoundarySegments[b]);

        BoundarySegments[b] = (int**) malloc(NBoundarySegments[b] * sizeof (int*));
        for (i = 0; i < NBoundarySegments[b]; i++) {
            BoundarySegments[b][i] = (int*) malloc(2 * sizeof (int));
            cdum = fgets(buff, bfsize, fp);
            sscanf(buff, "%d %d", &BoundarySegments[b][i][0], &BoundarySegments[b][i][1]);
            // Decrement for C Indexing Convention
            BoundarySegments[b][i][0]--;
            BoundarySegments[b][i][1]--;
        }
    }

    // Read in Number of Constant
    if (!feof(fp))
        cdum = fgets(buff, bfsize, fp);
    else return;
    if (!feof(fp))
        cdum = fgets(buff, bfsize, fp);
    else return;
    sscanf(buff, "%d", &NConstant);
    info("Number of Contants = %d", NConstant);

    // Allocate for Number of Constants
    if (NConstant > 0) {
        Constant = new double[NConstant];
        for (i = 0; i < NConstant; i++) {
            cdum = fgets(buff, bfsize, fp);
            sscanf(buff, "%lg", &Constant[i]);
        }
    }
    
    // Read in Number of Variables
    if (!feof(fp))
        cdum = fgets(buff, bfsize, fp);
    else return;
    if (!feof(fp))
        cdum = fgets(buff, bfsize, fp);
    if (!feof(fp))
    sscanf(buff, "%d", &NVariable);
    info("Number of Variables = %d", NVariable);

    if (NVariable > 0) {
        VariableName = new char*[NVariable];
#ifdef DEBUG
        if (VariableName == NULL)
            error("SLK_MeshReader: %s\n", "Error Allocating Memory 2");
#endif
        for (i = 0; i < NVariable; i++) {
            cdum = fgets(buff, bfsize, fp);
            b = strlen(buff);
            VariableName[i] = NULL;
            VariableName[i] = new char[b-1];
#ifdef DEBUG
            if (VariableName[i] == NULL)
                error("SLK_MeshReader: %s\n", "Error Allocating Memory 3");
#endif
            for (int j = 0; j < b-1; j++)
                VariableName[i][j] = buff[j];
        }
        
        // Allocate for Number of Variables
        Variable = new double*[NVariable];
#ifdef DEBUG
        if (Variable == NULL)
            error("SLK_MeshReader: %s\n", "Error Allocating Memory 4");
#endif
        for (i = 0; i < NVariable; i++) {
            Variable[i] = NULL;
            Variable[i] = new double[NNode];
#ifdef DEBUG
            if (Variable[i] == NULL)
                error("SLK_MeshReader: %s\n", "Error Allocating Memory 5");
#endif
        }
        for (i = 0; i < NNode; i++) {
            for (b = 0; b < NVariable; b++)
                idum = fscanf(fp, "%lg", &(Variable[b][i]));
        }
    }
    
    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::SLK_MeshWriter(const char* FileName) {
    int i, j, NQuad = 0;
    FILE *fp;

    // Open file for write
    if ((fp = fopen(FileName, "w")) == 0)
        error("SLK_MeshWriter: %s %s\n", "Couldn't Open Optimized Mesh File:", FileName);

    info("Writing Optimized Mesh File: %s", FileName);

    // Write out nodes
    fprintf(fp, "# Number of grid points\n");
    fprintf(fp, "%d\n", NNode);
    for (i = 0; i < NNode; i++)
        fprintf(fp, "%22.15e %22.15e\n", X[i], Y[i]);

    fprintf(fp, "# Number of blocks\n");
    fprintf(fp, "1\n");

    fprintf(fp, "# Number of triangular elements\n");
    fprintf(fp, "%d\n", NTri);
    for (i = 0; i < NTri; i++)
        fprintf(fp, "%d %d %d\n", Tri[i][0] + 1, Tri[i][1] + 1, Tri[i][2] + 1);

    fprintf(fp, "# Number of quadrilateral elements\n");
    fprintf(fp, "%d\n", NQuad);
    
    fprintf(fp, "# Number of boundaries\n");
    fprintf(fp, "%d\n", NBoundary);

    for (i = 0; i < NBoundary; i++) {
        fprintf(fp, "# Number of edges for boundary %d\n", i + 1);
        fprintf(fp, "%d\n", NBoundarySegments[i]);

        for (j = 0; j < NBoundarySegments[i]; j++)
            fprintf(fp, "%d %d\n", BoundarySegments[i][j][0] + 1, BoundarySegments[i][j][1] + 1);
    }

    // Write in Number of Constant
    fprintf(fp, "# Number of Constants\n");
    fprintf(fp, "%d\n", NConstant);
    if (NConstant > 0) {
        for (i = 0; i < NConstant; i++)
            fprintf(fp, "%lg\n", Constant[i]);
    }
    
    // Write in Number of Variables
    fprintf(fp, "# Number of Variables\n");
    fprintf(fp, "%d\n", NVariable);
    if (NVariable > 0) {
        for (i = 0; i < NVariable; i++)
            fprintf(fp, "%s\n", VariableName[i]);
        for (i = 0; i < NNode; i++) {
            for (j = 0; j < NVariable; j++)
                fprintf(fp, "%22.15e ", Variable[j][i]);
            fprintf(fp, "\n");
        }
    }
    
    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::SLK_GnuplotWriter(const char* FileName) {
    int i, n0, n1, n2;
    FILE *fp;

    // Open file for write
    if ((fp = fopen(FileName, "w")) == 0)
        error("SLK_GnuplotWriter: %s %s\n", "Couldn't Open Gnuplot File:", FileName);

    info("Writing Gnuplot File: %s", FileName);

    for (i = 0; i < NTri; i++) {
        n0 = Tri[i][0];
        n1 = Tri[i][1];
        n2 = Tri[i][2];
        fprintf(fp, "%22.15e %22.15e 0.0\n", X[n0], Y[n0]);
        fprintf(fp, "%22.15e %22.15e 0.0\n", X[n1], Y[n1]);
        fprintf(fp, "%22.15e %22.15e 0.0\n", X[n2], Y[n2]);
        fprintf(fp, "%22.15e %22.15e 0.0\n\n", X[n0], Y[n0]);
    }

    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::WKA_MeshReader(const char* FileName) {
    info("Not Implemented Yet");
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::WKA_MeshWriter(const char* FileName) {
    info("Not Implemented Yet");
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::WKA_GnuplotWriter(const char* FileName) {
    info("Not Implemented Yet");
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::Create_Interior_Boundary_Tag() {
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
void TriangleMeshDB::Delete_Interior_Boundary_Tag() {
    if (IBTag != NULL) {
        delete[] IBTag;
        IBTag = NULL;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::Create_Connectivity(int Sort) {
    Create_Node2Cell_Connectivity();
    Create_Cell2Cell_Connectivity();
    Create_Node2Node_Connectivity(Sort);
}

// *****************************************************************************
// *****************************************************************************
void TriangleMeshDB::Create_Node2Cell_Connectivity() {
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
void TriangleMeshDB::Create_Cell2Cell_Connectivity() {
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
void TriangleMeshDB::Create_Node2Node_Connectivity(int Sort) {
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

    // Check if Sorting is requested
    if (Sort == 0) {
        // No Sorting
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
    } else {
        // Perform sorting of connected nodes
        int iOrient, nextcell, tarcell;
        int n0, n1, n2, nc1;
        List *Sorter = new List();

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

            // Sort the Connected Nodes According to Counter Clock wise Orientations
            // This is done to Properly Contruct Virtual Control Volume
            if (Node2Node[iNode]->max > 0) {
                // Get the Right Most Triangle in Triangle Fan around the node
                if (Node2Cell[iNode]->max > 0) {
                    // Get the Starting Triangle
                    cellid = Node2Cell[iNode]->list[0];
                    n0 = n1 = n2 = -1;
                    // Triangle Search Loop Clock Wise
                    // max-1 because we dont want to find ourself where we started
                    for (nc1 = 0; nc1 < Node2Cell[iNode]->max - 1; nc1++) {
                        // Get the Edge => n1 to see right triangle
                        iOrient = -1;
                        for (i = 0; i < 3; i++) {
                            if (iNode == Tri[cellid][i])
                                iOrient = i;
                        }
                        switch (iOrient) {
                            case 0:
                                n0 = Tri[cellid][0];
                                n1 = Tri[cellid][1];
                                n2 = Tri[cellid][2];
                                break;
                            case 1:
                                n0 = Tri[cellid][1];
                                n1 = Tri[cellid][2];
                                n2 = Tri[cellid][0];
                                break;
                            case 2:
                                n0 = Tri[cellid][2];
                                n1 = Tri[cellid][0];
                                n2 = Tri[cellid][1];
                                break;
                            default:
                                break;
                        }
                        // Get the Triangle Common to Edge n0-n1 except current triangle
                        // For all Node0 connected Triangles
                        tarcell = -1;
                        for (nc = 0; nc < Node2Cell[n0]->max; nc++) {
                            // Get the triangle id of Node0 connected cell from list
                            nextcell = Node2Cell[n0]->list[nc];
                            // Check if Node0 and Node1 have common triangles
                            // Max of two cells can be common for an Edge and
                            // Min of One Cell if Edge is on Boundary
                            if ((nextcell != cellid) && (Node2Cell[n1]->Is_In_List(nextcell))) {
                                tarcell = nextcell;
                                nc = Node2Cell[n0]->max;
                            }
                        }
                        // Check if Already at the Right Mode Triangle in Fan
                        if (tarcell == -1)
                            // Exit Triangle Search
                            nc1 = Node2Cell[iNode]->max;
                        else
                            cellid = tarcell;
                    } // Triangle Search Loop Clock Wise

                    // Reset the Sorter
                    Sorter->Reset(0);
                    // Triangle Search Loop Anti-Clock Wise
                    for (nc1 = 0; nc1 < Node2Cell[iNode]->max; nc1++) {
                        // Get the Orientation
                        iOrient = -1;
                        for (i = 0; i < 3; i++) {
                            if (iNode == Tri[cellid][i])
                                iOrient = i;
                        }
                        switch (iOrient) {
                            case 0:
                                n0 = Tri[cellid][0];
                                n1 = Tri[cellid][1];
                                n2 = Tri[cellid][2];
                                break;
                            case 1:
                                n0 = Tri[cellid][1];
                                n1 = Tri[cellid][2];
                                n2 = Tri[cellid][0];
                                break;
                            case 2:
                                n0 = Tri[cellid][2];
                                n1 = Tri[cellid][0];
                                n2 = Tri[cellid][1];
                                break;
                            default:
                                break;
                        }

                        // Add Node1 to the Sorter List
                        Sorter->Check_List(n1);

                        // Get the Triangle Common to Edge n0-n2 except current triangle
                        // For all Node0 connected Triangles
                        tarcell = -1;
                        for (nc = 0; nc < Node2Cell[n0]->max; nc++) {
                            // Get the triangle id of Node0 connected cell from list
                            nextcell = Node2Cell[n0]->list[nc];
                            // Check if Node0 and Node1 have common triangles
                            // Max of two cells can be common for an Edge and
                            // Min of One Cell if Edge is on Boundary
                            if ((nextcell != cellid) && (Node2Cell[n2]->Is_In_List(nextcell))) {
                                tarcell = nextcell;
                                nc = Node2Cell[n0]->max;
                            }
                        }
                        // Check if Already at the Right Mode Triangle in Fan
                        if (tarcell == -1) {
                            // Exit Triangle Search
                            nc1 = Node2Cell[iNode]->max;
                            Sorter->Check_List(n2);
                        } else {
                            cellid = tarcell;
                        }
                    } // Triangle Search Loop Anti-Clock Wise

                    // Update the sorted List to Node2Node Connectivity
                    if (Sorter->max == Node2Node[iNode]->max) {
                        Node2Node[iNode]->Reset(0);
                        for (i = 0; i < Sorter->max; i++)
                            Node2Node[iNode]->Check_List(Sorter->list[i]);
                    } else {
                        error("Create_Node2Node_Connectivity: %s\n", "Sorting Problem");
                    }
                } // if Condition 2
            } // if Condition 1
        } // iNode Loop
    }
}

