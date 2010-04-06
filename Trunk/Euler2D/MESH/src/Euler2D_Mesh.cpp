/* 
 * File:   Euler2D_Mesh.cpp
 * Author: Ashish Gupta
 * 
 * Created on March 21, 2010, 7:50 PM
 */

#include <stdlib.h>

#include "List.h"
#include "Utils.h"
#include "Euler2D_Mesh.h"

// *****************************************************************************
// *****************************************************************************
Euler2D_Mesh::Euler2D_Mesh() {
    // Initialize
    Init();
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::Init() {
    restart      = 0;
    changeIndex  = 0;
    cell         = NULL;
    edge         = NULL;
    node         = NULL;
    boundaryEdge = NULL;
    boundaryNode = NULL;
    BNTag        = NULL;
}

// *****************************************************************************
// *****************************************************************************
Euler2D_Mesh::~Euler2D_Mesh() {
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
void Euler2D_Mesh::WKA_MeshReader(const char* FileName) {
    int i, icount;
    int node1, node2;
    int idum;
    FILE *inputMesh;
    char dumstring[100];
    int *tag = NULL;

    printf("============================================================\n");
    info("Reading Mesh File %s", FileName);

    // Open Mesh File
    if ((inputMesh = fopen(FileName, "r")) == (FILE *) NULL)
        error("Unable to Open Mesh File %s", FileName);

    // Read mesh sizes
    fscanf(inputMesh, "%d %d %d %d %d %d", &mesh.nnodes, &mesh.nedges, &mesh.ncells, &mesh.nbedges, &idum, &idum);
    info("NNodes = %d NEdges = %d NCells = %d NBEdge = %d", mesh.nnodes, mesh.nedges, mesh.ncells, mesh.nbedges);

    mesh.inside = mesh.ncells - mesh.nbedges;

    // Allocate Memory to store Connectivity
    edge = (EDGE *) calloc(mesh.nedges, sizeof (EDGE));
    cell = (CELL *) calloc(mesh.ncells, sizeof (CELL));
    boundaryEdge = (BOUNDARYEDGE *) calloc(mesh.nbedges, sizeof (BOUNDARYEDGE));

    icount = 0;
    for (i = 0; i < mesh.nedges; i++) {
        fscanf(inputMesh, "%d %d %d %d", &edge[i].node1, &edge[i].node2, &edge[i].cell1, &edge[i].cell2);
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
    fscanf(inputMesh, "\n");
    fgets(dumstring, 100, inputMesh);
    for (i = 0; i < mesh.nbedges; i++) {
        fscanf(inputMesh, "%d %d %lf %lf", &idum, &(boundaryEdge[i].bcType), &(boundaryEdge[i].c1), &(boundaryEdge[i].c2));
        if (changeIndex == 1)
            boundaryEdge[i].edgeNumber = idum - 1;
        else
            boundaryEdge[i].edgeNumber = idum;
    }

    /* Now read the mesh coordinates */
    fscanf(inputMesh, "\n");
    fgets(dumstring, 100, inputMesh);
    node = (NODE *) calloc(mesh.nnodes, sizeof (NODE));
    for (i = 0; i < mesh.nnodes; i++)
        fscanf(inputMesh, "%le %le", &node[i].x, &node[i].y);

    /* Fill boundaryNode array (used for Dirchlet boundary condition) */
    tag = (int *) calloc(mesh.nnodes, sizeof (int));
    for (i = 0; i < mesh.nnodes; i++) tag[i] = 0;

    /* Fill boundaryNode array */
    /* First find out how large to make the array */
    mesh.nbnodes = 0;
    for (i = 0; i < mesh.nbedges; i++) {
        node1 = edge[boundaryEdge[i].edgeNumber].node1;
        node2 = edge[boundaryEdge[i].edgeNumber].node2;
        if ((boundaryEdge[i].bcType >= 2000) && (boundaryEdge[i].bcType <= 2999)) {
            if (tag[node1] != 1) {
                mesh.nbnodes += 1;
                tag[node1] = 1;
            }
            if (tag[node2] != 1) {
                mesh.nbnodes += 1;
                tag[node2] = 1;
            }
        }
    }
    info("NBNodes = %d", mesh.nbnodes);

    boundaryNode = (BOUNDARYNODE *) calloc(mesh.nbnodes, sizeof (BOUNDARYNODE));
    mesh.nbnodes = 0;
    for (i = 0; i < mesh.nnodes; i++) tag[i] = 0;
    for (i = 0; i < mesh.nbedges; i++) {
        node1 = edge[boundaryEdge[i].edgeNumber].node1;
        node2 = edge[boundaryEdge[i].edgeNumber].node2;
        if ((boundaryEdge[i].bcType >= 2000) && (boundaryEdge[i].bcType <= 2999)) {
            if (tag[node1] != 1) {
                boundaryNode[mesh.nbnodes].nodeNumber = node1;
                if (boundaryEdge[i].bcType == 2000) {
                    boundaryNode[mesh.nbnodes].constant = node[node1].x;
                } else {
                    boundaryNode[mesh.nbnodes].constant = boundaryEdge[i].c1;
                }
                mesh.nbnodes += 1;
                tag[node1] = 1;
            }
            if (tag[node2] != 1) {
                boundaryNode[mesh.nbnodes].nodeNumber = node2;
                if (boundaryEdge[i].bcType == 2000) {
                    boundaryNode[mesh.nbnodes].constant = node[node2].x;
                } else {
                    boundaryNode[mesh.nbnodes].constant = boundaryEdge[i].c1;
                }
                mesh.nbnodes += 1;
                tag[node2] = 1;
            }
        }
    }

    // Free Memory
    if (tag != NULL)
        free(tag);

    printf("============================================================\n");
}

// *****************************************************************************
/*                                                                         */
/* Extract cell-to-node pointers from edge pointers                        */
/*                                                                         */
// *****************************************************************************
void Euler2D_Mesh::WKA_ExtractCells(void) {
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
    info("Cell Area: MAX = %lf MIN = %lf", areamax, areamin);
    info("Cell ID Min Area: %d Area = %le", imin, cell[imin].area);
    printf("============================================================\n");
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::Tag_Boundary_Nodes() {
    int i, n1, n2, bctype;
    
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
        // BCType 2: 2000-2999
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
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::Read_RestartFile(const char* FileName) {
    int iNode, var;
    FILE *fp;

    printf("============================================================\n");
    info("Reading Restart File %s", FileName);
    
    // Open Mesh File
    if ((fp = fopen(FileName, "rb")) == (FILE *) NULL)
        error("Read_RestartFile: Unable to Open Mesh File %s", FileName);
    
    var = 0;
    fread(&var, sizeof(int), 1, fp);
    if (var != mesh.nnodes)
        error("Read_RestartFile: Mismatched in Restart and Mesh File %s", FileName);
    var = 0;
    fread(&var, sizeof(int), 1, fp);
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        fread(&node[iNode].Q[0], sizeof (double), 1, fp);
        fread(&node[iNode].Q[1], sizeof (double), 1, fp);
        fread(&node[iNode].Q[2], sizeof (double), 1, fp);
        fread(&node[iNode].Q[3], sizeof (double), 1, fp);
    }

    fclose(fp);
    printf("============================================================\n");
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::Write_RestartFile(const char* FileName) {
    int iNode, var;
    FILE *fp;

    printf("============================================================\n");
    info("Writing Restart File %s", FileName);

    // Open Mesh File
    if ((fp = fopen(FileName, "wb")) == (FILE *) NULL)
        error("Write_RestartFile: Unable to Open Mesh File %s", FileName);
    
    fwrite(&mesh.nnodes, sizeof(int), 1, fp);
    var = 4;
    fwrite(&var, sizeof(int), 1, fp);
    for (iNode = 0; iNode < mesh.nnodes; iNode++) {
        fwrite(&node[iNode].Q[0], sizeof (double), 1, fp);
        fwrite(&node[iNode].Q[1], sizeof (double), 1, fp);
        fwrite(&node[iNode].Q[2], sizeof (double), 1, fp);
        fwrite(&node[iNode].Q[3], sizeof (double), 1, fp);
    }
    
    fclose(fp);
    printf("============================================================\n");
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::SLK_MeshWriter(const char* FileName) {
    int i, j, n1, n2, NQuad = 0;
    FILE *fp;

    // Open file for write
    if ((fp = fopen(FileName, "w")) == 0)
        error("SLK_MeshWriter: %s %s\n", "Couldn't Open Optimized Mesh File:", FileName);

    info("Writing Optimized Mesh File: %s", FileName);

    // Write out nodes
    fprintf(fp, "# Number of grid points\n");
    fprintf(fp, "%d\n", mesh.nnodes);
    for (i = 0; i < mesh.nnodes; i++)
        fprintf(fp, "%22.15e %22.15e\n", node[i].x, node[i].y);

    fprintf(fp, "# Number of blocks\n");
    fprintf(fp, "1\n");

    fprintf(fp, "# Number of triangular elements\n");
    fprintf(fp, "%d\n", mesh.inside);
    for (i = 0; i < mesh.inside; i++)
        fprintf(fp, "%d %d %d\n", cell[i].node1+1, cell[i].node2+1, cell[i].node3+1);

    fprintf(fp, "# Number of quadrilateral elements\n");
    fprintf(fp, "%d\n", NQuad);

    // Get the Number of Boundaries
    List bndy;
    for (i = 0; i < mesh.nbedges; i++)
        bndy.Check_List(boundaryEdge[i].bcType);

    fprintf(fp, "# Number of boundaries\n");
    fprintf(fp, "%d\n", bndy.max);

    if (bndy.max > 0) {
        int *nbs = NULL;
        nbs = (int*) malloc(bndy.max*sizeof(int));
        for (i = 0; i < bndy.max; i++)
            nbs[i] = 0;

        for (i = 0; i < mesh.nbedges; i++)
            nbs[bndy.Index(boundaryEdge[i].bcType)]++;
        
        for (i = 0; i < bndy.max; i++) {
            fprintf(fp, "# Number of edges for boundary %d\n", i + 1);
            fprintf(fp, "%d\n", nbs[i]);
            for (j = 0; j < mesh.nbedges; j++) {
                if (bndy.list[i] == boundaryEdge[j].bcType) {
                    n1 = edge[boundaryEdge[j].edgeNumber].node1;
                    n2 = edge[boundaryEdge[j].edgeNumber].node2;
                    fprintf(fp, "%d %d\n", n1+1, n2+1);
                }
            }
        }
    }

    // Constants
    NQuad = 4;
    fprintf(fp, "# Number of Constants\n");
    fprintf(fp, "%d\n", NQuad);
    fprintf(fp, "0.0\n");
    fprintf(fp, "0.0\n");
    fprintf(fp, "0.0\n");
    fprintf(fp, "1.4\n");
    // Variables
    fprintf(fp, "# Number of Variables\n");
    NQuad = 4;
    fprintf(fp, "%d\n", NQuad);
    fprintf(fp, "density\n");
    fprintf(fp, "x-momentum\n");
    fprintf(fp, "y-momentum\n");
    fprintf(fp, "energy\n");
    for (i = 0; i < mesh.nnodes; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e\n",
                node[i].Q[0], node[i].Q[1], node[i].Q[2], node[i].Q[3]);
    
    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::Write_Solution_GnuplotFile(const char* FileName) {
    int i;
    FILE *fp;

    // Open file for write
    if ((fp = fopen(FileName, "w")) == 0)
        error("Write_Solution_GnuplotFile: %s %s\n", "Couldn't Open Optimized Mesh File:", FileName);

    info("Writing Solution Gnuplot File: %s", FileName);
    
    // Write out nodes and Qs
    for (i = 0; i < mesh.nnodes; i++)
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
                node[i].x, node[i].y, node[i].Q[0], node[i].Q[1], node[i].Q[2], node[i].Q[3]);
    
    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::Write_TecplotFile(const char* FileName) {
    int i;
    double Pressure, Gamma;
    FILE *fp;

    // Open file for write
    if ((fp = fopen(FileName, "w")) == 0)
        error("Write_TecplotFile: %s %s\n", "Couldn't Open File:", FileName);

    info("Writing Tecplot File: %s", FileName);
    
    // Write out nodes
    fprintf(fp, "title = \"Triangular Mesh Field Data\"\n");
    fprintf(fp, "variables = \"x\", \"y\", \"Density\", \"U Velocity\", \"V Velocity\", \"Internal Energy\", \"Pressure\"\n");
    fprintf(fp, "zone n=%d, e=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n", mesh.nnodes, mesh.inside);

    // Write Grid Point and Solutions
    fprintf(fp, "\n");
    Gamma = 1.4;
    for (i = 0; i < mesh.nnodes; i++) {
        Pressure = (Gamma - 1.0)*(node[i].Q[3] - 0.5*node[i].Q[0]*((node[i].Q[1]/node[i].Q[0])*(node[i].Q[1]/node[i].Q[0])
                + (node[i].Q[2]/node[i].Q[0])*(node[i].Q[2]/node[i].Q[0])));
        fprintf(fp, "%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",
                node[i].x, node[i].y, node[i].Q[0], (node[i].Q[1]/node[i].Q[0]), (node[i].Q[2]/node[i].Q[0]), node[i].Q[3], Pressure);
    }

    // Write Triangle Connectivity
    fprintf(fp, "\n");
    for (i = 0; i < mesh.inside; i++)
        fprintf(fp, "%d %d %d\n", cell[i].node1+1, cell[i].node2+1, cell[i].node3+1);

    // Close File
    fclose(fp);
}

// *****************************************************************************
// *****************************************************************************
void Euler2D_Mesh::Compute_Geometric_Properties() {
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
        cell[icell].nx12   = 0.0;
        cell[icell].nx23   = 0.0;
        cell[icell].nx31   = 0.0;
        cell[icell].ny12   = 0.0;
        cell[icell].ny23   = 0.0;
        cell[icell].ny31   = 0.0;
        cell[icell].magc12 = 0.0;
        cell[icell].magc23 = 0.0;
        cell[icell].magc31 = 0.0;
        cell[icell].nxc12  = 0.0;
        cell[icell].nxc23  = 0.0;
        cell[icell].nxc31  = 0.0;
        cell[icell].nyc12  = 0.0;
        cell[icell].nyc23  = 0.0;
        cell[icell].nyc31  = 0.0;
    }
    
    for (inode = 0; inode < mesh.nnodes; inode++)
        node[inode].area = 0.0;

    for (iedge = 0; iedge < mesh.nedges; iedge++) {
        edge[iedge].nx  = 0.0;
        edge[iedge].ny  = 0.0;
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
        cell[icell].nxc12 = +(cell[icell].yc - 0.5*(node[n1].y + node[n2].y));
        cell[icell].nyc12 = -(cell[icell].xc - 0.5*(node[n1].x + node[n2].x));
        cell[icell].magc12 = sqrt(cell[icell].nxc12*cell[icell].nxc12 + cell[icell].nyc12*cell[icell].nyc12);
        cell[icell].nxc12 /= cell[icell].magc12;
        cell[icell].nyc12 /= cell[icell].magc12;

        // Edge 2-3
        cell[icell].nxc23 = +(cell[icell].yc - 0.5*(node[n2].y + node[n3].y));
        cell[icell].nyc23 = -(cell[icell].xc - 0.5*(node[n2].x + node[n3].x));
        cell[icell].magc23 = sqrt(cell[icell].nxc23*cell[icell].nxc23 + cell[icell].nyc23*cell[icell].nyc23);
        cell[icell].nxc23 /= cell[icell].magc23;
        cell[icell].nyc23 /= cell[icell].magc23;

        // Edge 3-1
        cell[icell].nxc31 = +(cell[icell].yc - 0.5*(node[n3].y + node[n1].y));
        cell[icell].nyc31 = -(cell[icell].xc - 0.5*(node[n3].x + node[n1].x));
        cell[icell].magc31 = sqrt(cell[icell].nxc31*cell[icell].nxc31 + cell[icell].nyc31*cell[icell].nyc31);
        cell[icell].nxc31 /= cell[icell].magc31;
        cell[icell].nyc31 /= cell[icell].magc31;

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
        cell[icell].nx12 = +(node[n2].y - node[n1].y);
        cell[icell].ny12 = -(node[n2].x - node[n1].x);
        cell[icell].mag12 = sqrt(cell[icell].nx12*cell[icell].nx12 + cell[icell].ny12*cell[icell].ny12);
        cell[icell].nx12 /= cell[icell].mag12;
        cell[icell].ny12 /= cell[icell].mag12;

        // Edge 2-3
        cell[icell].nx23 = +(node[n3].y - node[n2].y);
        cell[icell].ny23 = -(node[n3].x - node[n2].x);
        cell[icell].mag23 = sqrt(cell[icell].nx23*cell[icell].nx23 + cell[icell].ny23*cell[icell].ny23);
        cell[icell].nx23 /= cell[icell].mag23;
        cell[icell].ny23 /= cell[icell].mag23;

        // Edge 3-1
        cell[icell].nx31 = +(node[n1].y - node[n3].y);
        cell[icell].ny31 = -(node[n1].x - node[n3].x);
        cell[icell].mag31 = sqrt(cell[icell].nx31*cell[icell].nx31 + cell[icell].ny31*cell[icell].ny31);
        cell[icell].nx31 /= cell[icell].mag31;
        cell[icell].ny31 /= cell[icell].mag31;
    }

    for (iedge = 0; iedge < mesh.nedges; iedge++) {
        // Get the Nodes
        n1 = edge[iedge].node1;
        n2 = edge[iedge].node2;
        
        edge[iedge].nx  = +(node[n2].y - node[n1].y);
        edge[iedge].ny  = -(node[n2].x - node[n1].x);
        edge[iedge].mag = sqrt(edge[iedge].nx*edge[iedge].nx + edge[iedge].ny*edge[iedge].ny);
        edge[iedge].nx /= edge[iedge].mag;
        edge[iedge].ny /= edge[iedge].mag;
    }
}

