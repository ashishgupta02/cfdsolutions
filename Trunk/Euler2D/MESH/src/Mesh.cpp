/*
 * File:   Mesh.cpp
 * Author: Ashish Gupta
 *
 * Created on February 28, 2010, 9:13 PM
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Custom header files */
#include "Mesh.h"
#include "Utils.h"

// Global Mesh Variables
char *meshFile;
int changeIndex = 0;
MESH mesh;
CELL *cell;
EDGE *edge;
NODE *node;
BOUNDARYEDGE *boundaryEdge;
BOUNDARYNODE *boundaryNode;


/******************************** READEDGES ********************************/
/*                                                                         */
/* Reads mesh using edge pointers and initializes solution to linear       */
/*                                                                         */
/***************************************************************************/
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
/***************************************************************************/
void Simcenter_WKA_ReadEdges(void) {
    int i, icount;
    int node1, node2;
    int idum;
    FILE *inputMesh;
    char dumstring[100];
    int *tag = NULL;

    printf("============================================================\n");
    info("Reading Mesh File %s", meshFile);
    
    // Open Mesh File
    if ((inputMesh = fopen(meshFile, "r")) == (FILE *) NULL)
        error("Unable to Open Mesh File %s", meshFile);
    
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

/***************************** EXTRACTCELLS ********************************/
/*                                                                         */
/* Extract cell-to-node pointers from edge pointers                        */
/*                                                                         */
/***************************************************************************/
void Simcenter_WKA_ExtractCells(void) {
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
    for (i = 0; i < mesh.inside; i++) {
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


