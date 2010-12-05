/*******************************************************************************
 * File:        Mesh.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _MESH_H
#define	_MESH_H

typedef struct {
    int nnodes;
    int ncells; /* Includes ghost cells */
    int inside; /* Cells in interior of mesh */
    int nedges;
    int nbedges; /* Number of boundary edges */
    int nbnodes; /* Number of boundary nodes (used only for Dirchlet boundary conditions) */
} MESH;

typedef struct {
    int node1, node2, node3;
    double area;
    /* Median Dual Properties */
    double xc, yc;
    double unxc12, unyc12, magc12;
    double unxc23, unyc23, magc23;
    double unxc31, unyc31, magc31;
    /* Cell Edge Normals */
    double unx12, uny12, mag12;
    double unx23, uny23, mag23;
    double unx31, uny31, mag31;
} CELL;

typedef struct {
    int node1, node2;
    int cell1, cell2;
    double unx, uny, mag;
} EDGE;

typedef struct {
    double x, y;
    /* Median Dual Properties */
    double area;
    double Q[4];
    double Resi[4];
} NODE;

typedef struct {
    int edgeNumber;
    int bcType;
    double c1, c2;
} BOUNDARYEDGE;

typedef struct {
    int nodeNumber;
    int bcType;
    double constant;
} BOUNDARYNODE;

#endif	/* _MESH_H */

