/* 
 * File:   Mesh.h
 * Author: Ashish Gupta
 *
 * Created on February 28, 2010, 9:13 PM
 */

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
    double nxc12, nyc12, magc12;
    double nxc23, nyc23, magc23;
    double nxc31, nyc31, magc31;
    /* Cell Edge Normals */
    double nx12, ny12, mag12;
    double nx23, ny23, mag23;
    double nx31, ny31, mag31;
} CELL;

typedef struct {
    int node1, node2;
    int cell1, cell2;
    double nx, ny, mag;
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
    double constant;
} BOUNDARYNODE;

#endif	/* _MESH_H */

