/*******************************************************************************
 * File:        Commons.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef COMMONS_H
#define	COMMONS_H

#include "List.h"
#include "Vector3D.h"

#define PHY_DIM                         3
#define CELL_DIM                        3

// Cell
#define NUMBER_OF_ELEM_TYPES            6
#define TRI                             0
#define QUAD                            1
#define TETRA                           2
#define PYRA                            3
#define PRISM                           4
#define HEXA                            5

// Node
#define MAX_ELEM_NNODE                  8
#define NNODE_TRI                       3
#define NNODE_QUAD                      4
#define NNODE_TETRA                     4
#define NNODE_PYRA                      5
#define NNODE_PRISM                     6
#define NNODE_HEXA                      8

// Face
#define MAX_ELEM_NFACE                  6
#define NFACE_TRI                       1
#define NFACE_QUAD                      1
#define NFACE_TETRA                     4
#define NFACE_PYRA                      5
#define NFACE_PRISM                     5
#define NFACE_HEXA                      6

// Edge
#define MAX_ELEM_NEDGE                  12
#define NEDGE_TRI                       3
#define NEDGE_QUAD                      4
#define NEDGE_TETRA                     6
#define NEDGE_PYRA                      8
#define NEDGE_PRISM                     9
#define NEDGE_HEXA                      12

#define Invalid_Element                -999

#define _ASSERT(x) {assert(x);}

#define UNSIO_FLOWTYPE_UNKNOWN         -1
#define UNSIO_FLOWTYPE_INCOMPRESSIBLE   0
#define UNSIO_FLOWTYPE_COMPRESSIBLE     1
#define UNSIO_FLOWTYPE_VARMACHRUP       2

// Grid Connectivity Data Structure
extern int     nNode;
extern int     nCell;
extern int     nSurfCell;
extern int     nVolCell;
extern int     nEdge;
extern int     nBEdge;
extern int     nBNode;
extern int     nBC;
extern int     elemNode[NUMBER_OF_ELEM_TYPES];
extern int     elemFace[NUMBER_OF_ELEM_TYPES];
extern int     elemEdge[NUMBER_OF_ELEM_TYPES];
extern int     nElem[NUMBER_OF_ELEM_TYPES];
extern int     cellGlobalOffset[NUMBER_OF_ELEM_TYPES];
extern double *coordXYZ;
extern double *cVolume;
extern int    *cell2Node[NUMBER_OF_ELEM_TYPES];
extern int    *faceTag[2];
extern List  **node2Cell;
extern List  **node2Node;
extern List  **cell2Cell;
extern int    *edge2Node;
extern int    *crs_JA_Node2Node;
extern int    *crs_IA_Node2Node;
extern int    *crs_JA_Node2Cell;
extern int    *crs_IA_Node2Cell;
extern int    *crs_JA_Cell2Cell;
extern int    *crs_IA_Cell2Cell;

// Edge Data Structure
typedef struct edge_data {
    int      id;
    int      node[2];
    Vector3D areav;
} edge_data;

typedef struct bndry_edge_data {
    int      type;
    int      tag;
    int      node[2];
    Vector3D areav;
} bndry_edge_data;

extern bndry_edge_data *bndry_edge_info;
extern edge_data *int_edge_info;

// Basic Functions
void Commons_Init(void);
void Commons_Finalize(void);

// Create Connectivity Maps
void Create_Connectivity_Maps(int reOrder);
// Reordering of Graph to reduce Matrix-Computation Bandwidth
void Cuthill_Mckee_Reorder(void);
// Calculate Areas and Control Volume
void Calculate_Area_Volume();

void Initialize_Boundary_Condition();
void Create_Boundary_Condition();

#endif	/* COMMONS_H */
