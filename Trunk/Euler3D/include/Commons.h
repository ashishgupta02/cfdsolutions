/*******************************************************************************
 * File:        Commons.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef COMMONS_H
#define	COMMONS_H

#include <math.h>
#include "List.h"
#include "Vector3D.h"

#ifndef isnan
inline bool isnan(double x) {
    return x != x;
}
#endif

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
extern int    *bndType;
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
extern List  **node2Edge;
extern int    *edge2Node;
extern int    *crs_JA_Node2Node;
extern int    *crs_IA_Node2Node;
extern int    *crs_JA_Node2Cell;
extern int    *crs_IA_Node2Cell;
extern int    *crs_JA_Cell2Cell;
extern int    *crs_IA_Cell2Cell;
extern int    *crs_JA_Node2Edge;
extern int    *crs_IA_Node2Edge;

// Edge Data Structure
typedef struct Edge3D {
    int      type;
    int      tag;
    int      node[2];
    Vector3D areav;
} Edge3D;

extern Edge3D *bndEdge;
extern Edge3D *intEdge;

// Basic Functions
void Commons_Init(void);
void Commons_Finalize(void);

// Create Connectivity Maps
void Create_Connectivity_Maps(int reOrder);
// Free Excess Memory Used in Connectivity Creation
void Trim_Connectivity_Memory(void);
// Reordering of Graph to reduce Matrix-Computation Bandwidth
void Cuthill_Mckee_Reorder(void);
// Calculate Areas and Control Volume
void Calculate_Area_Volume(void);

#endif	/* COMMONS_H */

