/*******************************************************************************
 * File:        Commons.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"

int     nNode;
int     nCell;
int     nSurfCell;
int     nVolCell;
int     nEdge;
int     nBEdge;
int     nBNode;
int     nBC;
int    *bndType;
int     elemNode[NUMBER_OF_ELEM_TYPES];
int     elemFace[NUMBER_OF_ELEM_TYPES];
int     elemEdge[NUMBER_OF_ELEM_TYPES];
int     nElem[NUMBER_OF_ELEM_TYPES];
int     cellGlobalOffset[NUMBER_OF_ELEM_TYPES];
double *coordXYZ;
double *cVolume;
int    *cell2Node[NUMBER_OF_ELEM_TYPES];
int    *faceTag[2];
List  **node2Cell;
List  **node2Node;
List  **cell2Cell;
List  **node2Edge;
int    *edge2Node;
int    *crs_JA_Node2Node;
int    *crs_IA_Node2Node;
int    *crs_JA_Node2Cell;
int    *crs_IA_Node2Cell;
int    *crs_JA_Cell2Cell;
int    *crs_IA_Cell2Cell;
int    *crs_JA_Node2Edge;
int    *crs_IA_Node2Edge;
Edge3D *bndEdge;
Edge3D *intEdge;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Commons_Init(void) {
    nNode           = 0;
    nCell           = 0;
    nSurfCell       = 0;
    nVolCell        = 0;
    nEdge           = 0;
    nBEdge          = 0;
    nBNode          = 0;
    nBC             = 0;

    bndType         = NULL;
    
    elemNode[TRI]   = NNODE_TRI;
    elemNode[QUAD]  = NNODE_QUAD;
    elemNode[TETRA] = NNODE_TETRA;
    elemNode[PYRA]  = NNODE_PYRA;
    elemNode[PRISM] = NNODE_PRISM;
    elemNode[HEXA]  = NNODE_HEXA;

    elemFace[TRI]   = NFACE_TRI;
    elemFace[QUAD]  = NFACE_QUAD;
    elemFace[TETRA] = NFACE_TETRA;
    elemFace[PYRA]  = NFACE_PYRA;
    elemFace[PRISM] = NFACE_PRISM;
    elemFace[HEXA]  = NFACE_HEXA;
    
    elemEdge[TRI]   = NEDGE_TRI;
    elemEdge[QUAD]  = NEDGE_QUAD;
    elemEdge[TETRA] = NEDGE_TETRA;
    elemEdge[PYRA]  = NEDGE_PYRA;
    elemEdge[PRISM] = NEDGE_PRISM;
    elemEdge[HEXA]  = NEDGE_HEXA;
    
    nElem[TRI]      = 0;
    nElem[QUAD]     = 0;
    nElem[TETRA]    = 0;
    nElem[PYRA]     = 0;
    nElem[PRISM]    = 0;
    nElem[HEXA]     = 0;

    for (int i = TRI; i <= HEXA; i++)
        cellGlobalOffset[i] = 0;

    coordXYZ        = NULL;
    cVolume         = NULL;

    cell2Node[TRI]  = NULL;
    cell2Node[QUAD] = NULL;
    cell2Node[TETRA]= NULL;
    cell2Node[PYRA] = NULL;
    cell2Node[PRISM]= NULL;
    cell2Node[HEXA] = NULL;

    faceTag[TRI]    = NULL;
    faceTag[QUAD]   = NULL;

    node2Cell       = NULL;
    node2Node       = NULL;
    cell2Cell       = NULL;
    node2Edge       = NULL;
    edge2Node       = NULL;
    
    crs_JA_Node2Node= NULL;
    crs_IA_Node2Node= NULL;
    crs_JA_Node2Cell= NULL;
    crs_IA_Node2Cell= NULL;
    crs_JA_Cell2Cell= NULL;
    crs_IA_Cell2Cell= NULL;
    crs_JA_Node2Edge= NULL;
    crs_IA_Node2Edge= NULL;
    
    bndEdge = NULL;
    intEdge = NULL;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Commons_Finalize(void) {
    // Boundary
    if (bndType != NULL)
        delete[] bndType;
    bndType = NULL;
    
    // Coordinates
    if (coordXYZ != NULL)
        free(coordXYZ);
    coordXYZ = NULL;

    // Geometric Terms
    if (cVolume != NULL)
        free(cVolume);
    cVolume = NULL;

    // Connectivity Map Terms
    if (cell2Node[TRI] != NULL)
        free(cell2Node[TRI]);
    cell2Node[TRI] = NULL;

    if (cell2Node[QUAD] != NULL)
        free(cell2Node[QUAD]);
    cell2Node[QUAD] = NULL;

    if (cell2Node[TETRA] != NULL)
        free(cell2Node[TETRA]);
    cell2Node[TETRA] = NULL;

    if (cell2Node[PYRA] != NULL)
        free(cell2Node[PYRA]);
    cell2Node[PYRA] = NULL;

    if (cell2Node[PRISM] != NULL)
        free(cell2Node[PRISM]);
    cell2Node[PRISM] = NULL;

    if (cell2Node[HEXA] != NULL)
        free(cell2Node[HEXA]);
    cell2Node[HEXA] = NULL;

    if (faceTag[TRI] != NULL)
        free(faceTag[TRI]);
    faceTag[TRI] = NULL;

    if (faceTag[QUAD] != NULL)
        free(faceTag[QUAD]);
    faceTag[QUAD] = NULL;

    if (node2Cell != NULL) {
        for (int n = 0; n < nNode; n++)
            delete node2Cell[n];
        delete[] node2Cell;
    }
    node2Cell = NULL;

    if (node2Node != NULL) {
        for (int n = 0; n < nNode; n++)
            delete node2Node[n];
        delete[] node2Node;
    }
    node2Node = NULL;

    if (node2Edge != NULL) {
        for (int n = 0; n < nNode; n++)
            delete node2Edge[n];
        delete[] node2Edge;
    }
    node2Edge = NULL;

    if (cell2Cell != NULL) {
        for (int c = 0; c < nCell; c++)
            delete cell2Cell[c];
        delete[] cell2Cell;
    }
    cell2Cell = NULL;

    if (edge2Node != NULL)
        free(edge2Node);
    edge2Node = NULL;
    
    if (crs_JA_Node2Node != NULL)
        free(crs_JA_Node2Node);
    crs_JA_Node2Node = NULL;

    if (crs_IA_Node2Node != NULL)
        free(crs_IA_Node2Node);
    crs_IA_Node2Node = NULL;

    if (crs_JA_Node2Cell != NULL)
        free(crs_JA_Node2Cell);
    crs_JA_Node2Cell = NULL;

    if (crs_IA_Node2Cell != NULL)
        free(crs_IA_Node2Cell);
    crs_IA_Node2Cell = NULL;

    if (crs_JA_Cell2Cell != NULL)
        free(crs_JA_Cell2Cell);
    crs_JA_Cell2Cell = NULL;

    if (crs_IA_Cell2Cell != NULL)
        free(crs_IA_Cell2Cell);
    crs_IA_Cell2Cell = NULL;

    if (crs_JA_Node2Edge != NULL)
        free(crs_JA_Node2Edge);
    crs_JA_Node2Edge = NULL;

    if (crs_IA_Node2Edge != NULL)
        free(crs_IA_Node2Edge);
    crs_IA_Node2Edge = NULL;

    // Edge Data
    if (bndEdge != NULL)
        free(bndEdge);
    bndEdge = NULL;

    if (intEdge != NULL)
        free(intEdge);
    intEdge = NULL;
    
    printf("=============================================================================\n");
}

