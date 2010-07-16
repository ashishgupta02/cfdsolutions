/*******************************************************************************
 * File:        Commons.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"

int     nNode;
int     nCell;
int     nSurfCell;
int     nVolCell;
int     nEdge;
int     nBC;
int     elemNode[NUMBER_OF_ELEM_TYPES];
int     elemFace[NUMBER_OF_ELEM_TYPES];
int     nElem[NUMBER_OF_ELEM_TYPES];
int     cellGlobalOffset[NUMBER_OF_ELEM_TYPES];
double *coordXYZ;
int    *cell2Node[NUMBER_OF_ELEM_TYPES];
int    *faceTag[2];
List  **node2Cell;
List  **node2Node;
List  **cell2Cell;
int    *edge2Node;
int    *crs_JA_Node2Node;
int    *crs_IA_Node2Node;
int    *crs_JA_Node2Cell;
int    *crs_IA_Node2Cell;
int    *crs_JA_Cell2Cell;
int    *crs_IA_Cell2Cell;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Commons_Init(void) {
    nNode           = 0;
    nCell           = 0;
    nSurfCell       = 0;
    nVolCell        = 0;
    nEdge           = 0;
    nBC             = 0;
    
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
    
    nElem[TRI]      = 0;
    nElem[QUAD]     = 0;
    nElem[TETRA]    = 0;
    nElem[PYRA]     = 0;
    nElem[PRISM]    = 0;
    nElem[HEXA]     = 0;

    for (int i = TRI; i <= HEXA; i++)
        cellGlobalOffset[i] = 0;

    coordXYZ        = NULL;

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
    edge2Node       = NULL;
    
    crs_JA_Node2Node= NULL;
    crs_IA_Node2Node= NULL;
    crs_JA_Node2Cell= NULL;
    crs_IA_Node2Cell= NULL;
    crs_JA_Cell2Cell= NULL;
    crs_IA_Cell2Cell= NULL;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Commons_Finalize(void) {
    printf("=============================================================================\n");
}

