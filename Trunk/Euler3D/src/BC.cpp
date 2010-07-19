/*******************************************************************************
 * File:        BC.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Initialize_Boundary_Condition() {
    nBEdge = 0;
    nBNode = 0;

    // Calculate number of boundary edges
    for (int e = TRI; e <= QUAD; e++)
        for (int i = 0; i < nElem[e]; i++)
            nBEdge  += elemNode[e];
    nBNode = nBEdge;
    
    //Allocate  memory for boundary nodes and boundary edges
    bndry_edge_info = (bndry_edge_data*) malloc(nBEdge * sizeof (bndry_edge_data));
    // Allocate memory for edge
    int_edge_info = (edge_data*) malloc(nEdge * sizeof (edge_data));

    for (int i = 0; i < nBEdge; i++) {
        bndry_edge_info[i].type    = -1;
        bndry_edge_info[i].tag     = -1;
        bndry_edge_info[i].node[0] = -1;
        bndry_edge_info[i].node[1] = -1;
        bndry_edge_info[i].areav   = 0.0;
    }

    for (int i = 0; i < nEdge; i++) {
        int_edge_info[i].id      = -1;
        int_edge_info[i].node[0] = -1;
        int_edge_info[i].node[1] = -1;
        int_edge_info[i].areav   = 0.0;
    }
    
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Create_Boundary_Condition() {

}

