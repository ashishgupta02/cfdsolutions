/*******************************************************************************
 * File:        Area_Volume.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "List.h"
#include "Point3D.h"
#include "Vector3D.h"
#include "Commons.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
static void Get_Edge_Nodes_Faces_Connectivity(int cell_lid, int cell_type,
        int edge_num, int edge_node[], int faceL[], int faceR[]) {
    int node[MAX_ELEM_NNODE];
    int offset;

    edge_node[0] = -1;
    edge_node[1] = -1;
    for (int i = 0; i < 4; i++) {
        faceL[i] = -1;
        faceR[i] = -1;
    };
    
    if (cell_type == TRI || cell_type == QUAD)
        return;

    // Calculate the Offset for cell
    offset = cell_lid * elemNode[cell_type];

    // Get the cell nodes
    for (int i = 0; i < elemNode[cell_type]; i++)
        node[i] = cell2Node[cell_type][offset + i ];

    // Tetrahydra
    if (cell_type == TETRA) {
        switch (edge_num) {
            case 0:
                edge_node[0] = node[0];
                edge_node[1] = node[1];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceR[0] = 0;
                faceR[1] = 3;
                faceR[2] = 1;
                break;
            case 1:
                edge_node[0] = node[0];
                edge_node[1] = node[2];
                faceL[0] = 0;
                faceL[1] = 2;
                faceL[2] = 3;
                faceR[0] = 0;
                faceR[1] = 1;
                faceR[2] = 2;
                break;
            case 2:
                edge_node[0] = node[0];
                edge_node[1] = node[3];
                faceL[0] = 0;
                faceL[1] = 3;
                faceL[2] = 1;
                faceR[0] = 0;
                faceR[1] = 2;
                faceR[2] = 3;
                break;
            case 3:
                edge_node[0] = node[1];
                edge_node[1] = node[2];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceR[0] = 1;
                faceR[1] = 3;
                faceR[2] = 2;
                break;
            case 4:
                edge_node[0] = node[1];
                edge_node[1] = node[3];
                faceL[0] = 1;
                faceL[1] = 3;
                faceL[2] = 2;
                faceR[0] = 0;
                faceR[1] = 3;
                faceR[2] = 1;
                break;
            case 5:
                edge_node[0] = node[2];
                edge_node[1] = node[3];
                faceL[0] = 0;
                faceL[1] = 2;
                faceL[2] = 3;
                faceR[0] = 1;
                faceR[1] = 3;
                faceR[2] = 2;
                break;
        }
    }
    
    // Pyramid
    if (cell_type == PYRA) {
        switch (edge_num) {
            case 0:
                edge_node[0] = node[0];
                edge_node[1] = node[1];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceL[3] = 3;
                faceR[0] = 0;
                faceR[1] = 4;
                faceR[2] = 1;
                break;
            case 1:
                edge_node[0] = node[0];
                edge_node[1] = node[3];
                faceL[0] = 0;
                faceL[1] = 3;
                faceL[2] = 4;
                faceR[0] = 0;
                faceR[1] = 1;
                faceR[2] = 2;
                faceR[3] = 3;
                break;
            case 2:
                edge_node[0] = node[0];
                edge_node[1] = node[4];
                faceL[0] = 0;
                faceL[1] = 4;
                faceL[2] = 1;
                faceR[0] = 0;
                faceR[1] = 3;
                faceR[2] = 4;
                break;
            case 3:
                edge_node[0] = node[1];
                edge_node[1] = node[2];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceL[3] = 3;
                faceR[0] = 1;
                faceR[1] = 4;
                faceR[2] = 2;
                break;
            case 4:
                edge_node[0] = node[1];
                edge_node[1] = node[4];
                faceL[0] = 1;
                faceL[1] = 4;
                faceL[2] = 2;
                faceR[0] = 0;
                faceR[1] = 4;
                faceR[2] = 1;
                break;
            case 5:
                edge_node[0] = node[2];
                edge_node[1] = node[3];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceL[3] = 3;
                faceR[0] = 2;
                faceR[1] = 4;
                faceR[2] = 3;
                break;
            case 6:
                edge_node[0] = node[2];
                edge_node[1] = node[4];
                faceL[0] = 2;
                faceL[1] = 4;
                faceL[2] = 3;
                faceR[0] = 1;
                faceR[1] = 4;
                faceR[2] = 2;
                break;
            case 7:
                edge_node[0] = node[3];
                edge_node[1] = node[4];
                faceL[0] = 0;
                faceL[1] = 3;
                faceL[2] = 4;
                faceR[0] = 2;
                faceR[1] = 4;
                faceR[2] = 3;
                break;
        }
    }

    // Prism
    if (cell_type == PRISM) {
        switch (edge_num) {
            case 0:
                edge_node[0] = node[0];
                edge_node[1] = node[1];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceR[0] = 0;
                faceR[1] = 3;
                faceR[2] = 4;
                faceR[3] = 1;
                break;
            case 1:
                edge_node[0] = node[0];
                edge_node[1] = node[2];
                faceL[0] = 0;
                faceL[1] = 2;
                faceL[2] = 5;
                faceL[3] = 3;
                faceR[0] = 0;
                faceR[1] = 1;
                faceR[2] = 2;
                break;
            case 2:
                edge_node[0] = node[0];
                edge_node[1] = node[3];
                faceL[0] = 0;
                faceL[1] = 3;
                faceL[2] = 4;
                faceL[3] = 1;
                faceR[0] = 0;
                faceR[1] = 2;
                faceR[2] = 5;
                faceR[3] = 3;
                break;
            case 3:
                edge_node[0] = node[1];
                edge_node[1] = node[2];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceR[0] = 1;
                faceR[1] = 4;
                faceR[2] = 5;
                faceR[3] = 2;
                break;
            case 4:
                edge_node[0] = node[1];
                edge_node[1] = node[4];
                faceL[0] = 1;
                faceL[1] = 4;
                faceL[2] = 5;
                faceL[3] = 2;
                faceR[0] = 0;
                faceR[1] = 3;
                faceR[2] = 4;
                faceR[3] = 1;
                break;
            case 5:
                edge_node[0] = node[2];
                edge_node[1] = node[5];
                faceL[0] = 0;
                faceL[1] = 2;
                faceL[2] = 5;
                faceL[3] = 3;
                faceR[0] = 1;
                faceR[1] = 4;
                faceR[2] = 5;
                faceR[3] = 2;
                break;
            case 6:
                edge_node[0] = node[3];
                edge_node[1] = node[4];
                faceL[0] = 0;
                faceL[1] = 3;
                faceL[2] = 4;
                faceL[3] = 1;
                faceR[0] = 3;
                faceR[1] = 5;
                faceR[2] = 4;
                break;
            case 7:
                edge_node[0] = node[3];
                edge_node[1] = node[5];
                faceL[0] = 3;
                faceL[1] = 5;
                faceL[2] = 4;
                faceR[0] = 0;
                faceR[1] = 2;
                faceR[2] = 5;
                faceR[3] = 3;
                break;
            case 8:
                edge_node[0] = node[4];
                edge_node[1] = node[5];
                faceL[0] = 1;
                faceL[1] = 4;
                faceL[2] = 5;
                faceL[3] = 2;
                faceR[0] = 3;
                faceR[1] = 5;
                faceR[2] = 4;
                break;
        }
    }

    // Hexahedral
    if (cell_type == HEXA) {
        switch (edge_num) {
            case 0:
                edge_node[0] = node[0];
                edge_node[1] = node[1];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceL[3] = 3;
                faceR[0] = 0;
                faceR[1] = 4;
                faceR[2] = 5;
                faceR[3] = 1;
                break;
            case 1:
                edge_node[0] = node[0];
                edge_node[1] = node[3];
                faceL[0] = 0;
                faceL[1] = 3;
                faceL[2] = 7;
                faceL[3] = 4;
                faceR[0] = 0;
                faceR[1] = 1;
                faceR[2] = 2;
                faceR[3] = 3;
                break;
            case 2:
                edge_node[0] = node[0];
                edge_node[1] = node[4];
                faceL[0] = 0;
                faceL[1] = 4;
                faceL[2] = 5;
                faceL[3] = 1;
                faceR[0] = 0;
                faceR[1] = 3;
                faceR[2] = 7;
                faceR[3] = 4;
                break;
            case 3:
                edge_node[0] = node[1];
                edge_node[1] = node[2];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceL[3] = 3;
                faceR[0] = 1;
                faceR[1] = 5;
                faceR[2] = 6;
                faceR[3] = 2;
                break;
            case 4:
                edge_node[0] = node[1];
                edge_node[1] = node[5];
                faceL[0] = 1;
                faceL[1] = 5;
                faceL[2] = 6;
                faceL[3] = 2;
                faceR[0] = 0;
                faceR[1] = 4;
                faceR[2] = 5;
                faceR[3] = 1;
                break;
            case 5:
                edge_node[0] = node[2];
                edge_node[1] = node[3];
                faceL[0] = 0;
                faceL[1] = 1;
                faceL[2] = 2;
                faceL[3] = 3;
                faceR[0] = 2;
                faceR[1] = 6;
                faceR[2] = 7;
                faceR[3] = 3;
                break;
            case 6:
                edge_node[0] = node[2];
                edge_node[1] = node[6];
                faceL[0] = 2;
                faceL[1] = 6;
                faceL[2] = 7;
                faceL[3] = 3;
                faceR[0] = 1;
                faceR[1] = 5;
                faceR[2] = 6;
                faceR[3] = 2;
                break;
            case 7:
                edge_node[0] = node[3];
                edge_node[1] = node[7];
                faceL[0] = 0;
                faceL[1] = 3;
                faceL[2] = 7;
                faceL[3] = 4;
                faceR[0] = 2;
                faceR[1] = 6;
                faceR[2] = 7;
                faceR[3] = 3;
                break;
            case 8:
                edge_node[0] = node[4];
                edge_node[1] = node[5];
                faceL[0] = 0;
                faceL[1] = 4;
                faceL[2] = 5;
                faceL[3] = 1;
                faceR[0] = 4;
                faceR[1] = 7;
                faceR[2] = 6;
                faceR[3] = 5;
                break;
            case 9:
                edge_node[0] = node[4];
                edge_node[1] = node[7];
                faceL[0] = 4;
                faceL[1] = 7;
                faceL[2] = 6;
                faceL[3] = 5;
                faceR[0] = 0;
                faceR[1] = 3;
                faceR[2] = 7;
                faceR[3] = 4;
                break;
            case 10:
                edge_node[0] = node[5];
                edge_node[1] = node[6];

                faceL[0] = 1;
                faceL[1] = 5;
                faceL[2] = 6;
                faceL[3] = 2;
                faceR[0] = 4;
                faceR[1] = 7;
                faceR[2] = 6;
                faceR[3] = 5;
                break;
            case 11:
                edge_node[0] = node[6];
                edge_node[1] = node[7];
                faceL[0] = 2;
                faceL[1] = 6;
                faceL[2] = 7;
                faceL[3] = 3;
                faceR[0] = 4;
                faceR[1] = 7;
                faceR[2] = 6;
                faceR[3] = 5;
                break;
        }
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
static void Calculate_Boundary_Edge_Area() {
    int      n0, n1, n2, n3;
    Vector3D V1, V2, V3;
    Point3D  surcen, edgcen1, edgcen2;
    int bndry_node_count = 0;
    int bndry_edge_count = 0;
    
    // Triangle
    for (int i = 0; i < nElem[TRI]; i++) {
        bndry_edge_info[3 * i + 0].node[0] = cell2Node[TRI][3 * i + 0];
        bndry_edge_info[3 * i + 0].node[1] = nNode + bndry_node_count;
        bndry_edge_info[3 * i + 0].tag     = faceTag[TRI][i];

        bndry_edge_info[3 * i + 1].node[0] = cell2Node[TRI][3 * i + 1];
        bndry_edge_info[3 * i + 1].node[1] = nNode + bndry_node_count + 1;
        bndry_edge_info[3 * i + 1].tag     = faceTag[TRI][i];

        bndry_edge_info[3 * i + 2].node[0] = cell2Node[TRI][3 * i + 2];
        bndry_edge_info[3 * i + 2].node[1] = nNode + bndry_node_count + 2;
        bndry_edge_info[3 * i + 2].tag     = faceTag[TRI][i];

        bndry_node_count += 3;
        bndry_edge_count += 3;
    }

    // Quadrilateral
    for (int i = 0; i < nElem[QUAD]; i++) {
        bndry_edge_info[bndry_edge_count + 4 * i + 0].node[0] = cell2Node[QUAD][4 * i + 0];
        bndry_edge_info[bndry_edge_count + 4 * i + 0].node[1] = nNode + bndry_node_count;
        bndry_edge_info[bndry_edge_count + 4 * i + 0].tag     = faceTag[QUAD][i];

        bndry_edge_info[bndry_edge_count + 4 * i + 1].node[0] = cell2Node[QUAD][4 * i + 1];
        bndry_edge_info[bndry_edge_count + 4 * i + 1].node[1] = nNode + bndry_node_count + 1;
        bndry_edge_info[bndry_edge_count + 4 * i + 1].tag     = faceTag[QUAD][i];

        bndry_edge_info[bndry_edge_count + 4 * i + 2].node[0] = cell2Node[QUAD][4 * i + 2];
        bndry_edge_info[bndry_edge_count + 4 * i + 2].node[1] = nNode + bndry_node_count + 2;
        bndry_edge_info[bndry_edge_count + 4 * i + 2].tag     = faceTag[QUAD][i];

        bndry_edge_info[bndry_edge_count + 4 * i + 3].node[0] = cell2Node[QUAD][4 * i + 3];
        bndry_edge_info[bndry_edge_count + 4 * i + 3].node[1] = nNode + bndry_node_count + 3;
        bndry_edge_info[bndry_edge_count + 4 * i + 3].tag     = faceTag[QUAD][i];

        bndry_node_count += 4;
    }

    bndry_edge_count = 0;

    // Calculate the Area Associated with Boundary Edges (Node)
    // Triangle
    for (int i = 0; i < nElem[TRI]; i++) {
        // Get the Nodes
        n0 = cell2Node[TRI][3 * i + 0];
        n1 = cell2Node[TRI][3 * i + 1];
        n2 = cell2Node[TRI][3 * i + 2];

        // Calculate centroid
        surcen.pos[0] = (coordXYZ[3 * n0 + 0] + coordXYZ[3 * n1 + 0] + coordXYZ[3 * n2 + 0]) / 3.0;
        surcen.pos[1] = (coordXYZ[3 * n0 + 1] + coordXYZ[3 * n1 + 1] + coordXYZ[3 * n2 + 1]) / 3.0;
        surcen.pos[2] = (coordXYZ[3 * n0 + 2] + coordXYZ[3 * n1 + 2] + coordXYZ[3 * n2 + 2]) / 3.0;

        // -- Node 0
        // Calculate centroid of edges for node 0
        for (int j = 0; j < 3; j++) {
            edgcen1.pos[j] = (coordXYZ[3 * n0 + j] + coordXYZ[3 * n1 + j]) / 2.0;
            edgcen2.pos[j] = (coordXYZ[3 * n0 + j] + coordXYZ[3 * n2 + j]) / 2.0;
        }

        // Calculate actual vectors to be used in cross product
        for (int j = 0; j < 3; j++) {
            V1.vec[j] = edgcen1.pos[j] - coordXYZ[3 * n0 + j];
            V2.vec[j] = surcen.pos[j]  - coordXYZ[3 * n0 + j];
            V3.vec[j] = edgcen2.pos[j] - coordXYZ[3 * n0 + j];
        }
        // Add the two part area vector to the node
        bndry_edge_info[3 * i + 0].areav += 0.5 * (V1.cross(V2) + V2.cross(V3));

        // -- Node 1
        // Calculate centroid of edges for node 1
        for (int j = 0; j < 3; j++) {
            edgcen1.pos[j] = (coordXYZ[3 * n1 + j] + coordXYZ[3 * n2 + j]) / 2.0;
            edgcen2.pos[j] = (coordXYZ[3 * n1 + j] + coordXYZ[3 * n0 + j]) / 2.0;
        }

        // Calculate actual vectors to be used in cross product
        for (int j = 0; j < 3; j++) {
            V1.vec[j] = edgcen1.pos[j] - coordXYZ[3 * n1 + j];
            V2.vec[j] = surcen.pos[j]  - coordXYZ[3 * n1 + j];
            V3.vec[j] = edgcen2.pos[j] - coordXYZ[3 * n1 + j];
        }
        // Add the two part area vector to the node
        bndry_edge_info[3 * i + 1].areav += 0.5 * (V1.cross(V2) + V2.cross(V3));

        // -- Node 2
        //Calculate centroid of edges for node 1
        for (int j = 0; j < 3; j++) {
            edgcen1.pos[j] = (coordXYZ[3 * n2 + j] + coordXYZ[3 * n0 + j]) / 2.0;
            edgcen2.pos[j] = (coordXYZ[3 * n2 + j] + coordXYZ[3 * n1 + j]) / 2.0;
        }

        // Calculate actual vectors to be used in cross product
        for (int j = 0; j < 3; j++) {
            V1.vec[j] = edgcen1.pos[j] - coordXYZ[3 * n2 + j];
            V2.vec[j] = surcen.pos[j]  - coordXYZ[3 * n2 + j];
            V3.vec[j] = edgcen2.pos[j] - coordXYZ[3 * n2 + j];
        }
        // Add the two part area vector to the node
        bndry_edge_info[3 * i + 2].areav += 0.5 * (V1.cross(V2) + V2.cross(V3));

        // Increment the Processed Boundary Edges
        bndry_edge_count += 3;
    }

    // Quadrilateral
    for (int i = 0; i < nElem[QUAD]; i++) {
        // Get Nodes
        n0 = cell2Node[QUAD][4 * i + 0];
        n1 = cell2Node[QUAD][4 * i + 1];
        n2 = cell2Node[QUAD][4 * i + 2];
        n3 = cell2Node[QUAD][4 * i + 3];

        // Calculate centroid
        surcen.pos[0] = (coordXYZ[3 * n0 + 0] + coordXYZ[3 * n1 + 0]
                + coordXYZ[3 * n2 + 0] + coordXYZ[3 * n3 + 0]) / 4.0;
        surcen.pos[1] = (coordXYZ[3 * n0 + 1] + coordXYZ[3 * n1 + 1]
                + coordXYZ[3 * n2 + 1] + coordXYZ[3 * n3 + 1]) / 4.0;
        surcen.pos[2] = (coordXYZ[3 * n0 + 2] + coordXYZ[3 * n1 + 2]
                + coordXYZ[3 * n2 + 2] + coordXYZ[3 * n3 + 2]) / 4.0;

        // -- Node 0
        // Calculate centroid of edges for node 0
        for (int j = 0; j < 3; j++) {
            edgcen1.pos[j] = (coordXYZ[3 * n0 + j] + coordXYZ[3 * n1 + j]) / 2.0;
            edgcen2.pos[j] = (coordXYZ[3 * n0 + j] + coordXYZ[3 * n3 + j]) / 2.0;
        }

        // Calculate actual vectors to be used in cross product
        for (int j = 0; j < 3; j++) {
            V1.vec[j] = edgcen1.pos[j] - coordXYZ[3 * n0 + j];
            V2.vec[j] = surcen.pos[j]  - coordXYZ[3 * n0 + j];
            V3.vec[j] = edgcen2.pos[j] - coordXYZ[3 * n0 + j];
        }
        // Add the two part area vector to the node
        bndry_edge_info[4 * i + bndry_edge_count + 0].areav += 0.5 * (V1.cross(V2) + V2.cross(V3));

        // -- Node 1
        // Calculate centroid of edges for node 1
        for (int j = 0; j < 3; j++) {
            edgcen1.pos[j] = (coordXYZ[3 * n1 + j] + coordXYZ[3 * n2 + j]) / 2.0;
            edgcen2.pos[j] = (coordXYZ[3 * n1 + j] + coordXYZ[3 * n0 + j]) / 2.0;
        }

        // Calculate actual vectors to be used in cross product
        for (int j = 0; j < 3; j++) {
            V1.vec[j] = edgcen1.pos[j] - coordXYZ[3 * n1 + j];
            V2.vec[j] = surcen.pos[j]  - coordXYZ[3 * n1 + j];
            V3.vec[j] = edgcen2.pos[j] - coordXYZ[3 * n1 + j];
        }
        // Add the two part area vector to the node
        bndry_edge_info[4 * i + bndry_edge_count + 1].areav += 0.5 * (V1.cross(V2) + V2.cross(V3));

        // -- Node 2
        // Calculate centroid of edges for node 2
        for (int j = 0; j < 3; j++) {
            edgcen1.pos[j] = (coordXYZ[3 * n2 + j] + coordXYZ[3 * n3 + j]) / 2.0;
            edgcen2.pos[j] = (coordXYZ[3 * n2 + j] + coordXYZ[3 * n1 + j]) / 2.0;
        }

        // Calculate actual vectors to be used in cross product
        for (int j = 0; j < 3; j++) {
            V1.vec[j] = edgcen1.pos[j] - coordXYZ[3 * n2 + j];
            V2.vec[j] = surcen.pos[j]  - coordXYZ[3 * n2 + j];
            V3.vec[j] = edgcen2.pos[j] - coordXYZ[3 * n2 + j];
        }
        // Add the two part area vector to the node
        bndry_edge_info[4 * i + bndry_edge_count + 2].areav += 0.5 * (V1.cross(V2) + V2.cross(V3));

        // -- Node 3
        // Calculate centroid of edges for node 3
        for (int j = 0; j < 3; j++) {
            edgcen1.pos[j] = (coordXYZ[3 * n3 + j] + coordXYZ[3 * n0 + j]) / 2.0;
            edgcen2.pos[j] = (coordXYZ[3 * n3 + j] + coordXYZ[3 * n2 + j]) / 2.0;
        }

        // Calculate actual vectors to be used in cross product
        for (int j = 0; j < 3; j++) {
            V1.vec[j] = edgcen1.pos[j] - coordXYZ[3 * n3 + j];
            V2.vec[j] = surcen.pos[j]  - coordXYZ[3 * n3 + j];
            V3.vec[j] = edgcen2.pos[j] - coordXYZ[3 * n3 + j];
        }
        // Add the two part area vector to the node
        bndry_edge_info[4 * i + bndry_edge_count + 3].areav += 0.5 * (V1.cross(V2) + V2.cross(V3));
    }
} 

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
static void Calculate_Internal_Edge_Area_Volume() {
    int     offset;
    int     cell_lid;
    int     cell_type;
    int     faceL[4];
    int     faceR[4];
    int     edge_num;
    int     edge_node[2];
    int     test_edge;
    Vector3D V1, V2, V3, V4;
    Vector3D surav, areav;
    Point3D  volcen, edgcen, flcen, frcen;
    int      nnodes_faceL, nnodes_faceR;
    int      node[MAX_ELEM_NNODE];

    // Allocate memory to store control volume associate with node
    cVolume = (double*) malloc(nNode * sizeof (double));
    for (int i = 0; i < nNode; i++)
        cVolume[i] = 0.0;
    
    // Create Nodes to Edge Temporary Connectivity
    List **node2Edge;
    node2Edge = new List *[nNode];
    for (int i = 0; i < nNode; i++)
        node2Edge[i] = new List();

    for (int i = 0; i < nNode; i++) {
        for (int j = 0; j < nEdge; j++) {
            if (edge2Node[2 * j + 0] == i || edge2Node[2 * j + 1] == i)
                node2Edge[i]->Add_To_List(j);
        }
    }

    // Loop Over all Cells
    for (int i = 0; i < nCell; i++) {
        // Get local cell number and type
        cell_lid = i;
        for (int etype = TRI; etype <= HEXA; etype++) {
            if (i < cellGlobalOffset[etype]) {
                cell_type = etype;
                break;
            } else {
                cell_lid -= nElem[etype];
            }
        }
        offset = cell_lid * elemNode[cell_type];

        // Process Only Volume Cells
        if (cell_type != TRI && cell_type != QUAD) {
            // Get the cell nodes
            for (int j = 0; j < elemNode[cell_type]; j++)
                node[j] = cell2Node[cell_type][offset + j];

            // Determine centroid of the cell
            volcen = 0.0;
            for (int j = 0; j < elemNode[cell_type]; j++) {
                volcen.pos[0] += coordXYZ[PHY_DIM * node[j] + 0];
                volcen.pos[1] += coordXYZ[PHY_DIM * node[j] + 1];
                volcen.pos[2] += coordXYZ[PHY_DIM * node[j] + 2];
            }
            volcen /= elemNode[cell_type];

            // Loop over the edges of this cell
            for (int j = 0; j < elemEdge[cell_type]; j++) {
                edge_num = j;
                Get_Edge_Nodes_Faces_Connectivity(cell_lid, cell_type, edge_num, edge_node, faceL, faceR);

                // Decide how many nodes faceL and faceR have
                if (faceL[3] == -1)
                    nnodes_faceL = 3;
                else
                    nnodes_faceL = 4;

                if (faceR[3] == -1)
                    nnodes_faceR = 3;
                else
                    nnodes_faceR = 4;

                // Get the Actual nodes of Faces
                // Left Face
                for (int i = 0; i < nnodes_faceL; i++)
                    faceL[i] = node[faceL[i]];
                // Right Face
                for (int i = 0; i < nnodes_faceR; i++)
                    faceR[i] = node[faceR[i]];

                // Compute Left Face Centriod
                flcen = 0.0;
                for (int i = 0; i < nnodes_faceL; i++) {
                    flcen.pos[0] += coordXYZ[3 * faceL[i] + 0];
                    flcen.pos[1] += coordXYZ[3 * faceL[i] + 1];
                    flcen.pos[2] += coordXYZ[3 * faceL[i] + 2];
                }
                flcen /= (double) nnodes_faceL;

                // Compute Left Face Centriod
                frcen = 0.0;
                for (int i = 0; i < nnodes_faceR; i++) {
                    frcen.pos[0] += coordXYZ[3 * faceR[i] + 0];
                    frcen.pos[1] += coordXYZ[3 * faceR[i] + 1];
                    frcen.pos[2] += coordXYZ[3 * faceR[i] + 2];
                }
                frcen /= (double) nnodes_faceR;

                // Compute Edge Centroid
                edgcen = 0.0;
                for (int i = 0; i < 2; i++) {
                    edgcen.pos[0] += coordXYZ[3 * edge_node[i] + 0];
                    edgcen.pos[1] += coordXYZ[3 * edge_node[i] + 1];
                    edgcen.pos[2] += coordXYZ[3 * edge_node[i] + 2];
                }
                edgcen /= 2.0;

                // This is for faceL (face left)
                for (int i = 0; i < 3; i++) {
                    V1.vec[i] = volcen.pos[i] - edgcen.pos[i];
                    V2.vec[i] = flcen.pos[i] - edgcen.pos[i];
                    V3.vec[i] = edgcen.pos[i] - coordXYZ[3 * edge_node[1] + i];
                    V4.vec[i] = coordXYZ[3 * edge_node[0] + i] - edgcen.pos[i];
                }
                // Volume of small tetraheron = (1/6)*a.(bxc) => vector triple product
                cVolume[edge_node[1]] += (1.0 / 6.0)*(V3.dot(V1.cross(V2)));
                cVolume[edge_node[0]] += (1.0 / 6.0)*(V4.dot(V1.cross(V2)));
                surav = 0.5 * V1.cross(V2);
                // Store Area Vector
                areav = surav;

                // This is for faceR (face right)
                for (int i = 0; i < 3; i++) {
                    V1.vec[i] = frcen.pos[i] - edgcen.pos[i];
                    V2.vec[i] = volcen.pos[i] - edgcen.pos[i];
                    V3.vec[i] = edgcen.pos[i] - coordXYZ[3 * edge_node[1] + i];
                    V4.vec[i] = coordXYZ[3 * edge_node[0] + i] - edgcen.pos[i];
                }
                // Volume of small tetraheron = (1/6)*a.(bxc) => vector triple product
                cVolume[edge_node[1]] += (1.0 / 6.0)*(V3.dot(V1.cross(V2)));
                cVolume[edge_node[0]] += (1.0 / 6.0)*(V4.dot(V1.cross(V2)));
                surav = 0.5 * V1.cross(V2);
                // Store Add Area Vector
                areav += surav;

                // Reverse the area vector for this condition
                if (edge_node[0] > edge_node[1])
                    areav = -1.0 *areav;
                
                // Loop over edges having this node
                for (int k = 0; k < node2Edge[edge_node[0]]->max; k++) {
                    test_edge = node2Edge[edge_node[0]]->list[k];
                    // Look for edge having this nodes
                    if ((edge2Node[2 * test_edge + 0] == edge_node[0]
                            && edge2Node[2 * test_edge + 1] == edge_node[1])
                            || (edge2Node[2 * test_edge + 0] == edge_node[1]
                            && edge2Node[2 * test_edge + 1] == edge_node[0])) {
                        int_edge_info[test_edge].areav  += areav;
                        int_edge_info[test_edge].id      = test_edge;
                        int_edge_info[test_edge].node[0] = edge2Node[2 * test_edge + 0];
                        int_edge_info[test_edge].node[1] = edge2Node[2 * test_edge + 1];
                    }
                }
            }
        }
    }

    // Free the Memory
    for (int i = 0; i < nNode; i++)
        delete node2Edge[i];
    delete[] node2Edge;
}

void check_areas_closed() {
    int i;

    double *test_areas_closed;
    test_areas_closed = (double*) malloc(3 * nNode * sizeof (double));

    int bnode;
    int inode1, inode2;
    for (i = 0; i < 3 * nNode; i++)
        test_areas_closed[i] = 0.0;

    double add_over_nodes = 0.0;

    for (i = 0; i < nBEdge; i++) {
        bnode = bndry_edge_info[i].node[0];
        test_areas_closed[3 * bnode + 0] -= bndry_edge_info[i].areav.vec[0];
        test_areas_closed[3 * bnode + 1] -= bndry_edge_info[i].areav.vec[1];
        test_areas_closed[3 * bnode + 2] -= bndry_edge_info[i].areav.vec[2];
    }
    
    for (i = 0; i < nEdge; i++) {
        inode1 = int_edge_info[i].node[0];
        inode2 = int_edge_info[i].node[1];

        test_areas_closed[3 * inode1 + 0] -= int_edge_info[i].areav.vec[0];
        test_areas_closed[3 * inode1 + 1] -= int_edge_info[i].areav.vec[1];
        test_areas_closed[3 * inode1 + 2] -= int_edge_info[inode1].areav.vec[2];

        test_areas_closed[3 * inode2 + 0] += int_edge_info[i].areav.vec[0];
        test_areas_closed[3 * inode2 + 1] += int_edge_info[i].areav.vec[1];
        test_areas_closed[3 * inode2 + 2] += int_edge_info[i].areav.vec[2];
    }
    
    for (i = 0; i < nNode; i++)
        add_over_nodes = test_areas_closed[3 * i + 0]
                + test_areas_closed[3 * i + 1]
                + test_areas_closed[3 * i + 2];
    
    printf("Control Volume Leak Check -> %.15e\n", add_over_nodes);

    free(test_areas_closed);

    double temp_x = 0.0;
    double temp_y = 0.0;
    double temp_z = 0.0;

    for (i = 0; i < nBEdge; i++) {
        temp_x += bndry_edge_info[i].areav.vec[0];
        temp_y += bndry_edge_info[i].areav.vec[1];
        temp_z += bndry_edge_info[i].areav.vec[2];
    }

    temp_x = sqrt(temp_x * temp_x + temp_y * temp_y + temp_z * temp_z);
    printf("Surface leakage -> %.15e\n", temp_x);
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Calculate_Area_Volume() {
    printf("=============================================================================\n");
    info("Computing Edge Areas and Control Volumes");
    // Calculate Boundary Edge Area
    Calculate_Boundary_Edge_Area();
    // Calculate Internal Edge Area and Associate Control Volume
    Calculate_Internal_Edge_Area_Volume();
    
    // Reverse the sign of area vector
    for (int i = 0; i < nEdge; i++)
        int_edge_info[i].areav *= -1.0;
    for (int i = 0; i < nBEdge; i++)
        bndry_edge_info[i].areav *= -1.0;
    check_areas_closed();
}
