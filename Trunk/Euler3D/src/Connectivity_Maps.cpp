/*******************************************************************************
 * File:        Connectivity_Maps.cpp
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
//! Create Connectivity: Node2Cell
//------------------------------------------------------------------------------
static void Create_Connectivity_Node2Cell() {
    // Allocate memory to store connectivity
    node2Cell = new List *[nNode];
    for (int n = 0; n < nNode; n++)
        node2Cell[n] = new List();

    int gid;
    // Loop over all Cell Types
    for (int etype = TRI; etype <= HEXA ; etype++) {
        // Loop over cells of this type
        for (int c = 0; c < nElem[etype]; c++) {
            // Get the Global ID
            gid = c;
            for (int e = TRI; e < etype; e++)
                gid += nElem[e];
            // Loop over the nodes of this cell
            for (int n = 0; n < elemNode[etype]; n++)
                node2Cell[cell2Node[etype][elemNode[etype]*c + n]]->Add_To_List(gid);
        }
    }

    // Compress the Information in CRS Format
    int counter;
    crs_IA_Node2Cell = (int*) malloc((nNode + 1) * sizeof (int));
    crs_IA_Node2Cell[0] = 0;
    for (int n = 1; n < nNode + 1; n++) {
        crs_IA_Node2Cell[n] = node2Cell[n - 1]->max;
        crs_IA_Node2Cell[n] += crs_IA_Node2Cell[n - 1];
    }

    crs_JA_Node2Cell = (int*) malloc(crs_IA_Node2Cell[nNode] * sizeof (int));
    for (int n = 0; n < nNode; n++) {
        counter = 0;
        for (int i = crs_IA_Node2Cell[n]; i < crs_IA_Node2Cell[n + 1]; i++) {
            crs_JA_Node2Cell[i] = node2Cell[n]->list[counter];
            counter++;
        }
    }
}

//------------------------------------------------------------------------------
//! Create Connectivity: Node2Node
//------------------------------------------------------------------------------
static void Create_Connectivity_Node2Node() {
    int cell_gid  = -1;
    int cell_type = -1;
    int cell_lid  = -1;
    
    // Allocate memory to store connectivity
    node2Node = new List *[nNode];
    for (int n = 0; n < nNode; n++)
        node2Node[n] = new List();

    // Declare variables needed to build node to node hash table
    int node_list[8];
    for (int n = 0; n < 8; n++)
        node_list[n] = -1;
    
    int n0, n1, n2, n3, n4, n5, n6, n7;
    
    for (int n = 0; n < nNode; n++) {
        
        for (int i = 0; i < node2Cell[n]->max; i++) {
            for (int j = 0; j < 8; j++)
                node_list[j] = -1;
            cell_gid = node2Cell[n]->list[i];

            // Get local cell number and type
            cell_lid = cell_gid;
            for (int etype = TRI; etype <= HEXA; etype++) {
                if (cell_gid < cellGlobalOffset[etype]) {
                    cell_type = etype;
                    break;
                } else {
                    cell_lid -= nElem[etype];
                }
            }

            // Triangle
            if (cell_type == TRI) {
                n0 = cell2Node[TRI][NNODE_TRI * cell_lid + 0];
                n1 = cell2Node[TRI][NNODE_TRI * cell_lid + 1];
                n2 = cell2Node[TRI][NNODE_TRI * cell_lid + 2];

                if (n == n0) {
                    node_list[0] = n1;
                    node_list[1] = n2;
                }

                if (n == n1) {
                    node_list[0] = n0;
                    node_list[1] = n2;
                }

                if (n == n2) {
                    node_list[0] = n0;
                    node_list[1] = n1;
                }

                // Add the Nodes to list
                for (int j = 0; j < 2; j++) {
                    if (node_list[j] != -1)
                        node2Node[n]->Check_List(node_list[j]);
                }
            }

            // Quadrilateral
            if (cell_type == QUAD) {
                n0 = cell2Node[QUAD][NNODE_QUAD * cell_lid + 0];
                n1 = cell2Node[QUAD][NNODE_QUAD * cell_lid + 1];
                n2 = cell2Node[QUAD][NNODE_QUAD * cell_lid + 2];
                n3 = cell2Node[QUAD][NNODE_QUAD * cell_lid + 3];

                if (n == n0) {
                    node_list[0] = n1;
                    node_list[1] = n3;
                }

                if (n == n1) {
                    node_list[0] = n0;
                    node_list[1] = n2;
                }

                if (n == n2) {
                    node_list[0] = n1;
                    node_list[1] = n3;
                }

                if (n == n3) {
                    node_list[0] = n0;
                    node_list[1] = n2;
                }
                // Add the Nodes to list
                for (int j = 0; j < 2; j++) {
                    if (node_list[j] != -1)
                        node2Node[n]->Check_List(node_list[j]);
                }
            }

            // Tetrahedral
            if (cell_type == TETRA) {
                n0 = cell2Node[TETRA][NNODE_TETRA * cell_lid + 0];
                n1 = cell2Node[TETRA][NNODE_TETRA * cell_lid + 1];
                n2 = cell2Node[TETRA][NNODE_TETRA * cell_lid + 2];
                n3 = cell2Node[TETRA][NNODE_TETRA * cell_lid + 3];
                
                if (n == n0) {
                    node_list[0] = n1;
                    node_list[1] = n2;
                    node_list[2] = n3;
                }

                if (n == n1) {
                    node_list[0] = n0;
                    node_list[1] = n2;
                    node_list[2] = n3;
                }

                if (n == n2) {
                    node_list[0] = n1;
                    node_list[1] = n0;
                    node_list[2] = n3;
                }

                if (n == n3) {
                    node_list[0] = n0;
                    node_list[1] = n1;
                    node_list[2] = n2;
                }
                // Add the Nodes to list
                for (int j = 0; j < 3; j++) {
                    if (node_list[j] != -1)
                        node2Node[n]->Check_List(node_list[j]);
                }
            }

            // Pyramid
            if (cell_type == PYRA) {
                n0 = cell2Node[PYRA][NNODE_PYRA * cell_lid + 0];
                n1 = cell2Node[PYRA][NNODE_PYRA * cell_lid + 1];
                n2 = cell2Node[PYRA][NNODE_PYRA * cell_lid + 2];
                n3 = cell2Node[PYRA][NNODE_PYRA * cell_lid + 3];
                n4 = cell2Node[PYRA][NNODE_PYRA * cell_lid + 4];

                if (n == n0) {
                    node_list[0] = n1;
                    node_list[1] = n3;
                    node_list[2] = n4;
                }

                if (n == n1) {
                    node_list[0] = n2;
                    node_list[1] = n0;
                    node_list[2] = n4;
                }

                if (n == n2) {
                    node_list[0] = n1;
                    node_list[1] = n4;
                    node_list[2] = n3;
                }

                if (n == n3) {
                    node_list[0] = n4;
                    node_list[1] = n0;
                    node_list[2] = n2;
                }


                if (n == n4) {
                    node_list[0] = n0;
                    node_list[1] = n1;
                    node_list[2] = n2;
                    node_list[3] = n3;
                }
                // Add the Nodes to list
                for (int j = 0; j < 4; j++) {
                    if (node_list[j] != -1)
                        node2Node[n]->Check_List(node_list[j]);
                }
            }

            // Prism
            if (cell_type == PRISM) {
                n0 = cell2Node[PRISM][NNODE_PRISM * cell_lid + 0];
                n1 = cell2Node[PRISM][NNODE_PRISM * cell_lid + 1];
                n2 = cell2Node[PRISM][NNODE_PRISM * cell_lid + 2];
                n3 = cell2Node[PRISM][NNODE_PRISM * cell_lid + 3];
                n4 = cell2Node[PRISM][NNODE_PRISM * cell_lid + 4];
                n5 = cell2Node[PRISM][NNODE_PRISM * cell_lid + 5];

                if (n == n0) {
                    node_list[0] = n1;
                    node_list[1] = n3;
                    node_list[2] = n2;
                }

                if (n == n1) {
                    node_list[0] = n2;
                    node_list[1] = n0;
                    node_list[2] = n4;
                }

                if (n == n2) {
                    node_list[0] = n0;
                    node_list[1] = n1;
                    node_list[2] = n5;
                }

                if (n == n3) {
                    node_list[0] = n4;
                    node_list[1] = n0;
                    node_list[2] = n5;
                }

                if (n == n4) {
                    node_list[0] = n1;
                    node_list[1] = n3;
                    node_list[2] = n5;
                }

                if (n == n5) {
                    node_list[0] = n2;
                    node_list[1] = n4;
                    node_list[2] = n3;
                }

                // Add the Nodes to list
                for (int j = 0; j < 3; j++) {
                    if (node_list[j] != -1)
                        node2Node[n]->Check_List(node_list[j]);
                }
            }
            
            // Hexahedral
            if (cell_type == HEXA) {
                n0 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 0];
                n1 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 1];
                n2 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 2];
                n3 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 3];
                n4 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 4];
                n5 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 5];
                n6 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 6];
                n7 = cell2Node[HEXA][NNODE_HEXA * cell_lid + 7];
                
                if (n == n0) {
                    node_list[0] = n1;
                    node_list[1] = n4;
                    node_list[2] = n3;
                }

                if (n == n1) {
                    node_list[0] = n2;
                    node_list[1] = n0;
                    node_list[2] = n5;
                }

                if (n == n2) {
                    node_list[0] = n1;
                    node_list[1] = n3;
                    node_list[2] = n6;
                }

                if (n == n3) {
                    node_list[0] = n0;
                    node_list[1] = n2;
                    node_list[2] = n7;
                }

                if (n == n4) {
                    node_list[0] = n0;
                    node_list[1] = n7;
                    node_list[2] = n5;
                }

                if (n == n5) {
                    node_list[0] = n1;
                    node_list[1] = n6;
                    node_list[2] = n4;
                }

                if (n == n6) {
                    node_list[0] = n2;
                    node_list[1] = n5;
                    node_list[2] = n7;
                }

                if (n == n7) {
                    node_list[0] = n3;
                    node_list[1] = n6;
                    node_list[2] = n4;
                }
                // Add the Nodes to list
                for (int j = 0; j < 3; j++) {
                    if (node_list[j] != -1)
                        node2Node[n]->Check_List(node_list[j]);
                }
            }
        }
    }

    // Compress the Information in CRS Format
    int counter;
    crs_IA_Node2Node = (int*) malloc((nNode + 1) * sizeof (int));
    crs_IA_Node2Node[0] = 0;
    for (int n = 1; n < nNode + 1; n++) {
        crs_IA_Node2Node[n] = node2Node[n - 1]->max;
        crs_IA_Node2Node[n] += crs_IA_Node2Node[n - 1];
    }

    crs_JA_Node2Node = (int*) malloc(crs_IA_Node2Node[nNode] * sizeof (int));
    for (int n = 0; n < nNode; n++) {
        counter = 0;
        for (int i = crs_IA_Node2Node[n]; i < crs_IA_Node2Node[n + 1]; i++) {
            crs_JA_Node2Node[i] = node2Node[n]->list[counter];
            counter++;
        }
    }
}

//------------------------------------------------------------------------------
//! Create Connectivity: Cell2Cell
//------------------------------------------------------------------------------
static void Create_Connectivity_Cell2Cell() {
    int *faces[NUMBER_OF_ELEM_TYPES];
    int nnodes_faces[NUMBER_OF_ELEM_TYPES] = {3, 4, 12, 16, 18, 24};

    // This is the (number of nodes for all faces for each element type) X (number of elements)
    int total_face_nodes[NUMBER_OF_ELEM_TYPES];
    for (int i = TRI; i <= HEXA; i++) {
        total_face_nodes[i] = nnodes_faces[i] * nElem[i];
    }

    // Now allocate faces array. This will hold all the face information for each element type
    for (int i = TRI; i <= HEXA; i++) {
        faces[i] = (int*) malloc(total_face_nodes[i] * sizeof (int));
    }

    // Now we create the ia index array FOR EACH ELEMENT TYPE
    int *ia_faces[NUMBER_OF_ELEM_TYPES];
    
    for (int i = TRI; i <= HEXA; i++) {
        ia_faces[i] = (int*) malloc((elemFace[i] + 1) * sizeof (int));
    }

    // Define index (ia) array for each element type
    // Triangle
    ia_faces[TRI][0] = 0;
    ia_faces[TRI][1] = 3;

    // Quadrilateral
    ia_faces[QUAD][0] = 0;
    ia_faces[QUAD][1] = 4;

    // Tetrahedral
    ia_faces[TETRA][0] = 0;
    ia_faces[TETRA][1] = 3;
    ia_faces[TETRA][2] = 6;
    ia_faces[TETRA][3] = 9;
    ia_faces[TETRA][4] = 12;

    // Pyramid
    ia_faces[PYRA][0] = 0;
    ia_faces[PYRA][1] = 4;
    ia_faces[PYRA][2] = 7;
    ia_faces[PYRA][3] = 10;
    ia_faces[PYRA][4] = 13;
    ia_faces[PYRA][5] = 16;

    // Prism
    ia_faces[PRISM][0] = 0;
    ia_faces[PRISM][1] = 3;
    ia_faces[PRISM][2] = 7;
    ia_faces[PRISM][3] = 11;
    ia_faces[PRISM][4] = 15;
    ia_faces[PRISM][5] = 18;

    // Hexahedral
    ia_faces[HEXA][0] = 0;
    ia_faces[HEXA][1] = 4;
    ia_faces[HEXA][2] = 8;
    ia_faces[HEXA][3] = 12;
    ia_faces[HEXA][4] = 16;
    ia_faces[HEXA][5] = 20;
    ia_faces[HEXA][6] = 24;

    // Triangle - Get Faces
    for (int i = 0; i < nElem[TRI]; i++) {
        faces[TRI][0 + i * nnodes_faces[TRI]] = cell2Node[TRI][NNODE_TRI * i + 0];
        faces[TRI][1 + i * nnodes_faces[TRI]] = cell2Node[TRI][NNODE_TRI * i + 1];
        faces[TRI][2 + i * nnodes_faces[TRI]] = cell2Node[TRI][NNODE_TRI * i + 2];
    }

    // Quadrilateral - Get Faces
    for (int i = 0; i < nElem[QUAD]; i++) {
        faces[QUAD][0 + i * nnodes_faces[QUAD]] = cell2Node[QUAD][NNODE_QUAD * i + 0];
        faces[QUAD][1 + i * nnodes_faces[QUAD]] = cell2Node[QUAD][NNODE_QUAD * i + 1];
        faces[QUAD][2 + i * nnodes_faces[QUAD]] = cell2Node[QUAD][NNODE_QUAD * i + 2];
        faces[QUAD][3 + i * nnodes_faces[QUAD]] = cell2Node[QUAD][NNODE_QUAD * i + 3];
    }

    // Tetrahedral - - Get Faces
    for (int i = 0; i < nElem[TETRA]; i++) {
        //Face 1 = 1-2-3
        faces[TETRA][0 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 0];
        faces[TETRA][1 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 1];
        faces[TETRA][2 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 2];

        //Face 2: 1-4-2
        faces[TETRA][3 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 0];
        faces[TETRA][4 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 3];
        faces[TETRA][5 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 1];

        //Face 3: 1-3-4
        faces[TETRA][6 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 0];
        faces[TETRA][7 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 2];
        faces[TETRA][8 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 3];

        //Face 4: 2-4-3
        faces[TETRA][ 9 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 1];
        faces[TETRA][10 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 3];
        faces[TETRA][11 + i * nnodes_faces[TETRA]] = cell2Node[TETRA][NNODE_TETRA * i + 2];
    }

    // Pyramid - Get Faces
    for (int i = 0; i < nElem[PYRA]; i++) {
        //Face 1 = 1-2-3-4
        faces[PYRA][0 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 0];
        faces[PYRA][1 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 1];
        faces[PYRA][2 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 2];
        faces[PYRA][3 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 3];

        //Face 2: 1-5-2
        faces[PYRA][4 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 0];
        faces[PYRA][5 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 4];
        faces[PYRA][6 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 1];

        //Face 3: 1-4-5
        faces[PYRA][7 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 0];
        faces[PYRA][8 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 3];
        faces[PYRA][9 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 4];

        //Face 4: 2-5-3
        faces[PYRA][10 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 1];
        faces[PYRA][11 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 4];
        faces[PYRA][12 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 2];

        //Face 5: 2-4-3
        faces[PYRA][13 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 2];
        faces[PYRA][14 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 4];
        faces[PYRA][15 + i * nnodes_faces[PYRA]] = cell2Node[PYRA][NNODE_PYRA * i + 3];
    }

    // Prism - Get Faces
    for (int i = 0; i < nElem[PRISM]; i++) {
        //Face 1 = 1-2-3
        faces[PRISM][0 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 0];
        faces[PRISM][1 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 1];
        faces[PRISM][2 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 2];

        //Face 2: 1-4-5-2
        faces[PRISM][3 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 0];
        faces[PRISM][4 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 3];
        faces[PRISM][5 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 4];
        faces[PRISM][6 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 1];

        //Face 3: 1-3-6-4
        faces[PRISM][ 7 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 0];
        faces[PRISM][ 8 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 2];
        faces[PRISM][ 9 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 5];
        faces[PRISM][10 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 3];

        //Face 4: 2-5-6-3
        faces[PRISM][11 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 1];
        faces[PRISM][12 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 4];
        faces[PRISM][13 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 5];
        faces[PRISM][14 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 2];

        //Face 5: 4-6-5
        faces[PRISM][15 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 3];
        faces[PRISM][16 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 5];
        faces[PRISM][17 + i * nnodes_faces[PRISM]] = cell2Node[PRISM][NNODE_PRISM * i + 4];
    }

    // Hexahedral - Get Faces
    for (int i = 0; i < nElem[HEXA]; i++) {
        //Face 1 = 1-2-3-4
        faces[HEXA][0 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 0];
        faces[HEXA][1 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 1];
        faces[HEXA][2 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 2];
        faces[HEXA][3 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 3];

        //Face 2: 1-5-6-2
        faces[HEXA][4 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 0];
        faces[HEXA][5 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 4];
        faces[HEXA][6 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 5];
        faces[HEXA][7 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 1];

        //Face 3: 1-4-8-5
        faces[HEXA][ 8 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 0];
        faces[HEXA][ 9 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 3];
        faces[HEXA][10 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 7];
        faces[HEXA][11 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 4];

        //Face 4: 2-6-7-3
        faces[HEXA][12 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 1];
        faces[HEXA][13 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 5];
        faces[HEXA][14 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 6];
        faces[HEXA][15 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 2];

        //Face 5: 3-7-8-4
        faces[HEXA][16 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 2];
        faces[HEXA][17 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 6];
        faces[HEXA][18 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 7];
        faces[HEXA][19 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 3];

        //Face 5: 5-8-7-6
        faces[HEXA][20 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 4];
        faces[HEXA][21 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 7];
        faces[HEXA][22 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 6];
        faces[HEXA][23 + i * nnodes_faces[HEXA]] = cell2Node[HEXA][NNODE_HEXA * i + 5];
    }
    
    // Now create Cell2Cell Connectivity
    int cell_type = -1;
    int cell_lid  = -1;
    int cell_gid  = -1;
    
    // Allocate Memory to store connectivity
    cell2Cell = new List *[nCell];
    for (int c = 0; c < nCell; c++)
        cell2Cell[c] = new List();

    for (int c = 0; c < nCell; c++) {
        // Get local cell type
        cell_type = -1;
        for (int etype = TRI; etype <= HEXA; etype++) {
            if (c < cellGlobalOffset[etype]) {
                cell_type = etype;
                break;
            }
        }
        
        for (int i = 0; i < elemFace[cell_type]; i++)
            cell2Cell[c]->Add_To_List(-1);
    }

    int combo_3[6][3] ={
        {0, 1, 2},
        {0, 2, 1},
        {1, 0, 2},
        {1, 2, 0},
        {2, 0, 1},
        {2, 1, 0},
    };

    int combo_4[24][4] = {
        {0, 1, 2, 3},
        {0, 1, 3, 2},
        {0, 2, 1, 3},
        {0, 2, 3, 1},
        {0, 3, 1, 2},
        {0, 3, 2, 1},
        {1, 0, 2, 3},
        {1, 0, 3, 2},
        {1, 2, 0, 3},
        {1, 2, 3, 0},
        {1, 3, 0, 2},
        {1, 3, 2, 0},
        {2, 0, 1, 3},
        {2, 0, 3, 1},
        {2, 1, 0, 3},
        {2, 1, 3, 0},
        {2, 3, 0, 1},
        {2, 3, 1, 0},
        {3, 0, 1, 2},
        {3, 0, 2, 1},
        {3, 1, 0, 2},
        {3, 1, 2, 0},
        {3, 2, 0, 1},
        {3, 2, 1, 0},
    };

    int other_cell_gid  = -1;
    int other_cell_lid  = -1;
    int other_cell_type = -1;
    int index           = -1;
    int nf_cell_lid     = 0;
    int nf_other_cell_lid = 0;
    int node_c[4];
    int node_oc[4];
    int test;
    int count = 0;
    int ncombs;

    // Loop Over all Nodes
    for (int n = 0; n < nNode; n++) {
        // Loop over all cells connected to this node
        for (int c = 0; c < node2Cell[n]->max; c++) {
            cell_gid = node2Cell[n]->list[c];
            
            // Get local cell number and type
            cell_lid = cell_gid;
            for (int etype = TRI; etype <= HEXA; etype++) {
                if (cell_gid < cellGlobalOffset[etype]) {
                    cell_type = etype;
                    break;
                } else {
                    cell_lid -= nElem[etype];
                }
            }

            // Loop over other cells connected to this nodes except seed cell
            for (int oc = 0; oc < node2Cell[n]->max; oc++) {
                if (c != oc) {
                    other_cell_gid = node2Cell[n]->list[oc];

                    // Get local cell number and type
                    other_cell_lid = other_cell_gid;
                    for (int etype = TRI; etype <= HEXA; etype++) {
                        if (other_cell_gid < cellGlobalOffset[etype]) {
                            other_cell_type = etype;
                            break;
                        } else {
                            other_cell_lid -= nElem[etype];
                        }
                    }
                    
                    // Loop over faces of the seed cell
                    for (int cf = 0; cf < elemFace[cell_type]; cf++) {
                        
                        nf_cell_lid = ia_faces[cell_type][cf + 1] - ia_faces[cell_type][cf];

                        // Get the Nodes of the Seed Cell Face
                        node_c[0] = -1;
                        node_c[1] = -1;
                        node_c[2] = -1;
                        node_c[3] = -1;
                        count = 0;
                        for (int i = ia_faces[cell_type][cf]; i < ia_faces[cell_type][cf + 1]; i++) {
                            node_c[count] = faces[cell_type][cell_lid * nnodes_faces[cell_type] + i];
                            count++;
                        }
                        
                        // Loop over the faces of the other element
                        for (int ocf = 0; ocf < elemFace[other_cell_type]; ocf++) {
                            nf_other_cell_lid = ia_faces[other_cell_type][ocf + 1] - ia_faces[other_cell_type][ocf];

                            // Check if Seed Cell and Other Cell current face has same number of nodes
                            if (nf_other_cell_lid == nf_cell_lid) {
                                // Get the Nodes of the Other Cell Face
                                count = 0;
                                node_oc[0] = -2;
                                node_oc[1] = -2;
                                node_oc[2] = -2;
                                node_oc[3] = -2;
                                for (int i = ia_faces[other_cell_type][ocf]; i < ia_faces[other_cell_type][ocf + 1]; i++) {
                                    node_oc[count] = faces[other_cell_type][other_cell_lid * nnodes_faces[other_cell_type] + i];
                                    count++;
                                }

                                test = 0;
                                // For Triangle Face
                                if (nf_other_cell_lid == NNODE_TRI) {
                                    ncombs = 6;
                                    for (int i = 0; i < ncombs; i++) {
                                        if (node_c[0] == node_oc[combo_3[i][0]]
                                                && node_c[1] == node_oc[combo_3[i][1]]
                                                && node_c[2] == node_oc[combo_3[i][2]])
                                            test = 1;
                                    }
                                }

                                // For Quadrilateral Face
                                if (nf_other_cell_lid == NNODE_QUAD) {
                                    ncombs = 24;
                                    for (int i = 0; i < ncombs; i++) {
                                        if (node_c[0] == node_oc[combo_4[i][0]]
                                                && node_c[1] == node_oc[combo_4[i][1]]
                                                && node_c[2] == node_oc[combo_4[i][2]]
                                                && node_c[3] == node_oc[combo_4[i][3]])
                                            test = 1;
                                    }
                                }

                                if (test == 1)
                                    index = cf;
                            }
                        }
                    }
                    
                    if (index != -1)
                        cell2Cell[cell_gid]->list[index] = other_cell_gid;
                }
            }
        }
    }

    // Compress the Information in CRS Format
    int counter;
    crs_IA_Cell2Cell = (int*) malloc((nCell + 1) * sizeof (int));
    crs_IA_Cell2Cell[0] = 0;
    for (int c = 1; c < nCell + 1; c++) {
        crs_IA_Cell2Cell[c] = cell2Cell[c - 1]->max;
        crs_IA_Cell2Cell[c] += crs_IA_Cell2Cell[c - 1];
    }

    crs_JA_Cell2Cell = (int*) malloc(crs_IA_Cell2Cell[nCell] * sizeof (int));
    for (int c = 0; c < nCell; c++) {
        counter = 0;
        for (int i = crs_IA_Cell2Cell[c]; i < crs_IA_Cell2Cell[c + 1]; i++) {
            crs_JA_Cell2Cell[i] = cell2Cell[c]->list[counter];
            counter++;
        }
    }
    
    // Free Memory Used
    for (int i = 0; i < NUMBER_OF_ELEM_TYPES; i++)
        free(ia_faces[i]);

    for (int i = 0; i < NUMBER_OF_ELEM_TYPES; i++)
        free(faces[i]);
}

//------------------------------------------------------------------------------
//! Create Connectivity: Edge2Node
//------------------------------------------------------------------------------
static void Create_Connectivity_Edge2Node() {
    int counter;
    int nodeid;

    // Find the Number of Edges
    // Seed node less then connected node add edge
    nEdge = 0;
    for (int n = 0; n < nNode; n++) {
        for (int i = 0; i < node2Node[n]->max; i++) {
            nodeid = node2Node[n]->list[i];
            if (nodeid > n)
                nEdge++;
        }
    }

    // Now allocate memory to store connectivity
    edge2Node = (int*) malloc(2 * nEdge * sizeof (int));
    counter = 0;
    for (int n = 0; n < nNode; n++) {
        for (int i = 0; i < node2Node[n]->max; i++) {
            nodeid = node2Node[n]->list[i];
            if (nodeid > n) {
                edge2Node[2 * counter]     = n;
                edge2Node[2 * counter + 1] = nodeid;
                counter++;
            }
        }
    }
}

//------------------------------------------------------------------------------
//! Verify Surface Element Connectivity
//------------------------------------------------------------------------------
static void Verify_Surface_Connectivity() {
    int      scell_lid;
    int      scell_gid;
    int      scell_type;
    int      vcell_lid;
    int      vcell_gid;
    int      vcell_type;
    int      soffset;
    int      voffset;
    int      snode[4];
    int      vnode[8];
    double   dot_prod;
    Vector3D sv1;
    Vector3D sv2;
    Vector3D snv;
    Vector3D scvc;
    Point3D  scen;
    Point3D  vcen;
    
    // Loop Over Surface Cells Only
    for (int c = 0; c < nSurfCell; c++) {
        scell_gid = c;

        // Get local cell number and type
        scell_lid = scell_gid;
        for (int etype = TRI; etype <= HEXA; etype++) {
            if (scell_gid < cellGlobalOffset[etype]) {
                scell_type = etype;
                break;
            } else {
                scell_lid -= nElem[etype];
            }
        }
        
        // Get the Nodes of Surface Cell
        soffset = scell_lid * elemNode[scell_type];
        for (int i = 0; i < elemNode[scell_type]; i++)
            snode[i] = cell2Node[scell_type][soffset + i];
        
        // Vector n0->n1
        sv1.vec[0] = coordXYZ[PHY_DIM * snode[1] + 0] - coordXYZ[PHY_DIM * snode[0] + 0];
        sv1.vec[1] = coordXYZ[PHY_DIM * snode[1] + 1] - coordXYZ[PHY_DIM * snode[0] + 1];
        sv1.vec[2] = coordXYZ[PHY_DIM * snode[1] + 2] - coordXYZ[PHY_DIM * snode[0] + 2];
        
        switch (scell_type) {
            case TRI:
                // Vector n0->n2
                sv2.vec[0] = coordXYZ[PHY_DIM * snode[2] + 0] - coordXYZ[PHY_DIM * snode[0] + 0];
                sv2.vec[1] = coordXYZ[PHY_DIM * snode[2] + 1] - coordXYZ[PHY_DIM * snode[0] + 1];
                sv2.vec[2] = coordXYZ[PHY_DIM * snode[2] + 2] - coordXYZ[PHY_DIM * snode[0] + 2];
                break;
            case QUAD:
                // Vector n0->n3
                sv2.vec[0] = coordXYZ[PHY_DIM * snode[3] + 0] - coordXYZ[PHY_DIM * snode[0] + 0];
                sv2.vec[1] = coordXYZ[PHY_DIM * snode[3] + 1] - coordXYZ[PHY_DIM * snode[0] + 1];
                sv2.vec[2] = coordXYZ[PHY_DIM * snode[3] + 2] - coordXYZ[PHY_DIM * snode[0] + 2];
                break;
        }
        
        // Get the Unit Normal Vector of surface cell
        snv = sv1.cross(sv2);
        snv.normalize();
        
        // Find centroid of surface cell
        scen = 0.0;
        for (int i = 0; i < elemNode[scell_type]; i++) {
            scen.pos[0] += coordXYZ[PHY_DIM * snode[i] + 0];
            scen.pos[0] += coordXYZ[PHY_DIM * snode[i] + 1];
            scen.pos[0] += coordXYZ[PHY_DIM * snode[i] + 2];
        }
        scen /= (double) elemNode[scell_type];
        
        // Get the Volume Cell Connected to this Surface Cell
        vcell_gid = cell2Cell[c]->list[0];
        // Get local cell number and type
        vcell_lid = vcell_gid;
        for (int etype = TRI; etype <= HEXA; etype++) {
            if (vcell_gid < cellGlobalOffset[etype]) {
                vcell_type = etype;
                break;
            } else {
                vcell_lid -= nElem[etype];
            }
        }
        
        // Get the Nodes of Volume Cell
        voffset = vcell_lid * elemNode[vcell_type];
        for (int i = 0; i < elemNode[vcell_type]; i++)
            vnode[i] = cell2Node[vcell_type][voffset + i];

        // Find centroid of connected volume cell
        vcen = 0.0;
        for (int i = 0; i < elemNode[vcell_type]; i++) {
            vcen.pos[0] += coordXYZ[PHY_DIM * vnode[i] + 0];
            vcen.pos[1] += coordXYZ[PHY_DIM * vnode[i] + 1];
            vcen.pos[2] += coordXYZ[PHY_DIM * vnode[i] + 2];
        }
        vcen /= (double) elemNode[vcell_type];

        // Get Vector Connecting Surface Centriod to Volume Centriod
        scvc.vec[0] = vcen.pos[0] - scen.pos[0];
        scvc.vec[1] = vcen.pos[1] - scen.pos[1];
        scvc.vec[2] = vcen.pos[2] - scen.pos[2];
        scvc.normalize();

        // Check if surface cells require reverse order
        // Area vector should point inside => dot product > 0
        // Else reverse the nodes
        dot_prod = scvc.dot(snv);
        if (dot_prod < 0.0) {
            switch (scell_type) {
                case TRI:
                    cell2Node[TRI][soffset + 0] = snode[2];
                    cell2Node[TRI][soffset + 1] = snode[1];
                    cell2Node[TRI][soffset + 2] = snode[0];
                    break;
                case QUAD:
                    cell2Node[QUAD][soffset + 0] = snode[3];
                    cell2Node[QUAD][soffset + 1] = snode[2];
                    cell2Node[QUAD][soffset + 2] = snode[1];
                    cell2Node[QUAD][soffset + 3] = snode[0];
                    break;
            }
        }
    } 
}

//------------------------------------------------------------------------------
//! Create All Connectivity Maps
//------------------------------------------------------------------------------
void Create_Connectivity_Maps() {
    // Calculate total number of Cells
    nCell     = 0;
    nSurfCell = 0;
    nVolCell  = 0;
    for (int etype = TRI; etype <= HEXA; etype++) {
        nCell += nElem[etype];
        if (etype < TETRA)
            nSurfCell += nElem[etype];
        else
            nVolCell += nElem[etype];
    }

    // Compute the Cell Global ID Offsets
    cellGlobalOffset[TRI] = nElem[TRI];
    for (int i = QUAD; i <= HEXA; i++)
        cellGlobalOffset[i] = cellGlobalOffset[i-1] + nElem[i];
    
    printf("=============================================================================\n");
    // Create Node2Cell Connectivity
    info("Creating Node to Cell Connectivity");
    Create_Connectivity_Node2Cell();
    
    // Create Node2Node Connectivity
    info("Creating Node to Node Connectivity");
    Create_Connectivity_Node2Node();
    
    // Create Cell2Cell Connectivity
    info("Creating Cell to Cell Connectivity");
    Create_Connectivity_Cell2Cell();
    
    // Create Edge2Node Connectivity
    info("Creating Edge to Node Connectivity");
    Create_Connectivity_Edge2Node();

    info("Verifing Surface Connectivity");
    Verify_Surface_Connectivity();

    
    for (int n = 0; n < nNode; n++) {
        delete node2Cell[n];
        delete node2Node[n];
    }
    delete[] node2Cell;
    delete[] node2Node;

    for (int c = 0; c < nCell; c++)
        delete cell2Cell[c];
    delete [] cell2Cell;
}

