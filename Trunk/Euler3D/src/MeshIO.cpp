/*******************************************************************************
 * File:        MeshIO.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"
#include "MeshIO.h"

//------------------------------------------------------------------------------
//! UGRID Grid Connectivity Ordering
//------------------------------------------------------------------------------
static void UGrid_Translate_Winding(void) {
    static int translation[NUMBER_OF_ELEM_TYPES][MAX_ELEM_NNODE] = {
        {2, 1, 0},                  // Triangle
        {3, 2, 1, 0},               // Quadrilateral
        {0, 1, 2, 3},               // Tetrahedral
        {0, 3, 4, 1, 2},            // Pyramid
        {0, 1, 2, 3, 4, 5},         // Prism
        {0, 1, 2, 3, 4, 5, 6, 7}    // Hexahedral
    };

    static int oldorder[MAX_ELEM_NNODE];
    
    for (int etype = TRI; etype <= HEXA; etype++) {
        for (int c = 0; c < nElem[etype]; c++) {
            // Copy the Old Ordering
            for (int n = 0; n < elemNode[etype]; n++)
                oldorder[n] = cell2Node[etype][elemNode[etype] * c + n];
            // Update the New Ordering
            for (int n = 0; n < elemNode[etype]; n++)
                cell2Node[etype][elemNode[etype] * c + translation[etype][n]] = oldorder[n];
        }
    }
}


//------------------------------------------------------------------------------
//! Simcenter UGRID Reader
//------------------------------------------------------------------------------
void UGrid_Reader(const char* filename) {
    FILE *fp;
    int bdim = 256;
    char buff[bdim];
    
    if ((fp = fopen(filename, "r")) == NULL)
        error("Ugrid_Reader: Unable to Read Grid File - %s", filename);
    
    printf("=============================================================================\n");
    info("Reading Grid File: %s", filename);
    
    // Read in the first line, the number of each element type
    fgets(buff, bdim, fp);
    sscanf(buff, "%d %d %d %d %d %d %d", &nNode, &nElem[TRI],
            &nElem[QUAD], &nElem[TETRA], &nElem[PYRA], &nElem[PRISM], &nElem[HEXA]);

    // Allocate Memory for Coordinates
    coordXYZ = (double*) malloc(PHY_DIM * nNode * sizeof (double));
    for (int n = 0; n < PHY_DIM * nNode; n++)
        coordXYZ[n] = 0.0;
    
    // Allocate Cell to Node Pointers for All Element Types
    for (int c = 0; c < NUMBER_OF_ELEM_TYPES; c++) {
        cell2Node[c] = (int*) malloc(elemNode[c] * nElem[c] * sizeof (int));
        for (int n = 0; n < nElem[c]; n++)
            cell2Node[c][n] = -1;
    }
    
    // Allocate Surface tags for Triangles and Quads
    for (int c = TRI; c <= QUAD; c++) {
        faceTag[c] = (int*) malloc(nElem[c] * sizeof (int));
        for (int n = 0; n < nElem[c]; n++)
            faceTag[c][n] = 0;
    }
    
    // Read in x, y, z coordinates
    for (int n = 0; n < nNode; n++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%lg %lg %lg",
                &coordXYZ[PHY_DIM * n + 0], &coordXYZ[PHY_DIM * n + 1], &coordXYZ[PHY_DIM * n + 2]);
    }

    // Read in Triangle connectivity
    for (int c = 0; c < nElem[TRI]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d %d %d", &cell2Node[TRI][NNODE_TRI * c + 0],
                &cell2Node[TRI][NNODE_TRI * c + 1], &cell2Node[TRI][NNODE_TRI * c + 2]);
        cell2Node[TRI][NNODE_TRI * c + 0]--;
        cell2Node[TRI][NNODE_TRI * c + 1]--;
        cell2Node[TRI][NNODE_TRI * c + 2]--;
    }
    
    // Read in Quadrilatral connectivity
    for (int c = 0; c < nElem[QUAD]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d %d %d %d",
                &cell2Node[QUAD][NNODE_QUAD * c + 0], &cell2Node[QUAD][NNODE_QUAD * c + 1],
                &cell2Node[QUAD][NNODE_QUAD * c + 2], &cell2Node[QUAD][NNODE_QUAD * c + 3]);
        cell2Node[QUAD][NNODE_QUAD * c + 0]--;
        cell2Node[QUAD][NNODE_QUAD * c + 1]--;
        cell2Node[QUAD][NNODE_QUAD * c + 2]--;
        cell2Node[QUAD][NNODE_QUAD * c + 3]--;   
    }

    nBC = 0;
    // Read in Triangle Boundary tags
    for (int c = 0; c < nElem[TRI]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d", &faceTag[TRI][c]);
        nBC = MAX(nBC, faceTag[TRI][c]);
    }
    
    // Read in Quadrilateral Boundary tags
    for (int c = 0; c < nElem[QUAD]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d", &faceTag[QUAD][c]);
        nBC = MAX(nBC, faceTag[QUAD][c]);
    }
    
    // Read in Tetrahedral connectivity
    for (int c = 0; c < nElem[TETRA]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d %d %d %d",
                &cell2Node[TETRA][NNODE_TETRA * c + 0], &cell2Node[TETRA][NNODE_TETRA * c + 1],
                &cell2Node[TETRA][NNODE_TETRA * c + 2], &cell2Node[TETRA][NNODE_TETRA * c + 3]);
        cell2Node[TETRA][NNODE_TETRA * c + 0]--;
        cell2Node[TETRA][NNODE_TETRA * c + 1]--;
        cell2Node[TETRA][NNODE_TETRA * c + 2]--;
        cell2Node[TETRA][NNODE_TETRA * c + 3]--;
    }
    
    // Read in Pyramid connectivity
    for (int c = 0; c < nElem[PYRA]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d %d %d %d %d", &cell2Node[PYRA][NNODE_PYRA * c + 0],
                &cell2Node[PYRA][NNODE_PYRA * c + 1], &cell2Node[PYRA][NNODE_PYRA * c + 2],
                &cell2Node[PYRA][NNODE_PYRA * c + 3], &cell2Node[PYRA][NNODE_PYRA * c + 4]);
        cell2Node[PYRA][NNODE_PYRA * c + 0]--;
        cell2Node[PYRA][NNODE_PYRA * c + 1]--;
        cell2Node[PYRA][NNODE_PYRA * c + 2]--;
        cell2Node[PYRA][NNODE_PYRA * c + 3]--;
        cell2Node[PYRA][NNODE_PYRA * c + 4]--;
    }

    // Read in Prism connectivity
    for (int c = 0; c < nElem[PRISM]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d %d %d %d %d %d",
                &cell2Node[PRISM][6 * c + 0], &cell2Node[PRISM][NNODE_PRISM * c + 1],
                &cell2Node[PRISM][6 * c + 2], &cell2Node[PRISM][NNODE_PRISM * c + 3],
                &cell2Node[PRISM][6 * c + 4], &cell2Node[PRISM][NNODE_PRISM * c + 5]);
        cell2Node[PRISM][NNODE_PRISM * c + 0]--;
        cell2Node[PRISM][NNODE_PRISM * c + 1]--;
        cell2Node[PRISM][NNODE_PRISM * c + 2]--;
        cell2Node[PRISM][NNODE_PRISM * c + 3]--;
        cell2Node[PRISM][NNODE_PRISM * c + 4]--;
        cell2Node[PRISM][NNODE_PRISM * c + 5]--;
    }

    // Read in Hexahedral connectivity
    for (int c = 0; c < nElem[HEXA]; c++) {
        fgets(buff, bdim, fp);
        sscanf(buff, "%d %d %d %d %d %d %d %d",
                &cell2Node[HEXA][8 * c + 0], &cell2Node[HEXA][NNODE_HEXA * c + 1],
                &cell2Node[HEXA][8 * c + 2], &cell2Node[HEXA][NNODE_HEXA * c + 3],
                &cell2Node[HEXA][8 * c + 4], &cell2Node[HEXA][NNODE_HEXA * c + 5],
                &cell2Node[HEXA][8 * c + 6], &cell2Node[HEXA][NNODE_HEXA * c + 7]);
        cell2Node[HEXA][NNODE_HEXA * c + 0]--;
        cell2Node[HEXA][NNODE_HEXA * c + 1]--;
        cell2Node[HEXA][NNODE_HEXA * c + 2]--;
        cell2Node[HEXA][NNODE_HEXA * c + 3]--;
        cell2Node[HEXA][NNODE_HEXA * c + 4]--;
        cell2Node[HEXA][NNODE_HEXA * c + 5]--;
        cell2Node[HEXA][NNODE_HEXA * c + 6]--;
        cell2Node[HEXA][NNODE_HEXA * c + 7]--;
    }
    
    // Close file
    fclose(fp);

    // Bring Grid Connectivity to Proper Odering
    UGrid_Translate_Winding();
    
    // Print Out Mesh Output
    printf("=============================================================================\n");
    info("No of Nodes           : %d", nNode);
    info("No of Triangle        : %d", nElem[TRI]);
    info("No of Quadrilateral   : %d", nElem[QUAD]);
    info("No of Tetrahedral     : %d", nElem[TETRA]);
    info("No of Pyramid         : %d", nElem[PYRA]);
    info("No of Prism           : %d", nElem[PRISM]);
    info("No of Hexahedral      : %d", nElem[HEXA]);
    info("No of Surface Cells   : %d", nElem[TRI] + nElem[QUAD]);
    info("No of Volume Cells    : %d", nElem[TETRA] + nElem[PYRA] + nElem[PRISM] + nElem[HEXA]);
    info("No of Boundaries      : %d", nBC);
}

