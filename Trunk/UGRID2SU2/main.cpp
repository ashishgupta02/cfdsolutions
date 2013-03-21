/*******************************************************************************
 * File:        main.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"
#include "Trim_Utils.h"
#include "List.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <cstdio>


using namespace std;

#define PHY_DIM                         3
#define CELL_DIM                        3

// Cell
#define CELL_TYPE_NONE                 -1
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

#define UNSIO_FLOWTYPE_UNKNOWN         -1
#define UNSIO_FLOWTYPE_INCOMPRESSIBLE   0
#define UNSIO_FLOWTYPE_COMPRESSIBLE     1
#define UNSIO_FLOWTYPE_VARMACHRUP       2


// Grid Connectivity Data Structure
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

    printf("=============================================================================\n");
}

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
    char buff[256];
    char *dummy = NULL;
    List lnBC;
    
    if ((fp = fopen(filename, "r")) == NULL)
        error("Ugrid_Reader: Unable to Read Grid File - %s", filename);
    
    printf("=============================================================================\n");
    info("Reading Grid File: %s", filename);
    
    // Read in the first line, the number of each element type
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d %d %d %d %d %d %d", &nNode, &nElem[TRI],
            &nElem[QUAD], &nElem[TETRA], &nElem[PYRA], &nElem[PRISM], &nElem[HEXA]);

    // Allocate Memory for Coordinates
    coordXYZ = (double*) malloc(PHY_DIM * nNode * sizeof (double));
    for (int n = 0; n < PHY_DIM * nNode; n++)
        coordXYZ[n] = 0.0;
    
    // Allocate Cell to Node Pointers for All Element Types
    for (int c = 0; c < NUMBER_OF_ELEM_TYPES; c++) {
        // Check for No of Elements
        if (nElem[c] > 0) {
            cell2Node[c] = (int*) malloc(elemNode[c] * nElem[c] * sizeof (int));
            for (int n = 0; n < nElem[c]; n++)
                cell2Node[c][n] = -1;
        } else
            cell2Node[c] = NULL;
    }
    
    // Allocate Surface tags for Triangles and Quads
    for (int c = TRI; c <= QUAD; c++) {
        // Check for No of Elements
        if (nElem[c] > 0) {
            faceTag[c] = (int*) malloc(nElem[c] * sizeof (int));
            for (int n = 0; n < nElem[c]; n++)
                faceTag[c][n] = 0;
        } else
            faceTag[c] = NULL;
    }
    
    // Read in x, y, z coordinates
    for (int n = 0; n < nNode; n++) {
        dummy = fgets(buff, bdim, fp);
        sscanf(buff, "%lg %lg %lg",
                &coordXYZ[PHY_DIM * n + 0], &coordXYZ[PHY_DIM * n + 1], &coordXYZ[PHY_DIM * n + 2]);
    }

    // Read in Triangle connectivity
    for (int c = 0; c < nElem[TRI]; c++) {
        dummy = fgets(buff, bdim, fp);
        sscanf(buff, "%d %d %d", &cell2Node[TRI][NNODE_TRI * c + 0],
                &cell2Node[TRI][NNODE_TRI * c + 1], &cell2Node[TRI][NNODE_TRI * c + 2]);
        cell2Node[TRI][NNODE_TRI * c + 0]--;
        cell2Node[TRI][NNODE_TRI * c + 1]--;
        cell2Node[TRI][NNODE_TRI * c + 2]--;
    }
    
    // Read in Quadrilateral connectivity
    for (int c = 0; c < nElem[QUAD]; c++) {
        dummy = fgets(buff, bdim, fp);
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
        dummy = fgets(buff, bdim, fp);
        sscanf(buff, "%d", &faceTag[TRI][c]);
        lnBC.Check_List(faceTag[TRI][c]);
    }
    
    // Read in Quadrilateral Boundary tags
    for (int c = 0; c < nElem[QUAD]; c++) {
        dummy = fgets(buff, bdim, fp);
        sscanf(buff, "%d", &faceTag[QUAD][c]);
        lnBC.Check_List(faceTag[QUAD][c]);
    }
    nBC = lnBC.max;
    
    // Read in Tetrahedral connectivity
    for (int c = 0; c < nElem[TETRA]; c++) {
        dummy = fgets(buff, bdim, fp);
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
        dummy = fgets(buff, bdim, fp);
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
        dummy = fgets(buff, bdim, fp);
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
        dummy = fgets(buff, bdim, fp);
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

    // Bring Grid Connectivity to Proper Ordering
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

//------------------------------------------------------------------------------
//! Stanford SU2 Writer
//------------------------------------------------------------------------------
void SU2_Writer(const char* filename) {
    int i, j, ielem;
    FILE *fp;
    List lnBC;
    
    if ((fp = fopen(filename, "w")) == NULL)
        error("SU2_Writer: Unable to Write Grid File - %s", filename);
    
    printf("=============================================================================\n");
    info("Writing Grid File: %s", filename);
    
    // Write the Physical Dimensions
    fprintf(fp, "NDIME= %d\n", PHY_DIM);
    
    // Write the Number of Volumetric Elements
    fprintf(fp, "NELEM= %d\n", nElem[TETRA]+nElem[PYRA]+nElem[PRISM]+nElem[HEXA]);
    
    ielem = 0;
    // Write Tetrahedral Elements
    for (i = 0; i < nElem[TETRA]; i++) {
        fprintf(fp, "10\t%d\t%d\t%d\t%d\t%d\n", cell2Node[TETRA][NNODE_TETRA * i + 0],
                cell2Node[TETRA][NNODE_TETRA * i + 1], cell2Node[TETRA][NNODE_TETRA * i + 2],
                cell2Node[TETRA][NNODE_TETRA * i + 3], ielem);
        ielem++;
    }
    
    // Write Hexahedral Elements
    for (i = 0; i < nElem[HEXA]; i++) {
        fprintf(fp, "12\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", cell2Node[HEXA][NNODE_HEXA * i + 0],
                cell2Node[HEXA][NNODE_HEXA * i + 1], cell2Node[HEXA][NNODE_HEXA * i + 2],
                cell2Node[HEXA][NNODE_HEXA * i + 3], cell2Node[HEXA][NNODE_HEXA * i + 4],
                cell2Node[HEXA][NNODE_HEXA * i + 5], cell2Node[HEXA][NNODE_HEXA * i + 6],
                cell2Node[HEXA][NNODE_HEXA * i + 7], ielem);
        ielem++;
    }
    
    // Write Prism Elements
    for (i = 0; i < nElem[PRISM]; i++) {
        fprintf(fp, "13\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", cell2Node[PRISM][NNODE_PRISM * i + 0],
                cell2Node[PRISM][NNODE_PRISM * i + 1], cell2Node[PRISM][NNODE_PRISM * i + 2],
                cell2Node[PRISM][NNODE_PRISM * i + 3], cell2Node[PRISM][NNODE_PRISM * i + 4],
                cell2Node[PRISM][NNODE_PRISM * i + 5], ielem);
        ielem++;
    }
    
    // Write Pyramids Elements
    for (i = 0; i < nElem[PYRA]; i++) {
        fprintf(fp, "14\t%d\t%d\t%d\t%d\t%d\t%d\n", cell2Node[PYRA][NNODE_PYRA * i + 0],
                cell2Node[PYRA][NNODE_PYRA * i + 1], cell2Node[PYRA][NNODE_PYRA * i + 2],
                cell2Node[PYRA][NNODE_PYRA * i + 3], cell2Node[PYRA][NNODE_PYRA * i + 4],
                ielem);
        ielem++;
    }
    
    // Write the Coordinates
    fprintf(fp, "NPOIN= %d\n", nNode);
    for (i = 0 ; i < nNode; i++)
        fprintf(fp, "%22.15e\t%22.15e\t%22.15e\t%d\n", coordXYZ[3 * i + 0], coordXYZ[3 * i + 1], coordXYZ[3 * i + 2], i);
    
    // Write the Boundary Conditions
    fprintf(fp, "NMARK= %d\n", nBC);
    for (i = 0; i < nElem[TRI]; i++)
        lnBC.Check_List(faceTag[TRI][i]);
    for (i = 0; i < nElem[QUAD]; i++)
        lnBC.Check_List(faceTag[QUAD][i]);
    int *nElemBC = NULL;
    nElemBC = (int*) malloc(nBC*sizeof(int));
    for (i = 0; i < nBC; i++)
        nElemBC[i] = 0;
    
    for (i = 0; i < nElem[TRI]; i++) {
        for (j = 0; j < nBC; j++) {
            if (lnBC.list[j] == faceTag[TRI][i]) {
                nElemBC[j]++;
                break;
            }
        }
    }
    for (i = 0; i < nElem[QUAD]; i++) {
        for (j = 0; j < nBC; j++) {
            if (lnBC.list[j] == faceTag[QUAD][i]) {
                nElemBC[j]++;
                break;
            }
        }
    }
    
    for (i = 0; i < nBC; i++) {
        fprintf(fp, "MARKER_TAG= SURFACE_%d\n", lnBC.list[i]);
        fprintf(fp, "MARKER_ELEMS= %d\n", nElemBC[i]);
        for (ielem = 0; ielem < nElem[TRI]; ielem++) {
            if (lnBC.list[i] == faceTag[TRI][ielem]) {
                fprintf(fp, "5\t%d\t%d\t%d\n", cell2Node[TRI][NNODE_TRI * ielem + 0],
                        cell2Node[TRI][NNODE_TRI * ielem + 1], cell2Node[TRI][NNODE_TRI * ielem + 2]);
            }
        }
        for (ielem = 0; ielem < nElem[QUAD]; ielem++) {
            if (lnBC.list[i] == faceTag[QUAD][ielem]) {
                fprintf(fp, "9\t%d\t%d\t%d\t%d\n", cell2Node[QUAD][NNODE_QUAD * ielem + 0],
                        cell2Node[QUAD][NNODE_QUAD * ielem + 1], cell2Node[QUAD][NNODE_QUAD * ielem + 2],
                        cell2Node[QUAD][NNODE_QUAD * ielem + 3]);
            }
        }
    }
    
    // Close file
    fclose(fp);
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//! Main Code
//------------------------------------------------------------------------------
int main (int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "ERROR: Too few arguments\n");
        return EXIT_FAILURE;
    }
    
    string fugrid = argv[1];
    string fbmap, fsu2;
    istringstream var(fugrid);
    getline(var, fbmap, '.');
    fsu2 = fbmap;
    fbmap += ".bc";
    fsu2  += ".su2";
    cout << fugrid << " " << fbmap << " " << fsu2 << endl;
    char gridfile[256];    

    // Initialize the Common Data
    Commons_Init();
    
    // Read UGRID input Grid File
    strcpy(gridfile, fugrid.c_str());
    UGrid_Reader(gridfile);
    
    // Write the SU2 Grid File
    strcpy(gridfile, fsu2.c_str());
    SU2_Writer(gridfile);
    // Finalize the Common Data
    Commons_Finalize();
    return 0;
}
