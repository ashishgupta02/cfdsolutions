/*******************************************************************************
 * File:        MeshIO.cpp
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#include <string.h>

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"
#include "MeshIO.h"
#include "Solver.h"

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
    
    // Read in Quadrilatral connectivity
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
        nBC = MAX(nBC, faceTag[TRI][c]);
    }
    
    // Read in Quadrilateral Boundary tags
    for (int c = 0; c < nElem[QUAD]; c++) {
        dummy = fgets(buff, bdim, fp);
        sscanf(buff, "%d", &faceTag[QUAD][c]);
        nBC = MAX(nBC, faceTag[QUAD][c]);
    }
    
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

//------------------------------------------------------------------------------
//! VTK Solution Writer
//------------------------------------------------------------------------------
void VTK_Writer(const char* filename, int verbose) {
    int i, j;
    FILE *fp;

    if (verbose == 1)
        info("Writing VTK Solution File %s", filename);
    
    if ((fp = fopen(filename, "w")) == NULL)
        error("VTK_Writer: Unable to Write Solution File - %s", filename);
    
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Solver solution data\n");
    fprintf(fp, "ASCII \n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d double\n", nNode);
    for (i = 0; i < nNode; i++)
        fprintf(fp, "%f %f %f\n", coordXYZ[3 * i + 0], coordXYZ[3 * i + 1], coordXYZ[3 * i + 2]);


    /************PRINT CELLS************/
    int nelem_sum = 0;
    int nelem_nodes = 0;

    for (i = TRI; i <= HEXA; i++) nelem_sum += nElem[i];

    for (i = TRI; i <= HEXA; i++) nelem_nodes += nElem[i]*(elemNode[i] + 1);

    fprintf(fp, "CELLS %d %d\n", nelem_sum, nelem_nodes);

    for (i = TRI; i <= HEXA; i++) {
        for (j = 0; j < nElem[i]; j++) {
            if (i == TRI)fprintf(fp, "3 %d %d %d\n", cell2Node[TRI][3 * j + 0],
                    cell2Node[TRI][3 * j + 1], cell2Node[TRI][3 * j + 2]);
            if (i == QUAD)fprintf(fp, "4 %d %d %d %d\n",
                    cell2Node[QUAD][4 * j + 0], cell2Node[QUAD][4 * j + 1],
                    cell2Node[QUAD][4 * j + 2],cell2Node[QUAD][4 * j + 3]);
            if (i == TETRA)fprintf(fp, "4 %d %d %d %d\n",
                    cell2Node[TETRA][4 * j + 0], cell2Node[TETRA][4 * j + 1],
                    cell2Node[TETRA][4 * j + 2], cell2Node[TETRA][4 * j + 3]);
            if (i == PYRA)fprintf(fp, "5 %d %d %d %d %d\n",
                    cell2Node[PYRA][5 * j + 0], cell2Node[PYRA][5 * j + 1],
                    cell2Node[PYRA][5 * j + 2], cell2Node[PYRA][5 * j + 3],
                    cell2Node[PYRA][5 * j + 4]);
            if (i == PRISM)fprintf(fp, "6 %d %d %d %d %d %d\n",
                    cell2Node[PRISM][6 * j + 0], cell2Node[PRISM][6 * j + 1],
                    cell2Node[PRISM][6 * j + 2], cell2Node[PRISM][6 * j + 3],
                    cell2Node[PRISM][6 * j + 4], cell2Node[PRISM][6 * j + 5]);
            if (i == HEXA)fprintf(fp, "8 %d %d %d %d %d %d %d %d\n",
                    cell2Node[HEXA][8 * j + 0], cell2Node[HEXA][8 * j + 1],
                    cell2Node[HEXA][8 * j + 2], cell2Node[HEXA][8 * j + 3],
                    cell2Node[HEXA][8 * j + 4], cell2Node[HEXA][8 * j + 5],
                    cell2Node[HEXA][8 * j + 6], cell2Node[HEXA][8 * j + 7]);
        }
    }

    /************PRINT CELLS TYPES************/
    fprintf(fp, "CELL_TYPES %d\n", nelem_sum);
    for (i = TRI; i <= HEXA; i++) {
        for (j = 0; j < nElem[i]; j++) {
            if (i == TRI)fprintf(fp, "5\n");
            if (i == QUAD)fprintf(fp, "9\n");
            if (i == TETRA)fprintf(fp, "10\n");
            if (i == PYRA)fprintf(fp, "14\n");
            if (i == PRISM)fprintf(fp, "13\n");
            if (i == HEXA)fprintf(fp, "12\n");
        }
    }

    /***************DEFINE BOUNDARIES***************/
    fprintf(fp, "CELL_DATA %d\n", nelem_sum);
    fprintf(fp, "SCALARS BoundaryIds Int 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    int count = 0;

    for (i = TRI; i <= HEXA; i++) {
        if (i == TRI || i == QUAD) {
            for (j = 0; j < nElem[i]; j++) {
                fprintf(fp, "%d\n", faceTag[i][j]);
                count++;
            };
        }

        if (i != TRI && i != QUAD)
            for (j = 0; j < nElem[i]; j++)
                fprintf(fp, "0\n");

    }
    
    /*****PRINT OUT VARIABLES**************/
    fprintf(fp, "POINT_DATA %d\n", nNode);
    fprintf(fp, "SCALARS Density double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < nNode; i++)
        fprintf(fp, "%22.15e\n", Q1[i]);


    fprintf(fp, "SCALARS X_Velocity double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < nNode; i++)
        fprintf(fp, "%22.15e\n", Q2[i] / Q1[i]);


    fprintf(fp, "SCALARS Y_Velocity double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < nNode; i++)
        fprintf(fp, "%22.15e\n", Q3[i] / Q1[i]);

    fprintf(fp, "SCALARS Z_Velocity double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < nNode; i++)
        fprintf(fp, "%22.15e\n", Q4[i] / Q1[i]);

    double p, rho, et, u, v, w;

    fprintf(fp, "SCALARS Pressure double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < nNode; i++) {
        rho = Q1[i];
        u   = Q2[i] / rho;
        v   = Q3[i] / rho;
        w   = Q4[i] / rho;
        et  = Q5[i] / rho;
        p   = (Gamma - 1.0) * rho * (et - 0.5 * (u * u + v * v + w * w));

        fprintf(fp, "%22.15e\n", p);
    }

    fclose(fp);
}

