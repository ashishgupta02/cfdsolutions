/*
 * File:   TriangleMesh.cpp
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include "Utils.h"
#include "List.h"
#include "TriangleMesh.h"
#include "Grid_Utils.h"
#include "Bowyer_Watson.h"

int trimesh(int nn, int tdim, int nb, int nbs[], int*** bs,
        double x[], double y[], int tri[][3]) {
    int i, ntri;
    double Xmin, Ymin, Xmax, Ymax;
    double *BBoxX, *BBoxY;

    // Allocate Memory for Storing Cell2Cell Connectivity
    int (*cell2cell)[3];
    cell2cell = new int[tdim][3];
    
    // Data Stucture to store Node2Cell Connectivity
    List * node2cell[nn+4];
    for (i = 0; i < nn + 4; i++)
        node2cell[i] = new List();
    
    // Determine Extent of Domain
    Get_Domain_Extent(nn, x, y, &Xmin, &Ymin, &Xmax, &Ymax);
    
    // Create Bounding Box Points
    BBoxX = NULL;
    BBoxX = (double *) malloc((nn+4)*sizeof(double));
    if (BBoxX == NULL) {
        printf("ERROR: Triangle_Mesher - Memory Allocation Failure\n!");
        exit(EXIT_FAILURE);
    }
    BBoxY = NULL;
    BBoxY = (double *) malloc((nn+4)*sizeof(double));
    if (BBoxY == NULL) {
        printf("ERROR: Triangle_Mesher - Memory Allocation Failure\n!");
        exit(EXIT_FAILURE);
    }
    Create_BoundingBox_Domain(BBoxX, BBoxY, Xmin, Ymin, Xmax, Ymax);
    nn += 4;
    for (i = 4; i < nn; i++) {
        BBoxX[i] = x[i - 4];
        BBoxY[i] = y[i - 4];
    }

    // Add Bounding Box Triangle to Triangle List
    tri[0][0] = 0;
    tri[0][1] = 1;
    tri[0][2] = 2;
    tri[1][0] = 2;
    tri[1][1] = 3;
    tri[1][2] = 0;
    ntri = 2;
    
    // Create Node to Cell Connectivity - Brute Force Method
    Create_Node2Cell_Connectivity(4, ntri, tri, node2cell);

    // Remove Duplicate Entries in Node2Cell List
    for (i = 0; i < nn; i++)
        node2cell[i]->RemoveDuplicates();
    
    // Create Cell to Cell Connectivity - Brute Force Method
    Create_Cell2Cell_Connectivity(ntri, tri, cell2cell, node2cell);

    // Start Bowyer Watson Triangulation Algorithm
    int success = 0;
    success = Bowyer_Watson_Triangulation(nn, tdim, BBoxX, BBoxY,
            ntri, tri, cell2cell, node2cell);
    if (success == -1) {
        if (BBoxX != NULL)
            free(BBoxX);
        if (BBoxY != NULL)
            free(BBoxY);
        delete[] cell2cell;
        for (i = 0; i < nn; i++)
            delete node2cell[i];
        return success;
    }

    // Remove Duplicate Entries in Node2Cell List
    for (i = 0; i < nn; i++)
        node2cell[i]->RemoveDuplicates();
    
    // Increment the Boundary Node IDs to Take care of Bounding Box Nodes
    int bndry;
    for (bndry = 0; bndry < nb; bndry++) {
        for (i = 0; i < nbs[bndry]; i++) {
            bs[bndry][i][0] += 4;
            bs[bndry][i][1] += 4;
        }
    }
    
    // Start Boundary Edge Recovery
    Boundary_Edge_Recovery(nn, nb, nbs, bs, BBoxX, BBoxY,
            ntri, tri, cell2cell, node2cell);

    // Remove Duplicate Entries in Node2Cell List
    for (i = 0; i < nn; i++)
        node2cell[i]->RemoveDuplicates();

    // Start Flood Fill Algorithm to delete the external Triangles
    Flood_Fill_Algorithm(nb, nbs, bs, ntri, tri, cell2cell, node2cell);

    // Decrement the Node IDs to Take care of Bounding Box Nodes
    for (bndry = 0; bndry < nb; bndry++) {
        for (i = 0; i < nbs[bndry]; i++) {
            bs[bndry][i][0] -= 4;
            bs[bndry][i][1] -= 4;
        }
    }
    
    // Delete the Bounding Box
    Delete_BoundingBox_Domain(ntri,tri);

    // Finally Free the Memory
    for (i = 0; i < nn; i++)
        delete node2cell[i];
    delete[] cell2cell;
    if (BBoxX != NULL)
        free(BBoxX);
    if (BBoxY != NULL)
        free(BBoxY);
    return ntri;
}

