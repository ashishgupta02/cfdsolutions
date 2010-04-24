/*
 * File:   Bowyer_Watson.cpp
 * Author: Ashish Gupta
 *
 * Created on February 5, 2010, 7:22 PM
 */

#include "Grid_Utils.h"
#include "Bowyer_Watson.h"

// *****************************************************************************
// Bower Watson Algorithm to Create Triangulation into 2D Rectangular Domain
// Using Generic 2 Triangles generated and Seed Points.
// *****************************************************************************
int Bowyer_Watson_Triangulation(int nn, int tdim, double x[], double y[],
        int &ntri, int tri[][3], int cell2cell[][3], List* node2cell[]) {
    int i, j, k, t, n, n0, n1;
    int seed, edgeCount, nbrTri;

    // Intialize
    n0 = -1;
    n1 = -1;
    
    // Create Map Array use for Triangle Compression
    int* map = NULL;
    map = (int *) malloc(tdim*sizeof(int));
    if (map == NULL) {
        printf("ERROR: Boweyer_Watson_Triangulation - Memory Allocation Failure\n!");
        exit(EXIT_FAILURE);
    }
    
    // Create Intitial Memory to Store List of Voilated Triangles
    // Size of 60 is used for optimization purpose only
    List *violated_tri;
    violated_tri = new List();
    violated_tri->Redimension(60);

    // Create Intitial Memory to Store List of New Triangles
    // Size of 60 is used for optimization purpose only
    List *new_tri;
    new_tri = new List();
    new_tri->Redimension(60);
    
    // Create Initial Memory to Store Edges of Convex Hull Created
    // Size of 180 is used for optimization purpose only
    // Node: Size of SavedEdges Should be approximately 3 times the Number of Voilated Triangles
    int savedEdgesSize = 10;
    int (*savedEdges)[2];
    savedEdges = new int[savedEdgesSize][2];

#ifdef VERBOSE
    printf("START: Bowyer Watson Triangulation\n");
    printf("BOWYER_WATSON: Found (%d) Number of Points for Insertion\n", nn-4);
    printf("BOWYER_WATSON: ");
    fflush(stdout);
    int progress1, progress2 = -1;
#endif
    // Start inserting Non Bonding Box Seed Points
    for (n = 4; n < nn; n++) 
    {
#ifdef VERBOSE
        // Print Progress
        if (((n-4)%100) == 0) {
            progress1 = (int)(((double)n)/((double)nn)*100.0);
            if (progress1 != progress2) {
                progress2 = progress1;
                if ((progress1%10) == 0) {
                    printf(". ");
                    fflush(stdout);
                    //printf("BOWYER_WATSON: Triangulation (%d) Percent Complete\n", progress1);
                }
            }
        }
#endif
        // Find Triangle Containing Seed Point in the Triangle Mesh
        seed = 0;
        seed = Point_In_Triangle_Mesh(seed, x[n], y[n], x, y, tri, cell2cell);
        if (seed < 0) {
            printf("ERROR: Seed Point => %d\t = (%lf, %lf) Not Found in Triangle Mesh. Aborting !",
                    n-4, x[n], y[n]);
            exit(EXIT_FAILURE);
        }

        // Initialize List Maintaining the List of Violated Triangles, Reset for Memory Reuse
        violated_tri->Reset(0);
        
        // Seed Triangle: Add Neighbours Voilating the circle test and Add to Delete Triangle List
        Find_Neighbour_Triangles_Violating_CircleTest(seed, n, x, y, tri, cell2cell, violated_tri);

        // Create list of Saved Edges
        // Check Triangles to be deleted, save outer edges that will form convex hull
        if (savedEdgesSize < (3 * violated_tri->max)) {
            delete[] savedEdges;
            savedEdges = new int[3 * violated_tri->max][2];
            savedEdgesSize = 3 * violated_tri->max;
        }
        edgeCount = 0;
        for (i = 0; i < violated_tri->max; i++) {
            t = violated_tri->list[i];
            for (j = 0; j < 3; j++) {
                // Check if neighbor exists and is not in delete list, save that edge
                nbrTri = cell2cell[t][j];
                if ((nbrTri == -1) || (!(violated_tri->Is_In_List(nbrTri)))) {
                    switch (j) {
                        case 0:
                            n0 = tri[t][0];
                            n1 = tri[t][1];
                            break;
                        case 1:
                            n0 = tri[t][1];
                            n1 = tri[t][2];
                            break;
                        case 2:
                            n0 = tri[t][2];
                            n1 = tri[t][0];
                            break;
                    }
                    savedEdges[edgeCount][0] = n0;
                    savedEdges[edgeCount][1] = n1;
                    edgeCount++;
                }
            }
        }

        // Delete violated triangles from Node2Cell Connectivity
        // Create Original Triangle map
        int tri_d, node;
        for (t = 0; t < ntri; t++)
            map[t] = t;
        for (i = 0; i < violated_tri->max; i++) {
            tri_d = violated_tri->list[i];
            map[tri_d] = -1;
            // Remove Voilated Triangles from Node to Cell Connectivity
            for (j = 0; j < 3; j++) {
                node = tri[tri_d][j];
                node2cell[node]->Delete_From_List(tri_d);
            }
        }
        // Start Updatign the Connectivity Table removing Voilated triangles
        // Note: Only Node2Cell Connectivity is updated completly
        //       Cell2Cell Connectivity is updated separatly below
        j = 0;
        for (i = 0; i < ntri; i++) {
            if (map[i] != -1) {
                if (i != j) {
                    // Update Connectivity Table with new Location (j)
                    for (k = 0; k < 3; k++) {
                        node2cell[tri[i][k]]->Replace(i, j);
                        tri[j][k] = tri[i][k];
                        cell2cell[j][k] = cell2cell[i][k];
                    }
                    map[i] = j;
                }
                j++;
            }
        }

        // Update the Cell2Cell Connectivity Triangles IDs
        for (i = 0; i < ntri; i++) {
            for (k = 0; k < 3; k++) {
                if (cell2cell[i][k] != -1) {
                    cell2cell[i][k] = map[cell2cell[i][k]];
                }
            }
        }
        ntri = ntri - violated_tri->max;

        // Check the memory usage :
        // Total Initial triangles >= Created Triangles + New Triangles
        if (tdim < ntri + edgeCount) {
            if (map != NULL)
                free(map);
            delete[] savedEdges;
            delete violated_tri;
            delete new_tri;
            return -1;
        }

        // Initialize List Maintaining the List of New Triangles, Reset for Memory Reuse
        new_tri->Reset(0);
        
        // Create new triangles out of saved edges and seed node
        for (i = 0; i < edgeCount; i++) {
            tri[ntri][0] = n;
            tri[ntri][1] = savedEdges[i][0];
            tri[ntri][2] = savedEdges[i][1];
            new_tri->Add_To_List(ntri);
            node2cell[n]->Add_To_List(ntri);
            node2cell[savedEdges[i][0]]->Add_To_List(ntri);
            node2cell[savedEdges[i][1]]->Add_To_List(ntri);
            ntri++;
        }
        
        // Update Cell2Cell Connectivity with New Triangles
        Update_Cell2Cell_Connectivity(tri, cell2cell, node2cell, new_tri);
        
    }

#ifdef VERBOSE
    printf("\nBOWYER_WATSON: Created (%d) Triangles with (%d) Seed Points\n", ntri, n-4);
    printf("END: Bowyer Watson Triangulation\n");
#endif

    // Free Used Memory
    if (map != NULL)
        free(map);
    delete violated_tri;
    delete new_tri;
    delete[] savedEdges;

    return 0;
}
