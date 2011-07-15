/*******************************************************************************
 * File:        Cuthill_Mckee_Reorder.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "List.h"
#include "Commons.h"

//------------------------------------------------------------------------------
//! Reordering of Graph to reduce Matrix-Computation Bandwidth
//------------------------------------------------------------------------------
void Cuthill_Mckee_Reorder() {
    int *Q;
    int *R;
    int *S;
    int *newId;
    int Rdeg_max;
    int temp;
    int card_S_ini;
    
    // Allocate memory to store new node ids
    newId = (int*) malloc(nNode * sizeof (int));
    for (int i = 0; i < nNode; i++)
        newId[i] = -1;

    // Allocate memory and store old node ids
    Q = (int*) malloc(nNode * sizeof (int));
    for (int i = 0; i < nNode; i++)
        Q[i] = i;
    R = (int*) malloc(100 * sizeof (int));
    S = (int*) malloc(nNode * sizeof (int));

    // These variables hold number of elements in Q and R
    int card_Q, card_S, card_R;
    card_Q = nNode;
    card_S = 0;
    card_R = 0;
    
    /*Find seed and take it from Q*/
    //Here we are going to use as seed the coordinate with the minimum value of x
    int seed;
    seed = 2;
    double min_x = DBL_MAX;
    for (int i = 0; i < nNode; i++) {
        min_x = MIN(min_x, coordXYZ[PHY_DIM * i + 0]);
        if (fabs(min_x - coordXYZ[3 * i + 0]) < 10*DBL_ZERO)
            seed = i;
    }

    //This is the index which will be used to define out ordered list newId[i];
    int index;
    index = 0;

    // This defines initial newId
    newId[0] = seed;

    // Find the postition of our seed in the Q array and put a -1 in its place.
    // Decrement the cardnality of the Q array
    Q[seed] = -1;
    card_Q--;

    // While Q is not empty, reorder the vertices
    while (card_Q > 0) {
        //Step 1
        // This finds all unused vertices in the adj(z) intersect Q and puts them in R
        // R <- adj(z) intersect Q
        card_R = 0;
        for (int i = 0; i < node2Node[seed]->max; i++) {
            for (int j = 0; j < nNode; j++) {
                if (Q[j] == node2Node[seed]->list[i]) {
                    R[card_R] = node2Node[seed]->list[i];
                    card_R++;
                }
            }
        }

        // Step 2
        // Q <- Q - R
        for (int i = 0; i < nNode; i++) {
            if (Q[i] == seed)
                Q[i] = -1;

            for (int j = 0; j < card_R; j++) {
                if (Q[i] == R[j]) {
                    Q[i] = -1;
                    card_Q--;
                }
            }
        }

        // Step 3: Sort vertices in R by degree
        for (int i = 0; i < card_R; i++) {
            Rdeg_max = -1;
            for (int j = i; j < card_R; j++) {
                Rdeg_max = MAX(Rdeg_max, node2Node[R[j]]->max);
                if (Rdeg_max <= node2Node[R[j]]->max) {
                    temp = R[i];
                    R[i] = R[j];
                    R[j] = temp;
                }
            }
        }        
        card_S_ini = card_S;
        card_S += card_R;

        // Step 4: Fill S
        for (int i = 0; i < card_R; i++) {
            index++;
            newId[index] = R[i];
            S[card_S_ini + i] = R[i];
        }

        seed = S[0];
        
        // Move the 0th element off Queue and all nodes up
        for (int i = 1; i < card_S; i++)
            S[i - 1] = S[i];
        card_S--;
    }
    
    // Free Some memory
    free(Q);
    free(R);
    free(S);

    // Translate cell2Node Connectivity
    for (int e = TRI; e <= HEXA; e++) {
        for (int n = 0; n < elemNode[e]*nElem[e]; n++)
            cell2Node[e][n] = newId[cell2Node[e][n]];
    }
    
    // Translate coordXYZ
    int xyz_trans;
    double *xyz_copy;
    xyz_copy = (double*) calloc(3 * nNode, sizeof (double));
    for (int i = 0; i < PHY_DIM*nNode; i++)
        xyz_copy[i] = coordXYZ[i];
    
    for (int i = 0; i < nNode; i++) {
        xyz_trans = newId[ i ];
        coordXYZ[PHY_DIM * xyz_trans + 0] = xyz_copy[PHY_DIM * i + 0];
        coordXYZ[PHY_DIM * xyz_trans + 1] = xyz_copy[PHY_DIM * i + 1];
        coordXYZ[PHY_DIM * xyz_trans + 2] = xyz_copy[PHY_DIM * i + 2];
    }

    // Free Memory
    free(newId);
    free(xyz_copy);
}

