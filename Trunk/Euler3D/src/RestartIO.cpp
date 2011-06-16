/*******************************************************************************
 * File:        MeshIO.cpp
 * Author:      Ashish Gupta
 * Revision:    2
 ******************************************************************************/

#include <string.h>

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"
#include "RestartIO.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! Restart Solution Writer
//------------------------------------------------------------------------------
void Restart_Writer(const char* filename) {
    int i, var;
    FILE *fp;

    printf("=============================================================================\n");
    info("Writing Restart File %s", filename);

    if ((fp = fopen(filename, "wb")) == NULL)
        error("Restart_Writer: Unable to Write Solution File - %s", filename);

    // Restart Iteration
    fwrite(&RestartIteration, sizeof(int), 1, fp);

    // No of Physical Nodes
    fwrite(&nNode, sizeof(int), 1, fp);

    // No of Boundary/Ghost Nodes
    fwrite(&nBNode, sizeof(int), 1, fp);

    // No of Variables
    var = 5;
    fwrite(&var, sizeof(int), 1, fp);

    // Write the restart variables
    for (i = 0; i < (nNode + nBNode); i++) {
        fwrite(&Q1[i], sizeof (double), 1, fp);
        fwrite(&Q2[i], sizeof (double), 1, fp);
        fwrite(&Q3[i], sizeof (double), 1, fp);
        fwrite(&Q4[i], sizeof (double), 1, fp);
        fwrite(&Q5[i], sizeof (double), 1, fp);
    }

    fclose(fp);
}

//------------------------------------------------------------------------------
//! Restart Solution Reader
//------------------------------------------------------------------------------
void Restart_Reader(const char* filename) {
    int i, var;
    size_t sdum;
    FILE *fp;

    printf("=============================================================================\n");
    info("Reading Restart File %s", filename);

    if ((fp = fopen(filename, "rb")) == NULL)
        error("Restart_Reader: Unable to Open Solution File - %s", filename);

    // Restart Iteration
    sdum = fread(&RestartIteration, sizeof(int), 1, fp);

    // No of Physical Nodes
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    if (var != nNode)
        error("Restart_Reader: Mismatched in Restart and Mesh File %s - 1", filename);

    // No of Boundary/Ghost Nodes
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    if (var != nBNode)
        error("Restart_Reader: Mismatched in Restart and Mesh File %s - 2", filename);

    // No of Variables
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);

    // Read the restart variables
    for (i = 0; i < (nNode + nBNode); i++) {
        sdum = fread(&Q1[i], sizeof (double), 1, fp);
        sdum = fread(&Q2[i], sizeof (double), 1, fp);
        sdum = fread(&Q3[i], sizeof (double), 1, fp);
        sdum = fread(&Q4[i], sizeof (double), 1, fp);
        sdum = fread(&Q5[i], sizeof (double), 1, fp);
    }

    fclose(fp);
}


