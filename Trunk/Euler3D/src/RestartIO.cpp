/*******************************************************************************
 * File:        RestartIO.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"
#include "MeshIO.h"
#include "RestartIO.h"
#include "Solver.h"
#include "Material.h"

static int RestartCounter = 0;

//------------------------------------------------------------------------------
//! Check if Restart File is Requested
//------------------------------------------------------------------------------
void Check_Restart(int Iteration) {
    if (RestartCycle > 0 && RestartOutput != 0) {
        char filename[256];
        if (RestartCounter != (RestartCycle - 1)) {
            RestartCounter++;
            return;
        } else {
            RestartCounter = 0;
            str_blank(filename);
            sprintf(filename, "Restart_%d.srt", Iteration+1);
            // Now Write the file
            Restart_Writer(filename, 0);
            // Write VTK File
            str_blank(filename);
            sprintf(filename, "Solution_%d.vtk", Iteration+1);
            VTK_Writer(filename, 0);
        }
    }
}

//------------------------------------------------------------------------------
//! Restart Solution Writer
//------------------------------------------------------------------------------
void Restart_Writer(const char* filename, int verbose) {
    int i, var;
    FILE *fp;

    if (verbose == 1) {
        printf("=============================================================================\n");
        info("Writing Restart File %s", filename);
    }
    
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
    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CONSERVATIVE) {
        for (i = 0; i < (nNode + nBNode); i++) {
            fwrite(&Q1[i], sizeof (double), 1, fp);
            fwrite(&Q2[i], sizeof (double), 1, fp);
            fwrite(&Q3[i], sizeof (double), 1, fp);
            fwrite(&Q4[i], sizeof (double), 1, fp);
            fwrite(&Q5[i], sizeof (double), 1, fp);
        }
    }
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PRIMITIVE_PUT) {
        double p, T;
        for (i = 0; i < (nNode + nBNode); i++) {
            p = Q1[i] + Gauge_Pressure;
            T = Q5[i] + Gauge_Temperature;
            fwrite(&p, sizeof (double), 1, fp);
            fwrite(&Q2[i], sizeof (double), 1, fp);
            fwrite(&Q3[i], sizeof (double), 1, fp);
            fwrite(&Q4[i], sizeof (double), 1, fp);
            fwrite(&T, sizeof (double), 1, fp);
        }
    }
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_PRIMITIVE_RUP) {
        double p;
        for (i = 0; i < (nNode + nBNode); i++) {
            p = Q5[i] + Gauge_Pressure;
            fwrite(&Q1[i], sizeof (double), 1, fp);
            fwrite(&Q2[i], sizeof (double), 1, fp);
            fwrite(&Q3[i], sizeof (double), 1, fp);
            fwrite(&Q4[i], sizeof (double), 1, fp);
            fwrite(&p, sizeof (double), 1, fp);
        }
    }
    fclose(fp);
}

//------------------------------------------------------------------------------
//! Restart Solution Reader
//------------------------------------------------------------------------------
void Restart_Reader(const char* filename) {
    int i, var;
    size_t sdum;
    FILE *fp = NULL;

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
    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CONSERVATIVE) {
        for (i = 0; i < (nNode + nBNode); i++) {
            sdum = fread(&Q1[i], sizeof (double), 1, fp);
            sdum = fread(&Q2[i], sizeof (double), 1, fp);
            sdum = fread(&Q3[i], sizeof (double), 1, fp);
            sdum = fread(&Q4[i], sizeof (double), 1, fp);
            sdum = fread(&Q5[i], sizeof (double), 1, fp);
        }
    }
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PRIMITIVE_PUT) {
        for (i = 0; i < (nNode + nBNode); i++) {
            sdum = fread(&Q1[i], sizeof (double), 1, fp);
            Q1[i] -= Gauge_Pressure;
            sdum = fread(&Q2[i], sizeof (double), 1, fp);
            sdum = fread(&Q3[i], sizeof (double), 1, fp);
            sdum = fread(&Q4[i], sizeof (double), 1, fp);
            sdum = fread(&Q5[i], sizeof (double), 1, fp);
            Q5[i] -= Gauge_Temperature;
        }
    }
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_PRIMITIVE_RUP) {
        for (i = 0; i < (nNode + nBNode); i++) {
            sdum = fread(&Q1[i], sizeof (double), 1, fp);
            sdum = fread(&Q2[i], sizeof (double), 1, fp);
            sdum = fread(&Q3[i], sizeof (double), 1, fp);
            sdum = fread(&Q4[i], sizeof (double), 1, fp);
            sdum = fread(&Q5[i], sizeof (double), 1, fp);
            Q5[i] -= Gauge_Pressure;
        }
    }
    fclose(fp);
}


