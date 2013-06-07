/*******************************************************************************
 * File:        RestartIO.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
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
    int i, inode, var;
    FILE *fp = NULL;
    double IOConversion[NEQUATIONS];
    double Q[NEQUATIONS];
    
    // Initialize
    for (i = 0; i < NEQUATIONS; i++)
        IOConversion[i] = 1.0;
    
    if (verbose == 1) {
        printf("=============================================================================\n");
        info("Writing Restart File %s", filename);
    }
    
    if ((fp = fopen(filename, "wb")) == NULL)
        error("Restart_Writer:1: Unable to Write Solution File - %s", filename);

    // Restart Iteration
    fwrite(&RestartIteration, sizeof(int), 1, fp);

    // No of Physical Nodes
    fwrite(&nNode, sizeof(int), 1, fp);

    // No of Boundary/Ghost Nodes
    fwrite(&nBNode, sizeof(int), 1, fp);

    // No of Variables
    var = NEQUATIONS;
    fwrite(&var, sizeof(int), 1, fp);
    
    // Check the Material Computation Type and Get Appropriate Conversion
    switch (MaterialCompType) {
        case MATERIAL_COMP_TYPE_NDIM:
            if (RestartOutputIOType == MATERIAL_COMP_TYPE_DIM) {
                switch (VariableType) {
                    case VARIABLE_CON:
                        IOConversion[0] = Ref_Rho;
                        IOConversion[1] = Ref_Rho*Ref_Velocity;
                        IOConversion[2] = Ref_Rho*Ref_Velocity;
                        IOConversion[3] = Ref_Rho*Ref_Velocity;
                        IOConversion[4] = Ref_Rho*Ref_TotalEnergy;
                        break;
                    case VARIABLE_RUP:
                        IOConversion[0] = Ref_Rho;
                        IOConversion[1] = Ref_Velocity;
                        IOConversion[2] = Ref_Velocity;
                        IOConversion[3] = Ref_Velocity;
                        IOConversion[4] = Ref_Pressure;
                        break;
                    case VARIABLE_PUT:
                        IOConversion[0] = Ref_Pressure;
                        IOConversion[1] = Ref_Velocity;
                        IOConversion[2] = Ref_Velocity;
                        IOConversion[3] = Ref_Velocity;
                        IOConversion[4] = Ref_Temperature;
                        break;
                    case VARIABLE_RUT:
                        IOConversion[0] = Ref_Rho;
                        IOConversion[1] = Ref_Velocity;
                        IOConversion[2] = Ref_Velocity;
                        IOConversion[3] = Ref_Velocity;
                        IOConversion[4] = Ref_Temperature;
                        break;
                    default:
                        error("Restart_Writer:2: Undefined or Unsupported Variable Type - %d", VariableType);
                        break;
                }
            }
            break;
        case MATERIAL_COMP_TYPE_DIM:
            if (RestartOutputIOType == MATERIAL_COMP_TYPE_NDIM) {
                switch (VariableType) {
                    case VARIABLE_CON:
                        IOConversion[0] = 1.0/Ref_Rho;
                        IOConversion[1] = 1.0/(Ref_Rho*Ref_Velocity);
                        IOConversion[2] = 1.0/(Ref_Rho*Ref_Velocity);
                        IOConversion[3] = 1.0/(Ref_Rho*Ref_Velocity);
                        IOConversion[4] = 1.0/(Ref_Rho*Ref_TotalEnergy);
                        break;
                    case VARIABLE_RUP:
                        IOConversion[0] = 1.0/Ref_Rho;
                        IOConversion[1] = 1.0/Ref_Velocity;
                        IOConversion[2] = 1.0/Ref_Velocity;
                        IOConversion[3] = 1.0/Ref_Velocity;
                        IOConversion[4] = 1.0/Ref_Pressure;
                        break;
                    case VARIABLE_PUT:
                        IOConversion[0] = 1.0/Ref_Pressure;
                        IOConversion[1] = 1.0/Ref_Velocity;
                        IOConversion[2] = 1.0/Ref_Velocity;
                        IOConversion[3] = 1.0/Ref_Velocity;
                        IOConversion[4] = 1.0/Ref_Temperature;
                        break;
                    case VARIABLE_RUT:
                        IOConversion[0] = 1.0/Ref_Rho;
                        IOConversion[1] = 1.0/Ref_Velocity;
                        IOConversion[2] = 1.0/Ref_Velocity;
                        IOConversion[3] = 1.0/Ref_Velocity;
                        IOConversion[4] = 1.0/Ref_Temperature;
                        break;
                    default:
                        error("Restart_Writer:3: Undefined or Unsupported Variable Type - %d", VariableType);
                        break;
                }
            }
            break;
        default:
            error("Restart_Writer:4: Undefined Material Computation Type - %d", MaterialCompType);
            break;
    }
    
    // Write the restart variables
    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CON) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            Q[0] = Q1[inode]*IOConversion[0];
            Q[1] = Q2[inode]*IOConversion[1];
            Q[2] = Q3[inode]*IOConversion[2];
            Q[3] = Q4[inode]*IOConversion[3];
            Q[4] = Q5[inode]*IOConversion[4];
            for (i = 0; i < NEQUATIONS; i++)
                fwrite(&Q[i], sizeof (double), 1, fp);
        }
    }
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_RUP) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            Q[0] = Q1[inode]*IOConversion[0];
            Q[1] = Q2[inode]*IOConversion[1];
            Q[2] = Q3[inode]*IOConversion[2];
            Q[3] = Q4[inode]*IOConversion[3];
            Q[4] = (Q5[inode] + Gauge_Pressure)*IOConversion[4];
            for (i = 0; i < NEQUATIONS; i++)
                fwrite(&Q[i], sizeof (double), 1, fp);
        }
    }
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PUT) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            Q[0] = (Q1[inode] + Gauge_Pressure)*IOConversion[0];
            Q[1] = Q2[inode]*IOConversion[1];
            Q[2] = Q3[inode]*IOConversion[2];
            Q[3] = Q4[inode]*IOConversion[3];
            Q[4] = Q5[inode]*IOConversion[4];
            for (i = 0; i < NEQUATIONS; i++)
                fwrite(&Q[i], sizeof (double), 1, fp);
        }
    }
    // Primitive Variable Formulation Density Velocity Temperature
    if (VariableType == VARIABLE_RUT) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            Q[0] = Q1[inode]*IOConversion[0];
            Q[1] = Q2[inode]*IOConversion[1];
            Q[2] = Q3[inode]*IOConversion[2];
            Q[3] = Q4[inode]*IOConversion[3];
            Q[4] = Q5[inode]*IOConversion[4];
            for (i = 0; i < NEQUATIONS; i++)
                fwrite(&Q[i], sizeof (double), 1, fp);
        }
    }
    fclose(fp);
}

//------------------------------------------------------------------------------
//! Restart Solution Reader
//------------------------------------------------------------------------------
void Restart_Reader(const char* filename) {
    int i, inode, var;
    size_t sdum;
    FILE *fp = NULL;
    double IOConversion[NEQUATIONS];
    double Q[NEQUATIONS];
    
    // Initialize
    for (i = 0; i < NEQUATIONS; i++)
        IOConversion[i] = 1.0;
    
    printf("=============================================================================\n");
    info("Reading Restart File %s", filename);

    if ((fp = fopen(filename, "rb")) == NULL)
        error("Restart_Reader:1: Unable to Open Solution File - %s", filename);

    // Restart Iteration
    sdum = fread(&RestartIteration, sizeof(int), 1, fp);
    
    // No of Physical Nodes
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    if (var != nNode)
        error("Restart_Reader:2: Mismatched in Restart and Mesh File %s - 1", filename);

    // No of Boundary/Ghost Nodes
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    if (var != nBNode)
        error("Restart_Reader:3: Mismatched in Restart and Mesh File %s - 2", filename);

    // No of Variables
    var = 0;
    sdum = fread(&var, sizeof(int), 1, fp);
    
    // Check the Material Computation Type and Get Appropriate Conversion
    switch (MaterialCompType) {
        case MATERIAL_COMP_TYPE_NDIM:
            if (RestartInputIOType == MATERIAL_COMP_TYPE_DIM) {
                switch (VariableType) {
                    case VARIABLE_CON:
                        IOConversion[0] = 1.0/Ref_Rho;
                        IOConversion[1] = 1.0/(Ref_Rho*Ref_Velocity);
                        IOConversion[2] = 1.0/(Ref_Rho*Ref_Velocity);
                        IOConversion[3] = 1.0/(Ref_Rho*Ref_Velocity);
                        IOConversion[4] = 1.0/(Ref_Rho*Ref_TotalEnergy);
                        break;
                    case VARIABLE_RUP:
                        IOConversion[0] = 1.0/Ref_Rho;
                        IOConversion[1] = 1.0/Ref_Velocity;
                        IOConversion[2] = 1.0/Ref_Velocity;
                        IOConversion[3] = 1.0/Ref_Velocity;
                        IOConversion[4] = 1.0/Ref_Pressure;
                        break;
                    case VARIABLE_PUT:
                        IOConversion[0] = 1.0/Ref_Pressure;
                        IOConversion[1] = 1.0/Ref_Velocity;
                        IOConversion[2] = 1.0/Ref_Velocity;
                        IOConversion[3] = 1.0/Ref_Velocity;
                        IOConversion[4] = 1.0/Ref_Temperature;
                        break;
                    case VARIABLE_RUT:
                        IOConversion[0] = 1.0/Ref_Rho;
                        IOConversion[1] = 1.0/Ref_Velocity;
                        IOConversion[2] = 1.0/Ref_Velocity;
                        IOConversion[3] = 1.0/Ref_Velocity;
                        IOConversion[4] = 1.0/Ref_Temperature;
                        break;
                    default:
                        error("Restart_Reader:4: Undefined or Unsupported Variable Type - %d", VariableType);
                        break;
                }
            }
            break;
        case MATERIAL_COMP_TYPE_DIM:
            if (RestartInputIOType == MATERIAL_COMP_TYPE_NDIM) {
                switch (VariableType) {
                    case VARIABLE_CON:
                        IOConversion[0] = Ref_Rho;
                        IOConversion[1] = Ref_Rho*Ref_Velocity;
                        IOConversion[2] = Ref_Rho*Ref_Velocity;
                        IOConversion[3] = Ref_Rho*Ref_Velocity;
                        IOConversion[4] = Ref_Rho*Ref_TotalEnergy;
                        break;
                    case VARIABLE_RUP:
                        IOConversion[0] = Ref_Rho;
                        IOConversion[1] = Ref_Velocity;
                        IOConversion[2] = Ref_Velocity;
                        IOConversion[3] = Ref_Velocity;
                        IOConversion[4] = Ref_Pressure;
                        break;
                    case VARIABLE_PUT:
                        IOConversion[0] = Ref_Pressure;
                        IOConversion[1] = Ref_Velocity;
                        IOConversion[2] = Ref_Velocity;
                        IOConversion[3] = Ref_Velocity;
                        IOConversion[4] = Ref_Temperature;
                        break;
                    case VARIABLE_RUT:
                        IOConversion[0] = Ref_Rho;
                        IOConversion[1] = Ref_Velocity;
                        IOConversion[2] = Ref_Velocity;
                        IOConversion[3] = Ref_Velocity;
                        IOConversion[4] = Ref_Temperature;
                        break;
                    default:
                        error("Restart_Reader:5: Undefined or Unsupported Variable Type - %d", VariableType);
                        break;
                }
            }
            break;
        default:
            error("Restart_Reader:6: Undefined Material Computation Type - %d", MaterialCompType);
            break;
    }
  
    // Read the restart variables
    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CON) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            for (i = 0; i < NEQUATIONS; i++)
                sdum = fread(&Q[i], sizeof (double), 1, fp);
            
            Q1[inode] = Q[0]*IOConversion[0];
            Q2[inode] = Q[1]*IOConversion[1];
            Q3[inode] = Q[2]*IOConversion[2];
            Q4[inode] = Q[3]*IOConversion[3];
            Q5[inode] = Q[4]*IOConversion[4];
        }
    }
    // Primitive Variable Formulation Density Velocity Pressure
    warn("Restart_Reader: Genearalise This Fix");
    if (VariableType == VARIABLE_RUP) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            for (i = 0; i < NEQUATIONS; i++)
                sdum = fread(&Q[i], sizeof (double), 1, fp);
            
            Q1[inode] = Q[0]*IOConversion[0];
            Q2[inode] = Q[1]*IOConversion[1];
            Q3[inode] = Q[2]*IOConversion[2];
            Q4[inode] = Q[3]*IOConversion[3];
            Q[4] = Q[4]*Inf_Pressure;
            Q5[inode] = (Q[4] - Gauge_Pressure)*IOConversion[4];
        }
    }
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PUT) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            for (i = 0; i < NEQUATIONS; i++)
                sdum = fread(&Q[i], sizeof (double), 1, fp);
            
            Q1[inode] = (Q[0] - Gauge_Pressure)*IOConversion[0];
            Q2[inode] = Q[1]*IOConversion[1];
            Q3[inode] = Q[2]*IOConversion[2];
            Q4[inode] = Q[3]*IOConversion[3];
            Q5[inode] = Q[4]*IOConversion[4];
        }
    }
    // Primitive Variable Formulation Density Velocity Temperature
    if (VariableType == VARIABLE_RUT) {
        for (inode = 0; inode < (nNode + nBNode); inode++) {
            for (i = 0; i < NEQUATIONS; i++)
                sdum = fread(&Q[i], sizeof (double), 1, fp);
            
            Q1[inode] = Q[0]*IOConversion[0];
            Q2[inode] = Q[1]*IOConversion[1];
            Q3[inode] = Q[2]*IOConversion[2];
            Q4[inode] = Q[3]*IOConversion[3];
            Q5[inode] = Q[4]*IOConversion[4];
        }
    }
    fclose(fp);
    
//    if (MaterialType == MATERIAL_TYPE_NIST) {
//        if (MaterialCompType == MATERIAL_COMP_TYPE_NDIM) {
//            int node_L, node_R;
//            for (int iBEdge = 0; iBEdge < nBEdge; iBEdge++) {
//                node_L = bndEdge[iBEdge].node[0];
//                node_R = bndEdge[iBEdge].node[1];
//                Q1[node_R] = Q1[node_L];
//                Q2[node_R] = Q2[node_L];
//                Q3[node_R] = Q3[node_L];
//                Q4[node_R] = Q4[node_L];
//                Q5[node_R] = Q5[node_L];
//            }
//        }
//    }
    
    sdum++; // Done to remove compiler un-used warning
}

