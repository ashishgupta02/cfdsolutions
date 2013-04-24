/*******************************************************************************
 * File:        DebugSolver.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "List.h"
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void DebugNAN(void) {
    int i;
    double tmp;
    List DebugNode;
    List DebugBNode;
    FILE *fp;
    
    // Physical Nodes
    for (i = 0; i < nNode; i++) {
        tmp = Q1[i] + Q2[i] + Q3[i] + Q4[i] + Q5[i];
        if (!isfinite(tmp))
            DebugNode.Add_To_List(i);
        tmp = Res1_Conv[i] + Res2_Conv[i] + Res3_Conv[i] + Res4_Conv[i] + Res5_Conv[i];
        if (!isfinite(tmp))
            DebugNode.Add_To_List(i);
        tmp = Res1_Diss[i] + Res2_Diss[i] + Res3_Diss[i] + Res4_Diss[i] + Res5_Diss[i];
        if (!isfinite(tmp))
            DebugNode.Add_To_List(i);
    }
    DebugNode.RemoveDuplicates();

    // Ghost Nodes
    for (i = nNode; i < (nNode + nBNode); i++) {
        tmp = Q1[i] + Q2[i] + Q3[i] + Q4[i] + Q5[i];
        if (!isfinite(tmp))
            DebugBNode.Add_To_List(i);
    }
    DebugBNode.RemoveDuplicates();

    if ((fp = fopen("Debug.txt", "w")) == NULL)
        error("DebugNAN: Unable to Write Debug File - %s", "Debug.txt");
    DebugNode.print(fp);
    fprintf(fp, "\nBoundary Ghost Nodes \n");
    DebugBNode.print(fp);
    fprintf(fp, "\n");
    fclose(fp);
}