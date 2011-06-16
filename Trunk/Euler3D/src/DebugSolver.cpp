/*******************************************************************************
 * File:        DebugSolver.h
 * Author:      Ashish Gupta
 * Revision:    2
 ******************************************************************************/

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
    
    for (i = 0; i < nNode; i++) {
        tmp = Q1[i] + Q2[i] + Q3[i] + Q4[i] + Q5[i];
        if (!isfinite(tmp))
            DebugNode.Add_To_List(i);
        tmp = Res1[i] + Res2[i] + Res3[i] + Res4[i] + Res5[i];
        if (!isfinite(tmp))
            DebugNode.Add_To_List(i);
    }
    DebugNode.RemoveDuplicates();

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