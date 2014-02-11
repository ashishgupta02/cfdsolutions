/*[
 * Copyright 2006   Ashish Gupta
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the complete modified source code.  Modifications are to
 * be distributed as patches to the released version.  Permission to
 * distribute binaries produced by compiling modified sources is granted,
 * provided you
 *   1. distribute the corresponding source modifications from the
 *    released version in the form of a patch file along with the binaries,
 *   2. add special version identification to distinguish your version
 *    in addition to the base release version number,
 *   3. provide your name and address as the primary contact for the
 *    support of your modified version, and
 *   4. retain our contact information in regard to use of the base
 *    software.
 * Permission to distribute the released version of the source code along
 * with corresponding source modifications in the form of a patch file is
 * granted with same provisions 2 through 4 for binary distributions.
 *
 * This software is provided "as is" without express or implied warranty
 * to the extent permitted by applicable law.
]*/

/*
 * File		Gradient_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "PostProcessing_2D.h"
#include "Gradient_2D.h"
#include "Algebra.h"

/*--------------------------------------------------------------*/
static void FATALS(char *errmsg) {
    char cmd[129];

    sprintf(cmd, "error_exit {%s}", errmsg);
    printf("%s\n", cmd);
    exit(1);
}

/*---------------------------------------------------------------*/

/* Calculation of gradient at Cell Center using Weigthed Area and Normal Method
        Input Data: Cell-Center Value
        Output Data: Gradient at Cell Center 
 */
int CellGradientWANM_2D(Data_2D_Un *In, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    if (!(In->Size == NoCells2D))
        FATALS("CellGradientWANM_2D: Data Size Mismatched");

    if (OutX->Data == NULL)
        FATALS("CellGradientWANM_2D: OutX Memory Failed");

    if (OutY->Data == NULL)
        FATALS("CellGradientWANM_2D: OutY Memory Failed");


    return 0;
}

/*---------------------------------------------------------------*/

/* Calculation of gradient at Cell Center using Weigthed Area and Normal Method
        Input Data: Node Value
        Output Data: Gradient at Node
 */
int NodeGradientWANM_2D(Data_2D_Un *In, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    return 0;
}

/*---------------------------------------------------------------*/

/* Calculation of gradient at Cell Center using Least Square Method 
        Input Data: Cell-Center Value
        Output Data: Gradient at Cell Center
 */
int CellGradientLSM_2D(Data_2D_Un *In, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    int i, icell, ncells, CurCell;
    double *deltaX;
    double *deltaY;
    double *deltaData;
    double Var[2];
    double Xo, Yo;
    Matrix2 M, InvM;

    Var[0] = 0.0;
    Var[1] = 0.0;

    if (!(In->Size == NoCells2D))
        FATALS("CellGradientLSM_2D: Data Size Mismatched");

    if (OutX->Data == NULL)
        FATALS("CellGradientLSM_2D: OutX Memory Failed");

    if (OutY->Data == NULL)
        FATALS("CellGradientLSM_2D: OutY Memory Failed");

    for (icell = 0; icell < NoCells2D; icell++) {
        OutX->Data[icell] = 0.0;
        OutY->Data[icell] = 0.0;
    }

    for (icell = 0; icell < NoCells2D; icell++) {
        ncells = Cell2D[icell].NearCells;

        deltaX = (double *) malloc(ncells * sizeof (double));
        if (deltaX == NULL)
            FATALS("CellGradientLSM_2D: deltaX Memory Allocation Failed");

        deltaY = (double *) malloc(ncells * sizeof (double));
        if (deltaY == NULL)
            FATALS("CellGradientLSM_2D: deltaY Memory Allocation Failed");

        deltaData = (double *) malloc(ncells * sizeof (double));
        if (deltaData == NULL)
            FATALS("CellGradientLSM_2D: deltaData Memory Allocation Failed");

        Xo = Cell2D[icell].Centroid[0];
        Yo = Cell2D[icell].Centroid[1];

        for (i = 0; i < ncells; i++) {
            CurCell = Cell2D[icell].NearNeighbours[i];

            deltaX[i] = Cell2D[CurCell].Centroid[0] - Xo;
            deltaY[i] = Cell2D[CurCell].Centroid[1] - Yo;
            deltaData[i] = In->Data[CurCell] - In->Data[icell];
        }

        for (i = 0; i < 2; i++) {
            M.element[i][0] = 0.0;
            M.element[i][1] = 0.0;

            InvM.element[i][0] = 0.0;
            InvM.element[i][1] = 0.0;
        }

        for (i = 0; i < ncells; i++) {
            M.element[0][0] += (deltaX[i] * deltaX[i]);
            M.element[0][1] += (deltaX[i] * deltaY[i]);
            M.element[1][0] += (deltaY[i] * deltaX[i]);
            M.element[1][1] += (deltaY[i] * deltaY[i]);
        }

        if (Matrix2_Inverse(&M, &InvM) == 0.0)
            FATALS("CellGradient_2D: Matrix Determinant Value = 0.0");

        for (i = 0; i < ncells; i++) {
            Var[0] += (deltaData[i] * deltaX[i]);
            Var[1] += (deltaData[i] * deltaY[i]);
        }

        OutX->Data[icell] = (InvM.element[0][0] * Var[0])
                + (InvM.element[0][1] * Var[1]);
        OutY->Data[icell] = (InvM.element[1][0] * Var[0])
                + (InvM.element[1][1] * Var[1]);

        free(deltaX);
        free(deltaY);
        free(deltaData);
    }

    return 0;
}

/*---------------------------------------------------------------*/

/* Calculation of gradient at Node using Least Square Method 
        Input Data: Nodal Value
        Output Data: Gradient at Node
 */
int NodeGradientLSM_2D(Data_2D_Un *In, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    int i, inode, nnodes, CurNode;
    double *deltaX;
    double *deltaY;
    double *deltaData;
    double Var[2];
    double Xo, Yo;
    Matrix2 M, InvM;

    Var[0] = 0.0;
    Var[1] = 0.0;

    if (!(In->Size == NoNodes2D))
        FATALS("NodeGradientLSM_2D: Data Size Mismatched");

    if (OutX->Data == NULL)
        FATALS("NodeGradientLSM_2D: OutX Memory Failed");

    if (OutY->Data == NULL)
        FATALS("NodeGradientLSM_2D: OutY Memory Failed");

    for (inode = 0; inode < NoNodes2D; inode++) {
        OutX->Data[inode] = 0.0;
        OutY->Data[inode] = 0.0;
    }

    for (inode = 0; inode < NoNodes2D; inode++) {
        nnodes = NodeAdjacentNodes2D[inode].Size;

        deltaX = (double *) malloc(nnodes * sizeof (double));
        if (deltaX == NULL)
            FATALS("NodeGradientLSM_2D: deltaX Memory Allocation Failed");

        deltaY = (double *) malloc(nnodes * sizeof (double));
        if (deltaY == NULL)
            FATALS("NodeGradientLSM_2D: deltaY Memory Allocation Failed");

        deltaData = (double *) malloc(nnodes * sizeof (double));
        if (deltaData == NULL)
            FATALS("NodeGradientLSM_2D: deltaData Memory Allocation Failed");

        Xo = Node2D[inode].Coordinate[0];
        Yo = Node2D[inode].Coordinate[1];

        for (i = 0; i < nnodes; i++) {
            CurNode = NodeAdjacentNodes2D[inode].Data[i];

            deltaX[i] = Node2D[CurNode].Coordinate[0] - Xo;
            deltaY[i] = Node2D[CurNode].Coordinate[1] - Yo;
            deltaData[i] = In->Data[CurNode] - In->Data[inode];
        }

        for (i = 0; i < 2; i++) {
            M.element[i][0] = 0.0;
            M.element[i][1] = 0.0;

            InvM.element[i][0] = 0.0;
            InvM.element[i][1] = 0.0;
        }

        for (i = 0; i < nnodes; i++) {
            M.element[0][0] += (deltaX[i] * deltaX[i]);
            M.element[0][1] += (deltaX[i] * deltaY[i]);
            M.element[1][0] += (deltaY[i] * deltaX[i]);
            M.element[1][1] += (deltaY[i] * deltaY[i]);
        }

        if (Matrix2_Inverse(&M, &InvM) == 0.0)
            FATALS("CellGradient_2D: Matrix Determinant Value = 0.0");

        for (i = 0; i < nnodes; i++) {
            Var[0] += (deltaData[i] * deltaX[i]);
            Var[1] += (deltaData[i] * deltaY[i]);
        }

        OutX->Data[inode] = (InvM.element[0][0] * Var[0])
                + (InvM.element[0][1] * Var[1]);
        OutY->Data[inode] = (InvM.element[1][0] * Var[0])
                + (InvM.element[1][1] * Var[1]);

        free(deltaX);
        free(deltaY);
        free(deltaData);
    }

    return 0;
}

/*---------------------------------------------------------------*/

/* Calculation of gradient at Cell Center 
        Method: Line Integral Counter-Clock Wise
        Input Data: Nodal Value
        Output Data: Gradient at Cell Center
 */
int CellGradientLI_2D(Data_2D_Un *In, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    int nodesPerCell;
    int inode, icell, node[4];
    double x[4], y[4];
    double deltaY[4], deltaX[4], AvgData[4];

    if (!(In->Size == NoNodes2D))
        FATALS("CellGradientLI_2D: Data Size Mismatched");

    if (OutX->Data == NULL)
        FATALS("CellGradientLI_2D: OutX Memory Failed");

    if (OutY->Data == NULL)
        FATALS("CellGradientLI_2D: OutY Memory Failed");

    for (icell = 0; icell < NoCells2D; icell++) {
        OutX->Data[icell] = 0.0;
        OutY->Data[icell] = 0.0;
    }

    for (icell = 0; icell < NoCells2D; icell++) {
        nodesPerCell = Cell2D[icell].NodesPerCell;

        for (inode = 0; inode < nodesPerCell; inode++) {
            node[inode] = Cell2D[icell].ConnectNode[inode];
            x[inode] = Node2D[node[inode]].Coordinate[0];
            y[inode] = Node2D[node[inode]].Coordinate[1];
        }

        if (nodesPerCell == 3) {
            deltaY[0] = y[1] - y[0];
            deltaY[1] = y[2] - y[1];
            deltaY[2] = y[0] - y[2];

            deltaX[0] = x[1] - x[0];
            deltaX[1] = x[2] - x[1];
            deltaX[2] = x[0] - x[2];

            AvgData[0] = (In->Data[node[0]] + In->Data[node[1]]) / 2;
            AvgData[1] = (In->Data[node[1]] + In->Data[node[2]]) / 2;
            AvgData[2] = (In->Data[node[2]] + In->Data[node[0]]) / 2;
        } else {
            deltaY[0] = y[1] - y[0];
            deltaY[1] = y[2] - y[1];
            deltaY[2] = y[3] - y[2];
            deltaY[3] = y[0] - y[3];

            deltaX[0] = x[0] - x[1];
            deltaX[1] = x[1] - x[2];
            deltaX[2] = x[2] - x[3];
            deltaX[3] = x[3] - x[0];

            AvgData[0] = (In->Data[node[0]] + In->Data[node[1]]) / 2;
            AvgData[1] = (In->Data[node[1]] + In->Data[node[2]]) / 2;
            AvgData[2] = (In->Data[node[2]] + In->Data[node[3]]) / 2;
            AvgData[3] = (In->Data[node[3]] + In->Data[node[0]]) / 2;
        }

        for (inode = 0; inode < nodesPerCell; inode++) {
            OutX->Data[icell] += AvgData[inode] * deltaY[inode];
            OutY->Data[icell] += AvgData[inode] * deltaX[inode];
        }

        OutX->Data[icell] = OutX->Data[icell] / Cell2D[icell].Area;
        OutY->Data[icell] = OutY->Data[icell] / Cell2D[icell].Area;
    }

    return 0;
}

/*---------------------------------------------------------------*/

/* Calculation of gradient at Node (Weighted Area)
        Input Data: Cell Center Gradient
        Output Data: Gradient at Node
 */
int NodeGradientLI_2D(Data_2D_Un *InX, Data_2D_Un *InY, Data_2D_Un *OutX, Data_2D_Un *OutY) {
    int cellID;
    int inode, icell;

    if (!(InX->Size == NoCells2D))
        FATALS("NodeGradientLI_2D: Data Size Mismatched: 1");

    if (!(InY->Size == NoCells2D))
        FATALS("NodeGradientLI_2D: Data Size Mismatched: 2");

    if (OutX->Data == NULL)
        FATALS("NodeGradientLI_2D: OutX Memory Failed");

    if (OutY->Data == NULL)
        FATALS("NodeGradientLI_2D: OutY Memory Failed");

    for (inode = 0; inode < NoNodes2D; inode++) {
        Node2D[inode].WeightTotal = 0.0;

        OutX->Data[inode] = 0.0;
        OutY->Data[inode] = 0.0;

        for (icell = 0; icell < NodeAdjacentCells2D[inode].Size; icell++) {
            cellID = NodeAdjacentCells2D[inode].Data[icell];

            Node2D[inode].WeightTotal += Cell2D[cellID].Area;

            OutX->Data[inode] += InX->Data[icell] * Cell2D[cellID].Area;
            OutY->Data[inode] += InY->Data[icell] * Cell2D[cellID].Area;
        }

        OutX->Data[inode] = OutX->Data[inode] / Node2D[inode].WeightTotal;
        OutY->Data[inode] = OutY->Data[inode] / Node2D[inode].WeightTotal;
    }

    return 0;
}
