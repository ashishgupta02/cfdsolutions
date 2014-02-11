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
 * File		CellPointOperation_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "CellPointOperation_2D.h"

/*--------------------------------------------------------------*/
void PointInCell_2D(double x, double y, int cellID, int *flag) {
    int i, j, node1, node2;
    int result = 0;

    *flag = 0;

    for (i = 0, j = (Cell2D[cellID].NodesPerCell - 1); i < Cell2D[cellID].NodesPerCell; j = i, i++) {
        node1 = Cell2D[cellID].ConnectNode[j];
        node2 = Cell2D[cellID].ConnectNode[i];

        /* AntiClock Wise (cgns standard) */
        if (((y - Node2D[node1].Coordinate[1])*(Node2D[node2].Coordinate[0] - Node2D[node1].Coordinate[0]) -
                (x - Node2D[node1].Coordinate[0])*(Node2D[node2].Coordinate[1] - Node2D[node1].Coordinate[1])) >= 0.0)
            result++;
    }

    if (result == Cell2D[cellID].NodesPerCell)
        *flag = 1;
}

/*--------------------------------------------------------------*/
void CellContainPoint_2D(double x, double y, int *cellID) {
    int i, flag;

    *cellID = -1;

    for (i = 0; i < NoCells2D; i++) {
        PointInCell_2D(x, y, i, &flag);
        if (flag == 1) {
            *cellID = i;
            i = NoCells2D;
        }
    }
}
