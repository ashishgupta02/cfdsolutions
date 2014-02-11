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
 * File		VectorTopology_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include "Algebra.h"
#include "Vector.h"
#include "PostProcessing_2D.h"
#include "PPInitializeSolution_2D.h"
#include "CellPointOperation_2D.h"
#include "Interpolation_2D.h"
#include "Gradient_2D.h"
#include "Stream_2D.h"
#include "Streamline_2D.h"
#include "ScalarTopology_2D.h"
#include "VectorTopology_2D.h"

/* VectorTopology Components Arg1 = X, Arg2 = Y */
int VectorTopologyArg1 = -1;
int VectorTopologyArg2 = -1;

/* Used to check change in vector args */
static int CheckVecTopArg1 = -1;
static int CheckVecTopArg2 = -1;

/* No Critical Points */
static int NCriticalPoints = 0;
static int flagCP = 0;

/* Data Structure to Store Critical Points */
CriticalPoint_2D_Un CPoint2D;

/* Store Critical Points */
// static Point2 *CriticalPoint;

/* Topology Vector Data if CellCenter */
Data_2D_Un TopVecX;
Data_2D_Un TopVecY;

/* Gradient of Vector Data */
static Data_2D_Un TopGradVecXX;
static Data_2D_Un TopGradVecXY;
static Data_2D_Un TopGradVecYX;
static Data_2D_Un TopGradVecYY;

static int firstTop = 0;

/* Tolorence */
static double tol = 0.00001;

/*--------------------------------------------------------------*/
void SetTopologicalSeedPoints_2D(void) {
    int i, j, first, cellID;
    double *distance;
    double tmp;
    Point2 test;

    distance = (double *) malloc((CPoint2D.Size) * sizeof (double));
    if (distance == NULL) {
        Warn("SetTopologicalSeedPoints_2D: Memory Allocation Failed 2");
        return;
    }

    NSeedPoint = 0;

    for (i = 0; i < CPoint2D.Size; i++) {
        first = 0;
        for (j = 0; j < CPoint2D.Size; j++) {
            if (i == j)
                continue;

            tmp = pow((CPoint2D.CriticalPoint[i].x - CPoint2D.CriticalPoint[j].x), 2)
                    + pow((CPoint2D.CriticalPoint[i].y - CPoint2D.CriticalPoint[j].y), 2);
            tmp = sqrt(tmp);

            if (!first)
                distance[i] = tmp;

            if (distance[i] > tmp)
                distance[i] = tmp;

            first = 1;
        }

        /* Set Seed Points For Critical Points */
        switch (CPoint2D.Type[i]) {
            case Saddle:
                NSeedPoint = NSeedPoint + 12;
                break;
            case Source:
                NSeedPoint = NSeedPoint + 16;
                break;
            case Sink:
                NSeedPoint = NSeedPoint + 16;
                break;
            case Center:
                NSeedPoint = NSeedPoint + 5;
                break;
            case RepellingSpiral:
                NSeedPoint = NSeedPoint + 5;
                break;
            case AttractingSpiral:
                NSeedPoint = NSeedPoint + 5;
                break;
        }
    }

    if (CPoint2D.Size == 1)
        distance[0] = 10.0;

    /* Release Previous Allocate SeedPoint */
    if (SeedPoint != NULL) free(SeedPoint);
    if (SeedPoint != NULL) free(SeedCellID);

    /* Allocate Memory */
    SeedPoint = (Point2 *) malloc(NSeedPoint * sizeof (Point2));
    SeedCellID = (int *) malloc(NSeedPoint * sizeof (int));
    if ((SeedPoint == NULL) || (SeedCellID == NULL)) {
        Warn("SetTopologicalSeedPoints_2D: Memory Allocation Failed 1");
        return;
    }

    SeedID = 0;

    for (i = 0; i < CPoint2D.Size; i++) {
        if (distance[i] > 1.0)
            distance[i] = 10.0;

        /* Set Seed Points For Critical Points */
        switch (CPoint2D.Type[i]) {
            case Saddle:
                for (j = 0; j < 6; j++) {
                    test.x = CPoint2D.CriticalPoint[i].x + ((j - 3) * distance[i] / 6);
                    test.y = CPoint2D.CriticalPoint[i].y;

                    CellContainPoint_2D(test.x, test.y, &cellID);
                    if (cellID == -1)
                        continue;

                    SeedPoint[SeedID].x = test.x;
                    SeedPoint[SeedID].y = test.y;
                    SeedCellID[SeedID] = cellID;
                    SeedID++;
                }

                for (j = 0; j < 6; j++) {
                    test.x = CPoint2D.CriticalPoint[i].x;
                    test.y = CPoint2D.CriticalPoint[i].y + ((j - 3) * distance[i] / 6);

                    CellContainPoint_2D(test.x, test.y, &cellID);
                    if (cellID == -1)
                        continue;

                    SeedPoint[SeedID].x = test.x;
                    SeedPoint[SeedID].y = test.y;
                    SeedCellID[SeedID] = cellID;
                    SeedID++;
                }
                break;
            case Source:
                for (j = 0; j < 16; j++) {
                    test.x = CPoint2D.CriticalPoint[i].x + (cos(j * (PI / 16)) * distance[i] / 2);
                    test.y = CPoint2D.CriticalPoint[i].y + (sin(j * (PI / 16)) * distance[i] / 2);

                    CellContainPoint_2D(test.x, test.y, &cellID);
                    if (cellID == -1)
                        continue;

                    SeedPoint[SeedID].x = test.x;
                    SeedPoint[SeedID].y = test.y;
                    SeedCellID[SeedID] = cellID;
                    SeedID++;
                }
                break;
            case Sink:
                for (j = 0; j < 16; j++) {
                    test.x = CPoint2D.CriticalPoint[i].x + (cos(j * (PI / 16)) * distance[i] / 2);
                    test.y = CPoint2D.CriticalPoint[i].y + (sin(j * (PI / 16)) * distance[i] / 2);

                    CellContainPoint_2D(test.x, test.y, &cellID);
                    if (cellID == -1)
                        continue;

                    SeedPoint[SeedID].x = test.x;
                    SeedPoint[SeedID].y = test.y;
                    SeedCellID[SeedID] = cellID;
                    SeedID++;
                }
                break;
            case Center:
                for (j = 0; j < 5; j++) {
                    test.x = CPoint2D.CriticalPoint[i].x;
                    test.y = CPoint2D.CriticalPoint[i].y + ((5 - j) * distance[i] / 10);

                    CellContainPoint_2D(test.x, test.y, &cellID);
                    if (cellID == -1) {
                        test.x = CPoint2D.CriticalPoint[i].x;
                        test.y = CPoint2D.CriticalPoint[i].y - ((5 - j) * distance[i] / 10);

                        CellContainPoint_2D(test.x, test.y, &cellID);
                        if (cellID == -1)
                            continue;
                    }

                    SeedPoint[SeedID].x = test.x;
                    SeedPoint[SeedID].y = test.y;
                    SeedCellID[SeedID] = cellID;
                    SeedID++;
                }
                break;
            case RepellingSpiral:
                for (j = 0; j < 5; j++) {
                    test.x = CPoint2D.CriticalPoint[i].x;
                    test.y = CPoint2D.CriticalPoint[i].y + ((5 - j) * distance[i] / 10);

                    CellContainPoint_2D(test.x, test.y, &cellID);
                    if (cellID == -1) {
                        test.x = CPoint2D.CriticalPoint[i].x;
                        test.y = CPoint2D.CriticalPoint[i].y - ((5 - j) * distance[i] / 10);

                        CellContainPoint_2D(test.x, test.y, &cellID);
                        if (cellID == -1)
                            continue;
                    }

                    SeedPoint[SeedID].x = test.x;
                    SeedPoint[SeedID].y = test.y;
                    SeedCellID[SeedID] = cellID;
                    SeedID++;
                }
                break;
            case AttractingSpiral:
                for (j = 0; j < 5; j++) {
                    test.x = CPoint2D.CriticalPoint[i].x;
                    test.y = CPoint2D.CriticalPoint[i].y + ((5 - j) * distance[i] / 10);

                    CellContainPoint_2D(test.x, test.y, &cellID);
                    if (cellID == -1) {
                        test.x = CPoint2D.CriticalPoint[i].x;
                        test.y = CPoint2D.CriticalPoint[i].y - ((5 - j) * distance[i] / 10);

                        CellContainPoint_2D(test.x, test.y, &cellID);
                        if (cellID == -1)
                            continue;
                    }

                    SeedPoint[SeedID].x = test.x;
                    SeedPoint[SeedID].y = test.y;
                    SeedCellID[SeedID] = cellID;
                    SeedID++;
                }
                break;
        }
    }

    NSeedPoint = SeedID;

    if (distance != NULL) free(distance);
}

/*--------------------------------------------------------------*/
int AddCriticalPoint_2D(int cellID, Point2 CPoint, Point2 CEigenVector1, Point2 CEigenVector2, int CType) {
    int *tmp1, *tmp2;
    Point2 *tmp3, *tmp4, *tmp5;

    if (NCriticalPoints == 0)
        return 0;

    if (NCriticalPoints == 1) {
        CPoint2D.Size = NCriticalPoints;

        /* First Time Topological Operation only */
        if (!firstTop) {
            CPoint2D.Type = (int *) malloc(10 * sizeof (int));
            CPoint2D.Cell = (int *) malloc(10 * sizeof (int));
            if ((CPoint2D.Type == NULL) || (CPoint2D.Cell == NULL)) {
                Warn("AddCriticalPoint2D: Memory Allocation Failed 1");
                return 1;
            }

            CPoint2D.CriticalPoint = (Point2 *) malloc(10 * sizeof (Point2));
            if (CPoint2D.CriticalPoint == NULL) {
                Warn("AddCriticalPoint2D: Memory Allocation Failed 2");
                return 1;
            }

            CPoint2D.EigenVector1 = (Point2 *) malloc(10 * sizeof (Point2));
            CPoint2D.EigenVector2 = (Point2 *) malloc(10 * sizeof (Point2));
            if ((CPoint2D.EigenVector1 == NULL) || (CPoint2D.EigenVector2 == NULL)) {
                Warn("AddCriticalPoint2D: Memory Allocation Failed 3");
                return 1;
            }
        }

        CPoint2D.Type[0] = CType;
        CPoint2D.Cell[0] = cellID;
        CPoint2D.CriticalPoint[0].x = CPoint.x;
        CPoint2D.CriticalPoint[0].y = CPoint.y;
        CPoint2D.EigenVector1[0].x = CEigenVector1.x;
        CPoint2D.EigenVector1[0].y = CEigenVector1.y;
        CPoint2D.EigenVector2[0].x = CEigenVector2.x;
        CPoint2D.EigenVector2[0].y = CEigenVector2.y;

        flagCP = 1;
        return 0;
    }

    if ((NCriticalPoints != CPoint2D.Size) && (NCriticalPoints > CPoint2D.Size)) {
        if (NCriticalPoints > (flagCP * 10)) {
            flagCP++;

            /* Rellocate Memory */
            tmp1 = realloc(CPoint2D.Type, (flagCP * 10) * sizeof (int));
            tmp2 = realloc(CPoint2D.Cell, (flagCP * 10) * sizeof (int));
            if ((tmp1 == NULL) || (tmp2 == NULL)) {
                Warn("AddCriticalPoint2D: Memory Allocation Failed 4");
                return 1;
            }
            CPoint2D.Type = tmp1;
            CPoint2D.Cell = tmp2;

            tmp3 = realloc(CPoint2D.CriticalPoint, (flagCP * 10) * sizeof (Point2));
            if (tmp3 == NULL) {
                Warn("AddCriticalPoint2D: Memory Allocation Failed 5");
                return 1;
            }
            CPoint2D.CriticalPoint = tmp3;

            tmp4 = realloc(CPoint2D.EigenVector1, (flagCP * 10) * sizeof (Point2));
            if (tmp4 == NULL) {
                Warn("AddCriticalPoint2D: Memory Allocation Failed 6");
                return 1;
            }
            CPoint2D.EigenVector1 = tmp4;

            tmp5 = realloc(CPoint2D.EigenVector2, (flagCP * 10) * sizeof (Point2));
            if (tmp5 == NULL) {
                Warn("AddCriticalPoint2D: Memory Allocation Failed 7");
                return 1;
            }
            CPoint2D.EigenVector2 = tmp5;
        }

        printf("flagCP = %d\n", flagCP);

        /* Add Critical Points */
        CPoint2D.Size = NCriticalPoints;

        CPoint2D.Type[CPoint2D.Size - 1] = CType;
        CPoint2D.Cell[CPoint2D.Size - 1] = cellID;
        CPoint2D.CriticalPoint[CPoint2D.Size - 1].x = CPoint.x;
        CPoint2D.CriticalPoint[CPoint2D.Size - 1].y = CPoint.y;
        CPoint2D.EigenVector1[CPoint2D.Size - 1].x = CEigenVector1.x;
        CPoint2D.EigenVector1[CPoint2D.Size - 1].y = CEigenVector1.y;
        CPoint2D.EigenVector2[CPoint2D.Size - 1].x = CEigenVector2.x;
        CPoint2D.EigenVector2[CPoint2D.Size - 1].y = CEigenVector2.y;
    }

    return 0;
}

/*--------------------------------------------------------------*/
int GetCriticalPoint_2D(int cellID, Point2 *CPoint, Point2 *CEigenVector1, Point2 *CEigenVector2, int *CType) {
    int i, node[3], flag, count;
    double vectorX[3], vectorY[3];
    double A1, B1, C1, A2, B2, C2;
    double eigen1, eigen2, discri;

    switch (Solution2D.Location) {
        case Vertex:
            for (i = 0; i < 3; i++) {
                node[i] = Cell2D[cellID].ConnectNode[i];

                vectorX[i] = Vector2DX->Data[node[i]];
                vectorY[i] = Vector2DY->Data[node[i]];
            }

            /* Linear Interpolation */
            A1 = (vectorX[0] * (Node2D[node[1]].Coordinate[1] - Node2D[node[2]].Coordinate[1]))
                    + (vectorX[1] * (Node2D[node[2]].Coordinate[1] - Node2D[node[0]].Coordinate[1]))
                    + (vectorX[2] * (Node2D[node[0]].Coordinate[1] - Node2D[node[1]].Coordinate[1]));
            A2 = (vectorY[0] * (Node2D[node[1]].Coordinate[1] - Node2D[node[2]].Coordinate[1]))
                    + (vectorY[1] * (Node2D[node[2]].Coordinate[1] - Node2D[node[0]].Coordinate[1]))
                    + (vectorY[2] * (Node2D[node[0]].Coordinate[1] - Node2D[node[1]].Coordinate[1]));

            B1 = (vectorX[0] * (Node2D[node[2]].Coordinate[0] - Node2D[node[1]].Coordinate[0]))
                    + (vectorX[1] * (Node2D[node[0]].Coordinate[0] - Node2D[node[2]].Coordinate[0]))
                    + (vectorX[2] * (Node2D[node[1]].Coordinate[0] - Node2D[node[0]].Coordinate[0]));
            B2 = (vectorY[0] * (Node2D[node[2]].Coordinate[0] - Node2D[node[1]].Coordinate[0]))
                    + (vectorY[1] * (Node2D[node[0]].Coordinate[0] - Node2D[node[2]].Coordinate[0]))
                    + (vectorY[2] * (Node2D[node[1]].Coordinate[0] - Node2D[node[0]].Coordinate[0]));

            C1 = (vectorX[0] * ((Node2D[node[1]].Coordinate[0] * Node2D[node[2]].Coordinate[1])
                    - (Node2D[node[2]].Coordinate[0] * Node2D[node[1]].Coordinate[1])))
                    + (vectorX[0] * ((Node2D[node[2]].Coordinate[0] * Node2D[node[0]].Coordinate[1])
                    - (Node2D[node[0]].Coordinate[0] * Node2D[node[2]].Coordinate[1])))
                    + (vectorX[0] * ((Node2D[node[0]].Coordinate[0] * Node2D[node[1]].Coordinate[1])
                    - (Node2D[node[1]].Coordinate[0] * Node2D[node[0]].Coordinate[1])));
            C2 = (vectorY[0] * ((Node2D[node[1]].Coordinate[0] * Node2D[node[2]].Coordinate[1])
                    - (Node2D[node[2]].Coordinate[0] * Node2D[node[1]].Coordinate[1])))
                    + (vectorY[0] * ((Node2D[node[2]].Coordinate[0] * Node2D[node[0]].Coordinate[1])
                    - (Node2D[node[0]].Coordinate[0] * Node2D[node[2]].Coordinate[1])))
                    + (vectorY[0] * ((Node2D[node[0]].Coordinate[0] * Node2D[node[1]].Coordinate[1])
                    - (Node2D[node[1]].Coordinate[0] * Node2D[node[0]].Coordinate[1])));

            if (((B1 * A2) - (A1 * B2)) == 0.0)
                return 1;

            count = 0;
            /* Setting Tolerance and Identifying the Critical Point */
            for (i = 0; i < 9; i++) {
                switch (i) {
                    case 0:
                        break;
                    case 1:
                        C1 = C1 + tol;
                        break;
                    case 2:
                        C1 = C1 - tol;
                        break;
                    case 3:
                        C2 = C2 + tol;
                        break;
                    case 4:
                        C2 = C2 - tol;
                        break;
                    case 5:
                        C1 = C1 - tol;
                        C2 = C2 - tol;
                        break;
                    case 6:
                        C1 = C1 + tol;
                        C2 = C2 + tol;
                        break;
                    case 7:
                        C1 = C1 - tol;
                        C2 = C2 + tol;
                        break;
                    case 8:
                        C1 = C1 + tol;
                        C2 = C2 - tol;
                        break;
                }

                CPoint->x = ((C1 * B2) - (B1 * C2)) / ((B1 * A2) - (A1 * B2));
                CPoint->y = ((A1 * C2) - (C1 * A2)) / ((B1 * A2) - (A1 * B2));

                PointInCell_2D(CPoint->x, CPoint->y, cellID, &flag);
                if (flag == 1)
                    i = 9;
                else {
                    count++;
                    if (count == 9)
                        return 1;
                }
            }
            break;
        case CellCenter:
            /* Talior Series Expansion */
            A1 = TopGradVecXX.Data[cellID];
            A2 = TopGradVecYX.Data[cellID];

            B1 = TopGradVecXY.Data[cellID];
            B2 = TopGradVecYY.Data[cellID];

            C1 = Vector2DX->Data[cellID] - (A1 * Cell2D[cellID].Centroid[0]) - (B1 * Cell2D[cellID].Centroid[1]);
            C2 = Vector2DY->Data[cellID] - (A2 * Cell2D[cellID].Centroid[0]) - (B2 * Cell2D[cellID].Centroid[1]);

            if (((B1 * A2) - (A1 * B2)) == 0.0)
                return 1;

            CPoint->x = ((C1 * B2) - (B1 * C2)) / ((B1 * A2) - (A1 * B2));
            CPoint->y = ((A1 * C2) - (C1 * A2)) / ((B1 * A2) - (A1 * B2));

            PointInCell_2D(CPoint->x, CPoint->y, cellID, &flag);
            if (flag != 1)
                return 1;

            break;
    }

    /* Get Critical Point Classification Type */
    discri = ((A1 + B2) * (A1 + B2)) - (4 * ((A1 * B2) - (B1 * A2)));
    *CType = 0;

    /* Real Roots */
    if (discri >= 0.0) {
        eigen1 = ((A1 + B2) + sqrt(discri)) / 2;
        eigen2 = ((A1 + B2) - sqrt(discri)) / 2;

        if ((eigen1 * eigen2) < 0.0)
            *CType = Saddle;
        else {
            if ((eigen1 > 0.0) && (eigen2 > 0.0))
                *CType = Source;
            else {
                if ((eigen1 < 0.0) && (eigen2 < 0.0))
                    *CType = Sink;
            }
        }
    }
    /* Contains Imaginary Part */
    if (discri < 0.0) {
        eigen1 = (A1 + B2) / 2;
        eigen2 = (A1 + B2) / 2;

        if ((eigen1 == 0.0) && (eigen2 == 0.0))
            *CType = Center;
        else {
            if ((eigen1 > 0.0) && (eigen2 > 0.0))
                *CType = RepellingSpiral;
            else {
                if ((eigen1 < 0.0) && (eigen2 < 0.0))
                    *CType = AttractingSpiral;
            }
        }
    }

    /* Calculate Eigen Vectors */
    switch (*CType) {
        case Saddle:
            CEigenVector1->x = 1.0;
            CEigenVector1->y = (eigen1 - A1) / B1;

            CEigenVector2->x = 1.0;
            CEigenVector2->y = (eigen2 - A1) / B1;
            break;
        case Source:
            CEigenVector1->x = 1.0;
            CEigenVector1->y = (eigen1 - A1) / B1;

            CEigenVector2->x = 1.0;
            CEigenVector2->y = (eigen2 - A1) / B1;
            break;
        case Sink:
            CEigenVector1->x = 1.0;
            CEigenVector1->y = (eigen1 - A1) / B1;

            CEigenVector2->x = 1.0;
            CEigenVector2->y = (eigen2 - A1) / B1;
            break;
        case Center:
            CEigenVector1->x = 0.0;
            CEigenVector1->y = 0.0;

            CEigenVector2->x = 0.0;
            CEigenVector2->y = 0.0;
            break;
        case RepellingSpiral:
            CEigenVector1->x = 0.0;
            CEigenVector1->y = 0.0;

            CEigenVector2->x = 0.0;
            CEigenVector2->y = 0.0;
            break;
        case AttractingSpiral:
            CEigenVector1->x = 0.0;
            CEigenVector1->y = 0.0;

            CEigenVector2->x = 0.0;
            CEigenVector2->y = 0.0;
            break;
    }

    return 0;
}

/*--------------------------------------------------------------*/
void VectorTopology_2D(void) {
    int i, CType, flag = 0;
    Point2 CPoint, CEigenVector1, CEigenVector2;

    switch (StreamModule) {
        case 2:
            /* Vector Topology */
            if ((VectorTopologyArg1 == -1) || (VectorTopologyArg2 == -1))
                return;

            /* Should not go for first time */
            if (firstTop) {
                if ((VectorTopologyArg1 != CheckVecTopArg1) || (VectorTopologyArg2 != CheckVecTopArg2))
                    flag = 1;
            }

            break;
        case 3:
            /* Scalar Topology */
            flag = 1;
            break;
    }

    /* If First time */
    if ((!firstTop) || (flag == 1)) {
        NCriticalPoints = 0;
        CPoint2D.Size = 0;

        /* If Vector Topology */
        if (StreamModule == 2) {
            switch (Solution2D.Location) {
                case Vertex:
                    /* Take Data from Solution2D Structure */

                    /* Initialize Stream Pointers */
                    Vector2DX = Solution2D.Sols[VectorTopologyArg1];
                    Vector2DY = Solution2D.Sols[VectorTopologyArg2];

                    break;
                case CellCenter:
                    /* Initialize Stream Pointers */
                    Vector2DX = Solution2D.Sols[VectorTopologyArg1];
                    Vector2DY = Solution2D.Sols[VectorTopologyArg2];

                    /* For VectorX */
                    TopGradVecXX.Size = NoCells2D;
                    TopGradVecXY.Size = NoCells2D;
                    /* Memory Allocation */
                    TopGradVecXX.Data = (double *) malloc(NoCells2D * sizeof (double));
                    TopGradVecXY.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if ((TopGradVecXX.Data == NULL) || (TopGradVecXY.Data == NULL)) {
                        Warn("VectorTopology_2D: Memory Allocation Failed: 1");
                        return;
                    }
                    /* Calculate Gradient */
                    CellGradientLSM_2D(Vector2DX, &TopGradVecXX, &TopGradVecXY);

                    /* For VectorY */
                    TopGradVecYX.Size = NoCells2D;
                    TopGradVecYY.Size = NoCells2D;
                    /* Memory Allocation */
                    TopGradVecYX.Data = (double *) malloc(NoCells2D * sizeof (double));
                    TopGradVecYY.Data = (double *) malloc(NoCells2D * sizeof (double));
                    if ((TopGradVecYX.Data == NULL) || (TopGradVecYY.Data == NULL)) {
                        Warn("VectorTopology_2D: Memory Allocation Failed: 2");
                        return;
                    }
                    /* Calculate Gradient */
                    CellGradientLSM_2D(Vector2DY, &TopGradVecYX, &TopGradVecYY);

                    /* CellCenter Data to Vertex */
                    /*
                    TopVecX.Size = NoNodes2D;
                    TopVecY.Size = NoNodes2D;
                     */
                    /* Allocate Memory */
                    /*
                    TopVecX.Data = (double *) malloc(NoNodes2D * sizeof(double));
                    TopVecY.Data = (double *) malloc(NoNodes2D * sizeof(double));
                    if ((TopVecX.Data == NULL) || (TopVecY.Data == NULL)) {
                            Warn("VectorTopology_2D: Memory Allocation Failed");
                            return;
                    }
                     */
                    /* Interpolate */
                    /*
                    ArithmeticCC2Node_2D(Solution2D.Sols[VectorTopologyArg1], &TopVecX);
                    ArithmeticCC2Node_2D(Solution2D.Sols[VectorTopologyArg2], &TopVecY);
                     */
                    /* Initialize Stream Pointers */
                    /*
                    Vector2DX = &TopVecX;
                    Vector2DY = &TopVecY;
                     */
                    break;
            }
        }            /* If Scalar Topology */
        else {
            /* Initialize Stream Pointers */
            Vector2DX = &GraVecX;
            Vector2DY = &GraVecY;
        }

        for (i = 0; i < NoCells2D; i++) {
            /* If No Critical Point */
            if (GetCriticalPoint_2D(i, &CPoint, &CEigenVector1, &CEigenVector2, &CType))
                continue;

            NCriticalPoints++;
            AddCriticalPoint_2D(i, CPoint, CEigenVector1, CEigenVector2, CType);
        }

        if (NCriticalPoints == 0) {
            printf("No Topological Structure Found\n");
            return;
        } else {
            printf("Topological Structures:\n");
            for (i = 0; i < CPoint2D.Size; i++)
                printf("%d: Cell-ID = %d\tType = %d\tX = %lf\tY = %lf\n",
                    i + 1, CPoint2D.Cell[i], CPoint2D.Type[i], CPoint2D.CriticalPoint[i].x, CPoint2D.CriticalPoint[i].y);
        }
        firstTop = 1;

        /* Set Stream Seed Points According to Topology */
        SetTopologicalSeedPoints_2D();

        printf("NCritical Points = %d\n", NCriticalPoints);

        switch (Solution2D.Location) {
            case CellCenter:
                /* CellCenter Data to Vertex */
                TopVecX.Size = NoNodes2D;
                TopVecY.Size = NoNodes2D;

                /* Allocate Memory */
                TopVecX.Data = (double *) malloc(NoNodes2D * sizeof (double));
                TopVecY.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if ((TopVecX.Data == NULL) || (TopVecY.Data == NULL)) {
                    Warn("VectorTopology_2D: Memory Allocation Failed");
                    return;
                }

                /* Interpolate */
                ArithmeticCC2Node_2D(Solution2D.Sols[VectorTopologyArg1], &TopVecX);
                ArithmeticCC2Node_2D(Solution2D.Sols[VectorTopologyArg2], &TopVecY);
                /* Initialize Stream Pointers */
                Vector2DX = &TopVecX;
                Vector2DY = &TopVecY;

                break;
        }
    }
}
