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
 * File		PrePostProcessing_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include "Vector.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PrePostProcessing_2D.h"

/*--------------------------------------------------------------*/
int MinMaxEdge_2D(double *MinEdge, double *MaxEdge) {
    int i;

    *MinEdge = Edge2D[0].length;
    *MaxEdge = Edge2D[0].length;

    for (i = 1; i < NoEdges2D; i++) {
        if (*MinEdge > Edge2D[i].length)
            *MinEdge = Edge2D[i].length;

        if (*MaxEdge < Edge2D[i].length)
            *MaxEdge = Edge2D[i].length;
    }

    return 0;
}

/*********************************************************************
PolyCentroid_2D: Calculates the centroid (xCentroid, yCentroid) and area
of a polygon, given its vertices (x[0], y[0]) ... (x[n-1], y[n-1]). It
is assumed that the contour is closed, i.e., that the vertex following
(x[n-1], y[n-1]) is (x[0], y[0]).  The algebraic sign of the area is
positive for counterclockwise ordering of vertices in x-y plane;
otherwise negative.

Returned values:  0 for normal execution;  1 if the polygon is
degenerate (number of vertices < 3);  and 2 if area = 0 (and the
centroid is undefined).
 **********************************************************************/

/*---------------------------------------------------------------*/
int PolyCentroid_2D(double x[], double y[], int n,
        double *xCentroid, double *yCentroid, double *area) {

    register int i, j;
    double ai, atmp = 0, xtmp = 0, ytmp = 0;

    if (n < 3) return 1;
    for (i = n - 1, j = 0; j < n; i = j, j++) {
        ai = x[i] * y[j] - x[j] * y[i];
        atmp += ai;
        xtmp += (x[j] + x[i]) * ai;
        ytmp += (y[j] + y[i]) * ai;
    }

    *area = atmp / 2;
    if (atmp != 0) {
        *xCentroid = xtmp / (3 * atmp);
        *yCentroid = ytmp / (3 * atmp);
        return 0;
    }
    return 2;
}

/*---------------------------------------------------------------*/
int CounterClockCellNodes_2D(void) {
    int i, j, count, icell, nnodes;
    int Node[4];
    int Next[2];
    int edge[4];
    Point2 A, B, C, D;

    for (icell = 0; icell < NoCells2D; icell++) {
        nnodes = Cell2D[icell].NodesPerCell;
        count = 0;
        switch (nnodes) {
            case 3:
                /* Triangular Element */
                Node[0] = Cell2D[icell].ConnectNode[0];
                Node[1] = Cell2D[icell].ConnectNode[1];
                Node[2] = Cell2D[icell].ConnectNode[2];

                A.x = Node2D[Node[0]].Coordinate[0];
                A.y = Node2D[Node[0]].Coordinate[1];

                B.x = Node2D[Node[1]].Coordinate[0];
                B.y = Node2D[Node[1]].Coordinate[1];

                C.x = Node2D[Node[2]].Coordinate[0];
                C.y = Node2D[Node[2]].Coordinate[1];

                if (CrossO_2D(&A, &B, &C) < 0.0) {
                    Cell2D[icell].ConnectNode[1] = Node[2];
                    Cell2D[icell].ConnectNode[2] = Node[1];
                }

                break;
            case 4:
                /* Quadrilateral Element */
                Node[0] = Cell2D[icell].ConnectNode[0];
                Node[1] = Cell2D[icell].ConnectNode[1];
                Node[2] = Cell2D[icell].ConnectNode[2];
                Node[3] = Cell2D[icell].ConnectNode[3];

                edge[0] = Cell2D[icell].ConnectEdge[0];
                edge[1] = Cell2D[icell].ConnectEdge[1];
                edge[2] = Cell2D[icell].ConnectEdge[2];
                edge[3] = Cell2D[icell].ConnectEdge[3];

                for (i = 0; i < 4; i++) {
                    if ((Edge2D[edge[i]].ConnectedNodes[0] == Node[0]) ||
                            (Edge2D[edge[i]].ConnectedNodes[1] == Node[0])) {
                        for (j = 1; j < 4; j++) {
                            if (Edge2D[edge[i]].ConnectedNodes[0] == Node[j]) {
                                Next[count] = j;
                                count++;
                            } else if (Edge2D[edge[i]].ConnectedNodes[1] == Node[j]) {
                                Next[count] = j;
                                count++;
                            }
                        }
                    }
                }

                A.x = Node2D[Node[0]].Coordinate[0];
                A.y = Node2D[Node[0]].Coordinate[1];

                B.x = Node2D[Node[Next[0]]].Coordinate[0];
                B.y = Node2D[Node[Next[0]]].Coordinate[1];

                C.x = Node2D[Node[Next[1]]].Coordinate[0];
                C.y = Node2D[Node[Next[1]]].Coordinate[1];

                D.x = Cell2D[icell].Centroid[0];
                D.y = Cell2D[icell].Centroid[1];

                if (CrossO_2D(&D, &A, &B) < 0.0) {
                    if (CrossO_2D(&D, &A, &C) > 0.0) {
                        for (i = 1; i < 4; i++) {
                            if (Node[i] != Node[Next[0]]) {
                                if (Node[i] != Node[Next[1]])
                                    Cell2D[icell].ConnectNode[2] = Node[i];
                            }
                        }

                        Cell2D[icell].ConnectNode[1] = Node[Next[1]];
                        Cell2D[icell].ConnectNode[3] = Node[Next[0]];
                    } else {
                        Warn("CounterClockCellNodes_2D: Debug Algo:1");
                        return 1;
                    }
                } else {
                    if (CrossO_2D(&D, &A, &C) > 0.0) {
                        Warn("CounterClockCellNodes_2D: Debug Algo:2");
                        return 1;
                    } else {
                        for (i = 1; i < 4; i++) {
                            if (Node[i] != Node[Next[0]]) {
                                if (Node[i] != Node[Next[1]])
                                    Cell2D[icell].ConnectNode[2] = Node[i];
                            }
                        }

                        Cell2D[icell].ConnectNode[1] = Node[Next[0]];
                        Cell2D[icell].ConnectNode[3] = Node[Next[1]];
                    }
                }

                break;
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int CellCentroid_2D(void) {
    int icell, nodeID;
    int inode, nnodes;
    double sumx, sumy;

    for (icell = 0; icell < NoCells2D; icell++) {
        sumx = 0.0;
        sumy = 0.0;
        nnodes = Cell2D[icell].NodesPerCell;

        for (inode = 0; inode < nnodes; inode++) {
            nodeID = Cell2D[icell].ConnectNode[inode];
            sumx += Node2D[nodeID].Coordinate[0];
            sumy += Node2D[nodeID].Coordinate[1];
        }

        Cell2D[icell].Centroid[0] = sumx / (double) nnodes;
        Cell2D[icell].Centroid[1] = sumy / (double) nnodes;
    }
    return 0;
}

/*---------------------------------------------------------------*/
int CellArea_2D(void) {
    int icell, n0, n1, n2, n3;
    Point2 *vector;
    double area;
    int nnodes;

    for (icell = 0; icell < NoCells2D; icell++) {
        nnodes = Cell2D[icell].NodesPerCell;

        switch (nnodes) {
            case 3:
                vector = (Point2 *) malloc((nnodes - 1) * sizeof (Point2));
                if (vector == NULL) {
                    Warn("Memory Allocation Failed for vector");
                    return 1;
                }

                n0 = Cell2D[icell].ConnectNode[0];
                n1 = Cell2D[icell].ConnectNode[1];
                n2 = Cell2D[icell].ConnectNode[2];

                vector[0].x = Node2D[n1].Coordinate[0] - Node2D[n0].Coordinate[0];
                vector[0].y = Node2D[n1].Coordinate[1] - Node2D[n0].Coordinate[1];

                vector[1].x = Node2D[n2].Coordinate[0] - Node2D[n0].Coordinate[0];
                vector[1].y = Node2D[n2].Coordinate[1] - Node2D[n0].Coordinate[1];

                area = 0.5 * (vector[0].x * vector[1].y - vector[0].y * vector[1].x);

                if (area < 0.0) {
                    Cell2D[icell].ConnectNode[0] = n2;
                    Cell2D[icell].ConnectNode[2] = n0;
                }

                area = fabs(area);

                if (area < 1.0E-10) {
                    printf("Error: Triangle area=%lf, very small...\n", area);
                    exit(1);
                }

                Cell2D[icell].Area = fabs(area);
                free(vector);
                break;
            case 4:
                vector = (Point2 *) malloc((nnodes - 1) * sizeof (Point2));
                if (vector == NULL) {
                    Warn("Memory Allocation Failed for vector");
                    return 1;
                }

                n0 = Cell2D[icell].ConnectNode[0];
                n1 = Cell2D[icell].ConnectNode[1];
                n2 = Cell2D[icell].ConnectNode[2];
                n3 = Cell2D[icell].ConnectNode[3];

                vector[0].x = Node2D[n1].Coordinate[0] - Node2D[n0].Coordinate[0];
                vector[0].y = Node2D[n1].Coordinate[1] - Node2D[n0].Coordinate[1];

                vector[1].x = Node2D[n2].Coordinate[0] - Node2D[n0].Coordinate[0];
                vector[1].y = Node2D[n2].Coordinate[1] - Node2D[n0].Coordinate[1];

                vector[2].x = Node2D[n3].Coordinate[0] - Node2D[n0].Coordinate[0];
                vector[2].y = Node2D[n3].Coordinate[1] - Node2D[n0].Coordinate[1];

                area = 0.5 * ((vector[0].x * vector[1].y - vector[0].y * vector[1].x) +
                        (vector[1].x * vector[2].y - vector[1].y * vector[2].x));

                if (area < 0.0) {
                    printf("\n area of %d -ve ", icell);
                    Cell2D[icell].ConnectNode[0] = n2;
                    Cell2D[icell].ConnectNode[2] = n0;
                }

                area = fabs(area);

                if (area < 1.0E-10) {
                    printf("Error: Quadrilateral Area=%lf, very small\n", area);
                    exit(1);
                }

                Cell2D[icell].Area = fabs(area);
                free(vector);
                break;
        }
    }

    return 0;
}

/*---------------------------------------------------------------*/
int EdgeNormal_2D(void) {
    int iedge, tmp, n1p, n2p, n1, n2;
    int NodeList[8];
    int icell, edgeID, nedges;
    double x1, x2, y1, y2, isign, invLength;
    int *edgeNotVisited;

    if (NoEdges2D == 0) {
        Warn("No of Edges not known yet");
        return 1;
    }

    edgeNotVisited = (int *) malloc(NoEdges2D * sizeof (int));
    if (edgeNotVisited == NULL) {
        Warn("EdgeNormal_2D: Memory Failed");
        return 1;
    }

    for (iedge = 0; iedge < NoEdges2D; iedge++)
        edgeNotVisited[iedge] = 0;

    for (icell = 0; icell < NoCells2D; icell++) {
        nedges = Cell2D[icell].EdgesPerCell;

        switch (nedges) {
            case 3:
                NodeList[0] = Cell2D[icell].ConnectNode[0];
                NodeList[1] = Cell2D[icell].ConnectNode[1];
                NodeList[2] = Cell2D[icell].ConnectNode[1];
                NodeList[3] = Cell2D[icell].ConnectNode[2];
                NodeList[4] = Cell2D[icell].ConnectNode[2];
                NodeList[5] = Cell2D[icell].ConnectNode[0];
                break;
            case 4:
                NodeList[0] = Cell2D[icell].ConnectNode[0];
                NodeList[1] = Cell2D[icell].ConnectNode[1];
                NodeList[2] = Cell2D[icell].ConnectNode[1];
                NodeList[3] = Cell2D[icell].ConnectNode[2];
                NodeList[4] = Cell2D[icell].ConnectNode[2];
                NodeList[5] = Cell2D[icell].ConnectNode[3];
                NodeList[6] = Cell2D[icell].ConnectNode[3];
                NodeList[7] = Cell2D[icell].ConnectNode[0];
                break;
        }

        for (iedge = 0; iedge < nedges; iedge++) {
            n1 = NodeList[2 * iedge];
            n2 = NodeList[2 * iedge + 1];

            n1p = n1;
            n2p = n2;

            if (n2p < n1p) {
                tmp = n2p;
                n2p = n1p;
                n1p = tmp;
            }

            x1 = Node2D[n1].Coordinate[0];
            y1 = Node2D[n1].Coordinate[1];

            x2 = Node2D[n2].Coordinate[0];
            y2 = Node2D[n2].Coordinate[1];

            edgeID = GetEdgeNumber_2D(n1p, n2p);

            if (Edge2D[edgeID].Cell[1] == icell)
                isign = 1.0;
            else
                isign = -1.0;

            if (!edgeNotVisited[edgeID]) {
                Edge2D[edgeID].dx = x2 - x1;
                Edge2D[edgeID].dy = y2 - y1;
                Edge2D[edgeID].length = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

                if (Edge2D[edgeID].length == 0.0) {
                    Warn("Error : Edge length is zero");
                    return 1;
                }
                invLength = 1.0 / Edge2D[edgeID].length;
                Edge2D[edgeID].nx = 1.0 * isign * (y2 - y1) * invLength;
                Edge2D[edgeID].ny = -1.0 * isign * (x2 - x1) * invLength;
                edgeNotVisited[edgeID] = 1;
            }
        }
    }
    free(edgeNotVisited);
    return 0;
}

/*---------------------------------------------------------------*/
int PrePostProcessing_2D(void) {
    CellCentroid_2D();
    CellArea_2D();
    EdgeNormal_2D();
    CounterClockCellNodes_2D();

    return 0;
}
