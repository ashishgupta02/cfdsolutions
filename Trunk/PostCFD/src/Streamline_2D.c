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
 * File		Streamline_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include <GL/glut.h>
#include "Vector.h"
#include "Algebra.h"
#include "PostProcessing_2D.h"
#include "PPInitializeSolution_2D.h"
#include "Graphics_2D.h"
#include "DrawMesh_2D.h"
#include "ColorCode_2D.h"
#include "ColorPlot_2D.h"
#include "PPOptions_2D.h"
#include "CellPointOperation_2D.h"
#include "Interpolation_2D.h"
#include "Boundary_2D.h"
#include "Stream_2D.h"
#include "Streamline_2D.h"

double StreamTime = 0.0;
int NSeedPoint = 0;
int SeedRandom = 0;
int SeedLine = 0;
Point2 *SeedPoint;
int *SeedCellID;
int SeedID;
/* This Matrix Pointer is used to store Cononical Transformation Matrix of Cells */
/* Eta_Desh = T * X_Desh + K 
        T = CellConTran
        K = CellConConst 
 */
static Matrix2 *CellConTran;
static Point2 *CellConConst;
static double StreamTimeStep;


static Data_2D_Un StreamVelX;
static Data_2D_Un StreamVelY;

/* Streamline Equation */
/* d(Eta_Desh(t))/dt = C * Eta_Desh + e
        C = StreamCoef
        e = StreamConst
 */
static Matrix2 *StreamCoef;
static Point2 *StreamConst;

/* Use to Preprocess the GL Commands for faster display */
static GLuint StreamList;
static int notfirst = 0;
static double Time;

/* Streamline Arguments */
int StreamlineArg1 = -1;
int StreamlineArg2 = -1;

static int CheckStreamArg1;
static int CheckStreamArg2;

/* Min-Max Vector Value */
static double MinValueSol, MaxValueSol;

/*--------------------------------------------------------------*/
void DrawSeedPoint_2D(void) {
    int i;

    DrawMesh_2D();

    /* Display Seed Points */
    glColor3f(1.0, 0.0, 0.0);
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for (i = 0; i < NSeedPoint; i++) {
        glVertex2d((GLdouble) SeedPoint[i].x, (GLdouble) SeedPoint[i].y);
    }
    glEnd();
    glFlush();

}

/*--------------------------------------------------------------*/
void InitializeSeedPoint_2D(double x, double y) {
    int i, cellID, flag = 0;
    Point2 data;
    double tmp;

    if (!StreamlineFlag)
        return;

    /* For Stream Shade Line */
    if ((SeedLine == 1) && (SeedID == 2)) {
        if (SeedPoint[0].x != SeedPoint[1].x) {
            tmp = SeedPoint[1].x - SeedPoint[0].x;

            for (i = 0; i < (NSeedPoint - 2); i++) {
                if (flag == (NSeedPoint - 2))
                    continue;

                data.x = SeedPoint[0].x + ((i + 1)*(tmp / NSeedPoint));
                data.y = (((SeedPoint[1].y - SeedPoint[0].y) / (SeedPoint[1].x - SeedPoint[0].x)) * (data.x - SeedPoint[0].x))
                        + SeedPoint[0].y;

                CellContainPoint_2D(data.x, data.y, &cellID);
                if (cellID == -1) {
                    i--;
                    flag++;
                    continue;
                }

                SeedPoint[i + 2].x = data.x;
                SeedPoint[i + 2].y = data.y;
                SeedCellID[SeedID] = cellID;
                flag++;
                SeedID++;
            }
        } else {
            tmp = SeedPoint[1].y - SeedPoint[0].y;

            for (i = 0; i < (NSeedPoint - 2); i++) {
                if (flag == (NSeedPoint - 2))
                    continue;

                data.y = SeedPoint[0].y + ((i + 1)*(tmp / NSeedPoint));
                data.x = ((data.y - SeedPoint[0].y) * ((SeedPoint[1].x - SeedPoint[0].x) / (SeedPoint[1].y - SeedPoint[0].y)))
                        + SeedPoint[0].x;

                CellContainPoint_2D(data.x, data.y, &cellID);
                if (cellID == -1) {
                    i--;
                    flag++;
                    continue;
                }

                SeedPoint[i + 2].x = data.x;
                SeedPoint[i + 2].y = data.y;
                SeedCellID[SeedID] = cellID;
                flag++;
                SeedID++;
            }
        }

        NSeedPoint = SeedID;

        StreamlineFlag = 0;
        glutDisplayFunc(DrawSeedPoint_2D);
        return;
    }

    if (NSeedPoint == SeedID) {
        StreamlineFlag = 0;
        glutDisplayFunc(DrawSeedPoint_2D);
        return;
    }

    CellContainPoint_2D(x, y, &cellID);
    if (cellID == -1) {
        Warn("Invalid Point");
        return;
    }

    SeedPoint[SeedID].x = x;
    SeedPoint[SeedID].y = y;
    SeedCellID[SeedID] = cellID;
    SeedID++;
    /* Draw Point */
    glColor3f(1.0, 0.0, 0.0);
    glPointSize(2.0);
    glBegin(GL_POINTS);
    glVertex2d((GLdouble) x, (GLdouble) y);
    glEnd();
    glFlush();
    printf("Seed Point: %d\tX = %lf\tY = %lf\n", SeedID, x, y);
}

/*--------------------------------------------------------------*/
void GetSeedPoints_2D(void) {

    /* Release Previous Allocate SeedPoint */
    free(SeedPoint);
    free(SeedCellID);

    /* Allocate Memory */
    SeedPoint = (Point2 *) malloc(NSeedPoint * sizeof (Point2));
    SeedCellID = (int *) malloc(NSeedPoint * sizeof (int));
    if ((SeedPoint == NULL) || (SeedCellID == NULL)) {
        Warn("GetSeedPoint_2D: Memory Failed");
        StreamlineFlag = 0;
        return;
    }

    SeedID = 0;
    notfirst = 0;
}

/*--------------------------------------------------------------*/
int CellCononicalTransform_2D(void) {
    int i, node[3];
    double A;

    if (Cell2D[0].NodesPerCell != 3) {
        Warn("CellCononicalTransform_2D: Topology Not Supported");
        return 1;
    }

    /* Initialize Memory for Cell Cononical Transformation Matrix */
    CellConTran = (Matrix2 *) malloc(NoCells2D * sizeof (Matrix2));
    if (CellConTran == NULL) {
        Warn("CellCononicalTransform_2D: Memory Failed 1");
        free(CellConTran);
        return 1;
    }

    /* Initialize Memory for Cell Cononical Vector */
    CellConConst = (Point2 *) malloc(NoCells2D * sizeof (Point2));
    if (CellConConst == NULL) {
        Warn("CellCononicalTransform_2D: Memory Failed 2");
        free(CellConConst);
        free(CellConTran);
        return 1;
    }

    for (i = 0; i < NoCells2D; i++) {
        /* Get Cell Nodes */
        node[0] = Cell2D[i].ConnectNode[0];
        node[1] = Cell2D[i].ConnectNode[1];
        node[2] = Cell2D[i].ConnectNode[2];

        /* A = (x1-x0)(y2-y0) - (y1-y0)(x2-x0) */
        A = (((Node2D[node[1]].Coordinate[0] - Node2D[node[0]].Coordinate[0]) * (Node2D[node[2]].Coordinate[1] - Node2D[node[0]].Coordinate[1]))
                - ((Node2D[node[1]].Coordinate[1] - Node2D[node[0]].Coordinate[1]) * (Node2D[node[2]].Coordinate[0] - Node2D[node[0]].Coordinate[0])));

        CellConTran[i].element[0][0] = (Node2D[node[2]].Coordinate[1] - Node2D[node[0]].Coordinate[1]) / A;
        CellConTran[i].element[0][1] = (Node2D[node[0]].Coordinate[0] - Node2D[node[2]].Coordinate[0]) / A;
        CellConTran[i].element[1][0] = (Node2D[node[0]].Coordinate[1] - Node2D[node[1]].Coordinate[1]) / A;
        CellConTran[i].element[1][1] = (Node2D[node[1]].Coordinate[0] - Node2D[node[0]].Coordinate[0]) / A;

        CellConConst[i].x = ((CellConTran[i].element[0][0] * (-1.0 * Node2D[node[0]].Coordinate[0]))
                + (CellConTran[i].element[0][1] * (-1.0 * Node2D[node[0]].Coordinate[1])));

        CellConConst[i].y = ((CellConTran[i].element[1][0] * (-1.0 * Node2D[node[0]].Coordinate[0]))
                + (CellConTran[i].element[1][1] * (-1.0 * Node2D[node[0]].Coordinate[1])));
    }

    return 0;
}

/*--------------------------------------------------------------*/

/* If Velocity is Cell-Center StreamTimeStep_2D Interpolates to Nodes */
void StreamTimeStep_2D(void) {
    int i, j, VelX, VelY;
    double tmp0, tmp;
    int node[3];

    if (Cell2D[0].NodesPerCell != 3) {
        Warn("StreamTimeStep_2D: Topology Not Supported");
        return;
    }

    /* Streamline */
    if (StreamModule == 1) {
        VelX = StreamlineArg1;
        VelY = StreamlineArg2;

        switch (Solution2D.Location) {
            case Vertex:
                /* Get Initial Cell Node */
                node[0] = Cell2D[0].ConnectNode[0];

                tmp0 = pow(((CellConTran[0].element[0][0] * Solution2D.Sols[VelX]->Data[node[0]])
                        + (CellConTran[0].element[0][1] * Solution2D.Sols[VelY]->Data[node[0]])), 2)
                        + pow(((CellConTran[0].element[1][0] * Solution2D.Sols[VelX]->Data[node[0]])
                        + (CellConTran[0].element[1][1] * Solution2D.Sols[VelY]->Data[node[0]])), 2);

                for (i = 0; i < NoCells2D; i++) {
                    /* Get Cell Nodes */
                    node[0] = Cell2D[i].ConnectNode[0];
                    node[1] = Cell2D[i].ConnectNode[1];
                    node[2] = Cell2D[i].ConnectNode[2];

                    for (j = 0; j < 3; j++) {
                        tmp = pow(((CellConTran[i].element[0][0] * Solution2D.Sols[VelX]->Data[node[j]])
                                + (CellConTran[i].element[0][1] * Solution2D.Sols[VelY]->Data[node[j]])), 2)
                                + pow(((CellConTran[i].element[1][0] * Solution2D.Sols[VelX]->Data[node[j]])
                                + (CellConTran[i].element[1][1] * Solution2D.Sols[VelY]->Data[node[j]])), 2);

                        if (tmp > tmp0)
                            tmp0 = tmp;
                    }
                }
                break;
            case CellCenter:
                /* Initialize VelocityX and VelocityY for Nodes */
                StreamVelX.Size = NoNodes2D;
                StreamVelX.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (StreamVelX.Data == NULL) {
                    Warn("StreamTimeStep_2D: Memory Failed 1");
                    return;
                }

                StreamVelY.Size = NoNodes2D;
                StreamVelY.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (StreamVelY.Data == NULL) {
                    Warn("StreamTimeStep_2D: Memory Failed 2");
                    return;
                }

                /* Interpolate VelocityX and VelocityY to Nodes */
                ArithmeticCC2Node_2D(Solution2D.Sols[VelX], &StreamVelX);
                ArithmeticCC2Node_2D(Solution2D.Sols[VelY], &StreamVelY);

                /* Get Initial Cell Node */
                node[0] = Cell2D[0].ConnectNode[0];

                tmp0 = pow(((CellConTran[0].element[0][0] * StreamVelX.Data[node[0]])
                        + (CellConTran[0].element[0][1] * StreamVelY.Data[node[0]])), 2)
                        + pow(((CellConTran[0].element[1][0] * StreamVelX.Data[node[0]])
                        + (CellConTran[0].element[1][1] * StreamVelY.Data[node[0]])), 2);

                for (i = 0; i < NoCells2D; i++) {
                    /* Get Cell Nodes */
                    node[0] = Cell2D[i].ConnectNode[0];
                    node[1] = Cell2D[i].ConnectNode[1];
                    node[2] = Cell2D[i].ConnectNode[2];

                    for (j = 0; j < 3; j++) {
                        tmp = pow(((CellConTran[i].element[0][0] * StreamVelX.Data[node[j]])
                                + (CellConTran[i].element[0][1] * StreamVelY.Data[node[j]])), 2)
                                + pow(((CellConTran[i].element[1][0] * StreamVelX.Data[node[j]])
                                + (CellConTran[i].element[1][1] * StreamVelY.Data[node[j]])), 2);

                        if (tmp > tmp0)
                            tmp0 = tmp;
                    }
                }
                break;
        }
    }        /* Vector Topology */
    else {
        /* Get Initial Cell Node */
        node[0] = Cell2D[0].ConnectNode[0];

        tmp0 = pow(((CellConTran[0].element[0][0] * Vector2DX->Data[node[0]])
                + (CellConTran[0].element[0][1] * Vector2DY->Data[node[0]])), 2)
                + pow(((CellConTran[0].element[1][0] * Vector2DX->Data[node[0]])
                + (CellConTran[0].element[1][1] * Vector2DY->Data[node[0]])), 2);

        for (i = 0; i < NoCells2D; i++) {
            /* Get Cell Nodes */
            node[0] = Cell2D[i].ConnectNode[0];
            node[1] = Cell2D[i].ConnectNode[1];
            node[2] = Cell2D[i].ConnectNode[2];

            for (j = 0; j < 3; j++) {
                tmp = pow(((CellConTran[i].element[0][0] * Vector2DX->Data[node[j]])
                        + (CellConTran[i].element[0][1] * Vector2DY->Data[node[j]])), 2)
                        + pow(((CellConTran[i].element[1][0] * Vector2DX->Data[node[j]])
                        + (CellConTran[i].element[1][1] * Vector2DY->Data[node[j]])), 2);

                if (tmp > tmp0)
                    tmp0 = tmp;
            }
        }
    }

    tmp0 = sqrt(tmp0);

    StreamTimeStep = 1 / tmp0;
}

/*--------------------------------------------------------------*/
void StreamEqCoef_2D(void) {
    int i, j, VelX, VelY, node[3];
    double velocityX[3], velocityY[3];

    if (Cell2D[0].NodesPerCell != 3) {
        Warn("StreamEqCoeff_2D: Topology Not Supported");
        return;
    }

    /* Initialize Memory for Coefficients */
    StreamCoef = (Matrix2 *) malloc(NoCells2D * sizeof (Matrix2));
    StreamConst = (Point2 *) malloc(NoCells2D * sizeof (Point2));
    if ((StreamCoef == NULL) || (StreamConst == NULL)) {
        Warn("StreamEqCoef_2D: Memory Failed");
        return;
    }

    /* Streamline */
    if (StreamModule == 1) {
        VelX = StreamlineArg1;
        VelY = StreamlineArg2;

        switch (Solution2D.Location) {
            case Vertex:
                for (i = 0; i < NoCells2D; i++) {
                    for (j = 0; j < 3; j++) {
                        node[j] = Cell2D[i].ConnectNode[j];

                        velocityX[j] = Solution2D.Sols[VelX]->Data[node[j]];
                        velocityY[j] = Solution2D.Sols[VelY]->Data[node[j]];
                    }

                    StreamCoef[i].element[0][0] = (CellConTran[i].element[0][0] * (velocityX[1] - velocityX[0]))
                            + (CellConTran[i].element[0][1] * (velocityY[1] - velocityY[0]));
                    StreamCoef[i].element[0][1] = (CellConTran[i].element[0][0] * (velocityX[2] - velocityX[0]))
                            + (CellConTran[i].element[0][1] * (velocityY[2] - velocityY[0]));
                    StreamCoef[i].element[1][0] = (CellConTran[i].element[1][0] * (velocityX[1] - velocityX[0]))
                            + (CellConTran[i].element[1][1] * (velocityY[1] - velocityY[0]));
                    StreamCoef[i].element[1][1] = (CellConTran[i].element[1][0] * (velocityX[2] - velocityX[0]))
                            + (CellConTran[i].element[1][1] * (velocityY[2] - velocityY[0]));

                    StreamConst[i].x = (CellConTran[i].element[0][0] * velocityX[0]) + (CellConTran[i].element[0][1] * velocityY[0]);
                    StreamConst[i].y = (CellConTran[i].element[1][0] * velocityX[0]) + (CellConTran[i].element[1][1] * velocityY[0]);
                }
                break;
            case CellCenter:
                for (i = 0; i < NoCells2D; i++) {
                    for (j = 0; j < 3; j++) {
                        node[j] = Cell2D[i].ConnectNode[j];

                        velocityX[j] = StreamVelX.Data[node[j]];
                        velocityY[j] = StreamVelY.Data[node[j]];
                    }

                    StreamCoef[i].element[0][0] = (CellConTran[i].element[0][0] * (velocityX[1] - velocityX[0]))
                            + (CellConTran[i].element[0][1] * (velocityY[1] - velocityY[0]));
                    StreamCoef[i].element[0][1] = (CellConTran[i].element[0][0] * (velocityX[2] - velocityX[0]))
                            + (CellConTran[i].element[0][1] * (velocityY[2] - velocityY[0]));
                    StreamCoef[i].element[1][0] = (CellConTran[i].element[1][0] * (velocityX[1] - velocityX[0]))
                            + (CellConTran[i].element[1][1] * (velocityY[1] - velocityY[0]));
                    StreamCoef[i].element[1][1] = (CellConTran[i].element[1][0] * (velocityX[2] - velocityX[0]))
                            + (CellConTran[i].element[1][1] * (velocityY[2] - velocityY[0]));

                    StreamConst[i].x = (CellConTran[i].element[0][0] * velocityX[0]) + (CellConTran[i].element[0][1] * velocityY[0]);
                    StreamConst[i].y = (CellConTran[i].element[1][0] * velocityX[0]) + (CellConTran[i].element[1][1] * velocityY[0]);
                }
                break;
        }
    }        /* Vector Topology */
    else {
        for (i = 0; i < NoCells2D; i++) {
            for (j = 0; j < 3; j++) {
                node[j] = Cell2D[i].ConnectNode[j];

                velocityX[j] = Vector2DX->Data[node[j]];
                velocityY[j] = Vector2DY->Data[node[j]];
            }

            StreamCoef[i].element[0][0] = (CellConTran[i].element[0][0] * (velocityX[1] - velocityX[0]))
                    + (CellConTran[i].element[0][1] * (velocityY[1] - velocityY[0]));
            StreamCoef[i].element[0][1] = (CellConTran[i].element[0][0] * (velocityX[2] - velocityX[0]))
                    + (CellConTran[i].element[0][1] * (velocityY[2] - velocityY[0]));
            StreamCoef[i].element[1][0] = (CellConTran[i].element[1][0] * (velocityX[1] - velocityX[0]))
                    + (CellConTran[i].element[1][1] * (velocityY[1] - velocityY[0]));
            StreamCoef[i].element[1][1] = (CellConTran[i].element[1][0] * (velocityX[2] - velocityX[0]))
                    + (CellConTran[i].element[1][1] * (velocityY[2] - velocityY[0]));

            StreamConst[i].x = (CellConTran[i].element[0][0] * velocityX[0]) + (CellConTran[i].element[0][1] * velocityY[0]);
            StreamConst[i].y = (CellConTran[i].element[1][0] * velocityX[0]) + (CellConTran[i].element[1][1] * velocityY[0]);
        }
    }
}

/*--------------------------------------------------------------*/

/* Cell Cononical Coordinate from Physical Coordinates */
void CellPhy2ConCoordinate_2D(double x, double y, int cellID, Point2 *ConCoord) {
    ConCoord->x = (CellConTran[cellID].element[0][0] * x) + (CellConTran[cellID].element[0][1] * y) + CellConConst[cellID].x;
    ConCoord->y = (CellConTran[cellID].element[1][0] * x) + (CellConTran[cellID].element[1][1] * y) + CellConConst[cellID].y;
}

/*--------------------------------------------------------------*/

/* Physical Coordinate from Cononical Coordinates */
void CellCon2PhyCoordinate_2D(Point2 ConCoord, int cellID, Point2 *PhyCoord) {
    int node[3];

    node[0] = Cell2D[cellID].ConnectNode[0];
    node[1] = Cell2D[cellID].ConnectNode[1];
    node[2] = Cell2D[cellID].ConnectNode[2];

    PhyCoord->x = ((1 - ConCoord.x - ConCoord.y) * Node2D[node[0]].Coordinate[0])
            + (ConCoord.x * Node2D[node[1]].Coordinate[0])
            + (ConCoord.y * Node2D[node[2]].Coordinate[0]);
    PhyCoord->y = ((1 - ConCoord.x - ConCoord.y) * Node2D[node[0]].Coordinate[1])
            + (ConCoord.x * Node2D[node[1]].Coordinate[1])
            + (ConCoord.y * Node2D[node[2]].Coordinate[1]);
}

/*--------------------------------------------------------------*/
void StreamRungeKutta4_2D(int CurCellID, Point2 CurConCoord, Point2 *NextConCoord) {
    int i;
    double C[4];
    Matrix2 H;
    Matrix2 tmp;
    Point2 e;

    /* Calculate H Matrix */
    C[0] = StreamTimeStep;
    C[1] = pow(StreamTimeStep, 2) / 2;
    C[2] = pow(StreamTimeStep, 3) / 6;
    C[3] = pow(StreamTimeStep, 4) / 24;

    H.element[0][0] = 0.0;
    H.element[0][1] = 0.0;
    H.element[1][0] = 0.0;
    H.element[1][1] = 0.0;

    for (i = 4; i > 0; i--) {
        H.element[0][0] = H.element[0][0] + C[i - 1];
        H.element[1][1] = H.element[1][1] + C[i - 1];

        H.element[0][0] = (H.element[0][0] * StreamCoef[CurCellID].element[0][0])
                + (H.element[0][1] * StreamCoef[CurCellID].element[1][0]);

        H.element[0][1] = (H.element[0][0] * StreamCoef[CurCellID].element[0][1])
                + (H.element[0][1] * StreamCoef[CurCellID].element[1][1]);

        H.element[1][0] = (H.element[1][0] * StreamCoef[CurCellID].element[0][0])
                + (H.element[1][1] * StreamCoef[CurCellID].element[1][0]);

        H.element[1][1] = (H.element[1][0] * StreamCoef[CurCellID].element[0][1])
                + (H.element[1][1] * StreamCoef[CurCellID].element[1][1]);
    }

    H.element[0][0] = H.element[0][0] + 1.0;
    H.element[1][1] = H.element[1][1] + 1.0;

    /* Calculate e Vector */
    C[0] = StreamTimeStep / 2;
    C[1] = pow(StreamTimeStep, 2) / 6;
    C[3] = pow(StreamTimeStep, 3) / 24;

    tmp.element[0][0] = 0.0;
    tmp.element[0][1] = 0.0;
    tmp.element[1][0] = 0.0;
    tmp.element[1][1] = 0.0;

    for (i = 3; i > 0; i--) {
        tmp.element[0][0] = tmp.element[0][0] + C[i - 1];
        tmp.element[1][1] = tmp.element[1][1] + C[i - 1];

        tmp.element[0][0] = (tmp.element[0][0] * StreamCoef[CurCellID].element[0][0])
                + (tmp.element[0][1] * StreamCoef[CurCellID].element[1][0]);

        tmp.element[0][1] = (tmp.element[0][0] * StreamCoef[CurCellID].element[0][1])
                + (tmp.element[0][1] * StreamCoef[CurCellID].element[1][1]);

        tmp.element[1][0] = (tmp.element[1][0] * StreamCoef[CurCellID].element[0][0])
                + (tmp.element[1][1] * StreamCoef[CurCellID].element[1][0]);

        tmp.element[1][1] = (tmp.element[1][0] * StreamCoef[CurCellID].element[0][1])
                + (tmp.element[1][1] * StreamCoef[CurCellID].element[1][1]);
    }

    tmp.element[0][0] = tmp.element[0][0] + 1.0;
    tmp.element[1][1] = tmp.element[1][1] + 1.0;

    e.x = StreamTimeStep * ((tmp.element[0][0] * StreamConst[CurCellID].x) + (tmp.element[0][1] * StreamConst[CurCellID].y));
    e.y = StreamTimeStep * ((tmp.element[1][0] * StreamConst[CurCellID].x) + (tmp.element[1][1] * StreamConst[CurCellID].y));

    /* Calculate Next Cononical Coordinate */
    NextConCoord->x = (H.element[0][0] * CurConCoord.x) + (H.element[0][1] * CurConCoord.y) + e.x;
    NextConCoord->y = (H.element[1][0] * CurConCoord.x) + (H.element[1][1] * CurConCoord.y) + e.y;
}

/*--------------------------------------------------------------*/
void StreamNextCell_2D(int CurCellID, Point2 NextConCoord, int *NextCellID) {
    int i, flag, sEdge;
    double N[3];
    double MinN, MaxN;
    Point2 PhyCoord;

    /* Get Current Cell Physical Coordinates */
    CellCon2PhyCoordinate_2D(NextConCoord, CurCellID, &PhyCoord);

    for (;;) {
        /* Shape Functions */
        N[0] = 1 - NextConCoord.x - NextConCoord.y;
        N[1] = NextConCoord.x;
        N[2] = NextConCoord.y;

        MinN = N[0];
        MaxN = N[0];
        flag = 0;

        for (i = 1; i < 3; i++) {
            if (MaxN < N[i])
                MaxN = N[i];

            if (MinN > N[i]) {
                MinN = N[i];
                flag = i;
            }
        }

        /* If Next and Current Cells Are Same */
        if ((MaxN <= 1.0) && (MinN >= 0.0)) {
            *NextCellID = CurCellID;
            break;
        }

        switch (flag) {
            case 0:
                sEdge = Cell2D[CurCellID].ConnectEdge[1];
                break;
            case 1:
                sEdge = Cell2D[CurCellID].ConnectEdge[2];
                break;
            case 2:
                sEdge = Cell2D[CurCellID].ConnectEdge[0];
                break;
        }

        /* Boundary Edge */
        if (Edge2D[sEdge].NoOfNeighbourCells == 1) {
            *NextCellID = -1;
            break;
        }

        /* Get Current Search Cell ID */
        if (Edge2D[sEdge].Cell[0] != CurCellID)
            CurCellID = Edge2D[sEdge].Cell[0];
        else
            CurCellID = Edge2D[sEdge].Cell[1];

        /* Get Current Search Cell Cononical Coordinate */
        CellPhy2ConCoordinate_2D(PhyCoord.x, PhyCoord.y, CurCellID, &NextConCoord);
    }
}

/*--------------------------------------------------------------*/
void Streamline_2D(void) {
    int i, CurCellID, NextCellID, node[3], Color_Index;
    int flag = 0;
    Point2 CurConCoord, NextConCoord, PhyCoord, InterpVel;
    double timestep = 0.0;
    Data_2D_Un *VelX, *VelY;

    switch (StreamModule) {
        case 1:
            /* Streamlines */
            if ((StreamlineArg1 == -1) || (StreamlineArg2 == -1)) {
                Warn("Streamline_2D: Invalid Vector Initialization");
                return;
            }

            /* Should not go for first time */
            if (notfirst) {
                if ((StreamlineArg1 != CheckStreamArg1) || (StreamlineArg2 != CheckStreamArg2))
                    flag = 1;
            }

            /* Initialize Check Arguments */
            CheckStreamArg1 = StreamlineArg1;
            CheckStreamArg2 = StreamlineArg2;
            break;

        case 2:
            /* Vector Topology */
            flag = 1;
            break;
    }

    /* Check Seed Points */
    if (NSeedPoint != SeedID) {
        Warn("Streamline_2D: Seed Points Mismatched");
        return;
    }

    if ((!notfirst) || (flag == 1)) {
        /* Initialize Color */
        ColorEncode_2D();

        if (StreamModule == 1) {
            /* Get Min-Max Solution */
            SolutionMinMaxVector_2D(StreamlineArg1, StreamlineArg2, &MinValueSol, &MaxValueSol);
        } else {
            /* Get Min-Max Vector */
            MinMaxVector_2D(Vector2DX, Vector2DY, &MinValueSol, &MaxValueSol);
        }
    }

    /* Reset the Screen window */
    DrawInitial_2D();

    /* Display Mesh Name */
    DrawName_2D(-1);

    /* Draw the legend of color */
    if (LegendActivateFlag)
        DrawLegend_2D(MinValueSol, MaxValueSol, 1, 1);

    /* Initialize Drawing Area */
    InitializeDrawArea_2D();

    /* Draw Mesh */
    if (MeshActivateFlag)
        Mesh_2D();

    /* Display Preprocessed List */
    if ((notfirst) && (StreamTime == Time)) {
        glCallList(StreamList);

        /* Display Boundaries */
        if (BoundaryActivateFlag)
            DisplayBoundary_2D();

        return;
    }

    /* Check Stream Time */
    if (StreamTime == 0.0) {
        Warn("Streamline_2D: Stream Time = 0.0");
        return;
    }

    CellCononicalTransform_2D();
    StreamTimeStep_2D();
    printf("Stream Time Step = %lf\n", StreamTimeStep);

    i = (int) (StreamTime / StreamTimeStep);
    if (i > 100000) {
        Warn("Streamline_2D: StreamTime too large");
        return;
    }

    StreamEqCoef_2D();

    /* Get Velocity Data Pointers */
    if (StreamModule == 1) {
        switch (Solution2D.Location) {
            case Vertex:
                VelX = Solution2D.Sols[StreamlineArg1];
                VelY = Solution2D.Sols[StreamlineArg2];
                break;
            case CellCenter:
                VelX = &StreamVelX;
                VelY = &StreamVelY;
                break;
        }
    } else {
        /* Vector Topology Pointers */
        VelX = Vector2DX;
        VelY = Vector2DY;
    }

    notfirst = 1;
    Time = StreamTime;

    /* Delete Old List */
    glDeleteLists(StreamList, 1);

    /* Create New List */
    StreamList = glGenLists(1);

    glNewList(StreamList, GL_COMPILE_AND_EXECUTE);

    glShadeModel(GL_SMOOTH);

    for (i = 0; i < NSeedPoint; i++) {
        CurCellID = SeedCellID[i];

        PhyCoord.x = SeedPoint[i].x;
        PhyCoord.y = SeedPoint[i].y;

        timestep = 0.0;

        glBegin(GL_LINE_STRIP);
        do {
            /* Physical Coordinate of Current Point */
            CellPhy2ConCoordinate_2D(PhyCoord.x, PhyCoord.y, CurCellID, &CurConCoord);

            /* Get CurCell NodesID */
            node[0] = Cell2D[CurCellID].ConnectNode[0];
            node[1] = Cell2D[CurCellID].ConnectNode[1];
            node[2] = Cell2D[CurCellID].ConnectNode[2];

            /* Get Interpolation Velocity */
            InterpVel.x = ((1 - CurConCoord.x - CurConCoord.y) * VelX->Data[node[0]])
                    + (CurConCoord.x * VelX->Data[node[1]])
                    + (CurConCoord.y * VelX->Data[node[2]]);
            InterpVel.y = ((1 - CurConCoord.x - CurConCoord.y) * VelY->Data[node[0]])
                    + (CurConCoord.x * VelY->Data[node[1]])
                    + (CurConCoord.y * VelY->Data[node[2]]);

            /* Get Color Index */
            Color_Index = (int) ((sqrt(pow(InterpVel.x, 2) + pow(InterpVel.y, 2)) - MinValueSol) / (MaxValueSol - MinValueSol) * 255);

            /* Setting Drawing Color */
            glColor3f(Red[Color_Index], Green[Color_Index], Blue[Color_Index]);

            /* Draw Next Point */
            glVertex2d(PhyCoord.x, PhyCoord.y);

            /* Solve Runge-Kutta Method Order Four */
            StreamRungeKutta4_2D(CurCellID, CurConCoord, &NextConCoord);

            /* Physical Coordinate of Next Point */
            CellCon2PhyCoordinate_2D(NextConCoord, CurCellID, &PhyCoord);

            /* Find Next Cell Containing Stream Point */
            StreamNextCell_2D(CurCellID, NextConCoord, &NextCellID);

            /* Reached the Boundary */
            if (NextCellID == -1)
                break;

            /* Next Stream Point Cell */
            CurCellID = NextCellID;

            timestep = timestep + StreamTimeStep;
        } while (timestep < StreamTime);

        glEnd();
    }

    glEndList();

    glFlush();

    /* Display Boundaries */
    if (BoundaryActivateFlag)
        DisplayBoundary_2D();

    if (CellConTran != NULL) free(CellConTran);
    if (CellConConst != NULL) free(CellConConst);
    if (StreamCoef != NULL) free(StreamCoef);
    if (StreamConst != NULL) free(StreamConst);
    if (StreamModule == 1) {
        if (StreamVelX.Data != NULL) free(StreamVelX.Data);
        if (StreamVelY.Data != NULL) free(StreamVelY.Data);
    }
}
