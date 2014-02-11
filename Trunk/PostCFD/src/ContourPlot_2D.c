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
 * File		ContourPlot_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <GL/glut.h>
#include <cgnslib.h>
#include "PostProcessing_2D.h"
#include "PPInitializeSolution_2D.h"
#include "Interpolation_2D.h"
#include "PPOptions_2D.h"
#include "ColorPlot_2D.h"
#include "Graphics_2D.h"
#include "ColorCode_2D.h"
#include "DrawMesh_2D.h"
#include "Boundary_2D.h"

/* Use to Preprocess the GL Commands for faster display */
static GLuint ContourList;

/* -1 is because sol starts from 0 */
static int ContourPlotID = -1;

static int LevelNo = 0;

static int notfirst = 0;

/*--------------------------------------------------------------*/
int CreateContourPlot_2D(int SolID, int Level) {
    /* Used when Location is Cell-Center */
    Data_2D_Un tmp;
    int i, j, k, m, color_index;
    int Node1, Node2;
    /* Node Solution Values */
    double NSol1, NSol2;
    /* Interpolated Coordinate */
    double ConX, ConY;
    double MinValueSol, MaxValueSol, InterValueSol, DelValueSol;
    /* Cordinates for Flood Fill */
    double *FloodFillX1, *FloodFillX2, *FloodFillY1, *FloodFillY2;
    int FloodCount, FloodValue, NodeSum, l, Node3, Node[3];
    double NSol3;
    int *FloodColors, *FloodIndex, *FloodEdge1, *FloodEdge2;

    /* Disable Streamline Flag */
    if (StreamlineFlag == 1)
        StreamlineFlag = 0;

    /* Reset the Screen window */
    DrawInitial_2D();

    /* Initialize Drawing Area */
    InitializeDrawArea_2D();

    if (Cell2D[0].NodesPerCell != 3) {
        //MSG ("CreateContourPlot_2D: Only Trinagle for now");
        glFlush();
        return 1;
    }

    /* Initialize Color */
    ColorEncode_2D();

    /* Get Min-Max Solution */
    SolutionMinMax_2D(SolID, &MinValueSol, &MaxValueSol);

    /*Draw Mesh */
    if (MeshActivateFlag)
        Mesh_2D();

    if (Level < 10)
        Level = 10;

    /* Display Preprocessed List */
    if ((ContourPlotID == SolID) && (LevelNo == Level)) {
        glCallList(ContourList);
    } else {
        ContourPlotID = SolID;
        LevelNo = Level;

        if (notfirst) {
            /* Delete Old List */
            glDeleteLists(ContourList, 1);

            /* Create New List */
            ContourList = glGenLists(1);
        } else {
            /* Create New List */
            ContourList = glGenLists(1);
            notfirst = 1;
        }

        glNewList(ContourList, GL_COMPILE_AND_EXECUTE);

        DelValueSol = (MaxValueSol - MinValueSol) / Level;
        InterValueSol = MinValueSol;

        switch (Solution2D.Location) {
            case Vertex:
                /* 	i => Level No
                        j => Cell No
                        k, m => Node No
                 */
                /* For Each Level Value */
                //				for (i = 0; i < Level; i++) {
                //					/* Set Color for Level */
                //					color_index = (int)(((float)i/(float)Level) * 255);
                //					glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
                //					glLineWidth(1.0);
                //					
                //					for (j = 0; j < NoCells2D; j++) {
                //						glBegin(GL_LINE_STRIP);
                //						for (k = 0, m = (Cell2D[j].NodesPerCell - 1); k < Cell2D[j].NodesPerCell; m = k, k++) {
                //							/* Get Node IDs */
                //							Node1 = Cell2D[j].ConnectNode[m];
                //							Node2 = Cell2D[j].ConnectNode[k];
                //							
                //							/* Get Node Solution Value */
                //							NSol1 = Solution2D.Sols[SolID]->Data[Node1];
                //							NSol2 = Solution2D.Sols[SolID]->Data[Node2];
                //							
                //							/* Check if the Edge Contains a Contour point */
                //							if (((NSol1 <= InterValueSol) && (NSol2 > InterValueSol)) || ((NSol2 <= InterValueSol) && (NSol1 > InterValueSol))) {
                //								/* Linear Interpolation */
                //								/* Get X coordinate */
                //								ConX = ((Node2D[Node1].Coordinate[0]*(NSol2 - InterValueSol))/(NSol2 - NSol1)) 
                //									+ ((Node2D[Node2].Coordinate[0]*(InterValueSol - NSol1))/(NSol2 - NSol1));
                //								
                //								/* Get Y coordinate */
                //								ConY = ((Node2D[Node1].Coordinate[1]*(NSol2 - InterValueSol))/(NSol2 - NSol1)) 
                //									+ ((Node2D[Node2].Coordinate[1]*(InterValueSol - NSol1))/(NSol2 - NSol1));
                //									
                //								/* Display Point */
                //								glVertex2d((GLdouble)ConX, (GLdouble)ConY);
                //							}
                //						}
                //						glEnd();
                //					}
                //					/* Update InterValueSol for Next Level */
                //					InterValueSol = InterValueSol + DelValueSol;
                //				}
                /* Try to do flood fill */
                /* Allocate the Memory to hold the coordinates of Line Segments */
                FloodFillX1 = malloc(Level * sizeof (double));
                FloodFillX2 = malloc(Level * sizeof (double));
                FloodFillY1 = malloc(Level * sizeof (double));
                FloodFillY2 = malloc(Level * sizeof (double));
                /* Allocate the Memory to hold Flood Colors */
                FloodColors = malloc((Level + 1) * sizeof (int));
                /* Allocate the Memory to hold Valid Flood Index No */
                FloodIndex = malloc(Level * sizeof (int));
                /* Allocate the Memory to hold Edge Information */
                FloodEdge1 = malloc(Level * sizeof (int));
                FloodEdge2 = malloc(Level * sizeof (int));
                /* Assign the Color Index */
                for (i = 0; i < (Level + 1); i++) {
                    FloodColors[i] = (int) (((float) i / (float) (Level + 1)) * 255);
                }
                /* Initialise Level Index */
                for (i = 0; i < Level; i++) {
                    FloodIndex[i] = -1;
                }
                for (i = 0; i < NoCells2D; i++) {
                    InterValueSol = MinValueSol;
                    FloodCount = 0;
                    NodeSum = 0;
                    FloodValue = 0;
                    for (j = 0; j < Level; j++) {
                        l = 0;
                        for (k = 0, m = (Cell2D[i].NodesPerCell - 1); k < Cell2D[i].NodesPerCell; m = k, k++) {
                            /* Get Node IDs */
                            Node1 = Cell2D[i].ConnectNode[m];
                            Node2 = Cell2D[i].ConnectNode[k];

                            /* Get the Node Solution Value */
                            NSol1 = Solution2D.Sols[SolID]->Data[Node1];
                            NSol2 = Solution2D.Sols[SolID]->Data[Node2];
                            /* Check if the Edge Contains a Contour point */
                            if (((NSol1 <= InterValueSol) && (NSol2 > InterValueSol)) || ((NSol2 <= InterValueSol) && (NSol1 > InterValueSol))) {
                                /* Linear Interpolation */
                                /* First Coordinate */
                                if (l == 0) {
                                    /* Get X coordinate */
                                    FloodFillX1[FloodCount] = ((Node2D[Node1].Coordinate[0]*(NSol2 - InterValueSol)) / (NSol2 - NSol1))
                                            + ((Node2D[Node2].Coordinate[0]*(InterValueSol - NSol1)) / (NSol2 - NSol1));

                                    /* Get Y coordinate */
                                    FloodFillY1[FloodCount] = ((Node2D[Node1].Coordinate[1]*(NSol2 - InterValueSol)) / (NSol2 - NSol1))
                                            + ((Node2D[Node2].Coordinate[1]*(InterValueSol - NSol1)) / (NSol2 - NSol1));

                                    if (k == 0 && m == 2) {
                                        FloodEdge1[FloodCount] = 1;
                                    } else if (k == 1 && m == 0) {
                                        FloodEdge1[FloodCount] = 2;
                                    } else {
                                        FloodEdge1[FloodCount] = 4;
                                    }

                                    l = 1;
                                } else {
                                    /* Second Coordinate */
                                    /* Get X coordinate */
                                    FloodFillX2[FloodCount] = ((Node2D[Node1].Coordinate[0]*(NSol2 - InterValueSol)) / (NSol2 - NSol1))
                                            + ((Node2D[Node2].Coordinate[0]*(InterValueSol - NSol1)) / (NSol2 - NSol1));

                                    /* Get Y coordinate */
                                    FloodFillY2[FloodCount] = ((Node2D[Node1].Coordinate[1]*(NSol2 - InterValueSol)) / (NSol2 - NSol1))
                                            + ((Node2D[Node2].Coordinate[1]*(InterValueSol - NSol1)) / (NSol2 - NSol1));

                                    FloodIndex[FloodCount] = j;

                                    if (k == 0 && m == 2) {
                                        FloodEdge2[FloodCount] = 1;
                                    } else if (k == 1 && m == 0) {
                                        FloodEdge2[FloodCount] = 2;
                                    } else {
                                        FloodEdge2[FloodCount] = 4;
                                    }

                                    FloodCount++;
                                }

                            }
                        }
                        /* Update InterValueSol for Next Level */
                        InterValueSol = InterValueSol + DelValueSol;
                    }
                    /* Flood Fill Diagonistic Test */
                    /* No Contour Cut */
                    if (FloodCount == 0) {
                        Node1 = Cell2D[i].ConnectNode[0];
                        Node2 = Cell2D[i].ConnectNode[1];
                        Node3 = Cell2D[i].ConnectNode[2];
                        NSol1 = Solution2D.Sols[SolID]->Data[Node1];
                        NSol2 = Solution2D.Sols[SolID]->Data[Node2];
                        NSol3 = Solution2D.Sols[SolID]->Data[Node3];

                        InterValueSol = MinValueSol;
                        /* Get the Color Index */
                        for (j = 0; j < Level; j++) {
                            if (NSol1 < InterValueSol) {
                                color_index = FloodColors[j];
                                break;
                            }
                            InterValueSol = InterValueSol + DelValueSol;
                        }
                        /* Display the polygons */
                        glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
                        glBegin(GL_TRIANGLES);
                        glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                        glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                        glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                        glEnd();
                    }
                    /* Cell as Contour Segment */
                    if (FloodCount > 0) {
                        Node1 = Cell2D[i].ConnectNode[0];
                        Node2 = Cell2D[i].ConnectNode[1];
                        Node3 = Cell2D[i].ConnectNode[2];
                        NSol1 = Solution2D.Sols[SolID]->Data[Node1];
                        NSol2 = Solution2D.Sols[SolID]->Data[Node2];
                        NSol3 = Solution2D.Sols[SolID]->Data[Node3];
                        Node[0] = 0;
                        Node[1] = 0;
                        Node[2] = 0;
                        for (j = 0; j < FloodCount; j++) {
                            InterValueSol = MinValueSol + (FloodIndex[j] * DelValueSol);
                            /*
                            glColor3f(0.0, 0.0, 0.0);
                            glLineWidth(1.0);
                            glBegin(GL_LINES);
                               glVertex2d((GLdouble)FloodFillX1[j], (GLdouble)FloodFillY1[j]);
                               glVertex2d((GLdouble)FloodFillX2[j], (GLdouble)FloodFillY2[j]);
                            glEnd();
                             */
                            if (j == 0) {
                                /* Exploiting the Linear Interpolation of Triangular Element */
                                if (NSol1 < InterValueSol) {
                                    Node[0] = 1;
                                }
                                if (NSol2 < InterValueSol) {
                                    Node[1] = 1;
                                }
                                if (NSol3 < InterValueSol) {
                                    Node[2] = 1;
                                }
                                NodeSum = Node[0] + Node[1] + Node[2];
                                FloodValue = FloodEdge1[j] + FloodEdge2[j];
                                color_index = FloodColors[FloodIndex[j]];
                                glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
                                switch (FloodValue) {
                                    case 3:
                                        glBegin(GL_POLYGON);
                                        if (NodeSum == 2) {
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                            glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                            glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                        } else {
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                        }
                                        glEnd();
                                        break;
                                    case 5:
                                        glBegin(GL_POLYGON);
                                        if (NodeSum == 2) {
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                            glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                        } else {
                                            glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                        }
                                        glEnd();
                                        break;
                                    case 6:
                                        glBegin(GL_POLYGON);
                                        if (NodeSum == 2) {
                                            glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                            glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                        } else {
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                        }
                                        glEnd();
                                        break;
                                }
                            } else {
                                color_index = FloodColors[FloodIndex[j]];
                                glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
                                if (NodeSum == 2) {
                                    glBegin(GL_POLYGON);
                                    if (FloodValue == 3 || FloodValue == 6) {
                                        glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                        glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                        glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                        glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                    } else {
                                        glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                        glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                        glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                        glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                    }
                                    glEnd();
                                } else {
                                    /* Exploiting the Linear Interpolation of Triangular Element */
                                    if (NSol1 < InterValueSol) {
                                        Node[0] = 1;
                                    }
                                    if (NSol2 < InterValueSol) {
                                        Node[1] = 1;
                                    }
                                    if (NSol3 < InterValueSol) {
                                        Node[2] = 1;
                                    }
                                    NodeSum = Node[0] + Node[1] + Node[2];
                                    FloodValue = FloodEdge1[j] + FloodEdge2[j];
                                    if (NodeSum == 1) {
                                        glBegin(GL_POLYGON);
                                        if (FloodValue == 3 || FloodValue == 6) {
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                            glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                        } else {
                                            glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                            glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                            glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                            glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                        }
                                        glEnd();
                                    } else {
                                        if (Node[0] == 0) {
                                            switch (FloodEdge1[j - 1] + FloodEdge2[j - 1]) {
                                                case 5:
                                                    glBegin(GL_POLYGON);
                                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                                    glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                                    glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                                    glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                                    glEnd();
                                                    break;
                                                case 6:
                                                    glBegin(GL_POLYGON);
                                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                                    glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                                    glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                                    glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                                    glEnd();
                                                    break;
                                            }
                                        }
                                        if (Node[1] == 0) {
                                            switch (FloodEdge1[j - 1] + FloodEdge2[j - 1]) {
                                                case 3:
                                                    glBegin(GL_POLYGON);
                                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                                    glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                                    glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                                    glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                                    glEnd();
                                                    break;
                                                case 5:
                                                    glBegin(GL_POLYGON);
                                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                                    glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                                    glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                                    glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                                    glEnd();
                                                    break;
                                            }
                                        }
                                        if (Node[2] == 0) {
                                            switch (FloodEdge1[j - 1] + FloodEdge2[j - 1]) {
                                                case 3:
                                                    glBegin(GL_POLYGON);
                                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                                    glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                                    glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                                    glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                                    glEnd();
                                                    break;
                                                case 6:
                                                    glBegin(GL_POLYGON);
                                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                                    glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                                    glVertex2d((GLdouble) FloodFillX1[j - 1], (GLdouble) FloodFillY1[j - 1]);
                                                    glVertex2d((GLdouble) FloodFillX2[j - 1], (GLdouble) FloodFillY2[j - 1]);
                                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                                    glEnd();
                                                    break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        /* Fill the Last Polygon in Triangle if available */
                        j = j - 1;
                        color_index = FloodColors[FloodIndex[j] + 1];
                        glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
                        if (NodeSum == 2) {
                            switch (FloodEdge1[j] + FloodEdge2[j]) {
                                case 3:
                                    glBegin(GL_POLYGON);
                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                    glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                    glEnd();
                                    break;
                                case 5:
                                    glBegin(GL_POLYGON);
                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                    glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                    glEnd();
                                    break;
                                case 6:
                                    glBegin(GL_POLYGON);
                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                    glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                    glEnd();
                                    break;
                            }
                        }
                        if (NodeSum == 1) {
                            switch (FloodEdge1[j] + FloodEdge2[j]) {
                                case 3:
                                    glBegin(GL_POLYGON);
                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                    glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                    glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                    glEnd();
                                    break;
                                case 5:
                                    glBegin(GL_POLYGON);
                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                    glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                    glVertex2d((GLdouble) Node2D[Node2].Coordinate[0], (GLdouble) Node2D[Node2].Coordinate[1]);
                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                    glEnd();
                                    break;
                                case 6:
                                    glBegin(GL_POLYGON);
                                    glVertex2d((GLdouble) FloodFillX1[j], (GLdouble) FloodFillY1[j]);
                                    glVertex2d((GLdouble) FloodFillX2[j], (GLdouble) FloodFillY2[j]);
                                    glVertex2d((GLdouble) Node2D[Node3].Coordinate[0], (GLdouble) Node2D[Node3].Coordinate[1]);
                                    glVertex2d((GLdouble) Node2D[Node1].Coordinate[0], (GLdouble) Node2D[Node1].Coordinate[1]);
                                    glEnd();
                                    break;
                            }
                        }
                    }
                }

                /* Free the Memory */
                free(FloodFillX1);
                free(FloodFillX2);
                free(FloodFillY1);
                free(FloodFillY2);
                free(FloodColors);
                free(FloodIndex);
                free(FloodEdge1);
                free(FloodEdge2);

                break;
            case CellCenter:
                /* Initializing tmp variable */
                tmp.Size = NoNodes2D;
                /* Allocate Memory to tmp variable */
                tmp.Data = (double *) malloc(NoNodes2D * sizeof (double));
                if (tmp.Data == NULL) {
                    /* Insert Message Code */
                    break;
                }

                /* Interpolating Solution Data to Node */
                /* Arithmetic Interpolation for Time being */
                ArithmeticCC2Node_2D(Solution2D.Sols[SolID], &tmp);

                /* For Each Level Value */
                for (i = 0; i < Level; i++) {
                    /* Set Color for Level */
                    color_index = (int) (((float) i / (float) Level) * 255);
                    glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
                    glLineWidth(1.0);

                    for (j = 0; j < NoCells2D; j++) {
                        glBegin(GL_LINE_STRIP);
                        for (k = 0, m = (Cell2D[j].NodesPerCell - 1); k < Cell2D[j].NodesPerCell; m = k, k++) {
                            /* Get Node IDs */
                            Node1 = Cell2D[j].ConnectNode[m];
                            Node2 = Cell2D[j].ConnectNode[k];

                            /* Get Node Solution Value */
                            NSol1 = tmp.Data[Node1];
                            NSol2 = tmp.Data[Node2];

                            /* Check if the Edge Contains a Contour point */
                            if (((NSol1 <= InterValueSol) && (NSol2 > InterValueSol)) || ((NSol2 <= InterValueSol) && (NSol1 > InterValueSol))) {
                                /* Linear Interpolation */
                                /* Get X coordinate */
                                ConX = ((Node2D[Node1].Coordinate[0]*(NSol2 - InterValueSol)) / (NSol2 - NSol1))
                                        + ((Node2D[Node2].Coordinate[0]*(InterValueSol - NSol1)) / (NSol2 - NSol1));

                                /* Get Y coordinate */
                                ConY = ((Node2D[Node1].Coordinate[1]*(NSol2 - InterValueSol)) / (NSol2 - NSol1))
                                        + ((Node2D[Node2].Coordinate[1]*(InterValueSol - NSol1)) / (NSol2 - NSol1));

                                /* Display Point */
                                glVertex2d((GLdouble) ConX, (GLdouble) ConY);
                            }
                        }
                        glEnd();
                    }
                    /* Update InterValueSol for Next Level */
                    InterValueSol = InterValueSol + DelValueSol;
                }
                /* Free tmp variable memory */
                free(tmp.Data);
                break;
        }

        /* Draw Boundary */
        DisplayBoundary_2D();
        glEndList();
    }

    /* Draw the legend of color */
    if (LegendActivateFlag)
        DrawLegend_2D(MinValueSol, MaxValueSol, Level, -1);

    DrawName_2D(SolID);
    return 0;
}

/*--------------------------------------------------------------*/
void ContourPlot_2D(void) {
    CreateContourPlot_2D(ContourPlotArg, ContourLevel);
}
