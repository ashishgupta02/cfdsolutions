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
 * File		ColorPlot_2D.c
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
#include "Graphics_2D.h"
#include "DrawMesh_2D.h"
#include "ColorCode_2D.h"
#include "PPOptions_2D.h"
#include "ColorPlot_2D.h"
#include "Boundary_2D.h"

/* Use to Preprocess the GL Commands for faster display */
static GLuint ColorList;

/* Use to Preprocess the GL Commands for Legend display */
static GLuint LegendList;
static int LegendID = -1;

/* -1 is because sol starts from 0 */
static int ColorPlotID = -1;

/* for color */
static int notfirst = 0;
/* for legend */
static int notfirst1 = 0;

/*--------------------------------------------------------------*/
void DrawName_2D(int SolutionLocation) {
    char buff[33];

    glViewport(windowWidth / 10, 0, (9 * windowWidth) / 10, (windowHeight) / 10);

    /* Reset Projection Matrix Stack */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluOrtho2D(0.0, 100.0, 0.0, 10.0);

    /* Reset Model View Matrix Stack */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (SolutionLocation != -1) {
        /* Set Color */
        glColor3f(1.0f, 1.0f, 1.0f);

        if (SolutionLocation == -2)
            /* For Printing Mesh Type*/
            sprintf(buff, "%s", "Triangular Mesh: CGNS = TRI3");
        else
            /* Print Values of Max, Min Solution */
            sprintf(buff, "%s", Solution2D.Sols[SolutionLocation]->Name);

        /* Set the Raster Position */
        glRasterPos2d(20.0, 2.0);
        PrintString_2D(GLUT_BITMAP_HELVETICA_18, buff);
    }

    glFlush();
}

/*--------------------------------------------------------------*/
void DrawLegend_2D(double MinValue, double MaxValue, int levels, int SolLocation) {
    double left, right, bottom, top;
    double del_y; /* y step */
    double y1, y2; /* y coords */
    int color_index; /* color code */
    char buff[20]; /* buffer for printing */
    int i;

    /* Box Dimensions */
    left = 0.0;
    right = 10.0;
    bottom = 0.0;
    top = 100.0;

    glViewport(0, windowHeight / 10, (windowWidth) / 10, (9 * windowHeight) / 10);

    /* Reset Projection Matix */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluOrtho2D(left, right, bottom, top);

    /* Reset Model View Matrix Stack */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (levels <= 0)
        levels = 10;

    if (levels >= 20)
        levels = 20;

    /* For Color Plot */
    if (SolLocation != -1)
        levels = 20;

    if (notfirst1) {
        if (LegendID == SolLocation) {
            glCallList(LegendList);
        } else {
            LegendID = SolLocation;

            /* Delete Old List */
            glDeleteLists(LegendList, 1);

            /* Create New List */
            LegendList = glGenLists(1);

            glNewList(LegendList, GL_COMPILE_AND_EXECUTE);
            del_y = ((top - 5.0) - (bottom + 5.0)) / (double) levels;

            y1 = del_y + bottom + 5.0;
            y2 = bottom + 5.0;

            /* draw a horizontal divider (a tick) */
            glColor3f(0.0f, 0.7f, 0.7f);
            glBegin(GL_LINES);
            glVertex2d(3.0, y2);
            glVertex2d(4.0, y2);
            glEnd();

            for (i = 0; i < levels; i++) {
                color_index = (int) (((float) i / (float) levels) * 255);

                glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
                glBegin(GL_QUADS);
                glVertex2d(1.0, y1);
                glVertex2d(3.0, y1);
                glVertex2d(3.0, y2);
                glVertex2d(1.0, y2);
                glEnd();

                /* draw a horizontal divider (a tick) */
                glColor3f(0.0f, 0.7f, 0.7f);
                glBegin(GL_LINES);
                glVertex2d(3.0, y1);
                glVertex2d(4.0, y1);
                glEnd();

                y1 += del_y;
                y2 += del_y;

            }

            y2 = bottom + 5.0;

            /* Printing Values */
            glColor3f(1.0f, 1.0f, 1.0f);
            for (i = 0; i <= levels; i++) {
                /* Print Values of Max, Min Solution */
                glColor3f(1.0f, 1.0f, 1.0f);
                sprintf(buff, "%lf", MinValue + ((MaxValue - MinValue) / levels) * i);
                /* Set the Raster Position */
                glRasterPos2d(4.0, y2);
                PrintString_2D(GLUT_BITMAP_HELVETICA_10, buff);
                y2 += del_y;
            }
            glEndList();
        }
    } else {
        /* Create New List */
        LegendList = glGenLists(1);
        notfirst1 = 1;

        LegendID = SolLocation;

        glNewList(LegendList, GL_COMPILE_AND_EXECUTE);
        del_y = ((top - 5.0) - (bottom + 5.0)) / (double) levels;

        y1 = del_y + bottom + 5.0;
        y2 = bottom + 5.0;

        /* draw a horizontal divider (a tick) */
        glColor3f(0.0f, 0.7f, 0.7f);
        glBegin(GL_LINES);
        glVertex2d(3.0, y2);
        glVertex2d(4.0, y2);
        glEnd();

        for (i = 0; i < levels; i++) {
            color_index = (int) (((float) i / (float) levels) * 255);

            glColor3f(Red[color_index], Green[color_index], Blue[color_index]);
            glBegin(GL_QUADS);
            glVertex2d(1.0, y1);
            glVertex2d(3.0, y1);
            glVertex2d(3.0, y2);
            glVertex2d(1.0, y2);
            glEnd();

            /* draw a horizontal divider (a tick) */
            glColor3f(0.0f, 0.7f, 0.7f);
            glBegin(GL_LINES);
            glVertex2d(3.0, y1);
            glVertex2d(4.0, y1);
            glEnd();

            y1 += del_y;
            y2 += del_y;
        }

        y2 = bottom + 5.0;

        /* Printing Values */
        glColor3f(1.0f, 1.0f, 1.0f);
        for (i = 0; i <= levels; i++) {
            /* Print Values of Max, Min Solution */
            glColor3f(1.0f, 1.0f, 1.0f);
            sprintf(buff, "%lf", MinValue + ((MaxValue - MinValue) / levels) * i);
            /* Set the Raster Position */
            glRasterPos2d(4.0, y2);
            PrintString_2D(GLUT_BITMAP_HELVETICA_10, buff);
            y2 += del_y;
        }
        glEndList();
    }
}

/*--------------------------------------------------------------*/
void CreateColorPlot_2D(int SolLocation) {
    int i, j, Color_Index, inode;
    double MinValueSol, MaxValueSol;

    /* Disable Streamline Flag */
    if (StreamlineFlag == 1)
        StreamlineFlag = 0;

    /* Initialize Color */
    ColorEncode_2D();

    /* Get Min-Max Solution */
    SolutionMinMax_2D(SolLocation, &MinValueSol, &MaxValueSol);

    /* Reset the Screen window */
    DrawInitial_2D();

    /* Put Name of color plot */
    DrawName_2D(SolLocation);

    /* Initialize Drawing Area */
    InitializeDrawArea_2D();

    /* Display Preprocessed List */
    if (ColorPlotID == SolLocation) {
        glCallList(ColorList);
    } else {
        ColorPlotID = SolLocation;

        if (notfirst) {
            /* Delete Old List */
            glDeleteLists(ColorList, 1);

            /* Create New List */
            ColorList = glGenLists(1);
        } else {
            /* Create New List */
            ColorList = glGenLists(1);
            notfirst = 1;
        }

        glNewList(ColorList, GL_COMPILE_AND_EXECUTE);

        switch (Solution2D.Location) {
            case Vertex:
                for (i = 0; i < NoCells2D; i++) {
                    glShadeModel(GL_SMOOTH);
                    glBegin(GL_POLYGON);
                    for (j = 0; j < Cell2D[i].NodesPerCell; j++) {
                        /* Getting Node ID */
                        inode = Cell2D[i].ConnectNode[j];
                        /* Assigning Color for Solution Value */
                        Color_Index = (int) ((Solution2D.Sols[SolLocation]->Data[inode] - MinValueSol) / (MaxValueSol - MinValueSol) * 255);
                        /* Setting Drawing Color */
                        glColor3f(Red[Color_Index], Green[Color_Index], Blue[Color_Index]);
                        glVertex2d((GLdouble) Node2D[inode].Coordinate[0], (GLdouble) Node2D[inode].Coordinate[1]);
                    }
                    glEnd();
                }
                break;
            case CellCenter:
                for (i = 0; i < NoCells2D; i++) {
                    glShadeModel(GL_SMOOTH);
                    glBegin(GL_POLYGON);
                    /* Assigning Color for Solution Value */
                    Color_Index = (int) ((Solution2D.Sols[SolLocation]->Data[i] - MinValueSol) / (MaxValueSol - MinValueSol) * 255);
                    /* Setting Drawing Color */
                    glColor3f(Red[Color_Index], Green[Color_Index], Blue[Color_Index]);
                    for (j = 0; j < Cell2D[i].NodesPerCell; j++) {
                        /* Getting Node ID */
                        inode = Cell2D[i].ConnectNode[j];
                        glVertex2d((GLdouble) Node2D[inode].Coordinate[0], (GLdouble) Node2D[inode].Coordinate[1]);
                    }
                    glEnd();
                }
                break;
        }
        glEndList();
    }

    /* Draw Boundary */
    if (BoundaryActivateFlag)
        DisplayBoundary_2D();

    /* Draw the legend of color */
    if (LegendActivateFlag)
        DrawLegend_2D(MinValueSol, MaxValueSol, 1, SolLocation);

    glFlush();
}

/*--------------------------------------------------------------*/
void ColorPlot_2D(void) {
    CreateColorPlot_2D(ColorPlotArg);
}
