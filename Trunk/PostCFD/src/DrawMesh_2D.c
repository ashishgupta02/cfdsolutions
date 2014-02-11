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
 * File		DrawMesh_2D.c
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
#include "PPOptions_2D.h"
#include "Graphics_2D.h"
#include "ColorCode_2D.h"
#include "ColorPlot_2D.h"
#include "Boundary_2D.h"
#include "DrawMesh_2D.h"

double MinMaxX[2]; /* [O] = minvalue, [1] = maxvalue */
double MinMaxY[2]; /* [O] = minvalue, [1] = maxvalue */

/* Use to Preprocess the GL Commands for faster display */
static GLuint MeshList;
static GLuint MeshScalarList;

static int notfirst = 0;
static int notfirstscalar = 0;

static int MeshScalarID = -1;

/*--------------------------------------------------------------*/
void CoordinateMinMax_2D(void) {
    int inode;

    MinMaxX[0] = Node2D[0].Coordinate[0];
    MinMaxX[1] = Node2D[0].Coordinate[0];
    MinMaxY[0] = Node2D[0].Coordinate[1];
    MinMaxY[1] = Node2D[0].Coordinate[1];

    for (inode = 1; inode < NoNodes2D; inode++) {
        if (Node2D[inode].Coordinate[0] < MinMaxX[0])
            MinMaxX[0] = Node2D[inode].Coordinate[0];

        if (Node2D[inode].Coordinate[0] > MinMaxX[1])
            MinMaxX[1] = Node2D[inode].Coordinate[0];

        if (Node2D[inode].Coordinate[1] < MinMaxY[0])
            MinMaxY[0] = Node2D[inode].Coordinate[1];

        if (Node2D[inode].Coordinate[1] > MinMaxY[1])
            MinMaxY[1] = Node2D[inode].Coordinate[1];
    }
}

/*--------------------------------------------------------------*/
void DrawMesh_2D(void) {
    int icell, jnode;
    int node[3];

    /* Reset the Screen window */
    DrawInitial_2D();

    /* -2 is to Display Mesh Name */
    if (!StreamlineFlag)
        DrawName_2D(-2);

    /* Initialize Drawing Area */
    InitializeDrawArea_2D();

    /* Display Preprocessed List */
    if (notfirst) {
        glCallList(MeshList);
    } else {
        notfirst = 1;
        /* Create New List */
        MeshList = glGenLists(1);

        glNewList(MeshList, GL_COMPILE_AND_EXECUTE);
        glLineWidth(1.0);
        glColor3f(1.0f, 1.0f, 1.0f);

        for (icell = 0; icell < NoCells2D; icell++) {
            for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++)
                node[jnode] = Cell2D[icell].ConnectNode[jnode];

            glBegin(GL_LINE_LOOP);
            glVertex2d((GLdouble) Node2D[node[0]].Coordinate[0], (GLdouble) Node2D[node[0]].Coordinate[1]);
            glVertex2d((GLdouble) Node2D[node[1]].Coordinate[0], (GLdouble) Node2D[node[1]].Coordinate[1]);
            glVertex2d((GLdouble) Node2D[node[2]].Coordinate[0], (GLdouble) Node2D[node[2]].Coordinate[1]);
            glEnd();
        }

        glEndList();
    }

    /* Display Boundaries */
    if (BoundaryActivateFlag)
        DisplayBoundary_2D();

    glFlush();
}

/*--------------------------------------------------------------*/
void MeshScalarPlot_2D(int SolID) {
    int Color_Index;
    int icell, jnode;
    int node[3];
    double MinValueSol, MaxValueSol;

    /* Initialize Color */
    ColorEncode_2D();

    /* Get Min-Max Solution */
    SolutionMinMax_2D(SolID, &MinValueSol, &MaxValueSol);

    /* Reset the Screen Window */
    DrawInitial_2D();

    /* Display Scalar Name */
    DrawName_2D(SolID);

    /* Draw the legend of color */
    if (LegendActivateFlag)
        DrawLegend_2D(MinValueSol, MaxValueSol, 1, SolID);

    /* Initialize Drawing Area */
    InitializeDrawArea_2D();

    /* Display Preprocessed List */
    if (MeshScalarID == SolID) {
        glCallList(MeshScalarList);
    } else {
        MeshScalarID = SolID;

        if (notfirstscalar) {
            /* Delete Old List */
            glDeleteLists(MeshScalarList, 1);

            /* Create New List */
            MeshScalarList = glGenLists(1);
        } else {
            /* Create New List */
            MeshScalarList = glGenLists(1);
            notfirstscalar = 1;
        }

        glNewList(MeshScalarList, GL_COMPILE_AND_EXECUTE);

        switch (Solution2D.Location) {
            case Vertex:
                for (icell = 0; icell < NoCells2D; icell++) {
                    for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++)
                        node[jnode] = Cell2D[icell].ConnectNode[jnode];

                    glShadeModel(GL_SMOOTH);
                    glBegin(GL_LINE_LOOP);
                    for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++) {
                        /* Assigning Color for Solution Value */
                        Color_Index = (int) ((Solution2D.Sols[SolID]->Data[node[jnode]] - MinValueSol) / (MaxValueSol - MinValueSol) * 255);

                        /* Setting Drawing Color */
                        glColor3f(Red[Color_Index], Green[Color_Index], Blue[Color_Index]);
                        glVertex2d((GLdouble) Node2D[node[jnode]].Coordinate[0], (GLdouble) Node2D[node[jnode]].Coordinate[1]);
                    }
                    glEnd();
                }
                break;
            case CellCenter:
                for (icell = 0; icell < NoCells2D; icell++) {
                    for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++)
                        node[jnode] = Cell2D[icell].ConnectNode[jnode];

                    glShadeModel(GL_SMOOTH);
                    /* Assigning Color for Solution Value */
                    Color_Index = (int) ((Solution2D.Sols[SolID]->Data[icell] - MinValueSol) / (MaxValueSol - MinValueSol) * 255);
                    /* Setting Drawing Color */
                    glColor3f(Red[Color_Index], Green[Color_Index], Blue[Color_Index]);

                    glBegin(GL_LINE_LOOP);
                    for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++)
                        glVertex2d((GLdouble) Node2D[node[jnode]].Coordinate[0], (GLdouble) Node2D[node[jnode]].Coordinate[1]);
                    glEnd();
                }
                break;
        }
        glEndList();
    }

    /* Display Boundaries */
    if (BoundaryActivateFlag)
        DisplayBoundary_2D();

    glFlush();
}

/*--------------------------------------------------------------*/
void MeshScalar_2D(void) {
    MeshScalarPlot_2D(MeshScalarArg);
}

/*--------------------------------------------------------------*/
void Mesh_2D(void) {
    int icell, jnode;
    int node[3];

    glLineWidth(1.0);
    for (icell = 0; icell < NoCells2D; icell++) {
        for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++)
            node[jnode] = Cell2D[icell].ConnectNode[jnode];

        glColor3f(0.0, 0.2, 0.0);

        glBegin(GL_LINE_LOOP);
        glVertex2d((GLdouble) Node2D[node[0]].Coordinate[0], (GLdouble) Node2D[node[0]].Coordinate[1]);
        glVertex2d((GLdouble) Node2D[node[1]].Coordinate[0], (GLdouble) Node2D[node[1]].Coordinate[1]);
        glVertex2d((GLdouble) Node2D[node[2]].Coordinate[0], (GLdouble) Node2D[node[2]].Coordinate[1]);
        glEnd();
    }
}
