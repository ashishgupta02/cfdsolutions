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
 * File		Graphics_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */

/* User Defined */
#include "Global.h"
#include <GL/glut.h>
#include "CellPointOperation_2D.h"
#include "DrawMesh_2D.h"
#include "Boundary_2D.h"
#include "PPOptions_2D.h"
#include "Graphics_2D.h"
#include "ColorPlot_2D.h"
#include "Streamline_2D.h"

/* Extern from Graphics_2D.h */
int MainWindowID;

/* Size of Window */
int windowHeight = 600;
int windowWidth = 600;

/* Position of Window */
int winPositionX = 100;
int winPositionY = 150;

/* World Aspect Ratio */
double WorldAspect;

/*--------------------------------------------------------------*/
static void MouseClick_2D(int button, int state, int x, int y) {
    GLint viewport[4];
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint tmp_y;
    /* World Coordinate */
    GLdouble wx, wy, wz;

    /* If in Streamline Mode */
    if (StreamlineFlag) {
        switch (button) {
            case GLUT_LEFT_BUTTON:
                if (state == GLUT_DOWN) {
                    glGetIntegerv(GL_VIEWPORT, viewport);
                    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
                    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);

                    tmp_y = windowHeight - (GLint) y;

                    /* Get World Coordinates */
                    gluUnProject((GLdouble) x, (GLdouble) tmp_y, 0.0, modelMatrix, projMatrix, viewport, &wx, &wy, &wz);

                    /* Initialize the SeedPoint Array */
                    InitializeSeedPoint_2D(wx, wy);
                }
                break;
            case GLUT_MIDDLE_BUTTON:
                break;
        }
    }
}

/*--------------------------------------------------------------*/
static void MouseMove_2D(int x, int y) {
    x = 0;
    y = 0;
}

/*--------------------------------------------------------------*/

/* font rendering-- use stroked fonts */
void PrintString_2D(void *font, char *str) {
    int i, l = strlen(str);

    for (i = 0; i < l; i++)
        glutBitmapCharacter(font, *str++);
}

/*--------------------------------------------------------------*/
static void DisplayInit_2D(void) {
    glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glFlush();
}

/*--------------------------------------------------------------*/
static void Reshape_2D(GLsizei width, GLsizei height) {
    windowWidth = width;
    windowHeight = height;
}

/*--------------------------------------------------------------*/
void DrawAxis_2D(void) {
    glViewport(0, 0, (windowWidth) / 10, (windowHeight) / 10);

    /* Reset Projection Matrix Stack */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    /* 1.05 is multiplied to give display tolerance */
    gluOrtho2D(0, 10, 0, 10);

    /* Reset Model View Matrix Stack */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glLineWidth(1.0);
    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_LINES);
    glVertex2d(1.0f, 1.0f);
    glVertex2d(9.0f, 1.0f);
    glEnd();
    glColor3f(0.0f, 0.0f, 1.0f);
    glBegin(GL_LINES);
    glVertex2d(1.0f, 1.0f);
    glVertex2d(1.0f, 9.0f);
    glEnd();
    glFlush();
}

/*--------------------------------------------------------------*/
void DrawInitial_2D(void) {
    /* Reset Screen Window */
    DisplayInit_2D();

    if (!StreamlineFlag) {
        /* Draw Axis */
        DrawAxis_2D();

        /* Initialize Name Area */
        DrawName_2D(-1);
    }
}

/*--------------------------------------------------------------*/
void InitializeDrawArea_2D(void) {
    double RangeX, RangeY;
    GLint x, y, width, height;

    RangeX = MinMaxX[1] - MinMaxX[0];
    RangeY = MinMaxY[1] - MinMaxY[0];

    //switch (StreamlineFlag) {
    //case 0:
    //	/* StreamlineFlag Disable */
    if ((RangeX < RangeY) && (windowWidth < windowHeight)) {
        x = (GLint) (windowWidth * (11 - (9 * WorldAspect)) / 20);
        y = (GLint) (windowHeight / 10);
        width = (GLint) ((9 * windowWidth * WorldAspect) / 10);
        height = (GLint) ((9 * windowHeight) / 10);
    } else {
        if (WorldAspect > windowWidth / windowHeight) {
            x = (GLint) (windowWidth / 10);
            y = (GLint) (windowHeight * (11 - (9 / WorldAspect)) / 20);
            width = (GLint) ((9 * windowWidth) / 10);
            height = (GLint) ((9 * (windowHeight / WorldAspect)) / 10);
        } else {
            x = (GLint) (((11 * windowWidth) - 9 * (windowHeight * WorldAspect)) / 20);
            y = (GLint) (windowHeight / 10);
            width = (GLint) ((9 * windowHeight * WorldAspect) / 10);
            height = (GLint) ((9 * windowHeight) / 10);
        }
    }
    //	break;
    //case 1:
    //	/* StreamlineFlag Enable */
    //	if ((RangeX < RangeY) && (windowWidth < windowHeight)) {
    //		x	 = (GLint)((windowWidth*(1 - WorldAspect))/2);
    //		y	 = (GLint)(0);
    //		width	 = (GLint)(windowWidth * WorldAspect);
    //		height	 = (GLint)(windowHeight);
    //	}
    //	else {
    //		if (WorldAspect > windowWidth/windowHeight) {
    //			x	 = (GLint)(0);
    //			y	 = (GLint)((windowHeight*(1 - (1/WorldAspect)))/2);
    //			width	 = (GLint)(windowWidth);
    //			height	 = (GLint)(windowHeight/WorldAspect);
    //		}
    //		else {
    //			x	 = (GLint)((windowWidth - (windowHeight * WorldAspect))/2);
    //			y	 = (GLint)(0);
    //			width	 = (GLint)(windowHeight * WorldAspect);
    //			height	 = (GLint)(windowHeight);
    //		}
    //	}
    //	break;
    //}

    glViewport(x, y, width, height);

    /* Reset the projection matrix stack */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    /* 1.05 is multiplied to give display tolerance */
    gluOrtho2D((MinMaxX[0]*1.05), (MinMaxX[1]*1.05), (MinMaxY[0]*1.05), (MinMaxY[1]*1.05));

    /* Reset the model view matrix stack */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/*--------------------------------------------------------------*/
void Graphics_2D(void) {
    char window_title [ 80 ];
    strcpy(window_title, "PostCFD");

    /* Get Min-Max of World Window */
    CoordinateMinMax_2D();
    printf("MinX = %lf, MaxX = %lf, MinY = %lf, MaxY = %lf\n", MinMaxX[0], MinMaxX[1], MinMaxY[0], MinMaxY[1]);

    printf("Building Display List.......\n");

    /* Setting World Aspect Ratio */
    WorldAspect = (MinMaxX[1] - MinMaxX[0]) / (MinMaxY[1] - MinMaxY[0]);

    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(winPositionX, winPositionY);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInit(pargc, pargv);
    MainWindowID = glutCreateWindow(window_title);
    MainMenu_2D();
    glutReshapeFunc(Reshape_2D);
    glutDisplayFunc(DrawMesh_2D);
    glutMouseFunc(MouseClick_2D);
    glutMotionFunc(MouseMove_2D);
    glutMainLoop();
}
