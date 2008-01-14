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
 * File		VectorPlot_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include <GL/glut.h>
#include <cgnslib.h>
#include "PostProcessing_2D.h"
#include "PrePostProcessing_2D.h"
#include "PPInitializeSolution_2D.h"
#include "PPOptions_2D.h"
#include "Graphics_2D.h"
#include "ColorCode_2D.h"
#include "DrawMesh_2D.h"
#include "Boundary_2D.h"
#include "ColorPlot_2D.h"
#include "VectorPlot_2D.h"

/* Use to Preprocess the GL Commands for faster display */
static GLuint VectorList;

/* -1 is because sol starts from 0 */
static int VectorPlotIDX = -1;
static int VectorPlotIDY = -1;

static int notfirst = 0;

/*--------------------------------------------------------------*/
int CreateVectorPlot_2D (int SolIDX, int SolIDY) {
	double MinEdge, MaxEdge;
	double MinValueSol, MaxValueSol;
	double SolX, SolY;
	/* Arrow Variables */
	double X1, Y1, X2, Y2, X3, Y3, X4, Y4;
	double VectorLength, VectorAngle, ArrowAngle = 0.5, temp;
	int i, Color_Index;
	
	/* Disable Streamline Flag */
	if (StreamlineFlag == 1)
		StreamlineFlag = 0;

	//ColorPlot_2D();
	/* Reset the Screen window */
	DrawInitial_2D();
	
	/* Initialize Drawing Area */
	InitializeDrawArea_2D();
	
	/* Solution Vector Magnitude Min-Max */
	SolutionMinMaxVector_2D(SolIDX, SolIDY, &MinValueSol, &MaxValueSol);
	
	/* Initialize Color */
	ColorEncode_2D();
	
	/* Display boundaries */
	if(BoundaryActivateFlag)
		DisplayBoundary_2D();

	/* Draw Mesh */		
	if (MeshActivateFlag)
		Mesh_2D();
				
	/* Get Min-Max Edge Value */
	MinMaxEdge_2D(&MinEdge, &MaxEdge);
	
	if (MaxEdge > (2*MinEdge)) {
		if (MaxEdge > (5*MinEdge)) {
			temp = MinEdge;
			MinEdge = (MinEdge + MaxEdge)/5;
			MaxEdge = (temp + MaxEdge)/5 + temp;
		}
		else {
		 	MaxEdge = 2 * MinEdge;
		}
	}
	
	/* Display Preprocessed List */
	if((VectorPlotIDX == SolIDX) && (VectorPlotIDY == SolIDY)) {
		glCallList(VectorList);
	}
	else {
		VectorPlotIDX = SolIDX;
		VectorPlotIDY = SolIDY;
		
		if (notfirst) {
			/* Delete Old List */ 
			glDeleteLists(VectorList, 1);
			
			/* Create New List */
			VectorList = glGenLists(1);
		}
		else {
			/* Create New List */
			VectorList = glGenLists(1);
			notfirst = 1;
		}
				
		glNewList(VectorList, GL_COMPILE_AND_EXECUTE);		
		switch (Solution2D.Location) {
			case Vertex:
				for (i = 0; i < NoNodes2D; i++) {
					if (Node2D[i].Flag)
						continue;
					
					SolX = Solution2D.Sols[SolIDX]->Data[i];
					SolY = Solution2D.Sols[SolIDY]->Data[i];
					
					/* Get Color_Index */
					Color_Index = (int)(((sqrt(SolX*SolX + SolY*SolY) - MinValueSol)/(MaxValueSol - MinValueSol)) * 255);
					
					/* Get Vector Arrow Length */
					VectorLength = MinEdge + (((sqrt(SolX*SolX + SolY*SolY) - MinValueSol)/(MaxValueSol - MinValueSol)) * (MaxEdge - MinEdge));
								
					VectorAngle = atan(fabs(SolY/SolX));
					
					X1 = Node2D[i].Coordinate[0];
					Y1 = Node2D[i].Coordinate[1];
					
					X2 = VectorLength * cos(VectorAngle);
					Y2 = VectorLength * sin(VectorAngle);
					
					X3 = (0.2 * VectorLength) * cos(VectorAngle - ArrowAngle);
					Y3 = (0.2 * VectorLength) * sin(VectorAngle - ArrowAngle);
					
					X4 = (0.2 * VectorLength) * cos(VectorAngle + ArrowAngle);
					Y4 = (0.2 * VectorLength) * sin(VectorAngle + ArrowAngle);
					
					if ((SolX < 0) && (SolY >= 0)) {
						X2 = -X2;
						X3 = -X3;
						X4 = -X4;
					}
					if ((SolX < 0) && (SolY < 0)) {
						X2 = -X2;
						X3 = -X3;
						X4 = -X4;
						Y2 = -Y2;
						Y3 = -Y3;
						Y4 = -Y4;
					}
					if ((SolX >= 0) && (SolY < 0)) {
						Y2 = -Y2;
						Y3 = -Y3;
						Y4 = -Y4;
					}
					
					/* Setting Drawing Color */
					glColor3f(Red[Color_Index], Green[Color_Index], Blue[Color_Index]);
					//glColor3f(0.0f, 0.0f, 0.0f);
					
					/* Draw Arrow */
					glLineWidth(0.05);
					glBegin(GL_LINE_STRIP);
						glVertex2d((GLdouble)X1, (GLdouble)Y1);
						glVertex2d((GLdouble)(X1 + X2), (GLdouble)(Y1 + Y2));
						glVertex2d((GLdouble)(X1 + X2), (GLdouble)(Y1 + Y2));
						glVertex2d((GLdouble)(X1 + X2 - X3), (GLdouble)(Y1 + Y2 - Y3));
						glVertex2d((GLdouble)(X1 + X2), (GLdouble)(Y1 + Y2));
						glVertex2d((GLdouble)(X1 + X2 - X4), (GLdouble)(Y1 + Y2 - Y4));
					glEnd();
				}
				break;
			case CellCenter:
				for (i = 0; i < NoCells2D; i++) {
					SolX = Solution2D.Sols[SolIDX]->Data[i];
					SolY = Solution2D.Sols[SolIDY]->Data[i];
					
					/* Get Color_Index */
					Color_Index = (int)(((sqrt(SolX*SolX + SolY*SolY) - MinValueSol)/(MaxValueSol - MinValueSol)) * 255);
					
					/* Setting Drawing Color */
					glColor3f(Red[Color_Index], Green[Color_Index], Blue[Color_Index]);
	
					/* Get Vector Arrow Length */
					VectorLength = MinEdge + (((sqrt(SolX*SolX + SolY*SolY) - MinValueSol)/(MaxValueSol - MinValueSol)) * (MaxEdge - MinEdge));
					
					VectorAngle = atan(fabs(SolY/SolX));
					
					X1 = Cell2D[i].Centroid[0];
					Y1 = Cell2D[i].Centroid[1];
					
					X2 = VectorLength * cos(VectorAngle);
					Y2 = VectorLength * sin(VectorAngle);
					
					X3 = (0.2 * VectorLength) * cos(VectorAngle - ArrowAngle);
					Y3 = (0.2 * VectorLength) * sin(VectorAngle - ArrowAngle);
					
					X4 = (0.2 * VectorLength) * cos(VectorAngle + ArrowAngle);
					Y4 = (0.2 * VectorLength) * sin(VectorAngle + ArrowAngle);
					
					if ((SolX < 0) && (SolY >= 0)) {
						X2 = -X2;
						X3 = -X3;
						X4 = -X4;
					}
					if ((SolX < 0) && (SolY < 0)) {
						X2 = -X2;
						X3 = -X3;
						X4 = -X4;
						Y2 = -Y2;
						Y3 = -Y3;
						Y4 = -Y4;
					}
					if ((SolX >= 0) && (SolY < 0)) {
						Y2 = -Y2;
						Y3 = -Y3;
						Y4 = -Y4;
					}
					
					/* Draw Arrow */
					glLineWidth(0.05);
					glBegin(GL_LINES);
						glVertex2d(X1, Y1);
						glVertex2d((X1 + X2), (Y1 + Y2));
						glVertex2d((X1 + X2), (Y1 + Y2));
						glVertex2d((X1 + X2 - X3), (Y1 + Y2 - Y3));
						glVertex2d((X1 + X2), (Y1 + Y2));
						glVertex2d((X1 + X2 - X4), (Y1 + Y2 - Y4));
					glEnd();				
				}
				break;
		}
		glEndList();
	}
	
	/* Put Name of color plot */
	DrawName_2D(SolIDX);
		
	/* Draw the legend of color */
	if(LegendActivateFlag)
		DrawLegend_2D(MinValueSol, MaxValueSol, 1, 1);
	
	return 0;
}

/*--------------------------------------------------------------*/
void VectorPlot_2D (void) {
	CreateVectorPlot_2D(VectorPlotArg1, VectorPlotArg2);
}
