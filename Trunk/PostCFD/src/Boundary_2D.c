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
 * File		Boundary_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "Error.h"
#include <cgnslib.h>
#include <GL/glut.h>
#include "CGNSIO.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PPOptions_2D.h"
#include "Boundary_2D.h"

BoundaryCondition_2D_Un Boundary2D;

/* Boundary Name */
char BoCoName[26][33] = {
	"BCTypeNull", "BCTypeUserDefined",
	"BCAxisymmetricWedge", "BCDegenerateLine", "BCDegeneratePoint",
	"BCDirichlet", "BCExtrapolate", "BCFarfield", "BCGeneral", "BCInflow",
	"BCInflowSubsonic", "BCInflowSupersonic", "BCNeumann", "BCOutflow",
	"BCOutflowSubsonic", "BCOutflowSupersonic", "BCSymmetryPlane",
	"BCSymmetryPolar", "BCTunnelInflow", "BCTunnelOutflow", "BCWall",
	"BCWallInviscid", "BCWallViscous", "BCWallViscousHeatFlux",
	"BCWallViscousIsothermal", "FamilySpecified"
};

/*---------------------------------------------------------------*/
int InitializeBoundaryCondition_2D (ZONE *P) {
	int i, bc, node1, node2, edgeID;
	BOCO *boco;
	
	Boundary2D.NoBoCo = P->nbocos;
	Boundary2D.BoCo = (Boundary_2D_Un *) malloc(Boundary2D.NoBoCo * sizeof(Boundary_2D_Un));
	if (Boundary2D.BoCo == NULL) {
		MSG("InitializeBoundaryCondition_2D: Failed 1");
		return 1;
	}
	
	for (i = 0; i < NoNodes2D; i++)
		Node2D[i].Flag = 0;
	
	for (boco = P->bocos, bc = 1; bc <= P->nbocos; bc++, boco++) {
		if (boco->ptype != PointList) {
			MSG("InitializeBoundaryCondition: Error - Expecting Point List");
			return 1;
		}
		
		if((boco->type >= 0) || (boco->type <= 25)) {
			strcpy(Boundary2D.BoCo[bc-1].Name, boco->name);
			strcpy(Boundary2D.BoCo[bc-1].TypeName, BoCoName[boco->type]);
			Boundary2D.BoCo[bc-1].Type = boco->type;
			printf("BoCo = %s\tType = %d\tType Name = %s\n", Boundary2D.BoCo[bc-1].Name, Boundary2D.BoCo[bc-1].Type, Boundary2D.BoCo[bc-1].TypeName);
		}
		
		Boundary2D.BoCo[bc-1].Size = boco->npnts;
		Boundary2D.BoCo[bc-1].BPoints = (int *) malloc (boco->npnts * sizeof(int));
		if (Boundary2D.BoCo[bc-1].BPoints == NULL) {
			MSG("InitializeBoundaryCondition_2D: Failed 2");
			return 1;
		}
		
		for (i = 0; i < (boco->npnts); i++)
			Node2D[boco->pnts[i]].Flag = 1;
		
		
		for (i = 0; i < (boco->npnts - 1); i++) {
			node1 = boco->pnts[i];
			node2 = boco->pnts[i+1];
			Boundary2D.BoCo[bc-1].BPoints[i] = boco->pnts[i];
			Boundary2D.BoCo[bc-1].BPoints[i+1] = boco->pnts[i+1];
			i++;
			edgeID = GetEdgeNumber_2D(node1, node2); 
			Edge2D[edgeID].Flag = boco->type;
		}
	}
	
	return 0;
}

/*---------------------------------------------------------------*/
static void BoundaryColor_2D (float* BColor, int Type) {
	BColor[0] = (float) (Type)/26;
	BColor[1] = (float) (26 - Type)/(Type + 26);
	BColor[2] = (float)(26 - Type)/26;
}

/*---------------------------------------------------------------*/
void DisplayBoundary_2D (void) {
	unsigned int i, j;
	int nodeid[2];
	float BColor[3];
	
	for (i = 0; i < Boundary2D.NoBoCo; i++) {
		BoundaryColor_2D(BColor, Boundary2D.BoCo[i].Type);
		
		glColor3f(BColor[0], BColor[1], BColor[2]);
		glLineWidth(2.0);
		glBegin(GL_LINES);
		
		for (j = 0; j < (Boundary2D.BoCo[i].Size - 1); j++) {
			nodeid[0] = Boundary2D.BoCo[i].BPoints[j];
			nodeid[1] = Boundary2D.BoCo[i].BPoints[j+1];
			j++;
			glVertex2d(Node2D[nodeid[0]].Coordinate[0], Node2D[nodeid[0]].Coordinate[1]);
			glVertex2d(Node2D[nodeid[1]].Coordinate[0], Node2D[nodeid[1]].Coordinate[1]);
		}
		glEnd();
	}
}
