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
 * File		Interpolation_2D.c
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

/* User Defined */
#include "Global.h"
#include "Error.h"
#include "PostProcessing_2D.h"
#include "PostProcessingInitialize_2D.h"
#include "PPInitializeSolution_2D.h"
#include "Interpolation_2D.h"

/* Global use to convert Cell-Center Solution to Vertex */
static int InterpolationOptionCC2Node = 0;

/* Interpolation Option For Cell-Center to Node */
static char *usgmsg[] = {
	"usage  : Interpolation [options]",
	"options:",
	"	1	= Weighted Area Averaging",
	"	2	= Inverse Distance Averaging",
	"	3	= Laplacian Averaging",
	"	4	= Arithmetic Averaging",
	NULL
};

/*---------- usage --------------------------------------------------
 * Display usage message and exit
 *-------------------------------------------------------------------*/

static int usage (char **usgmsg) {
	int n;

	for (n = 0; NULL != usgmsg[n]; n++)
		fprintf (stderr, "%s\n", usgmsg[n]);

	return 0;
}

/*---------------------------------------------------------------*/
int SetInterpolationOptionCC2Node_2D (void) {
	/* Getting Interpolation Option */
	printf("Cell Center To Node Interpolation:\n");
	usage(usgmsg);
	printf("Give Interpolation Module Option: ");
	scanf("%d", &InterpolationOptionCC2Node);
	
	return 0;
}

/*---------------------------------------------------------------*/
int WeightedAreaCC2Node_2D (Data_2D_Un *In, Data_2D_Un *Out) {
	int cellID;
	int inode, jcell;
	double weightCell;
	double data;
	
	for (inode = 0; inode < NoNodes2D; inode++ ) {
		Node2D[inode].WeightTotal = 0.0;
		
		Out->Data[inode] = 0.0;
		
		for (jcell = 0; jcell < NodeAdjacentCells2D[inode].Size; jcell++) {
			cellID = NodeAdjacentCells2D[inode].Data[jcell];
			
			weightCell = Cell2D[cellID].Area;
			
			Node2D[inode].WeightTotal += weightCell;
			
			data = In->Data[cellID];
			
			Out->Data[inode] += data * weightCell;	
		}
		
		Out->Data[inode] = Out->Data[inode]/Node2D[inode].WeightTotal;
	}
	
	return 0;
} 

/*---------------------------------------------------------------*/
int InverseDistanceCC2Node_2D(Data_2D_Un *In, Data_2D_Un *Out) {
	int cellID;
	int inode, jcell;
	double xNode, yNode, xCell, yCell, distance, invDistance;
	double weightCell;
	double data;
	
	for (inode = 0; inode < NoNodes2D; inode++ ) {
		Node2D[inode].WeightTotal = 0.0;
		
		Out->Data[inode] = 0.0;
		
		for (jcell = 0; jcell < NodeAdjacentCells2D[inode].Size; jcell++) {
			cellID = NodeAdjacentCells2D[inode].Data[jcell];
			
			xNode = Node2D[inode].Coordinate[0];
			yNode = Node2D[inode].Coordinate[1];
			xCell = Cell2D[cellID].Centroid[0];
			yCell = Cell2D[cellID].Centroid[1];
			distance = sqrt(pow(xNode - xCell, 2) + pow(yNode - yCell, 2));
			invDistance = 1.0/distance;
			
			weightCell = invDistance;
			
			Node2D[inode].WeightTotal += weightCell;
			
			data = In->Data[cellID];
			
			Out->Data[inode] += data * weightCell;
		}
		
		Out->Data[inode] = Out->Data[inode]/Node2D[inode].WeightTotal;
	}

	return 0;
}

/*---------------------------------------------------------------*/
int LaplacianCC2Node_2D (Data_2D_Un *In, Data_2D_Un *Out) {
	int cellID;
	int inode, jcell;
	double xNode, yNode, xCell, yCell, distance;
	double weightCell;
	
	for (inode = 0; inode < NoNodes2D; inode++ ) {
		Node2D[inode].WeightTotal = 0.0;
		
		Out->Data[inode] = 0.0;
		
		xNode = Node2D[inode].Coordinate[0];
		yNode = Node2D[inode].Coordinate[1];
		
		for (jcell = 0; jcell < NodeAdjacentCells2D[inode].Size; jcell++) {
			cellID = NodeAdjacentCells2D[inode].Data[jcell];
						
			xCell = Cell2D[cellID].Centroid[0];
			yCell = Cell2D[cellID].Centroid[1];
			distance = sqrt(pow(xNode - xCell, 2) + pow(yNode - yCell, 2));
			
			weightCell = distance;
			
			Node2D[inode].WeightTotal += weightCell;
			
			Out->Data[inode] += (In->Data[cellID] * weightCell);
		}
		
		Out->Data[inode] = Out->Data[inode]/Node2D[inode].WeightTotal;
	}

	return 0;
}

/*---------------------------------------------------------------*/
int ArithmeticCC2Node_2D (Data_2D_Un *In, Data_2D_Un *Out) {
	int cellID;
	int inode, jcell;
	
	for (inode = 0; inode < NoNodes2D; inode++ ) {
		
		Out->Data[inode] = 0.0;
		
		for (jcell = 0; jcell < NodeAdjacentCells2D[inode].Size; jcell++) {
			cellID = NodeAdjacentCells2D[inode].Data[jcell];
			
			Out->Data[inode] += (In->Data[cellID]);
		}
		
		Out->Data[inode] = Out->Data[inode]/NodeAdjacentCells2D[inode].Size;
	}
	
	return 0;
}

/*---------------------------------------------------------------*/
int CellCenter2Node_2D (Data_2D_Un *In, Data_2D_Un *Out) {
	if (InterpolationOptionCC2Node == 0) {
		/* Getting Interpolation Option */
		printf("Cell Center To Node Interpolation:\n");
		usage(usgmsg);
		printf("Give Interpolation Module Option: ");
		scanf("%d", &InterpolationOptionCC2Node);
	}
	
	switch (InterpolationOptionCC2Node) {
		case 1:
			/* Weighted Area Averaging */
			WeightedAreaCC2Node_2D(In, Out);
			break;
		case 2:
			/* Inverse Distance Averaging */
			InverseDistanceCC2Node_2D(In, Out);
			break;
		case 3:
			/* Laplacian Averaging */
			LaplacianCC2Node_2D(In, Out);
			break;
		case 4:
			/* Arithmetic Averaging */
			ArithmeticCC2Node_2D(In, Out);
			break;
		default:
			/* Loop */
			MSG("CellCenter2Node_2D: Invalid Interpolation Option");
			return 1;
	}
	return 0;
}

/*****************************************************************/
/*
	NOTE: Below Interpolation Options Works only for triangles
*/ 
/*---------------------------------------------------------------*/

/* Interpolation Option For Node to Cell Center*/
static char *usgmsg1[] = {
	"usage  : Interpolation [options]",
	"options:",
	"	1	= Linear",
	NULL
};

/*---------------------------------------------------------------*/
int LinearNode2CC_2D(Data_2D_Un *In, Data_2D_Un *Out) {
	int i, icell, jnode;
	int node[3];
	double sai[3], alpha[3], beta[3], gamma[3];
	
	for (icell = 0; icell < NoCells2D; icell++ ) {
		for (i = 0; i < 3; i++) {
			sai[i] = 0.0;
			alpha[i] = 0.0;
			beta[i] = 0.0;
			gamma[i] = 0.0;
		}
		
		for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++)
			node[i] = Cell2D[icell].ConnectNode[i];
		
		alpha[0] = (Node2D[node[1]].Coordinate[0] * Node2D[node[2]].Coordinate[1]) -
				(Node2D[node[2]].Coordinate[0] * Node2D[node[1]].Coordinate[1]); 
		alpha[1] = (Node2D[node[2]].Coordinate[0] * Node2D[node[0]].Coordinate[1]) -
				(Node2D[node[0]].Coordinate[0] * Node2D[node[2]].Coordinate[1]); 
		alpha[2] = (Node2D[node[0]].Coordinate[0] * Node2D[node[1]].Coordinate[1]) -
				(Node2D[node[1]].Coordinate[0] * Node2D[node[0]].Coordinate[1]); 
		
		beta[0] = Node2D[node[1]].Coordinate[1] - Node2D[node[2]].Coordinate[1];
		beta[1] = Node2D[node[2]].Coordinate[1] - Node2D[node[0]].Coordinate[1];
		beta[2] = Node2D[node[0]].Coordinate[1] - Node2D[node[1]].Coordinate[1];
		
		gamma[0] = Node2D[node[2]].Coordinate[0] - Node2D[node[1]].Coordinate[0];
		gamma[1] = Node2D[node[0]].Coordinate[0] - Node2D[node[2]].Coordinate[0];
		gamma[2] = Node2D[node[1]].Coordinate[0] - Node2D[node[0]].Coordinate[0];
		
		Out->Data[icell] = 0.0;
		
		for (jnode = 0; jnode < Cell2D[icell].NodesPerCell; jnode++) {
			sai[jnode] = (alpha[jnode] + 
				beta[jnode]*Cell2D[icell].Centroid[0] + 
				gamma[jnode]*Cell2D[icell].Centroid[1])/(2*Cell2D[icell].Area);
			
			Out->Data[icell]  += (In->Data[node[jnode]] * sai[jnode]);
		}
	}
	
	return 0;
}

/*---------------------------------------------------------------*/
int Node2CellCenter_2D (Data_2D_Un *In, Data_2D_Un *Out) {
	int select;
	
	/* Getting Interpolation Options */
	printf("Node To Cell Center Interpolation:\n");
	usage(usgmsg1);
	printf("Give Interpolation Module Option: ");
	scanf("%d", &select);
	
	switch (select) {
		case 1:
			/* Linear Interpolation */
			LinearNode2CC_2D(In, Out);
			break;
		default:
			/* Loop */
			MSG("Node2CellCenter_2D: Invalid Interpolation Option");
			return 1;
	}

	return 0;
}
