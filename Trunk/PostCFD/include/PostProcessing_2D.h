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
 * File		PostProcessing_2D.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

#ifndef __POSTPROCESSING_2D_H__
#define __POSTPROCESSING_2D_H__

#define INTERNAL	0

#define EXTERNAL	1

/* Post-Processing Data Structure */
/* In array: X = *[0], Y = *[1], Z = *[2] */

/* 2D Post-Processing Data Structure */
typedef struct PP_2D_Node {
	/* Flag Denote Bounday or Internal Node */
	int	Flag;
	double	Coordinate[2];
	double	WeightTotal;
} Node_2D;

typedef struct PP_2D_Edge {
	int	Flag, NoOfNeighbourCells;
	double	MidCoordinate[2];
	int	ConnectedNodes[2];
	double	length;
	double	dx, dy, nx, ny;
	/* Cell[0] and Cell[1] are adjacent cells of edge. For boundary cell with
	NoOfNeighbouringCells = 1, Cell[0] = 0 */
	int	Cell[2];
	int	VisitFlag;
} Edge_2D;

typedef struct PP_2D_Cell {
	int	Flag, NodesPerCell, EdgesPerCell, NearCells;
	double	Area;
	double	Centroid[2];
	int	*ConnectNode, *ConnectEdge, *NearNeighbours;
} Cell_2D;

typedef struct PP_2D_Vector {
	int	Size;
	int	Capacity;
	int	*Data;
} Vector_2D;

/* Bigger Picture */
typedef struct PP_2D_Block {
	int	id;
	char	Name[33];
	Node_2D	*Node2D;
	Cell_2D *Cell2D;
	Edge_2D	*Edge2D;
	Vector_2D *NodeAdjacentNodes2D;
	Vector_2D *NodeAdjacentCells2D;
} Block_2D;

extern int	NoCells2D, NoEdges2D, NoNodes2D, NoOfBoundaryEdges2D, MaxEdgeNum2D;
extern Node_2D	*Node2D;
extern Cell_2D	*Cell2D;
extern Edge_2D	*Edge2D;
extern Vector_2D *NodeAdjacentNodes2D;
extern Vector_2D *NodeAdjacentCells2D;
extern int	*IA_2D;
extern int	*JA_2D;

/* Function Declearation */
extern int PostProcessing (const char *);

#endif
