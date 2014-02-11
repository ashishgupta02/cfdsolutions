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
 * File		PostProcessingInitialize_2D.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
 */



#ifndef __PostProcessingInitialize_2D_H__
#define __PostProcessingInitialize_2D_H__

#include "CGNSIO.h"
#include "PostProcessing_2D.h"

int compress_sort(int *, int, int);
int sort_bubble(int *, int);
int GetCellAdjacentCellsQuad_2D(int);
int GetCellAdjacentCellsTri_2D(int);
int GetCellAdjacentCells_2D(void);
int GetNodeAdjacentCells_2D(int *, int *, int);
int CountTotalNonZeroEntries_2D(void);
int CreateSymmetricCSRFormat_2D(void);
int FV_InsertSort_2D(int *, int);
int ResizeObject_2D(Vector_2D *);
int Push_Back_2D(Vector_2D *, int);
int CreateNodeAdjacentCells_2D(void);
int GetEdgeNumber_2D(int, int);
int CreateQuadEdgeList_2D(int);
int CreateTriangleEdgeList_2D(int);
int CreateCellEdgeList_2D(void);
int GetQuadrilateralEdges_2D(int, int *);
int GetTriangleEdges_2D(int, int *);
int GetCellEdgeNodes_2D(int, int, int *);
int CreateEdgeCellList_2D(void);
int SearchItem_2D(int, int);
int SearchInsert_2D(int, int);
int CreateNodeAdjacentNodes_2D(void);
int CreatEdgeStructure_2D(void);
int InitializeUnstructuredGrid_2D(ZONE *);
int InitializeStructuredGrid_2D(ZONE *);
int InitializeGrid_2D(ZONE *);

#endif
