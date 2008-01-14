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
 * File		VectorTopology_2D.h
 * Author	Ashish Gupta
 * Date		26/06/2006
 * Version	0.1
*/

#ifndef __VectorTopology_2D_H__
#define __VectorTopolopy_2D_H__

#include "Vector.h"

#define Saddle			1
#define Source			2
#define Sink			3
#define Center			4
#define RepellingSpiral		5
#define AttractingSpiral	6

typedef struct PP_2D_CriticalPoint_Un {
	int Size;
	int *Type;
	int *Cell;
	Point2 *CriticalPoint;
	Point2 *EigenVector1;
	Point2 *EigenVector2;
} CriticalPoint_2D_Un;

/* Vector Components Arg1 = X, Arg2 = Y */
extern int VectorTopologyArg1, VectorTopologyArg2;

extern CriticalPoint_2D_Un CPoint2D;

/* Topology Vector Data if CellCenter */
extern Data_2D_Un TopVecX;
extern Data_2D_Un TopVecY;

void SetTopologicalSeedPoints_2D (void);
int AddCriticalPoint_2D (int, Point2, Point2, Point2, int);
int GetCriticalPoint_2D (int, Point2 *, Point2 *, Point2 *, int *);
void VectorTopology_2D (void);

#endif
