/*******************************************************************************
 * File:        MeshIO.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _MESHIO_H
#define	_MESHIO_H

void UGrid_Reader(const char* filename);
void VTK_Writer(const char* filename, int verbose);

#endif	/* _MESHIO_H */

