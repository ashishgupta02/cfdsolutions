/*******************************************************************************
 * File:        MeshIO.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef MESHIO_H
#define	MESHIO_H

void UGrid_Reader(const char* filename);
void VTK_Writer(const char* filename, int verbose);

#endif	/* MESHIO_H */

