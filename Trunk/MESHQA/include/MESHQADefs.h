/*******************************************************************************
 * File:        MESHQADefs.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _MESHQADEFS_H
#define	_MESHQADEFS_H

#include "MESHQACommon.h"

#define MESHQA_DBL_MIN 1.0E-30
#define MESHQA_DBL_MAX 1.0E+30

//! Define the number of all cell quality types (corresponding to MESHQAEnums.h).
#define MESHQA_NQUALITYTYPES 31

//! Define the number of Quality Parameter for Edge.
#define MESHQA_NQUALITY_EDGE2 1

//! Define the number of Quality Parameter for Triangle.
#define MESHQA_NQUALITY_TRI3 13

//! Define the number of Quality Parameter for Quadralateral.
#define MESHQA_NQUALITY_QUAD4 23

//! Define the number of Quality Parameter for Tetra.
#define MESHQA_NQUALITY_TETRA4 16

//! Define the number of Quality Parameter for Pyramid.
#define MESHQA_NQUALITY_PYRA5 1

//! Define the number of Quality Parameter for Prism.
#define MESHQA_NQUALITY_PRISM6 1

//! Define the number of Quality Parameter for Hexa.
#define MESHQA_NQUALITY_HEXA8 20

//! Define the MIXED Element connectivity Type
#define MESHQA_MIXED 20

//! Define the number of all cell types (corresponding to MESHQAEnums.h).
#define MESHQA_NCELLTYPES 8

//! Define the maximum number of a cell type (corresponding to MESHQAEnums.h).
#define MESHQA_CELLTYPE_MAX 8

//! Define the number of unstructured cell types (corresponding to MESHQAEnums.h).
#define MESHQA_NUNSTRUCTCELLTYPES 8

//! Define the number of structured cell types (corresponding to MESHQAEnums.h).
#define MESHQA_NSTRUCTCELLTYPES 5

//! Define the maximum number of cell nodes.
#define MESHQA_NMAXCELLNODES 8

//! Define the maximum number of cell edges.
#define MESHQA_NMAXCELLEDGES 12

//! Define the maximum number of cell faces.
#define MESHQA_NMAXCELLFACES 6

//! Define the maximum number of nodes on an edge of a cell.
#define MESHQA_NMAXEDGENODES 2

//! Define the maximum number of nodes on a face of a cell.
#define MESHQA_NMAXFACENODES 4

//! Define the number of nodes forming an edge.
#define MESHQA_NNODES_EDGE2 2

//! Define the number of nodes forming a triangle.
#define MESHQA_NNODES_TRI3 3

//! Define the number of nodes forming a quadrilateral.
#define MESHQA_NNODES_QUAD4 4

//! Define the number of nodes forming a tetrahedron.
#define MESHQA_NNODES_TETRA4 4

//! Define the number of nodes forming a pyramid.
#define MESHQA_NNODES_PYRA5 5

//! Define the number of nodes forming a prism.
#define MESHQA_NNODES_PRISM6 6

//! Define the number of nodes forming a hexahedron.
#define MESHQA_NNODES_HEXA8 8

//! Define the number of edges forming a triangle.
#define MESHQA_NEDGES_TRI3 3

//! Define the number of edges forming a quadrilateral.
#define MESHQA_NEDGES_QUAD4 4

//! Define the number of edges forming a tetrahedron.
#define MESHQA_NEDGES_TETRA4 6

//! Define the number of edges forming a pyramid.
#define MESHQA_NEDGES_PYRA5 8

//! Define the number of edges forming a prism.
#define MESHQA_NEDGES_PRISM6 9

//! Define the number of edges forming a hexahedron.
#define MESHQA_NEDGES_HEXA8 12

//! Define the number of faces forming a Triangle.
#define MESHQA_NFACES_TRI3 1

//! Define the number of faces forming a quadrilateral.
#define MESHQA_NFACES_QUAD4 1

//! Define the number of faces forming a tetrahedron.
#define MESHQA_NFACES_TETRA4 4

//! Define the number of faces forming a pyramid.
#define MESHQA_NFACES_PYRA5 5

//! Define the number of faces forming a prism.
#define MESHQA_NFACES_PRISM6 5

//! Define the number of faces forming a hexahedron.
#define MESHQA_NFACES_HEXA8 6

#endif	/* _MESHQADEFS_H */

