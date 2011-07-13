/*******************************************************************************
 * File:        MESHQACellInfo.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "MESHQACellInfo.h"

//-----------------------------------------------------------------------------
//
//  Class Variables
//

const int MESHQACellInfo::cNNodes[ MESHQA_NCELLTYPES+1 ] = {
  0, // CT_Undefined
  1, // CT_Node
  2, // CT_Edge2
  3, // CT_Tri3
  4, // CT_Quad4
  4, // CT_Tetra4
  5, // CT_Pyra5
  6, // CT_Prism6
  8  // CT_Hexa8
};


const int MESHQACellInfo::cNEdges[ MESHQA_NCELLTYPES+1 ] = {
  0, // CT_Undefined
  0, // CT_Node
  0, // CT_Edge2
  3, // CT_Tri3
  4, // CT_Quad4
  6, // CT_Tetra4
  8, // CT_Pyra5
  9, // CT_Prism6
  12 // CT_Hexa8
};


const int MESHQACellInfo::cNFaces[ MESHQA_NCELLTYPES+1 ] = {
  0, // CT_Undefined
  0, // CT_Node
  0, // CT_Edge2
  0, // CT_Tri3
  0, // CT_Quad4
  4, // CT_Tetra4
  5, // CT_Pyra5
  5, // CT_Prism6
  6  // CT_Hexa8
};


const int MESHQACellInfo::cEdgeNodes[ MESHQA_NCELLTYPES+1 ][ MESHQA_NMAXCELLEDGES ][ MESHQA_NMAXEDGENODES ] = {
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}}, // CT_Undefined
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}}, // CT_Node
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}}, // CT_Edge2
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{ 1, 2},{ 2, 0},{ 0, 1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}}, // CT_Tri3
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{ 0, 1},{ 1, 2},{ 2, 3},{ 3, 0},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}}, // CT_Quad4
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{ 1, 2},{ 2, 0},{ 0, 1},{ 0, 3},{ 1, 3},{ 2, 3},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}}, // CT_Tetra4
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{ 0, 1},{ 1, 2},{ 2, 3},{ 3, 0},{ 0, 4},{ 1, 4},{ 2, 4},{ 3, 4},{-1,-1},{-1,-1},{-1,-1},{-1,-1}}, // CT_Pyra5
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{ 1, 2},{ 2, 0},{ 0, 1},{ 4, 5},{ 5, 3},{ 3, 4},{ 0, 3},{ 1, 4},{ 2, 5},{-1,-1},{-1,-1},{-1,-1}}, // CT_Prism6
  //  0       1       2       3       4       5       6       7       8       9      10      11
  {{ 0, 1},{ 1, 2},{ 2, 3},{ 3, 0},{ 4, 5},{ 5, 6},{ 6, 7},{ 7, 4},{ 0, 4},{ 1, 5},{ 2, 6},{ 3, 7}}  // CT_Hexa8
};


const int MESHQACellInfo::cFaceNodes[ MESHQA_NCELLTYPES+1 ][ MESHQA_NMAXCELLFACES ][ MESHQA_NMAXFACENODES ] = {
  //     0             1             2             3             4             5
  {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}, // CT_Undefined
  //     0             1             2             3             4             5
  {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}, // CT_Node
  //     0             1             2             3             4             5
  {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}, // CT_Edge2
  //     0             1             2             3             4             5
  {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}, // CT_Tri3
  //     0             1             2             3             4             5
  {{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}, // CT_Quad4
  //     0             1             2             3             4             5
  {{ 1, 0, 2,-1},{ 1, 3, 0,-1},{ 2, 3, 1,-1},{ 0, 3, 2,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}, // CT_Tetra4
  //     0             1             2             3             4             5
  {{ 1, 0, 3, 2},{ 1, 4, 0,-1},{ 2, 4, 1,-1},{ 3, 4, 2,-1},{ 0, 4, 3,-1},{-1,-1,-1,-1}}, // CT_Pyra5
  //     0             1             2             3             4             5
  {{ 1, 0, 2,-1},{ 3, 4, 5,-1},{ 1, 4, 3, 0},{ 1, 2, 5, 4},{ 2, 0, 3, 5},{-1,-1,-1,-1}}, // CT_Prism6
  //     0             1             2             3             4             5
  {{ 1, 0, 3, 2},{ 4, 5, 6, 7},{ 1, 5, 4, 0},{ 2, 6, 5, 1},{ 3, 7, 6, 2},{ 0, 4, 7, 3}}  // CT_Hexa8
};

const MESHQAEnums::CellType MESHQACellInfo::cFaceTypes[ MESHQA_NCELLTYPES+1 ][ MESHQA_NMAXCELLFACES ] = {
  {MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined}, // CT_Undefined
  {MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined}, // CT_Node
  {MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined}, // CT_Edge2
  {MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined}, // CT_Tri3
  {MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined}, // CT_Quad4
  {MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Undefined, MESHQAEnums::CT_Undefined}, // CT_Tetra4
  {MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Undefined}, // CT_Pyra5
  {MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Tri3,      MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Undefined}, // CT_Prism6
  {MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Quad4,     MESHQAEnums::CT_Quad4}      // CT_Hexa8
};


//-----------------------------------------------------------------------------
//
//  GetNNodes
//

int MESHQACellInfo::GetNNodes(MESHQAEnums::CellType t)
{
  return NNodes(t);
}


//-----------------------------------------------------------------------------
//
//  GetNEdges
//

int MESHQACellInfo::GetNEdges(MESHQAEnums::CellType t)
{
  return NEdges(t);
}


//-----------------------------------------------------------------------------
//
//  GetNFaces
//

int MESHQACellInfo::GetNFaces(MESHQAEnums::CellType t)
{
  return NFaces(t);
}


//-----------------------------------------------------------------------------
//
//  GetEdgeNode
//

int MESHQACellInfo::GetEdgeNode(MESHQAEnums::CellType t, int edge, int node)
{
  return EdgeNode(t, edge, node);
}


//-----------------------------------------------------------------------------
//
//  GetFaceNode
//

int MESHQACellInfo::GetFaceNode(MESHQAEnums::CellType t, int face, int node)
{
  return FaceNode(t, face, node);
}


//-----------------------------------------------------------------------------
//
//  GetFaceNode
//

MESHQAEnums::CellType MESHQACellInfo::GetFaceType(MESHQAEnums::CellType t, int face)
{
  return FaceType(t, face);
}

