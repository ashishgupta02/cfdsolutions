/*******************************************************************************
 * File:        MESHQACellInfo.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _MESHQACELLINFO_H
#define	_MESHQACELLINFO_H

#include "MESHQACommon.h"
#include "MESHQADefs.h"
#include "MESHQAEnums.h"

#ifdef MESHQA_SAFETYCHECKS
#define MESHQACELLINFO_SAFETYCHECKS
#endif

//-----------------------------------------------------------------------------
//
//  MESHQACellInfo
//

//! This class contains information about the different cell types.

/*!
  It contains the number of nodes, edges, faces, etc. for each
  cell type and the definition of edges and faces forming a cell.
 */

class MESHQACellInfo {
public:

    // --- Class Methods ---

    // ---------- INLINE METHODS ----------

    //! Return the number of nodes forming a cell.
    /*!
      \param t the cell type to be queried.
      \return The number of nodes forming a cell.
     */
    static inline int NNodes(MESHQAEnums::CellType t);

    //! Return the number of edges forming a cell.
    /*!
      \param t the cell type to be queried.
      \return The number of edges forming a cell.
     */
    static inline int NEdges(MESHQAEnums::CellType t);

    //! Return the number of faces forming a cell.
    /*!
      \param t the cell type to be queried.
      \return The number of faces forming a cell.
     */
    static inline int NFaces(MESHQAEnums::CellType t);

    //! Return the node of an edge of a cell for the given cell type.
    /*!
      \param t the cell type to be queried.
      \param edge the cell edge to be queried.
      \param node the node of the edge to be queried.
      \return The node of an edge of a cell for the given cell type.
     */
    static inline int EdgeNode(MESHQAEnums::CellType t, int edge, int node);

    //! Return the node of a face of a cell for the given cell type.
    /*!
      \param t the cell type to be queried.
      \param face the cell face to be queried.
      \param node the node of the face to be queried.
      \return The node of a face of a cell for the given cell type.
     */
    static inline int FaceNode(MESHQAEnums::CellType t, int face, int node);

    //! Return the type of a face of a cell for the given cell type.
    /*!
      \param t the cell type to be queried.
      \param face the cell face to be queried.
      \return The type of a face of a cell for the given cell type.
     */
    static inline MESHQAEnums::CellType FaceType(MESHQAEnums::CellType t, int face);


    // ---------- METHODS ----------

    //! Return the number of nodes forming a cell.
    /*!
      \param t the cell type to be queried.
      \return The number of nodes forming a cell.
     */
    static int GetNNodes(MESHQAEnums::CellType t);

    //! Return the number of edges forming a cell.
    /*!
      \param t the cell type to be queried.
      \return The number of edges forming a cell.
     */
    static int GetNEdges(MESHQAEnums::CellType t);

    //! Return the number of faces forming a cell.
    /*!
      \param t the cell type to be queried.
      \return The number of faces forming a cell.
     */
    static int GetNFaces(MESHQAEnums::CellType t);

    //! Return the node of an edge of a cell for the given cell type.
    /*!
      \param t the cell type to be queried.
      \param edge the cell edge to be queried.
      \param node the node of the edge to be queried.
      \return The node of an edge of a cell for the given cell type.
     */
    static int GetEdgeNode(MESHQAEnums::CellType t, int edge, int node);

    //! Return the node of a face of a cell for the given cell type.
    /*!
      \param t the cell type to be queried.
      \param face the cell face to be queried.
      \param node the node of the face to be queried.
      \return The node of a face of a cell for the given cell type.
     */
    static int GetFaceNode(MESHQAEnums::CellType t, int face, int node);

    //! Return the type of a face of a cell for the given cell type.
    /*!
      \param t the cell type to be queried.
      \param face the cell face to be queried.
      \return The type of a face of a cell for the given cell type.
     */
    static MESHQAEnums::CellType GetFaceType(MESHQAEnums::CellType t, int face);


    // --- Class Variables ---

    //! Array containing the number of nodes forming a cell for each cell type.
    static const int cNNodes[ MESHQA_NCELLTYPES + 1 ];

    //! Array containing the number of edges forming a cell for each cell type.
    static const int cNEdges[ MESHQA_NCELLTYPES + 1 ];

    //! Array containing the number of faces forming a cell for each cell type.
    static const int cNFaces[ MESHQA_NCELLTYPES + 1 ];

    //! Array containing the nodes of the edges forming a cell for each cell type.
    static const int cEdgeNodes[ MESHQA_NCELLTYPES + 1 ][ MESHQA_NMAXCELLEDGES ][ MESHQA_NMAXEDGENODES ];

    //! Array containing the nodes of the faces forming a cell for each cell type.
    static const int cFaceNodes[ MESHQA_NCELLTYPES + 1 ][ MESHQA_NMAXCELLFACES ][ MESHQA_NMAXFACENODES ];

    //! Array containing the cell type of the faces forming a cell for each cell type.
    static const MESHQAEnums::CellType cFaceTypes[ MESHQA_NCELLTYPES + 1 ][ MESHQA_NMAXCELLFACES ];


protected:

    // --- Methods ---

    //! Ensure the cell type to be valid, otherwise stop.
    /*!
      param t the cell type to be checked.
     */
    static inline void CheckCellType(MESHQAEnums::CellType t);
};


//-----------------------------------------------------------------------------
//
//  Abbreviations
//

#define MESHQANNodes(t) (MESHQACellInfo::NNodes(t))
#define MESHQANEdges(t) (MESHQACellInfo::NEdges(t))
#define MESHQANFaces(t) (MESHQACellInfo::NFaces(t))


//-----------------------------------------------------------------------------
//
//  NNodes
//

inline int MESHQACellInfo::NNodes(MESHQAEnums::CellType t) {
#ifdef MESHQACELLINFO_SAFETYCHECKS
    CheckCellType(t);
#endif

    return cNNodes[t];
}


//-----------------------------------------------------------------------------
//
//  NEdges
//

inline int MESHQACellInfo::NEdges(MESHQAEnums::CellType t) {
#ifdef MESHQACELLINFO_SAFETYCHECKS
    CheckCellType(t);
#endif

    return cNEdges[t];
}


//-----------------------------------------------------------------------------
//
//  NFaces
//

inline int MESHQACellInfo::NFaces(MESHQAEnums::CellType t) {
#ifdef MESHQACELLINFO_SAFETYCHECKS
    CheckCellType(t);
#endif

    return cNFaces[t];
}


//-----------------------------------------------------------------------------
//
//  EdgeNode
//

inline int MESHQACellInfo::EdgeNode(MESHQAEnums::CellType t, int edge, int node) {
#ifdef MESHQACELLINFO_SAFETYCHECKS
    // --- check cell type ---
    CheckCellType(t);
    // --- check edge number ---
    if (edge < 0 || edge >= MESHQA_NMAXCELLEDGES) {
        MESHQAErrorMacro("MESHQACellInfo::EdgeNode() illegal edge number!");
    }
    // --- check edge node number ---
    if (node < 0 || node >= MESHQA_NMAXEDGENODES) {
        MESHQAErrorMacro("MESHQACellInfo::EdgeNode() illegal edge node number!");
    }
#endif

    return cEdgeNodes[t][edge][node];
}


//-----------------------------------------------------------------------------
//
//  FaceNode
//

inline int MESHQACellInfo::FaceNode(MESHQAEnums::CellType t, int face, int node) {
#ifdef MESHQACELLINFO_SAFETYCHECKS
    // --- check cell type ---
    CheckCellType(t);
    // --- check face number ---
    if (face < 0 || face >= MESHQA_NMAXCELLFACES) {
        MESHQAErrorMacro("MESHQACellInfo::FaceNode() illegal face number!");
    }
    // --- check face node number ---
    if (node < 0 || node >= MESHQA_NMAXFACENODES) {
        MESHQAErrorMacro("MESHQACellInfo::FaceNode() illegal face node number!");
    }
#endif

    return cFaceNodes[t][face][node];
}


//-----------------------------------------------------------------------------
//
//  FaceType
//

inline MESHQAEnums::CellType MESHQACellInfo::FaceType(MESHQAEnums::CellType t, int face) {
#ifdef MESHQACELLINFO_SAFETYCHECKS
    // --- check cell type ---
    CheckCellType(t);
    // --- check face number ---
    if (face < 0 || face >= MESHQA_NMAXCELLFACES) {
        MESHQAErrorMacro("MESHQACellInfo::FaceType() illegal face number!");
    }
#endif

    return cFaceTypes[t][face];
}


//-----------------------------------------------------------------------------
//
//  CheckCellType
//

inline void MESHQACellInfo::CheckCellType(MESHQAEnums::CellType t) {
    if (t != MESHQAEnums::CT_Undefined &&
            t != MESHQAEnums::CT_Node &&
            t != MESHQAEnums::CT_Edge2 &&
            t != MESHQAEnums::CT_Tri3 &&
            t != MESHQAEnums::CT_Quad4 &&
            t != MESHQAEnums::CT_Tetra4 &&
            t != MESHQAEnums::CT_Pyra5 &&
            t != MESHQAEnums::CT_Prism6 &&
            t != MESHQAEnums::CT_Hexa8) {
        MESHQAErrorMacro("MESHQACellInfo::CheckCellType() illegal cell type!");
    }
}


#endif	/* _MESHQACELLINFO_H */

