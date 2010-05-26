/*******************************************************************************
 * File:        MESHQAEnums.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _MESHQAENUMS_H
#define	_MESHQAENUMS_H

#include "MESHQACommon.h"
#include "MESHQADefs.h"

//-----------------------------------------------------------------------------
//
//  MESHQAEnums
//

//! In this class common enumeration types are defined.

class MESHQAEnums {
public:

    //! The mesh cell types.

    /*!
      Explicit values are assigned for supporting the use of the cell type as index.
     */
    enum CellType {
        CT_Undefined = 0,
        // the unstructured cell types
        CT_Node = 1,
        CT_Edge2 = 2,
        CT_Tri3 = 3,
        CT_Quad4 = 4,
        CT_Tetra4 = 5,
        CT_Pyra5 = 6,
        CT_Prism6 = 7,
        CT_Hexa8 = 8
    };

    //! The Cell Quality types.
    /*!
      Explicit values are assigned for supporting the use of the quality type as index.
     */
    enum QualityType {
        QT_Undefined            = 0,
        // the cell Quality types
        QT_Length               = 1,
        QT_Aspect               = 2,
        QT_AspectGamma          = 3,
        QT_Area                 = 4,
        QT_SmallestAngle        = 5,
        QT_LargestAngle         = 6,
        QT_Condition            = 7,
        QT_Jacobian             = 8,
        QT_NormalizedJacobian   = 9,
        QT_Shear                = 10,
        QT_Shape                = 11,
        QT_RelativeSize         = 12,
        QT_ShapeSize            = 13,
        QT_Skew                 = 14,
        QT_Taper                = 15,
        QT_Warpage              = 16,
        QT_Stretch              = 17,
        QT_Oddy                 = 18,
        QT_Volume               = 19,
        QT_Diagonal             = 20,
        QT_Dimension            = 21,
        QT_EdgeRatio            = 22,
        QT_MaxEdgeRatio         = 23,
        QT_AspectFrobenius      = 24,
        QT_MedAspectFrobenius   = 25,
        QT_MaxAspectFrobenius   = 26,
        QT_Distortion           = 27,
        QT_ShearSize            = 28,
        QT_RadiusRatio          = 29,
        QT_AspectBeta           = 30,
        QT_CollapseRatio        = 31
    };

    //! Check whether the given cell type is an unstructured cell type.
    /*!
      \return True if the given cell type is an unstructured cell type, otherwise false.
     */
    static bool IsUnstructCellType(CellType t);
    
    //! All cell types in a list.
    static const CellType cAllCellTypes[MESHQA_NCELLTYPES];

    //! All Quality types in a list.
    static const QualityType cAllQualityTypes[MESHQA_NQUALITYTYPES];
    
    //! All unstructured cell types in a list.
    static const CellType cUnstructCellTypes[MESHQA_NUNSTRUCTCELLTYPES];
    
    //! Return the corresponding cell type.
    /*!
      \param t integer defining the cell type.
      \return The corresponding cell type.
     */
    static CellType Int2CellType(int t);

    //! Return the mesh cell type name.
    /*!
     \param t the cell type for which a name is returned.
     \return The pointer to a C-string containing the name of the cell type.
    */
    static const char* CellTypeToString(CellType t);

    //! Returns the Index of Quality according the Cell Type
    static int GetQualityIndex(CellType ct, QualityType qt);

    //! Returns the Quality Type from Integer according to the Cell Type
    static QualityType Int2QualityType(CellType ct, int t);

    //! Returns the Number of Quality Associated with Cell Type
    static int GetNumberQualityType(CellType ct);

    //! Returns the Quality Type Name
    static const char* QualityTypeToString(QualityType qt);
};

#endif	/* _MESHQAENUMS_H */

