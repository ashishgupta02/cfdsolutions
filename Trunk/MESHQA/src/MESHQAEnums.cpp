/*******************************************************************************
 * File:        MESHQAEnums.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "MESHQAEnums.h"

//-----------------------------------------------------------------------------
//
//  Class variables
//

const MESHQAEnums::CellType MESHQAEnums::cAllCellTypes[MESHQA_NCELLTYPES] = {
    CT_Node,
    CT_Edge2,
    CT_Tri3,
    CT_Quad4,
    CT_Tetra4,
    CT_Pyra5,
    CT_Prism6,
    CT_Hexa8
};


const MESHQAEnums::CellType MESHQAEnums::cUnstructCellTypes[MESHQA_NUNSTRUCTCELLTYPES] = {
    CT_Node,
    CT_Edge2,
    CT_Tri3,
    CT_Quad4,
    CT_Tetra4,
    CT_Pyra5,
    CT_Prism6,
    CT_Hexa8
};

const MESHQAEnums::QualityType MESHQAEnums::cAllQualityTypes[MESHQA_NQUALITYTYPES] = {
    QT_Length,
    QT_Aspect,
    QT_AspectGamma,
    QT_Area,
    QT_SmallestAngle,
    QT_LargestAngle,
    QT_Condition,
    QT_Jacobian,
    QT_NormalizedJacobian,
    QT_Shear,
    QT_Shape,
    QT_RelativeSize,
    QT_ShapeSize,
    QT_Skew,
    QT_Taper,
    QT_Warpage,
    QT_Stretch,
    QT_Oddy,
    QT_Volume,
    QT_Diagonal,
    QT_Dimension,
    QT_EdgeRatio,
    QT_MaxEdgeRatio,
    QT_AspectFrobenius,
    QT_MedAspectFrobenius,
    QT_MaxAspectFrobenius,
    QT_Distortion,
    QT_ShearSize,
    QT_RadiusRatio,
    QT_AspectBeta,
    QT_CollapseRatio
};

//-----------------------------------------------------------------------------
//
//  IsUnstructCellType
//

bool MESHQAEnums::IsUnstructCellType(CellType t) {
    for (int i = 0; i < MESHQA_NUNSTRUCTCELLTYPES; i++) {
        if (t == cUnstructCellTypes[i]) {
            return true;
        }
    }
    return false;
}


//-----------------------------------------------------------------------------
//
//  Int2CellType
//

MESHQAEnums::CellType MESHQAEnums::Int2CellType(int t) {
    for (int i = 0; i < MESHQA_NCELLTYPES; i++) {
        int ct = cAllCellTypes[i];
        if (ct == t) {
            return cAllCellTypes[i];
        }
    }

    return CT_Undefined;
}

//-----------------------------------------------------------------------------
//
//  CellTypeToString
//

const char* MESHQAEnums::CellTypeToString(CellType t) {
    switch (t) {
        case CT_Undefined: return "Undefined";
        case CT_Node: return "Node";
        case CT_Edge2: return "Edge2";
        case CT_Tri3: return "Tri3";
        case CT_Quad4: return "Quad4";
        case CT_Tetra4: return "Tetra4";
        case CT_Pyra5: return "Pyra5";
        case CT_Prism6: return "Prism6";
        case CT_Hexa8: return "Hexa8";
        default:
            break;
    }

    return "unknown cell type";
}

//-----------------------------------------------------------------------------
//
//  GetQualityIndex
//

int MESHQAEnums::GetQualityIndex(CellType ct, QualityType qt) {
    if (ct == CT_Undefined)
        return -1;
    if (qt == QT_Undefined)
        return -1;

    int qindex = -1;

    switch (ct) {
        case CT_Node:
            qindex = -1;
            break;
        case CT_Edge2:
            switch (qt) {
                case QT_Length:
                    qindex = 0;
                    break;
                default:
                    qindex = -1;
            }
            break;
        case CT_Tri3:
            switch (qt) {
                case QT_Aspect:
                    qindex = 0;
                    break;
                case QT_Area:
                    qindex = 1;
                    break;
                case QT_SmallestAngle:
                    qindex = 2;
                    break;
                case QT_LargestAngle:
                    qindex = 3;
                    break;
                case QT_Condition:
                    qindex = 4;
                    break;
                case QT_NormalizedJacobian:
                    qindex = 5;
                    break;
                case QT_Shape:
                    qindex = 6;
                    break;
                case QT_RelativeSize:
                    qindex = 7;
                    break;
                case QT_ShapeSize:
                    qindex = 8;
                    break;
                case QT_EdgeRatio:
                    qindex = 9;
                    break;
                case QT_AspectFrobenius:
                    qindex = 10;
                    break;
                case QT_Distortion:
                    qindex = 11;
                    break;
                case QT_RadiusRatio:
                    qindex = 12;
                    break;
                default:
                    qindex = -1;
            }
            break;
        case CT_Quad4:
            switch (qt) {
                case QT_Aspect:
                    qindex = 0;
                    break;
                case QT_Area:
                    qindex = 1;
                    break;
                case QT_SmallestAngle:
                    qindex = 2;
                    break;
                case QT_LargestAngle:
                    qindex = 3;
                    break;
                case QT_Condition:
                    qindex = 4;
                    break;
                case QT_Jacobian:
                    qindex = 5;
                    break;
                case QT_NormalizedJacobian:
                    qindex = 6;
                    break;
                case QT_Shear:
                    qindex = 7;
                    break;
                case QT_Shape:
                    qindex = 8;
                    break;
                case QT_RelativeSize:
                    qindex = 9;
                    break;
                case QT_ShapeSize:
                    qindex = 10;
                    break;
                case QT_Skew:
                    qindex = 11;
                    break;
                case QT_Taper:
                    qindex = 12;
                    break;
                case QT_Warpage:
                    qindex = 13;
                    break;
                case QT_Stretch:
                    qindex = 14;
                    break;
                case QT_Oddy:
                    qindex = 15;
                    break;
                case QT_EdgeRatio:
                    qindex = 16;
                    break;
                case QT_MaxEdgeRatio:
                    qindex = 17;
                    break;
                case QT_MedAspectFrobenius:
                    qindex = 18;
                    break;
                case QT_MaxAspectFrobenius:
                    qindex = 19;
                    break;
                case QT_Distortion:
                    qindex = 20;
                    break;
                case QT_ShearSize:
                    qindex = 21;
                    break;
                case QT_RadiusRatio:
                    qindex = 22;
                    break;
                default:
                    qindex = -1;
            }
            break;
        case CT_Tetra4:
            switch (qt) {
                case QT_Aspect:
                    qindex = 0;
                    break;
                case QT_AspectGamma:
                    qindex = 1;
                    break;
                case QT_SmallestAngle:
                    qindex = 2;
                    break;
                case QT_Condition:
                    qindex = 3;
                    break;
                case QT_Jacobian:
                    qindex = 4;
                    break;
                case QT_NormalizedJacobian:
                    qindex = 5;
                    break;
                case QT_Shape:
                    qindex = 6;
                    break;
                case QT_RelativeSize:
                    qindex = 7;
                    break;
                case QT_ShapeSize:
                    qindex = 8;
                    break;
                case QT_Volume:
                    qindex = 9;
                    break;
                case QT_EdgeRatio:
                    qindex = 10;
                    break;
                case QT_AspectFrobenius:
                    qindex = 11;
                    break;
                case QT_Distortion:
                    qindex = 12;
                    break;
                case QT_RadiusRatio:
                    qindex = 13;
                    break;
                case QT_AspectBeta:
                    qindex = 14;
                    break;
                case QT_CollapseRatio:
                    qindex = 15;
                    break;
                default:
                    qindex = -1;
            }
            break;
        case CT_Pyra5:
            switch (qt) {
                case QT_Volume:
                    qindex = 0;
                    break;
                default:
                    qindex = -1;
            }
            break;
        case CT_Prism6:
            switch (qt) {
                case QT_Volume:
                    qindex = 0;
                    break;
                default:
                    qindex = -1;
            }
            break;
        case CT_Hexa8:
            switch (qt) {
                case QT_Condition:
                    qindex = 0;
                    break;
                case QT_Jacobian:
                    qindex = 1;
                    break;
                case QT_NormalizedJacobian:
                    qindex = 2;
                    break;
                case QT_Shear:
                    qindex = 3;
                    break;
                case QT_Shape:
                    qindex = 4;
                    break;
                case QT_RelativeSize:
                    qindex = 5;
                    break;
                case QT_ShapeSize:
                    qindex = 6;
                    break;
                case QT_Skew:
                    qindex = 7;
                    break;
                case QT_Taper:
                    qindex = 8;
                    break;
                case QT_Stretch:
                    qindex = 9;
                    break;
                case QT_Oddy:
                    qindex = 10;
                    break;
                case QT_Volume:
                    qindex = 11;
                    break;
                case QT_Diagonal:
                    qindex = 12;
                    break;
                case QT_Dimension:
                    qindex = 13;
                    break;
                case QT_EdgeRatio:
                    qindex = 14;
                    break;
                case QT_MaxEdgeRatio:
                    qindex = 15;
                    break;
                case QT_MedAspectFrobenius:
                    qindex = 16;
                    break;
                case QT_MaxAspectFrobenius:
                    qindex = 17;
                    break;
                case QT_Distortion:
                    qindex = 18;
                    break;
                case QT_ShearSize:
                    qindex = 19;
                    break;
                default:
                    qindex = -1;
            }
            break;
        case CT_Undefined:
            qindex = -1;
            break;
    }

    return qindex;
}

//-----------------------------------------------------------------------------
//
//  Int2QualityType
//
//! Returns the Quality Type from Integer according to the Cell Type

MESHQAEnums::QualityType MESHQAEnums::Int2QualityType(CellType ct, int t) {
    for (int i = 0; i < MESHQA_NQUALITYTYPES; i++) {
        int qt = GetQualityIndex(ct, cAllQualityTypes[i]);
        if (qt == t) {
            return cAllQualityTypes[i];
        }
    }
    return QT_Undefined;
}

//-----------------------------------------------------------------------------
//
//  GetQualityIndex
//

int MESHQAEnums::GetNumberQualityType(CellType ct) {
    int nquality = 0;

    switch (ct) {
        case CT_Undefined:
            nquality = 0;
            break;
        case CT_Node:
            nquality = 0;
            break;
        case CT_Edge2:
            nquality = MESHQA_NQUALITY_EDGE2;
            break;
        case CT_Tri3:
            nquality = MESHQA_NQUALITY_TRI3;
            break;
        case CT_Quad4:
            nquality = MESHQA_NQUALITY_QUAD4;
            break;
        case CT_Tetra4:
            nquality = MESHQA_NQUALITY_TETRA4;
            break;
        case CT_Pyra5:
            nquality = MESHQA_NQUALITY_PYRA5;
            break;
        case CT_Prism6:
            nquality = MESHQA_NQUALITY_PRISM6;
            break;
        case CT_Hexa8:
            nquality = MESHQA_NQUALITY_HEXA8;
            break;
    }

    return nquality;
}

//-----------------------------------------------------------------------------
//
//  QualityTypeToString
//

const char* MESHQAEnums::QualityTypeToString(QualityType qt) {

    switch (qt) {
        case QT_Undefined:          return "Undefined";
        case QT_Length:             return "Length";
        case QT_Aspect:             return "Aspect";
        case QT_AspectGamma:        return "AspectGamma";
        case QT_Area:               return "Area";
        case QT_SmallestAngle:      return "SmallestAngle";
        case QT_LargestAngle:       return "LargestAngle";
        case QT_Condition:          return "Condition";
        case QT_Jacobian:           return "Jacobian";
        case QT_NormalizedJacobian: return "NormalizedJacobian";
        case QT_Shear:              return "Shear";
        case QT_Shape:              return "Shape";
        case QT_RelativeSize:       return "RelativeSize";
        case QT_ShapeSize:          return "ShapeSize";
        case QT_Skew:               return "Skew";
        case QT_Taper:              return "Taper";
        case QT_Warpage:            return "Warpage";
        case QT_Stretch:            return "Stretch";
        case QT_Oddy:               return "Oddy";
        case QT_Volume:             return "Volume";
        case QT_Diagonal:           return "Diagonal";
        case QT_Dimension:          return "Dimension";
        case QT_EdgeRatio:          return "EdgeRatio";
        case QT_MaxEdgeRatio:       return "MaxEdgeRatio";
        case QT_AspectFrobenius:    return "AspectFrobenius";
        case QT_MedAspectFrobenius: return "MedAspectFrobenius";
        case QT_MaxAspectFrobenius: return "MaxAspectFrobenius";
        case QT_Distortion:         return "Distortion";
        case QT_ShearSize:          return "ShearSize";
        case QT_RadiusRatio:        return "RadiusRatio";
        case QT_AspectBeta:         return "AspectBeta";
        case QT_CollapseRatio:      return "CollapseRatio";
        default:
            break;
    }
    
    return "unknown quality type";
}

