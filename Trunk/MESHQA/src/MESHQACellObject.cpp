/*******************************************************************************
 * File:        MESHQACellObject.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "MESHQACellObject.h"
#include "MESHQAQualityMetrics.h"

//! Default Constructor
MESHQACellObject::MESHQACellObject()
{
    state       = 0;
    normalstate = 0;
    id          = 0;
    cellType    = MESHQAEnums::CT_Undefined;
}

//! Overloaded Constructor
//! Just sets the Cell Type at the time of Object creation
MESHQACellObject::MESHQACellObject(MESHQAEnums::CellType ctype)
{
    state       = 0;
    normalstate = 0;
    id          = 0;
    cellType    = ctype;
}

//! Copy Constructor
//! Provide Way to have Deep Copy
//MESHQACellObject::MESHQACellObject(const MESHQACellObject& orig)
//{
//    cout << "Not Implemented: MESHQACellObject::MESHQACellObject(const MESHQACellObject& orig)" << endl;
// TODO
//}

//! Distructor for MESHQACellObject Object
//! Calls reset to free all the used resources
MESHQACellObject::~MESHQACellObject()
{
    reset();
}

// Basic Quality Operations
//! This function ensures that all memory is released
void MESHQACellObject::reset()
{
    if (!state) {
        normalstate = 0;
        id          = 0;
        cellType    = MESHQAEnums::CT_Undefined;
        return;
    }

    // Free the Vectors
    ver.clear();
    quality.clear();
    normal.clear();
    id          = 0;
    cellType    = MESHQAEnums::CT_Undefined;
    normalstate = 0;
    state       = 0;
    
    return;
}

//! Adds the new vertex to Cell Object Vector
void MESHQACellObject::addVertex(VERTEX* pver)
{
    // Check if Vertex exists
    if (pver == NULL) return;

    // Add vertex at the end
    ver.push_back(pver);
    
    state = 1;
    
    return;
}

//! Just the cell type to the Cell Object
void MESHQACellObject::setCellType(MESHQAEnums::CellType ctype)
{
    cellType = ctype;
    return;
}

//! Just return the Cell Type of Cell Object
MESHQAEnums::CellType MESHQACellObject::getCellType() const
{
    return cellType;
}

//! Global Initialize
//! Does Not Initialize if Nodes are not set Properly
//! Initialize only those quality which are available
void MESHQACellObject::intialize()
{
    // Check the state
    if (!state) return;
    
    // Verify the Cell Type
    // Check if vectors are appropriatly created
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            // TODO Not Implemented Yet
            return;
            break;
        case MESHQAEnums::CT_Tri3:
            if (ver.size() != MESHQA_NNODES_TRI3) return;
            break;
        case MESHQAEnums::CT_Quad4:
            if (ver.size() != MESHQA_NNODES_QUAD4) return;
            break;
        case MESHQAEnums::CT_Tetra4:
            if (ver.size() != MESHQA_NNODES_TETRA4) return;
            break;
        case MESHQAEnums::CT_Pyra5:
            if (ver.size() != MESHQA_NNODES_PYRA5) return;
            break;
        case MESHQAEnums::CT_Prism6:
            if (ver.size() != MESHQA_NNODES_PRISM6) return;
            break;
        case MESHQAEnums::CT_Hexa8:
            if (ver.size() != MESHQA_NNODES_HEXA8) return;
            break;
        default:
            return;
    }

    // Create the iterator
    map<int, double>::iterator iter;
    iter = quality.begin();
    size_t size = quality.size();

    // Initialize the already added Quality
    for (int i = 0; i < (int)size; i++) {
        iter->second = 0.0;
        iter++;
    }
    
    return;
}

//! Intialize Perticular Quality Type
//! Does Not Initialize if Nodes are not set Properly
//! Initialize only those quality which are available
void MESHQACellObject::intialize_quality(MESHQAEnums::QualityType qualityType)
{
// Check the state
    if (!state) return;

    // Verify the Cell Type
    // Check if vectors are appropriatly created
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            // TODO Not Implemented Yet
            return;
            break;
        case MESHQAEnums::CT_Tri3:
            if (ver.size() != MESHQA_NNODES_TRI3) return;
            break;
        case MESHQAEnums::CT_Quad4:
            if (ver.size() != MESHQA_NNODES_QUAD4) return;
            break;
        case MESHQAEnums::CT_Tetra4:
            if (ver.size() != MESHQA_NNODES_TETRA4) return;
            break;
        case MESHQAEnums::CT_Pyra5:
            if (ver.size() != MESHQA_NNODES_PYRA5) return;
            break;
        case MESHQAEnums::CT_Prism6:
            if (ver.size() != MESHQA_NNODES_PRISM6) return;
            break;
        case MESHQAEnums::CT_Hexa8:
            if (ver.size() != MESHQA_NNODES_HEXA8) return;
            break;
        default:
            return;
    }

    // Create the iterator
    map<int, double>::iterator iter;
    
    // Get the Quality Key
    int qualId = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // Quality Key Not Valid
    if (qualId == -1) return;
    
    // Find the Quality Key in Quality Map
    // Else Add new Quality into the map
    iter = quality.find(qualId);

    if (iter != quality.end())
        iter->second = 0.0;
    
    return;
}

//! Initialize Normals
void MESHQACellObject::intialize_normal()
{
    // Check the state
    if (!state) return;

    // Verify the Cell Type
    // Check if vectors are appropriatly created
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            // TODO Not Implemented Yet
            return;
            break;
        case MESHQAEnums::CT_Tri3:
            if (ver.size() != MESHQA_NNODES_TRI3) return;
            if (normal.size() != MESHQA_NFACES_TRI3*3)
                normal.clear();
            else {
                for (int i = 0; i < MESHQA_NFACES_TRI3*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Quad4:
            if (ver.size() != MESHQA_NNODES_QUAD4) return;
            if (normal.size() != MESHQA_NFACES_QUAD4*3)
                normal.clear();
            else {
                for (int i = 0; i < MESHQA_NFACES_QUAD4*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Tetra4:
            if (ver.size() != MESHQA_NNODES_TETRA4) return;
            if (normal.size() != MESHQA_NFACES_TETRA4*3)
                normal.clear();
            else {
                for (int i = 0; i < MESHQA_NFACES_TETRA4*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Pyra5:
            if (ver.size() != MESHQA_NNODES_PYRA5) return;
            if (normal.size() != MESHQA_NFACES_PYRA5*3)
                normal.clear();
            else {
                for (int i = 0; i < MESHQA_NFACES_PYRA5*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Prism6:
            if (ver.size() != MESHQA_NNODES_PRISM6) return;
            if (normal.size() != MESHQA_NFACES_PRISM6*3)
                normal.clear();
            else {
                for (int i = 0; i < MESHQA_NFACES_PRISM6*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Hexa8:
            if (ver.size() != MESHQA_NNODES_HEXA8) return;
            if (normal.size() != MESHQA_NFACES_HEXA8*3)
                normal.clear();
            else {
                for (int i = 0; i < MESHQA_NFACES_HEXA8*3; i++)
                    normal[i] = 0.0;
            }
            break;
        default:
            return;
    }

    // Set Normal State
    normalstate = 0;
    return;
}

//! Global Analyze
void MESHQACellObject::analyze()
{
    // Check the state
    if (!state) return;
    
    size_t size = ver.size();
    static double coordinates[MESHQA_NNODES_HEXA8][3];
    
    // Verify the Cell Type
    // Check if vectors are appropriatly created
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            // TODO Not Implemented Yet
            return;
            break;
        case MESHQAEnums::CT_Tri3:
            if (size != MESHQA_NNODES_TRI3) return;
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            break;
        case MESHQAEnums::CT_Quad4:
            if (size != MESHQA_NNODES_QUAD4) return;
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            break;
        case MESHQAEnums::CT_Tetra4:
            if (size != MESHQA_NNODES_TETRA4) return;
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            break;
        case MESHQAEnums::CT_Pyra5:
            if (size != MESHQA_NNODES_PYRA5) return;
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            coordinates[4][0] = ver[4]->x;
            coordinates[4][1] = ver[4]->y;
            coordinates[4][2] = ver[4]->z;
            break;
        case MESHQAEnums::CT_Prism6:
            if (size != MESHQA_NNODES_PRISM6) return;
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            coordinates[4][0] = ver[4]->x;
            coordinates[4][1] = ver[4]->y;
            coordinates[4][2] = ver[4]->z;
            coordinates[5][0] = ver[5]->x;
            coordinates[5][1] = ver[5]->y;
            coordinates[5][2] = ver[5]->z;
            break;
        case MESHQAEnums::CT_Hexa8:
            if (size != MESHQA_NNODES_HEXA8) return;
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            coordinates[4][0] = ver[4]->x;
            coordinates[4][1] = ver[4]->y;
            coordinates[4][2] = ver[4]->z;
            coordinates[5][0] = ver[5]->x;
            coordinates[5][1] = ver[5]->y;
            coordinates[5][2] = ver[5]->z;
            coordinates[6][0] = ver[6]->x;
            coordinates[6][1] = ver[6]->y;
            coordinates[6][2] = ver[6]->z;
            coordinates[7][0] = ver[7]->x;
            coordinates[7][1] = ver[7]->y;
            coordinates[7][2] = ver[7]->z;
            break;
        default:
            return;
    }
    
    // Create the iterator
    map<int, double>::iterator iter;

    // Get the Quality Key
    int nqual = MESHQAEnums::GetNumberQualityType(cellType);
    
    double qvalue = 0.0;
    MESHQAEnums::QualityType qualityType;
    for (int i = 0; i < nqual; i++) {
        // Get the Quality Type from ID
        qualityType = MESHQAEnums::Int2QualityType(cellType, i);
        qvalue = 0.0;
        // Update the Coordinates and Compute the Quality
        switch (cellType) {
            case MESHQAEnums::CT_Tri3 :
                qvalue = MESHQAQualityMetrics::analyze_Tri3(qualityType, coordinates);
                break;
            case MESHQAEnums::CT_Quad4 :
                qvalue = MESHQAQualityMetrics::analyze_Quad4(qualityType, coordinates);
                break;
            case MESHQAEnums::CT_Tetra4 :
                qvalue = MESHQAQualityMetrics::analyze_Tetra4(qualityType, coordinates);
                break;
            case MESHQAEnums::CT_Pyra5 :
                qvalue = MESHQAQualityMetrics::analyze_Pyra5(qualityType, coordinates);
                break;
            case MESHQAEnums::CT_Prism6 :
                qvalue = MESHQAQualityMetrics::analyze_Prism6(qualityType, coordinates);
                break;
            case MESHQAEnums::CT_Hexa8 :
                qvalue = MESHQAQualityMetrics::analyze_Hexa8(qualityType, coordinates);
                break;
            default:
                return;
        }

        // Find the Quality Key in Quality Map
        // Else Add new Quality into the map
        iter = quality.find(i);

        if (iter != quality.end())
            iter->second = qvalue;
        else
            quality.insert(pair<int, double>(i, qvalue));
    }

    return;
}

//! Analyze Perticular Quality Type
void MESHQACellObject::analyze_quality(MESHQAEnums::QualityType qualityType)
{
    // Check the state
    if (!state) return;
    
    size_t size = ver.size();
    static double coordinates[MESHQA_NNODES_HEXA8][3];

    // Verify the Cell Type
    // Check if vectors are appropriatly created
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            // TODO Not Implemented Yet
            return;
            break;
        case MESHQAEnums::CT_Tri3:
            if (size != MESHQA_NNODES_TRI3) return;
            break;
        case MESHQAEnums::CT_Quad4:
            if (size != MESHQA_NNODES_QUAD4) return;
            break;
        case MESHQAEnums::CT_Tetra4:
            if (size != MESHQA_NNODES_TETRA4) return;
            break;
        case MESHQAEnums::CT_Pyra5:
            if (size != MESHQA_NNODES_PYRA5) return;
            break;
        case MESHQAEnums::CT_Prism6:
            if (size != MESHQA_NNODES_PRISM6) return;
            break;
        case MESHQAEnums::CT_Hexa8:
            if (size != MESHQA_NNODES_HEXA8) return;
            break;
        default:
            return;
    }

    // Create the iterator
    map<int, double>::iterator iter;

    // Get the Quality Key
    int qualId = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    
    // Quality Key Not Valid
    if (qualId == -1) return;

    double qvalue = 0.0;
    
    // Update the Coordinates and Compute the Quality
    switch (cellType) {
        case MESHQAEnums::CT_Tri3:
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            qvalue = MESHQAQualityMetrics::analyze_Tri3(qualityType, coordinates);
            break;
        case MESHQAEnums::CT_Quad4:
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            qvalue = MESHQAQualityMetrics::analyze_Quad4(qualityType, coordinates);
            break;
        case MESHQAEnums::CT_Tetra4:
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            qvalue = MESHQAQualityMetrics::analyze_Tetra4(qualityType, coordinates);
            break;
        case MESHQAEnums::CT_Pyra5:
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            coordinates[4][0] = ver[4]->x;
            coordinates[4][1] = ver[4]->y;
            coordinates[4][2] = ver[4]->z;
            qvalue = MESHQAQualityMetrics::analyze_Pyra5(qualityType, coordinates);
            break;
        case MESHQAEnums::CT_Prism6:
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            coordinates[4][0] = ver[4]->x;
            coordinates[4][1] = ver[4]->y;
            coordinates[4][2] = ver[4]->z;
            coordinates[5][0] = ver[5]->x;
            coordinates[5][1] = ver[5]->y;
            coordinates[5][2] = ver[5]->z;
            qvalue = MESHQAQualityMetrics::analyze_Prism6(qualityType, coordinates);
            break;
        case MESHQAEnums::CT_Hexa8:
            coordinates[0][0] = ver[0]->x;
            coordinates[0][1] = ver[0]->y;
            coordinates[0][2] = ver[0]->z;
            coordinates[1][0] = ver[1]->x;
            coordinates[1][1] = ver[1]->y;
            coordinates[1][2] = ver[1]->z;
            coordinates[2][0] = ver[2]->x;
            coordinates[2][1] = ver[2]->y;
            coordinates[2][2] = ver[2]->z;
            coordinates[3][0] = ver[3]->x;
            coordinates[3][1] = ver[3]->y;
            coordinates[3][2] = ver[3]->z;
            coordinates[4][0] = ver[4]->x;
            coordinates[4][1] = ver[4]->y;
            coordinates[4][2] = ver[4]->z;
            coordinates[5][0] = ver[5]->x;
            coordinates[5][1] = ver[5]->y;
            coordinates[5][2] = ver[5]->z;
            coordinates[6][0] = ver[6]->x;
            coordinates[6][1] = ver[6]->y;
            coordinates[6][2] = ver[6]->z;
            coordinates[7][0] = ver[7]->x;
            coordinates[7][1] = ver[7]->y;
            coordinates[7][2] = ver[7]->z;
            qvalue = MESHQAQualityMetrics::analyze_Hexa8(qualityType, coordinates);
            break;
        default:
            return;
    }

    // Find the Quality Key in Quality Map
    // Else Add new Quality into the map
    iter = quality.find(qualId);
    
    if (iter != quality.end())
        iter->second = qvalue;
    else
        quality.insert(pair<int, double>(qualId, qvalue));

    return;
}

//! Analyze the Normal of the Cells
void MESHQACellObject::analyze_normal()
{
    // TODO Compute Normals
    // Check the state
    if (!state) return;

    double qvalue = 0.0;

    // Verify the Cell Type
    // Check if vectors are appropriatly created
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            // TODO Not Implemented Yet
            return;
            break;
        case MESHQAEnums::CT_Tri3:
            if (ver.size() != MESHQA_NNODES_TRI3) return;
            if (normal.size() != MESHQA_NFACES_TRI3*3) {
                normal.clear();
                for (int i = 0; i < MESHQA_NFACES_TRI3*3; i++)
                    normal.push_back(qvalue);
            } else {
                for (int i = 0; i < MESHQA_NFACES_TRI3*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Quad4:
            if (ver.size() != MESHQA_NNODES_QUAD4) return;
            if (normal.size() != MESHQA_NFACES_QUAD4*3) {
                normal.clear();
                for (int i = 0; i < MESHQA_NFACES_QUAD4*3; i++)
                    normal.push_back(qvalue);
            } else {
                for (int i = 0; i < MESHQA_NFACES_QUAD4*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Tetra4:
            if (ver.size() != MESHQA_NNODES_TETRA4) return;
            if (normal.size() != MESHQA_NFACES_TETRA4*3) {
                normal.clear();
                for (int i = 0; i < MESHQA_NFACES_TETRA4*3; i++)
                    normal.push_back(qvalue);
            } else {
                for (int i = 0; i < MESHQA_NFACES_TETRA4*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Pyra5:
            if (ver.size() != MESHQA_NNODES_PYRA5) return;
            if (normal.size() != MESHQA_NFACES_PYRA5*3) {
                normal.clear();
                for (int i = 0; i < MESHQA_NFACES_PYRA5*3; i++)
                    normal.push_back(qvalue);
            } else {
                for (int i = 0; i < MESHQA_NFACES_PYRA5*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Prism6:
            if (ver.size() != MESHQA_NNODES_PRISM6) return;
            if (normal.size() != MESHQA_NFACES_PRISM6*3) {
                normal.clear();
                for (int i = 0; i < MESHQA_NFACES_PRISM6*3; i++)
                    normal.push_back(qvalue);
            } else {
                for (int i = 0; i < MESHQA_NFACES_PRISM6*3; i++)
                    normal[i] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Hexa8:
            if (ver.size() != MESHQA_NNODES_HEXA8) return;
            if (normal.size() != MESHQA_NFACES_HEXA8*3) {
                normal.clear();
                for (int i = 0; i < MESHQA_NFACES_HEXA8*3; i++)
                    normal.push_back(qvalue);
            } else {
                for (int i = 0; i < MESHQA_NFACES_HEXA8*3; i++)
                    normal[i] = 0.0;
            }
            break;
        default:
            return;
    }

    // Set Normal State
    normalstate = 1;
    return;
}

//! Returns the Quality Value
double MESHQACellObject::get_Quality(MESHQAEnums::QualityType qualityType)
{
    // Check if Cell Object is Initialized or Analyzed
    if (!state)
        return 0.0;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0.0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0.0;

    // Get the Quality ID for Cell Type and Return the value
    int qualId = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    map<int, double>::iterator iter;
    // Find the Quality in the map
    iter = quality.find(qualId);
    if (iter != quality.end())
        return iter->second;

    return 0.0;
}

//! Returns the Normal Vector for Cell Type
//! Returns normal computed vector if Object is Analyzed
//! Returns the empty vector if not Analyzed
vector<double> MESHQACellObject::get_normal() const
{
    vector<double> vec;
    if (normalstate)
        return normal;

    return vec;
}

//! Analyze the Cell Node Orientation
//! Return True if Orientation are according to SIDS conventions
//! Return False if Orientation fails
bool MESHQACellObject::analyze_orientation()
{
    // TODO
    return true;
}

