/*******************************************************************************
 * File:        MESHQAElemset.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include "MESHQAElemset.h"
#include "MESHQACellObject.h"

//! Default Constructors : Initailize the Data Members
MESHQAElemset::MESHQAElemset()
{
    CFDDatabaseElemset = NULL;
    CFDDatabaseVertex  = NULL;
    state       = 0;
    type        = 0;
    id          = 0;
    zoneid      = 0;
    baseid      = 0;
    name.clear();
}

//! Overloaded Constructor
//! Just sets the CFD DatabaseElemset pointer at the time of Object creation
//! Just sets the CFD DatabaseVertex pointer at the time of Object creation
MESHQAElemset::MESHQAElemset(VERTEX* CFDDBVertex, ELEMSET* CFDDBElemset)
{
    CFDDatabaseElemset = NULL;
    CFDDatabaseVertex  = NULL;
    if (CFDDBElemset != NULL)
        CFDDatabaseElemset = CFDDBElemset;
    if (CFDDBVertex != NULL)
        CFDDatabaseVertex = CFDDBVertex;
    state  = 0;
    type        = 0;
    id          = 0;
    zoneid      = 0;
    baseid      = 0;
    name.clear();
}

//! Copy Constructor
//! Provide Way to have Deep Copy
MESHQAElemset::MESHQAElemset(const MESHQAElemset& orig)
{
    cout << "Not Implemented: MESHQAElemset::MESHQAElemset(const MESHQAElemset& orig)" << endl;
// TODO
}

//! Distructor for MESHQAElemset Object
//! Calls reset to free all the used resources
MESHQAElemset::~MESHQAElemset()
{
    reset();
}

// Operations related to Database

//! Returns the MESHQAElemset Object State
//! state = 0 => Object is Empty with No Data
//! state = 1 => Object contains Data
int MESHQAElemset::get_state() const
{
    return state;
}

//! Returns the Type of Elemset
//! Type = Tri3, Quad4, Tetra4, Pyra5, Prism6, Hexa8
//! Type = MIXED => contains more then one Type
int MESHQAElemset::get_type() const
{
    if (state)
        return type;
    else
        return 0;
}

//! Sets the CFDDBElemset in MESHQAElemset Object
//! Returns True if Success
//! Returns False if Fails => Already Other Database is Avaliable
bool MESHQAElemset::set_CFDDatabaseElemset(ELEMSET* CFDDBElemset)
{
    if (state)
        return false;
    else {
        if (CFDDBElemset != NULL)
            CFDDatabaseElemset = CFDDBElemset;
        else
            return false;
    }
    return true;
}

//! Sets the CFDDBVertex in MESHQAElemset Object
//! Returns True if Success
//! Returns False if Fails => Already Other Database is Avaliable
bool MESHQAElemset::set_CFDDatabaseVertex(VERTEX* CFDDBVertex)
{
    if (state)
        return false;
    else {
        if (CFDDBVertex != NULL)
            CFDDatabaseVertex = CFDDBVertex;
        else
            return false;
    }
    return true;
}

//! Checks for Availablity of CFDDBElemset Database and recursively fills the MESHQAElemset
//! returns none => Check the state if successfull
bool MESHQAElemset::update()
{
    int ncells = 0;
    int count  = 0;
    int cellid;
    MESHQACellObject CellObject;
    
    // if Root already contains Data return
    if (state) return true;

    // Check if data base is available
    if (CFDDatabaseElemset == NULL) return false;

    // Check if Database has any data
    if (CFDDatabaseElemset->csize <= 0) return false;

    // Check if Database actually have any data
    if (CFDDatabaseElemset->conn == NULL) return false;

    // Check if Coordinates are available
    if (CFDDatabaseVertex == NULL) return false;

    // Create the Database
    type = CFDDBType2MESHQAType(CFDDatabaseElemset->type);
    cellid = CFDDatabaseElemset->start;
    switch (type) {
        case MESHQAEnums::CT_Edge2:
            // TODO : Not Implemented Yet
            return false;
            break;
        case MESHQAEnums::CT_Tri3:
            // Remove already existing Triangles if any
            vTri3.clear();
            // Remove all vector related data from Cell Object
            CellObject.reset();
            ncells = CFDDatabaseElemset->csize/3;
            count = 0;
            for (int i = 0; i < ncells; i++) {
                CellObject.id = cellid;
                cellid++;
                CellObject.setCellType(MESHQAEnums::CT_Tri3);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count  ]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                count = count + 3;
                vTri3.add(CellObject);
                CellObject.reset();
            }
            break;
        case MESHQAEnums::CT_Quad4:
            // Remove already existing Quadrilaterals if any
            vQuad4.clear();
            // Remove all vector related data from Cell Object
            CellObject.reset();
            ncells = CFDDatabaseElemset->csize/4;
            count = 0;
            for (int i = 0; i < ncells; i++) {
                CellObject.id = cellid;
                cellid++;
                CellObject.setCellType(MESHQAEnums::CT_Quad4);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count  ]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                count = count + 4;
                vQuad4.add(CellObject);
                CellObject.reset();
            }
            break;
        case MESHQAEnums::CT_Tetra4:
            // Remove already existing Tetra if any
            vTetra4.clear();
            // Remove all vector related data from Cell Object
            CellObject.reset();
            ncells = CFDDatabaseElemset->csize/4;
            count = 0;
            for (int i = 0; i < ncells; i++) {
                CellObject.id = cellid;
                cellid++;
                CellObject.setCellType(MESHQAEnums::CT_Tetra4);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count  ]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                count = count + 4;
                vTetra4.add(CellObject);
                CellObject.reset();
            }
            break;
        case MESHQAEnums::CT_Pyra5:
            // Remove already existing Pyramid if any
            vPyra5.clear();
            // Remove all vector related data from Cell Object
            CellObject.reset();
            ncells = CFDDatabaseElemset->csize/5;
            count = 0;
            for (int i = 0; i < ncells; i++) {
                CellObject.id = cellid;
                cellid++;
                CellObject.setCellType(MESHQAEnums::CT_Pyra5);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count  ]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                count = count + 5;
                vPyra5.add(CellObject);
                CellObject.reset();
            }
            break;
        case MESHQAEnums::CT_Prism6:
            // Remove already existing Prism if any
            vPrism6.clear();
            // Remove all vector related data from Cell Object
            CellObject.reset();
            ncells = CFDDatabaseElemset->csize/6;
            count = 0;
            for (int i = 0; i < ncells; i++) {
                CellObject.id = cellid;
                cellid++;
                CellObject.setCellType(MESHQAEnums::CT_Prism6);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count  ]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+5]]);
                count = count + 6;
                vPrism6.add(CellObject);
                CellObject.reset();
            }
            break;
        case MESHQAEnums::CT_Hexa8:
            // Remove already existing Hexa if any
            vHexa8.clear();
            // Remove all vector related data from Cell Object
            CellObject.reset();
            ncells = CFDDatabaseElemset->csize/8;
            count = 0;
            for (int i = 0; i < ncells; i++) {
                CellObject.id = cellid;
                cellid++;
                CellObject.setCellType(MESHQAEnums::CT_Hexa8);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count  ]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+5]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+6]]);
                CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+7]]);
                count = count + 8;
                vHexa8.add(CellObject);
                CellObject.reset();
            }
            break;
        case MESHQA_MIXED:
            int curtype;
            int elemsize;
            // Remove already existing Cells if any
            vEdge2.clear();
            vTri3.clear();
            vQuad4.clear();
            vTetra4.clear();
            vPyra5.clear();
            vPrism6.clear();
            vHexa8.clear();
            // Remove all vector related data from Cell Object
            CellObject.reset();
            curtype = 0;
            count = 0;
            elemsize = CFDDatabaseElemset->end - CFDDatabaseElemset->start + 1;
            for (int i = 0; i < elemsize; i++) {
                curtype = CFDDBType2MESHQAType(CFDDatabaseElemset->conn[count]);
                CellObject.id = cellid;
                cellid++;
                switch(curtype) {
                    case MESHQAEnums::CT_Edge2:
                        // TODO : Not Implemented Yet
                        return false;
                        break;
                    case MESHQAEnums::CT_Tri3:
                        CellObject.setCellType(MESHQAEnums::CT_Tri3);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                        count = count + 4;
                        vTri3.add(CellObject);
                        break;
                    case MESHQAEnums::CT_Quad4:
                        CellObject.setCellType(MESHQAEnums::CT_Quad4);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                        count = count + 5;
                        vQuad4.add(CellObject);
                        break;
                    case MESHQAEnums::CT_Tetra4:
                        CellObject.setCellType(MESHQAEnums::CT_Tetra4);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                        count = count + 5;
                        vTetra4.add(CellObject);
                        break;
                    case MESHQAEnums::CT_Pyra5:
                        CellObject.setCellType(MESHQAEnums::CT_Pyra5);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+5]]);
                        count = count + 6;
                        vPyra5.add(CellObject);
                        break;
                    case MESHQAEnums::CT_Prism6:
                        CellObject.setCellType(MESHQAEnums::CT_Prism6);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+5]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+6]]);
                        count = count + 7;
                        vPrism6.add(CellObject);
                        break;
                    case MESHQAEnums::CT_Hexa8:
                        CellObject.setCellType(MESHQAEnums::CT_Hexa8);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+1]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+2]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+3]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+4]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+5]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+6]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+7]]);
                        CellObject.addVertex(&CFDDatabaseVertex[CFDDatabaseElemset->conn[count+8]]);
                        count = count + 9;
                        vHexa8.add(CellObject);
                        break;
                    case MESHQAEnums::CT_Undefined:
                        return false;
                        break;
                }
                // Remove all vector related data from Cell Object
                CellObject.reset();
            }
            break;
        case MESHQAEnums::CT_Undefined:
            return false;
            break;
    }
    
    // Now Set the State
    state = 1;
    return true;
}

//! Recursively frees the data acquired by MESHQAElemset Object
//! This function ensures that all memory is released
void MESHQAElemset::reset()
{
    // No data simply return
    if (!state) {
        CFDDatabaseElemset = NULL;
        CFDDatabaseVertex  = NULL;
        return;
    }

    // Edge2
    if (!vEdge2.empty())
        vEdge2.clear();

    // Tri3
    if (!vTri3.empty())
        vTri3.clear();

    // Quad4
    if (!vQuad4.empty())
        vQuad4.clear();

    // Tetra4
    if (!vTetra4.empty())
        vTetra4.clear();

    // Pyra5
    if (!vPyra5.empty())
        vPyra5.clear();

    // Prism6
    if (!vPrism6.empty())
        vPrism6.clear();

    // Hexa8
    if (!vHexa8.empty())
        vHexa8.clear();
    
    type = 0;
    CFDDatabaseElemset = NULL;
    CFDDatabaseVertex  = NULL;
    state  = 0;
    return;
}

//! Converts the CFDDB Cell Types to MESHQA Cell Type
//! Special Case for MIXED Cell Type Element Connectivity
int MESHQAElemset::CFDDBType2MESHQAType(int CFDDBCellType) {
    int ctype;

    switch(CFDDBCellType) {
        // BAR_2 => Edge2
        case 2:
            ctype = MESHQAEnums::CT_Edge2;
            break;
        // TRI_3
        case 5:
            ctype = MESHQAEnums::CT_Tri3;
            break;
        // QUAD_4
        case 7:
            ctype = MESHQAEnums::CT_Quad4;
            break;
        // TETRA_4
        case 10:
            ctype = MESHQAEnums::CT_Tetra4;
            break;
        // PYRA_5
        case 12:
            ctype = MESHQAEnums::CT_Pyra5;
            break;
        // PENTA_6
        case 14:
            ctype = MESHQAEnums::CT_Prism6;
            break;
        // HEXA_8
        case 17:
            ctype = MESHQAEnums::CT_Hexa8;
            break;
        // MIXED
        case 20:
            ctype = MESHQA_MIXED;
            break;
        default:
            ctype = MESHQAEnums::CT_Undefined;
    }
    return ctype;
}

// Basic Operations

//! Global Initialize the Quality Data
void MESHQAElemset::initialize()
{
    // Nothing to Initialize
    if (!state) return;

    // Edge2
    if (!vEdge2.empty())
        vEdge2.initialize();

    // Tri3
    if (!vTri3.empty())
        vTri3.initialize();

    // Quad4
    if (!vQuad4.empty())
        vQuad4.initialize();

    // Tetra4
    if (!vTetra4.empty())
        vTetra4.initialize();

    // Pyra5
    if (!vPyra5.empty())
        vPyra5.initialize();

    // Prism6
    if (!vPrism6.empty())
        vPrism6.initialize();

    // Hexa8
    if (!vHexa8.empty())
        vHexa8.initialize();

    return;
}

//! Intialize Perticular Cell Type
void MESHQAElemset::initialize_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Initialize
    if (!state) return;

    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            vEdge2.initialize();
            break;
        case MESHQAEnums::CT_Tri3:
            vTri3.initialize();
            break;
        case MESHQAEnums::CT_Quad4:
            vQuad4.initialize();
            break;
        case MESHQAEnums::CT_Tetra4:
            vTetra4.initialize();
            break;
        case MESHQAEnums::CT_Pyra5:
            vPyra5.initialize();
            break;
        case MESHQAEnums::CT_Prism6:
            vPrism6.initialize();
            break;
        case MESHQAEnums::CT_Hexa8:
            vHexa8.initialize();
            break;
        default:
            cout <<"MESHQAElemset::initialize_celltype: Exception" << endl;
    }
}

//! Intialize Perticular Quality Type
void MESHQAElemset::initialize_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;
    
    // Edge2
    if (!vEdge2.empty())
        vEdge2.initialize_quality(qualityType);

    // Tri3
    if (!vTri3.empty())
        vTri3.initialize_quality(qualityType);

    // Quad4
    if (!vQuad4.empty())
        vQuad4.initialize_quality(qualityType);

    // Tetra4
    if (!vTetra4.empty())
        vTetra4.initialize_quality(qualityType);

    // Pyra5
    if (!vPyra5.empty())
        vPyra5.initialize_quality(qualityType);

    // Prism6
    if (!vPrism6.empty())
        vPrism6.initialize_quality(qualityType);

    // Hexa8
    if (!vHexa8.empty())
        vHexa8.initialize_quality(qualityType);
}

//! Intialize Perticular Cell Type and Quality Type
void MESHQAElemset::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;

    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            vEdge2.initialize_quality(qualityType);
            break;
        case MESHQAEnums::CT_Tri3:
            vTri3.initialize_quality(qualityType);
            break;
        case MESHQAEnums::CT_Quad4:
            vQuad4.initialize_quality(qualityType);
            break;
        case MESHQAEnums::CT_Tetra4:
            vTetra4.initialize_quality(qualityType);
            break;
        case MESHQAEnums::CT_Pyra5:
            vPyra5.initialize_quality(qualityType);
            break;
        case MESHQAEnums::CT_Prism6:
            vPrism6.initialize_quality(qualityType);
            break;
        case MESHQAEnums::CT_Hexa8:
            vHexa8.initialize_quality(qualityType);
            break;
        default:
            cout <<"MESHQAElemset::initialize_celltype_quality: Exception" << endl;
    }
}

//! Global Analyze the Quality Data
void MESHQAElemset::analyze()
{
    // Nothing to Analyze
    if (!state) return;

    print_statistics();
    
    // Edge2
    if (!vEdge2.empty())
        vEdge2.analyze();

    // Tri3
    if (!vTri3.empty())
        vTri3.analyze();

    // Quad4
    if (!vQuad4.empty())
        vQuad4.analyze();

    // Tetra4
    if (!vTetra4.empty())
        vTetra4.analyze();

    // Pyra5
    if (!vPyra5.empty())
        vPyra5.analyze();

    // Prism6
    if (!vPrism6.empty())
        vPrism6.analyze();

    // Hexa8
    if (!vHexa8.empty())
        vHexa8.analyze();
}

//! Analyze Perticular Cell Type
void MESHQAElemset::analyze_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Analyze
    if (!state) return;
    
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            vEdge2.analyze();
            break;
        case MESHQAEnums::CT_Tri3:
            vTri3.analyze();
            break;
        case MESHQAEnums::CT_Quad4:
            vQuad4.analyze();
            break;
        case MESHQAEnums::CT_Tetra4:
            vTetra4.analyze();
            break;
        case MESHQAEnums::CT_Pyra5:
            vPyra5.analyze();
            break;
        case MESHQAEnums::CT_Prism6:
            vPrism6.analyze();
            break;
        case MESHQAEnums::CT_Hexa8:
            vHexa8.analyze();
            break;
        default:
            cout <<"MESHQAElemset::analyze_celltype: Exception" << endl;
    }
}

//! Analyze Perticular Quality Type
void MESHQAElemset::analyze_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;
    
    // Edge2
    if (!vEdge2.empty())
        vEdge2.analyze_quality(qualityType);

    // Tri3
    if (!vTri3.empty())
        vTri3.analyze_quality(qualityType);

    // Quad4
    if (!vQuad4.empty())
        vQuad4.analyze_quality(qualityType);

    // Tetra4
    if (!vTetra4.empty())
        vTetra4.analyze_quality(qualityType);

    // Pyra5
    if (!vPyra5.empty())
        vPyra5.analyze_quality(qualityType);

    // Prism6
    if (!vPrism6.empty())
        vPrism6.analyze_quality(qualityType);

    // Hexa8
    if (!vHexa8.empty())
        vHexa8.analyze_quality(qualityType);
}

//! Analyze Perticular Cell Type and Quality Type
void MESHQAElemset::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;

    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            vEdge2.analyze_quality(qualityType);
            break;
        case MESHQAEnums::CT_Tri3:
            vTri3.analyze_quality(qualityType);
            break;
        case MESHQAEnums::CT_Quad4:
            vQuad4.analyze_quality(qualityType);
            break;
        case MESHQAEnums::CT_Tetra4:
            vTetra4.analyze_quality(qualityType);
            break;
        case MESHQAEnums::CT_Pyra5:
            vPyra5.analyze_quality(qualityType);
            break;
        case MESHQAEnums::CT_Prism6:
            vPrism6.analyze_quality(qualityType);
            break;
        case MESHQAEnums::CT_Hexa8:
            vHexa8.analyze_quality(qualityType);
            break;
        default:
            cout <<"MESHQAElemset::analyze_celltype_quality: Exception" << endl;
    }
}

// Print Quality Statistics
void MESHQAElemset::print_statistics()
{
    // Nothing to Print
    if (!state) return;

    cout.width(70);
    cout.fill('=');
    cout << "\n";
    cout << "\t \t \t Elemset ID = " << id+1;
    cout << "\n" ;
    cout << "\t \t \t Elemset Name = " << name;
    cout << "\n";
}

// Print Quality Statistics of Perticular CellType
void MESHQAElemset::print_statistics_celltype(MESHQAEnums::CellType cellType)
{

}

// Print Quality Statistics of Perticular Quality Type
void MESHQAElemset::print_statistics_quality(MESHQAEnums::QualityType qualityType)
{

}

// Print Quality Statistics of Perticular Cell Type and Quality Type
void MESHQAElemset::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{

}

