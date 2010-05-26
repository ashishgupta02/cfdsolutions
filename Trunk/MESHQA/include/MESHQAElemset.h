/*******************************************************************************
 * File:        MESHQAElemset.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _MESHQAELEMSET_H
#define	_MESHQAELEMSET_H

#include "MESHQACommon.h"
#include "MESHQAEnums.h"
#include "MESHQACellQualityVector.h"
#include "corestruct.h"

using namespace std;

class MESHQAElemset {
// Private Data - Attributes
private:
    int state;
    int type;
    ELEMSET* CFDDatabaseElemset;
    VERTEX*  CFDDatabaseVertex;
// Protected Data - Attributes
protected:
    MESHQACellQualityVector vEdge2;
    MESHQACellQualityVector vTri3;
    MESHQACellQualityVector vQuad4;
    MESHQACellQualityVector vTetra4;
    MESHQACellQualityVector vPyra5;
    MESHQACellQualityVector vPrism6;
    MESHQACellQualityVector vHexa8;
// Public Data - Attributes
public:
    int id;
    int zoneid;
    int baseid;
    string name;
// Public Member Functions
public:
    // Constructors, Distructor and Operators
    MESHQAElemset();
    MESHQAElemset(VERTEX* CFDDBVertex, ELEMSET* CFDDBElemset);
    MESHQAElemset(const MESHQAElemset& orig);
    virtual ~MESHQAElemset();

    // Operations related to Database
    int  get_state() const;
    int  get_type() const;
    bool set_CFDDatabaseVertex(VERTEX* CFDDBVertex);
    bool set_CFDDatabaseElemset(ELEMSET* CFDDBElemset);
    bool update();
    void reset();
    int CFDDBType2MESHQAType(int CFDDBCellType);
    
    // Basic Operations
    // Global Initialize
    void initialize();
    // Intialize Perticular Cell Type
    void initialize_celltype(MESHQAEnums::CellType cellType);
    // Intialize Perticular Quality Type
    void initialize_quality(MESHQAEnums::QualityType qualityType);
    // Intialize Perticular Cell Type and Quality Type
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    // Global Analyze
    void analyze();
    // Analyze Perticular Cell Type
    void analyze_celltype(MESHQAEnums::CellType cellType);
    // Analyze Perticular Quality Type
    void analyze_quality(MESHQAEnums::QualityType qualityType);
    // Analyze Perticular Cell Type and Quality Type
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
protected:
    // Print Quality Statistics
    void print_statistics();
    // Print Quality Statistics of Perticular CellType
    void print_statistics_celltype(MESHQAEnums::CellType cellType);
    // Print Quality Statistics of Perticular Quality Type
    void print_statistics_quality(MESHQAEnums::QualityType qualityType);
    // Print Quality Statistics of Perticular Cell Type and Quality Type
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
};

#endif	/* _MESHQAELEMSET_H */

