/*******************************************************************************
 * File:        MESHQABase.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _MESHQABASE_H
#define	_MESHQABASE_H

#include "MESHQACommon.h"
#include "MESHQAEnums.h"
#include "corestruct.h"

using namespace std;

// Forward Declaration
class MESHQAZone;

class MESHQABase {
// Private Data - Attributes
private:
    int   state;
    BASE* CFDDatabaseBase;
// Protected Data - Attributes
protected:
    int celldim;
    int phydim;
    int nzones;
    MESHQAZone* zones;
// Public Data - Attributes
public:
    int id;
    string name;
// Public Member Functions
public:
    // Constructors, Distructor and Operators
    MESHQABase();
    MESHQABase(BASE* CFDDBBase);
    MESHQABase(const MESHQABase& orig);
    virtual ~MESHQABase();

    // Operations related to Database
    int  get_celldim() const;
    int  get_phydim() const;
    int  get_nzones() const;
    int  get_state() const;
    bool set_CFDDatabaseBase(BASE* CFDDBBase);
    bool update();
    void reset();

    // Basic Operations
    // Global Initialize
    void initialize();
    void initialize(int zoneid);
    void initialize(int zoneid, int elemsetid);
    // Intialize Perticular Cell Type
    void initialize_celltype(MESHQAEnums::CellType cellType);
    void initialize_celltype(MESHQAEnums::CellType cellType, int zoneid);
    void initialize_celltype(MESHQAEnums::CellType cellType, int zoneid, int elemsetid);
    // Intialize Perticular Quality Type
    void initialize_quality(MESHQAEnums::QualityType qualityType);
    void initialize_quality(MESHQAEnums::QualityType qualityType, int zoneid);
    void initialize_quality(MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid);
    // Intialize Perticular Cell Type and Quality Type
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid);
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid);
    // Global Analyze
    void analyze();
    void analyze(int zoneid);
    void analyze(int zoneid, int elemsetid);
    // Analyze Perticular Cell Type
    void analyze_celltype(MESHQAEnums::CellType cellType);
    void analyze_celltype(MESHQAEnums::CellType cellType, int zoneid);
    void analyze_celltype(MESHQAEnums::CellType cellType, int zoneid, int elemsetid);
    // Analyze Perticular Quality Type
    void analyze_quality(MESHQAEnums::QualityType qualityType);
    void analyze_quality(MESHQAEnums::QualityType qualityType, int zoneid);
    void analyze_quality(MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid);
    // Analyze Perticular Cell Type and Quality Type
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid);
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid);
protected:
    // Print Quality Statistics
    void print_statistics();
    void print_statistics(int zoneid);
    void print_statistics(int zoneid, int elemsetid);
    // Print Quality Statistics of Perticular CellType
    void print_statistics_celltype(MESHQAEnums::CellType cellType);
    void print_statistics_celltype(MESHQAEnums::CellType cellType, int zoneid);
    void print_statistics_celltype(MESHQAEnums::CellType cellType, int zoneid, int elemsetid);
    // Print Quality Statistics of Perticular Quality Type
    void print_statistics_quality(MESHQAEnums::QualityType qualityType);
    void print_statistics_quality(MESHQAEnums::QualityType qualityType, int zoneid);
    void print_statistics_quality(MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid);
    // Print Quality Statistics of Perticular Cell Type and Quality Type
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid);
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid);
};

#endif	/* _MESHQABASE_H */

