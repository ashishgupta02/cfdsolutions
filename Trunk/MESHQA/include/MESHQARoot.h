/*******************************************************************************
 * File:        MESHQARoot.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _MESHQAROOT_H
#define	_MESHQAROOT_H

#include "MESHQACommon.h"
#include "MESHQAEnums.h"
#include "corestruct.h"

// Forward Declaration
class MESHQABase;

class MESHQARoot {
// Private Data - Attributes
private:
    int   state;
    ROOT* CFDDatabase;
// Protected Data - Attributes
protected:
    int nbases;
    MESHQABase* bases;
// Public Member Functions
public:
    // Constructors, Distructor and Operators
    MESHQARoot();
    MESHQARoot(ROOT* CFDDB);
    MESHQARoot(const MESHQARoot& orig);
    virtual ~MESHQARoot();

    // Operations related to Database
    int  get_nbases() const;
    int  get_state() const;
    bool set_CFDDatabase(ROOT* CFDDB);
    bool update();
    void reset();

    // Basic Quality Operations
    // Global Initialize
    void initialize();
    void initialize(int baseid);
    void initialize(int baseid, int zoneid);
    void initialize(int baseid, int zoneid, int elemsetid);
    // Intialize Perticular Cell Type
    void initialize_celltype(MESHQAEnums::CellType cellType);
    void initialize_celltype(MESHQAEnums::CellType cellType, int baseid);
    void initialize_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid);
    void initialize_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid, int elemsetid);
    // Intialize Perticular Quality Type
    void initialize_quality(MESHQAEnums::QualityType qualityType);
    void initialize_quality(MESHQAEnums::QualityType qualityType, int baseid);
    void initialize_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid);
    void initialize_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid);
    // Intialize Perticular Cell Type and Quality Type
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid);
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid);
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid);
    // Global Analyze
    void analyze();
    void analyze(int baseid);
    void analyze(int baseid, int zoneid);
    void analyze(int baseid, int zoneid, int elemsetid);
    // Analyze Perticular Cell Type
    void analyze_celltype(MESHQAEnums::CellType cellType);
    void analyze_celltype(MESHQAEnums::CellType cellType, int baseid);
    void analyze_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid);
    void analyze_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid, int elemsetid);
    // Analyze Perticular Quality Type
    void analyze_quality(MESHQAEnums::QualityType qualityType);
    void analyze_quality(MESHQAEnums::QualityType qualityType, int baseid);
    void analyze_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid);
    void analyze_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid);
    // Analyze Perticular Cell Type and Quality Type
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid);
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid);
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid);
protected:
    // Print Quality Statistics
    void print_statistics();
    void print_statistics(int baseid);
    void print_statistics(int baseid, int zoneid);
    void print_statistics(int baseid, int zoneid, int elemsetid);
    // Print Quality Statistics of Perticular CellType
    void print_statistics_celltype(MESHQAEnums::CellType cellType);
    void print_statistics_celltype(MESHQAEnums::CellType cellType, int baseid);
    void print_statistics_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid);
    void print_statistics_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid, int elemsetid);
    // Print Quality Statistics of Perticular Quality Type
    void print_statistics_quality(MESHQAEnums::QualityType qualityType);
    void print_statistics_quality(MESHQAEnums::QualityType qualityType, int baseid);
    void print_statistics_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid);
    void print_statistics_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid);
    // Print Quality Statistics of Perticular Cell Type and Quality Type
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid);
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid);
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid);
};

#endif	/* _MESHQAROOT_H */

