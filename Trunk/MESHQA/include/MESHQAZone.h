/*******************************************************************************
 * File:        MESHQAZone.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _MESHQAZONE_H
#define	_MESHQAZONE_H

#include "MESHQACommon.h"
#include "MESHQAEnums.h"
#include "corestruct.h"

using namespace std;

// Forward Declaration
class MESHQAElemset;

class MESHQAZone {
// Private Data - Attributes
private:
    int   state;
    ZONE* CFDDatabaseZone;
// Protected Data - Attributes
protected:
    int nverts;
    int ncells;
    int nesets;
    MESHQAElemset* esets;
// Public Data - Attributes
public:
    int id;
    int baseid;
    string name;
// Public Member Functions
public:
    // Constructors, Distructor and Operators
    MESHQAZone();
    MESHQAZone(ZONE* CFDDBZone);
    MESHQAZone(const MESHQAZone& orig);
    virtual ~MESHQAZone();

    // Operations related to Database
    int  get_nverts() const;
    int  get_ncells() const;
    int  get_nesets() const;
    int  get_state() const;
    bool set_CFDDatabaseZone(ZONE* CFDDBZone);
    bool update();
    void reset();

    // Basic Quality Operations
    // Global Initialize
    void initialize();
    void initialize(int elemsetid);
    // Intialize Perticular Cell Type
    void initialize_celltype(MESHQAEnums::CellType cellType);
    void initialize_celltype(MESHQAEnums::CellType cellType, int elemsetid);
    // Intialize Perticular Quality Type
    void initialize_quality(MESHQAEnums::QualityType qualityType);
    void initialize_quality(MESHQAEnums::QualityType qualityType, int elemsetid);
    // Intialize Perticular Cell Type and Quality Type
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int elemsetid);
    // Global Analyze
    void analyze();
    void analyze(int elemsetid);
    // Analyze Perticular Cell Type
    void analyze_celltype(MESHQAEnums::CellType cellType);
    void analyze_celltype(MESHQAEnums::CellType cellType, int elemsetid);
    // Analyze Perticular Quality Type
    void analyze_quality(MESHQAEnums::QualityType qualityType);
    void analyze_quality(MESHQAEnums::QualityType qualityType, int elemsetid);
    // Analyze Perticular Cell Type and Quality Type
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int elemsetid);
protected:
    // Print Quality Statistics
    void print_statistics();
    void print_statistics(int elemsetid);
    // Print Quality Statistics of Perticular CellType
    void print_statistics_celltype(MESHQAEnums::CellType cellType);
    void print_statistics_celltype(MESHQAEnums::CellType cellType, int elemsetid);
    // Print Quality Statistics of Perticular Quality Type
    void print_statistics_quality(MESHQAEnums::QualityType qualityType);
    void print_statistics_quality(MESHQAEnums::QualityType qualityType, int elemsetid);
    // Print Quality Statistics of Perticular Cell Type and Quality Type
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    void print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int elemsetid);
};

#endif	/* _MESHQAZONE_H */

