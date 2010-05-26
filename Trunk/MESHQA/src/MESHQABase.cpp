/*******************************************************************************
 * File:        MESHQABase.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include "MESHQABase.h"
#include "MESHQAZone.h"

//! Default Constructors : Initailize the Data Members
MESHQABase::MESHQABase()
{
    CFDDatabaseBase = NULL;
    state       = 0;
    celldim     = 0;
    phydim      = 0;
    nzones      = 0;
    id          = 0;
    name.clear();
    zones       = NULL;
}

//! Overloaded Constructor
//! Just sets the CFD DatabaseBase pointer at the time of Object creation
MESHQABase::MESHQABase(BASE* CFDDBBase)
{
    CFDDatabaseBase = NULL;
    if (CFDDBBase != NULL)
        CFDDatabaseBase = CFDDBBase;
    state  = 0;
    celldim     = 0;
    phydim      = 0;
    nzones      = 0;
    id          = 0;
    name.clear();
    zones       = NULL;
}

//! Copy Constructor
//! Provide Way to have Deep Copy
MESHQABase::MESHQABase(const MESHQABase& orig)
{
    cout << "Not Implemented: MESHQABase::MESHQABase(const MESHQABase& orig)" << endl;
// TODO
}

//! Distructor for MESHQABase Object
//! Calls reset to free all the used resources
MESHQABase::~MESHQABase()
{
    reset();
}

// Operations related to Database
//! Returns the number of celldim in MESHQABase Object
//! celldim = 0 => if state = 0
//! celldim = 2 => if state = 1 => Face Cells
//! celldim = 3 => if state = 1 => Volume Cells
int MESHQABase::get_celldim() const
{
    if (state)
        return celldim;
    else
        return 0;
}

//! Returns the number of phydim in MESHQABase Object
//! phydim = 0 => if state = 0
//! phydim = 1 => if state = 1 => 1D
//! phydim = 2 => if state = 1 => 2D
//! phydim = 3 => if state = 1 => 3D
int MESHQABase::get_phydim() const
{
    if (state)
        return phydim;
    else
        return 0;
}

//! Returns the number of zones in MESHQABase Object
int MESHQABase::get_nzones() const
{
    if (state)
        return nzones;
    else
        return 0;
}

//! Returns the MESHQABase Object State
//! state = 0 => Object is Empty with No Data
//! state = 1 => Object contains Data
int MESHQABase::get_state() const
{
    return state;
}

//! Sets the CFDDBBase in MESHQABase Object
//! Returns True if Success
//! Returns False if Fails => Already Other Database is Avaliable
bool MESHQABase::set_CFDDatabaseBase(BASE* CFDDBBase)
{
    if (state)
        return false;
    else {
        if (CFDDBBase != NULL)
            CFDDatabaseBase = CFDDBBase;
        else
            return false;
    }
    return true;
}

//! Checks for Availablity of CFDDBBase Database and recursively fills the MESHQABase
//! returns none => Check the state if successfull
bool MESHQABase::update()
{
    // if Root already contains Data return
    if (state) return true;

    // Check if data base is available
    if (CFDDatabaseBase == NULL) return false;

    // Check if Database has any data
    if (CFDDatabaseBase->nzones <= 0) return false;

    // Check if Database actually have any data
    if (CFDDatabaseBase->zones == NULL) return false;

    // Create the Database
    celldim = CFDDatabaseBase->celldim;
    phydim  = CFDDatabaseBase->phydim;
    nzones  = CFDDatabaseBase->nzones;
    zones  = new MESHQAZone[nzones];
    for (int i = 0; i < nzones; i++) {
        zones[i].id     = i;
        zones[i].baseid = id;
        zones[i].name   = CFDDatabaseBase->zones[i].name;
        // Set CFDDatabaseZone Data
        if (!zones[i].set_CFDDatabaseZone(&CFDDatabaseBase->zones[i])) {
            delete[] zones;
            nzones  = 0;
            phydim  = 0;
            celldim = 0;
            return false;
        }
        if (!zones[i].update()) {
            // Reset all the Zones which where already successfully updated
            for (int j = 0; j < i; j++)
                zones[j].reset();

            delete[] zones;
            nzones  = 0;
            phydim  = 0;
            celldim = 0;
            return false;
        }
    }
    
    // Now Set the State
    state = 1;
    return true;
}

//! Recursively frees the data acquired by MESHQABase Object
//! This function ensures that all memory is released
void MESHQABase::reset()
{
    // No data simply return
    if(!state) {
        CFDDatabaseBase = NULL;
        return;
    }

    for (int i = 0; i < nzones; i++)
        zones[i].reset();

    delete[] zones;
    nzones  = 0;
    celldim = 0;
    phydim  = 0;
    CFDDatabaseBase = NULL;
    state = 0;
    return;
}

// Basic Quality Operations

//! Global Initialize the Quality Data
void MESHQABase::initialize()
{
    // Nothing to Initialize
    if (!state) return;
    
    for (int i = 0; i < nzones; i++)
        zones[i].initialize();
}

//! Initialize perticular Zone Quality Data
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::initialize(int zoneid)
{
    // Nothing to Initialize
    if (!state) return;

    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].initialize();
}

//! Initialize perticular Zone and Perticular Elemset in the zone
//! Zoneid range: 0 to (nzones-1)
//! Elemsetid range: 0 to (nesets-1)
void MESHQABase::initialize(int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].initialize(elemsetid);
    }
}

//! Intialize Perticular Cell Type
void MESHQABase::initialize_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    for (int i = 0; i < nzones; i++)
        zones[i].initialize_celltype(cellType);
}

//! Initialize perticular Zone Quality Data and Perticular Cell Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::initialize_celltype(MESHQAEnums::CellType cellType, int zoneid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].initialize_celltype(cellType);
}

//! Initialize perticular Zone, Elemset and Cell Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::initialize_celltype(MESHQAEnums::CellType cellType, int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].initialize_celltype(cellType, elemsetid);
    }
}

//! Intialize Perticular Quality Type
void MESHQABase::initialize_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nzones; i++)
        zones[i].initialize_quality(qualityType);
}

//! Initialize perticular Zone and Perticular Quality Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::initialize_quality(MESHQAEnums::QualityType qualityType, int zoneid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].initialize_quality(qualityType);
}

//! Initialize perticular Zone, Elemset and Quality Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::initialize_quality(MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].initialize_quality(qualityType, elemsetid);
    }
}

//! Intialize Perticular Cell Type and Quality Type
void MESHQABase::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nzones; i++)
        zones[i].initialize_celltype_quality(cellType, qualityType);
}

//! Initialize perticular Zone and Perticular Quality Type and Cell Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].initialize_celltype_quality(cellType, qualityType);
}

//! Initialize perticular Zone, Elemset and Quality Type and Cell Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].initialize_celltype_quality(cellType, qualityType, elemsetid);
    }
}

//! Global Analyze the Quality Data
void MESHQABase::analyze()
{
    // Nothing to Analyze
    if (!state) return;

    print_statistics();
    
    for (int i = 0; i < nzones; i++)
        zones[i].analyze();
}

//! Analyze perticular Zone Quality Data
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::analyze(int zoneid)
{
    // Nothing to Analyze
    if (!state) return;
    
    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].analyze();
}

//! Analyze perticular Zone and Perticular Elemset in the zone
//! Zoneid range: 0 to (nzones-1)
//! Elemsetid range: 0 to (nesets-1)
void MESHQABase::analyze(int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].analyze(elemsetid);
    }
}

//! Analyze Perticular Cell Type
void MESHQABase::analyze_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    for (int i = 0; i < nzones; i++)
        zones[i].analyze_celltype(cellType);
}

//! Analyze perticular Zone Quality Data and Perticular Cell Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::analyze_celltype(MESHQAEnums::CellType cellType, int zoneid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].analyze_celltype(cellType);
}

//! Analyze perticular Zone, Elemset and Cell Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::analyze_celltype(MESHQAEnums::CellType cellType, int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].analyze_celltype(cellType, elemsetid);
    }
}

//! Analyze Perticular Quality Type
void MESHQABase::analyze_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;
    
    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nzones; i++)
        zones[i].analyze_quality(qualityType);
}

//! Analyze perticular Zone and Perticular Quality Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::analyze_quality(MESHQAEnums::QualityType qualityType, int zoneid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].analyze_quality(qualityType);
}

//! Analyze perticular Zone, Elemset and Quality Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::analyze_quality(MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].analyze_quality(qualityType, elemsetid);
    }
}

//! Analyze Perticular Cell Type and Quality Type
void MESHQABase::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nzones; i++)
        zones[i].analyze_celltype_quality(cellType, qualityType);
}

//! Analyze perticular Zone and Perticular Quality Type and Cell Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones))
        zones[zoneid].analyze_celltype_quality(cellType, qualityType);
}

//! Analyze perticular Zone, Elemset and Quality Type and Cell Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((zoneid >= 0) && (zoneid < nzones)) {
        int nesets = zones[zoneid].get_nesets();
        if ((elemsetid >= 0) && (elemsetid < nesets))
            zones[zoneid].analyze_celltype_quality(cellType, qualityType, elemsetid);
    }
}

//! Print Quality Statistics
void MESHQABase::print_statistics()
{
    // Nothing to Print
    if (!state) return;

    cout.width(70);
    cout.fill('+');
    cout << "\n";
    cout << "\t Base ID = " << id+1;
    cout << "\n" ;
    cout << "\t Base Name = " << name;
    cout << "\n";
    cout << "\t Number of Zones = " << nzones;
    cout << "\n";
}

//! Print Quality Statistics
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::print_statistics(int zoneid)
{

}

//! Print Quality Statistics
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::print_statistics(int zoneid, int elemsetid)
{

}

//! Print Quality Statistics of Perticular CellType
void MESHQABase::print_statistics_celltype(MESHQAEnums::CellType cellType)
{

}

//! Print Quality Statistics of Perticular CellType
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::print_statistics_celltype(MESHQAEnums::CellType cellType, int zoneid)
{

}

//! Print Quality Statistics of Perticular CellType
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::print_statistics_celltype(MESHQAEnums::CellType cellType, int zoneid, int elemsetid)
{

}

//! Print Quality Statistics of Perticular Quality Type
void MESHQABase::print_statistics_quality(MESHQAEnums::QualityType qualityType)
{

}

//! Print Quality Statistics of Perticular Quality Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::print_statistics_quality(MESHQAEnums::QualityType qualityType, int zoneid)
{

}

//! Print Quality Statistics of Perticular Quality Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::print_statistics_quality(MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid)
{

}

//! Print Quality Statistics of Perticular Cell Type and Quality Type
void MESHQABase::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{

}

//! Print Quality Statistics of Perticular Cell Type and Quality Type
//! Zoneid range: 0 to (nzones-1)
void MESHQABase::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid)
{

}

//! Print Quality Statistics of Perticular Cell Type and Quality Type
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQABase::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int zoneid, int elemsetid)
{

}

