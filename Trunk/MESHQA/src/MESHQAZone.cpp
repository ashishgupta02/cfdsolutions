/*******************************************************************************
 * File:        MESHQAZone.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "MESHQAZone.h"
#include "MESHQAElemset.h"

//! Default Constructors : Initailize the Data Members
MESHQAZone::MESHQAZone()
{
    CFDDatabaseZone = NULL;
    state       = 0;
    nverts      = 0;
    ncells      = 0;
    nesets      = 0;
    id          = 0;
    baseid      = 0;
    name.clear();
    esets       = NULL;
}

//! Overloaded Constructor
//! Just sets the CFD DatabaseZone pointer at the time of Object creation
MESHQAZone::MESHQAZone(ZONE* CFDDBZone)
{
    CFDDatabaseZone = NULL;
    if (CFDDBZone != NULL)
        CFDDatabaseZone = CFDDBZone;
    state  = 0;
    nverts      = 0;
    ncells      = 0;
    nesets      = 0;
    id          = 0;
    baseid      = 0;
    name.clear();
    esets       = NULL;
}

//! Copy Constructor
//! Provide Way to have Deep Copy
MESHQAZone::MESHQAZone(const MESHQAZone& orig)
{
    cout << "Not Implemented: MESHQAZone::MESHQAZone(const MESHQAZone& orig)" << endl;
// TODO
}

//! Distructor for MESHQAZone Object
//! Calls reset to free all the used resources
MESHQAZone::~MESHQAZone()
{
    reset();
}

// Operations related to Database
//! Returns the number of Vertex in Perticular Zone
int MESHQAZone::get_nverts() const
{
    if (state)
        return nverts;
    else
        return 0;
}

//! Returns the number of Cells in Perticular Zone
int MESHQAZone::get_ncells() const
{
    if (state)
        return ncells;
    else
        return 0;
}

//! Returns the number of Element Sets in MESHQAZone Object
int MESHQAZone::get_nesets() const
{
    if (state)
        return nesets;
    else
        return 0;
}

//! Returns the MESHQAZone Object State
//! state = 0 => Object is Empty with No Data
//! state = 1 => Object contains Data
int MESHQAZone::get_state() const
{
    return state;
}

//! Sets the CFDDBZone in MESHQAZone Object
//! Returns True if Success
//! Returns False if Fails => Already Other Database is Avaliable
bool MESHQAZone::set_CFDDatabaseZone(ZONE* CFDDBZone)
{
    if (state)
        return false;
    else {
        if (CFDDBZone != NULL)
            CFDDatabaseZone = CFDDBZone;
        else
            return false;
    }
    return true;
}

//! Checks for Availablity of CFDDBZone Database and recursively fills the MESHQAZone
//! returns none => Check the state if successfull
bool MESHQAZone::update()
{
    // if Root already contains Data return
    if (state) return true;

    // Check if data base is available
    if (CFDDatabaseZone == NULL) return false;

    // Check if Database has any data
    if (CFDDatabaseZone->nesets <= 0) return false;

    // Check if Database actually have any data
    if (CFDDatabaseZone->esets == NULL) return false;

    // Check if Database has Coordinates
    if (CFDDatabaseZone->nverts <= 0) return false;

    // Check if Database actually have Coordinates
    if (CFDDatabaseZone->verts == NULL) return false;
    
    // Create the Database
    nverts  = CFDDatabaseZone->nverts;
    ncells  = CFDDatabaseZone->dim[1];
    nesets  = CFDDatabaseZone->nesets;
    esets   = new MESHQAElemset[nesets];
    
    for (int i = 0; i < nesets; i++) {
        esets[i].id     = i;
        esets[i].zoneid = id;
        esets[i].baseid = baseid;
        esets[i].name   = CFDDatabaseZone->esets[i].name;
        
        // Set CFDDatabaseElemset Data
        if ((!esets[i].set_CFDDatabaseElemset(&CFDDatabaseZone->esets[i])) ||
                (!esets[i].set_CFDDatabaseVertex(CFDDatabaseZone->verts)))
        {
            delete[] esets;
            nverts  = 0;
            ncells  = 0;
            nesets  = 0;
            return false;
        }
        if (!esets[i].update()) {
            // Reset all the Elemsets which where already successfully updated
            for (int j = 0; j < i; j++)
                esets[j].reset();

            delete[] esets;
            nesets  = 0;
            ncells  = 0;
            nverts  = 0;
            return false;
        }
    }
    
    // Now Set the State
    state = 1;
    return true;
}

//! Recursively frees the data acquired by MESHQABase Object
//! This function ensures that all memory is released
void MESHQAZone::reset()
{
    // No data simply return
    if(!state) {
        CFDDatabaseZone = NULL;
        return;
    }

    for (int i = 0; i < nesets; i++)
        esets[i].reset();

    delete[] esets;
    nverts      = 0;
    ncells      = 0;
    nesets      = 0;
    CFDDatabaseZone = NULL;
    state  = 0;
    return;
}

// Basic Quality Operations

//! Global Initialize the Quality Data
void MESHQAZone::initialize()
{
    // Nothing to Initialize
    if (!state) return;
    
    for (int i = 0; i < nesets; i++)
        esets[i].initialize();
}

//! Initialize Perticular Elemset in the zone
//! Elemsetid range: 0 to (nesets-1)
void MESHQAZone::initialize(int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].initialize();
}

//! Intialize Perticular Cell Type
void MESHQAZone::initialize_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    for (int i = 0; i < nesets; i++)
        esets[i].initialize_celltype(cellType);
}

//! Initialize perticular Elemset and Cell Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::initialize_celltype(MESHQAEnums::CellType cellType, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].initialize_celltype(cellType);
}

//! Intialize Perticular Quality Type
void MESHQAZone::initialize_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nesets; i++)
        esets[i].initialize_quality(qualityType);
}

//! Initialize perticular Elemset and Quality Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::initialize_quality(MESHQAEnums::QualityType qualityType, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].initialize_quality(qualityType);
}

//! Intialize Perticular Cell Type and Quality Type
void MESHQAZone::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nesets; i++)
        esets[i].initialize_celltype_quality(cellType, qualityType);
}

//! Initialize perticular Elemset, Cell Type and Quality Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].initialize_celltype_quality(cellType, qualityType);
}

//! Global Analyze the Quality Data
void MESHQAZone::analyze()
{
    // Nothing to Initialize
    if (!state) return;

    print_statistics();
    
    for (int i = 0; i < nesets; i++)
        esets[i].analyze();
}

//! Analyze Perticular Elemset in the zone
//! Elemsetid range: 0 to (nesets-1)
void MESHQAZone::analyze(int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;
    
    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].analyze();
}

//! Analyze Perticular Cell Type
void MESHQAZone::analyze_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    for (int i = 0; i < nesets; i++)
        esets[i].analyze_celltype(cellType);
}

//! Analyze perticular Elemset and Cell Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::analyze_celltype(MESHQAEnums::CellType cellType, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].analyze_celltype(cellType);
}

//! Analyze Perticular Quality Type
void MESHQAZone::analyze_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;
    
    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;
    
    for (int i = 0; i < nesets; i++)
        esets[i].analyze_quality(qualityType);
}

//! Analyze perticular Elemset and Quality Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::analyze_quality(MESHQAEnums::QualityType qualityType, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].analyze_quality(qualityType);
}

//! Analyze Perticular Cell Type and Quality Type
void MESHQAZone::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nesets; i++)
        esets[i].analyze_celltype_quality(cellType, qualityType);
}

//! Analyze perticular Elemset, Cell Type and Quality Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((elemsetid >= 0) && (elemsetid < nesets))
        esets[elemsetid].analyze_celltype_quality(cellType, qualityType);
}

//! Print Quality Statistics
void MESHQAZone::print_statistics()
{
    // Nothing to Print
    if (!state) return;

    cout.width(70);
    cout.fill('-');
    cout << "\n";
    cout << "\t \t Zone ID = " << id+1;
    cout << "\n" ;
    cout << "\t \t Zone Name = " << name;
    cout << "\n";
    cout << "\t \t Number of Elemsets = " << nesets;
    cout << "\n";
}

//! Print Quality Statistics
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::print_statistics(int elemsetid)
{

}

//! Print Quality Statistics of Perticular CellType
void MESHQAZone::print_statistics_celltype(MESHQAEnums::CellType cellType)
{

}

//! Print Quality Statistics of Perticular CellType
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::print_statistics_celltype(MESHQAEnums::CellType cellType, int elemsetid)
{

}

//! Print Quality Statistics of Perticular Quality Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::print_statistics_quality(MESHQAEnums::QualityType qualityType)
{

}

//! Print Quality Statistics of Perticular Quality Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::print_statistics_quality(MESHQAEnums::QualityType qualityType, int elemsetid)
{

}

//! Print Quality Statistics of Perticular Cell Type and Quality Type
void MESHQAZone::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{

}

//! Print Quality Statistics of Perticular Cell Type and Quality Type
//! Elmsetid range: 0 to (nesets-1)
void MESHQAZone::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int elemsetid)
{

}

