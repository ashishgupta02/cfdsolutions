/*******************************************************************************
 * File:        MESHQARoot.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "MESHQARoot.h"
#include "MESHQABase.h"

//! Default Constructors : Initailize the Data Members
MESHQARoot::MESHQARoot()
{
    CFDDatabase = NULL;
    state       = 0;
    nbases      = 0;
    bases       = NULL;
}

//! Overloaded Constructor
//! Just sets the CFD Database pointer at the time of Object creation
MESHQARoot::MESHQARoot(ROOT* CFDDB)
{
    CFDDatabase = NULL;
    if (CFDDB != NULL)
        CFDDatabase = CFDDB;
    state  = 0;
    nbases = 0;
    bases  = NULL;
}

//! Copy Constructor
//! Provide Way to have Deep Copy
MESHQARoot::MESHQARoot(const MESHQARoot& orig)
{
    cout << "Not Implemented: MESHQARoot::MESHQARoot(const MESHQARoot& orig)" << endl;
// TODO
}

//! Distructor for MESHQARoot Object
//! Calls reset to free all the used resources
MESHQARoot::~MESHQARoot()
{
    reset();
}

// Operations related to Database
//! Returns the number of bases in MESHQARoot Object
//! nbases = 0 => if state = 0
//! nbases = number => if state = 1
int MESHQARoot::get_nbases() const
{
    if (state)
        return nbases;
    else
        return 0;
}

//! Returns the MESHQARoot Object State
//! state = 0 => Object is Empty with No Data
//! state = 1 => Object contains Data
int MESHQARoot::get_state() const
{
    return state;
}

//! Sets the CFDDB in MESHQARoot Object
//! Returns True if Success
//! Returns False if Fails => Already Other Database is Avaliable
bool MESHQARoot::set_CFDDatabase(ROOT* CFDDB)
{
    if (state)
        return false;
    else {
        if (CFDDB != NULL)
            CFDDatabase = CFDDB;
        else
            return false;
    }
    return true;
}

//! Checks for Availablity of CFDDB Database and recursively fills the MESHQA
//! returns none => Check the state if successfull
bool MESHQARoot::update()
{
    // if Root already contains Data return
    if (state) return true;

    // Check if data base is available
    if (CFDDatabase == NULL) return false;

    // Check if Database has any data
    if (CFDDatabase->nbases <= 0) return false;

    // Check if Database actually have any data
    if (CFDDatabase->bases == NULL) return false;
    
    // Create the Database
    nbases = CFDDatabase->nbases;
    bases  = new MESHQABase[nbases];
    for (int i = 0; i < nbases; i++) {
        bases[i].id = i;
        bases[i].name = CFDDatabase->bases[i].name;
        // Set CFDDatabaseBase Data
        if (!bases[i].set_CFDDatabaseBase(&CFDDatabase->bases[i])) {
            delete[] bases;
            nbases = 0;
            return false;
        }
        if (!bases[i].update()) {
            // Reset all the Bases which where already successfully updated
            for (int j = 0; j < i; j++)
                bases[j].reset();

            delete[] bases;
            nbases = 0;
            return false;
        }
    }
    
    // Now Set the State
    state = 1;
    return true;
}

//! Recursively frees the data acquired by MESHQARoot Object
//! This function ensures that all memory is released
void MESHQARoot::reset()
{
    // No data simply return
    if(!state) {
        CFDDatabase = NULL;
        return;
    }

    for (int i = 0; i < nbases; i++)
        bases[i].reset();

    delete[] bases;
    nbases = 0;
    CFDDatabase = NULL;
    state = 0;
    return;
}

// Basic Quality Operations

//! Global Initialize the Quality Data
void MESHQARoot::initialize()
{
    // Nothing to Initialize
    if (!state) return;
    
    for (int i = 0; i < nbases; i++)
        bases[i].initialize();
}

//! Initialize perticular Base Quality Data
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::initialize(int baseid)
{
    // Nothing to Initialize
    if (!state) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].initialize();
}

//! Initialize perticular Base and Perticular Zone in the base
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::initialize(int baseid, int zoneid)
{
    // Nothing to Initialize
    if (!state) return;
    
    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].initialize(zoneid);
    }
}

//! Initialize perticular Base, Zone and Elemset
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::initialize(int baseid, int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].initialize(zoneid, elemsetid);
    }
}

//! Intialize Perticular Cell Type
void MESHQARoot::initialize_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    for (int i = 0; i < nbases; i++)
        bases[i].initialize_celltype(cellType);
}

//! Initialize perticular Base Quality Data and Perticular Cell Type
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::initialize_celltype(MESHQAEnums::CellType cellType, int baseid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].initialize_celltype(cellType);
}

//! Initialize perticular Base and Perticular Zone in the base and Cell Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::initialize_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].initialize_celltype(cellType, zoneid);
    }
}

//! Initialize perticular Base, Zone, Elemset and Cell Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::initialize_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].initialize_celltype(cellType, zoneid, elemsetid);
    }
}

//! Intialize Perticular Quality Type
void MESHQARoot::initialize_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nbases; i++)
        bases[i].initialize_quality(qualityType);
}

//! Initialize perticular Base Quality Data and Perticular Quality Type
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::initialize_quality(MESHQAEnums::QualityType qualityType, int baseid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].initialize_quality(qualityType);
}

//! Initialize perticular Base and Perticular Zone in the base and Quality Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::initialize_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].initialize_quality(qualityType, zoneid);
    }
}

//! Initialize perticular Base, Zone, Elemset and Quality Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::initialize_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].initialize_quality(qualityType, zoneid, elemsetid);
    }
}

//! Intialize Perticular Cell Type and Quality Type
void MESHQARoot::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;
    
    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nbases; i++)
        bases[i].initialize_celltype_quality(cellType, qualityType);
}

//! Initialize perticular Base Quality Data and Perticular Quality Type of CellType
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].initialize_celltype_quality(cellType, qualityType, baseid);
}

//! Initialize perticular Base and Perticular Zone in the base and Quality Type of CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].initialize_celltype_quality(cellType, qualityType, zoneid);
    }
}

//! Initialize perticular Base, Zone, Elemset and Quality Type of CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::initialize_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid)
{
    // Nothing to Initialize
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].initialize_celltype_quality(cellType, qualityType, zoneid, elemsetid);
    }
}

//! Global Analyze of Quality
void MESHQARoot::analyze()
{
    // Nothing to Analyze
    if (!state) return;
    
    print_statistics();

    for (int i = 0; i < nbases; i++)
        bases[i].analyze();
}

//! Analyze perticular Base Quality Data
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::analyze(int baseid)
{
    // Nothing to Analyze
    if (!state) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].analyze();
}

//! Analyze perticular Base and Perticular Zone in the base
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::analyze(int baseid, int zoneid)
{
    // Nothing to Analyze
    if (!state) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].analyze(zoneid);
    }
}

//! Analyze perticular Base, Zone and Elemset
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::analyze(int baseid, int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].analyze(zoneid, elemsetid);
    }
}

//! Analyze Perticular Cell Type
void MESHQARoot::analyze_celltype(MESHQAEnums::CellType cellType)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    for (int i = 0; i < nbases; i++)
        bases[i].analyze_celltype(cellType);
}

//! Analyze perticular Base Quality Data and Perticular Cell Type
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::analyze_celltype(MESHQAEnums::CellType cellType, int baseid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].analyze_celltype(cellType);
}

//! Analyze perticular Base and Perticular Zone in the base and Cell Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::analyze_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].analyze_celltype(cellType, zoneid);
    }
}

//! Initialize perticular Base, Zone, Elemset and Cell Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::analyze_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].analyze_celltype(cellType, zoneid, elemsetid);
    }
}

//! Analyze Perticular Quality Type
void MESHQARoot::analyze_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;
    
    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nbases; i++)
        bases[i].analyze_quality(qualityType);
}

//! Analyze perticular Base Quality Data and Perticular Quality Type
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::analyze_quality(MESHQAEnums::QualityType qualityType, int baseid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].analyze_quality(qualityType);
}

//! Analyze perticular Base and Perticular Zone in the base and Quality Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::analyze_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].analyze_quality(qualityType, zoneid);
    }
}

//! Analyze perticular Base, Zone, Elemset and Quality Type
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::analyze_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].analyze_quality(qualityType, zoneid, elemsetid);
    }
}

//! Analyze Perticular Cell Type and Quality Type
void MESHQARoot::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    for (int i = 0; i < nbases; i++)
        bases[i].analyze_celltype_quality(cellType, qualityType);
}

//! Analyze perticular Base Quality Data and Perticular Quality Type of CellType
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases))
        bases[baseid].analyze_celltype_quality(cellType, qualityType, baseid);
}

//! Analyze perticular Base and Perticular Zone in the base and Quality Type of CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            bases[baseid].analyze_celltype_quality(cellType, qualityType, zoneid);
    }
}

//! Analyze perticular Base, Zone, Elemset and Quality Type of CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::analyze_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid)
{
    // Nothing to Analyze
    if (!state) return;

    // Check the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined) return;

    // Check the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined) return;

    if ((baseid >= 0) && (baseid < nbases)) {
        int nzones = bases[baseid].get_nzones();
        if ((zoneid >= 0) && (zoneid < nzones))
            if (elemsetid >= 0)
                bases[baseid].analyze_celltype_quality(cellType, qualityType, zoneid, elemsetid);
    }
}

//! Print Quality Statistics
void MESHQARoot::print_statistics()
{
    // Nothing to Print
    if (!state) return;

    // This Method of printing should be replaced
    cout.width(110);
    cout.fill('*');
    cout << "\n";
    cout << "Analyzing Quality the Mesh:" ;
    cout << "\n" ;
    cout << "\t Number of Bases = " << nbases ;
    cout << "\n" ;
}

//! Print Quality Statistics
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::print_statistics(int baseid)
{
    // Nothing to Print
    if (!state) return;
}

//! Print Quality Statistics
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::print_statistics(int baseid, int zoneid)
{
    // Nothing to Print
    if (!state) return;
}

//! Print Quality Statistics
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::print_statistics(int baseid, int zoneid, int elemsetid)
{
    // Nothing to Print
    if (!state) return;
}

//! Print Quality Statistics of Perticular CellType
void MESHQARoot::print_statistics_celltype(MESHQAEnums::CellType cellType)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::print_statistics_celltype(MESHQAEnums::CellType cellType, int baseid)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::print_statistics_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::print_statistics_celltype(MESHQAEnums::CellType cellType, int baseid, int zoneid, int elemsetid)
{

}

//! Print Quality Statistics of Perticular Quality Type
void MESHQARoot::print_statistics_quality(MESHQAEnums::QualityType qualityType)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::print_statistics_quality(MESHQAEnums::QualityType qualityType, int baseid)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::print_statistics_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::print_statistics_quality(MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid)
{

}

//! Print Quality Statistics of Perticular Cell Type and Quality Type
void MESHQARoot::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
void MESHQARoot::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
void MESHQARoot::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid)
{

}

//! Print Quality Statistics of Perticular CellType
//! Baseid range: 0 to (nbases-1)
//! Zoneid range: 0 to (nzones-1)
//! Elmsetid range: 0 to (nesets-1)
void MESHQARoot::print_statistics_celltype_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, int baseid, int zoneid, int elemsetid)
{

}

