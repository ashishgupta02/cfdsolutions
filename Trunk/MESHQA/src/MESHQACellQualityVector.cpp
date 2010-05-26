/*******************************************************************************
 * File:        MESHQACellQualityVector.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include "MESHQACellQualityVector.h"
#include "MESHQACellObject.h"

//! Default Constructor
MESHQACellQualityVector::MESHQACellQualityVector()
{
    state  = 0;
    nCells = 0;
    for (int i = 0; i < MESHQA_NQUALITYTYPES; i++) {
        nCellsAboveMaxQuality[i]    = 0;
        nCellsBelowMinQuality[i]    = 0;
        MaxQualityCellId[i]         = 0;
        MinQualityCellId[i]         = 0;
        MaxQuality[i]               = -MESHQA_DBL_MAX;
        MinQuality[i]               = MESHQA_DBL_MAX;
        StdDevQuality[i]            = 0.0;
        AvgQuality[i]               = 0.0;
        MaxPassQuality[i]           = 0.0;
        MinPassQuality[i]           = 0.0;
        QualityState[i]             = 0;
    }
}

//! Copy Constructor
//! Provide Way to have Deep Copy
MESHQACellQualityVector::MESHQACellQualityVector(const MESHQACellQualityVector& orig)
{
    cout << "Not Implemented: MESHQACellQualityVector::MESHQACellQualityVector(const MESHQACellQualityVector& orig)" << endl;
// TODO
}

//! Distructor
//! Calls clear to free all the used resources
MESHQACellQualityVector::~MESHQACellQualityVector()
{
    clear();
}

//! This function ensures that all memory is released
void MESHQACellQualityVector::clear()
{
    // No data simply return
    if (!state) {
        CellsVector.clear();
        nCells = 0;
    }

    for (int i = 0; i < nCells; i++) {
        CellsVector[i].reset();
    }
    CellsVector.clear();

    for (int i = 0; i < MESHQA_NQUALITYTYPES; i++) {
        nCellsAboveMaxQuality[i]    = 0;
        nCellsBelowMinQuality[i]    = 0;
        MaxQualityCellId[i]         = 0;
        MinQualityCellId[i]         = 0;
        MaxQuality[i]               = -MESHQA_DBL_MAX;
        MinQuality[i]               = MESHQA_DBL_MAX;
        StdDevQuality[i]            = 0.0;
        AvgQuality[i]               = 0.0;
        MaxPassQuality[i]           = 0.0;
        MinPassQuality[i]           = 0.0;
        QualityState[i]             = 0;
    }
    
    nCells = 0;
    state  = 0;
    return;
}

//! Checks if CellsVector has CellObjects
bool MESHQACellQualityVector::empty() const
{
    if (!state)
        return false;

    if(!CellsVector.empty())
        return false;
    
    return true;
}

//! Returns of the Size of CellsVector
size_t MESHQACellQualityVector::size() const
{
    if (!state)
        return 0;

    if (nCells != (int)CellsVector.size())
        cout << "MESHQACellQualityVector::size: Exception" << endl;

    return (size_t) nCells;
}

//! Adds CellObject to CellsVector
void MESHQACellQualityVector::add(const MESHQACellObject& cellObject)
{
    // Add an Cell with the value specified to the end of the vector
    CellsVector.push_back(cellObject);
    
    nCells++;
    state = 1;
    return;
}

//! Global Initialize the Quality Data
void MESHQACellQualityVector::initialize()
{
    // Nothing to Initialize
    if (!state)
        return;
    
    // Check if Cell Type is set for Cell Objects
    MESHQAEnums::CellType cellType;
    cellType = CellsVector[0].getCellType();
    if (cellType == MESHQAEnums::CT_Undefined)
        return;
    
    for (int i = 0; i < nCells; i++)
        CellsVector[i].intialize();
    
    // Set the Quality State and Values
    for (int i = 0; i < MESHQAEnums::GetNumberQualityType(cellType); i++) {
        nCellsAboveMaxQuality[i]    = 0;
        nCellsBelowMinQuality[i]    = 0;
        MaxQualityCellId[i]         = 0;
        MinQualityCellId[i]         = 0;
        MaxQuality[i]               = -MESHQA_DBL_MAX;
        MinQuality[i]               = MESHQA_DBL_MAX;
        StdDevQuality[i]            = 0.0;
        AvgQuality[i]               = 0.0;
        QualityState[i]             = 0;
    }
    
    return;
}

//! Intialize Perticular Quality Type
void MESHQACellQualityVector::initialize_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Initialize
    if (!state)
        return;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return;

    // Check if Cell Type is set for Cell Objects
    MESHQAEnums::CellType cellType;
    cellType = CellsVector[0].getCellType();
    if (cellType == MESHQAEnums::CT_Undefined)
        return;
    
    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return;
    
    for (int i = 0; i < nCells; i++) {
        CellsVector[i].intialize_quality(qualityType);
    }

    // Set the Quality State and Values
    nCellsAboveMaxQuality[qindex]   = 0;
    nCellsBelowMinQuality[qindex]   = 0;
    MaxQualityCellId[qindex]        = 0;
    MinQualityCellId[qindex]        = 0;
    MaxQuality[qindex]              = -MESHQA_DBL_MAX;
    MinQuality[qindex]              = MESHQA_DBL_MAX;
    StdDevQuality[qindex]           = 0.0;
    AvgQuality[qindex]              = 0.0;
    QualityState[qindex]            = 0;
    
    return;
}

//! Global Analyze the Quality Data
void MESHQACellQualityVector::analyze()
{
    // Nothing to Analyze
    if (!state)
        return;

    // Check if Cell Type is set for Cell Objects
    MESHQAEnums::CellType cellType;
    cellType = CellsVector[0].getCellType();
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    for (int i = 0; i < nCells; i++) {
        CellsVector[i].analyze();
    }

    // Update Statistical Data
    update_statistics(cellType);

    // Print Statistical Data
    print_statistics(cellType);
    
    // Set the Quality State and Values
    for (int i = 0; i < MESHQAEnums::GetNumberQualityType(cellType); i++)
        QualityState[i] = 1;
    
    return;
}

//! Analyze Perticular Quality Type
void MESHQACellQualityVector::analyze_quality(MESHQAEnums::QualityType qualityType)
{
    // Nothing to Analyze
    if (!state)
        return;
    
    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return;
    
    // Check if Cell Type is set for Cell Objects
    MESHQAEnums::CellType cellType;
    cellType = CellsVector[0].getCellType();
    if (cellType == MESHQAEnums::CT_Undefined)
        return;
    
    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return;

    // Check if Analyze is actually required
    if (QualityState[qindex] != 1) {

        for (int i = 0; i < nCells; i++) {
            CellsVector[i].analyze_quality(qualityType);
        }

        // Update Statistical Data
        update_statistics_quality(cellType, qualityType);
    }
    
    // Print the Quality Statistics
    print_statistics_quality(cellType, qualityType);
    
    // Set the Quality State
    QualityState[qindex] = 1;
    
    return;
}

//! Sets the Maximum Pass Value of the Quality associated with perticular cell type
void MESHQACellQualityVector::setMaxPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, double maxValue)
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return;

    MaxPassQuality[qindex] = maxValue;
    return;
}

//! Sets the Minimum Pass Value of the Quality associated with perticular cell type
void MESHQACellQualityVector::setMinPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, double minValue)
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return;

    MinPassQuality[qindex] = minValue;
    return;
}


//! Returns the Maximum Pass Value of the Quality associated with perticular cell type
double MESHQACellQualityVector::getMaxPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0.0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0.0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0.0;

    return MaxPassQuality[qindex];
}

//! Returns the Minimum Pass Value of the Quality associated with perticular cell type
double MESHQACellQualityVector::getMinPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0.0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0.0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0.0;

    return MinPassQuality[qindex];
}

//! Returns the Maximum Value of the Quality associated with perticular cell type
double MESHQACellQualityVector::getMaxQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0.0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0.0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0.0;

    return MaxQuality[qindex];
}

//! Returns the Minimum Value of the Quality associated with perticular cell type
double MESHQACellQualityVector::getMinQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0.0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0.0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0.0;

    return MinQuality[qindex];
}

//! Returns the Average Value of the Quality associated with perticular cell type
double MESHQACellQualityVector::getAvgQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0.0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0.0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0.0;

    return AvgQuality[qindex];
}

//! Returns the Standarad Deviation Value of the Quality associated with perticular cell type
double MESHQACellQualityVector::getStdDevQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0.0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0.0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0.0;

    return StdDevQuality[qindex];
}

//! Returns the Number of Cells above Maximum Value of the Quality associated with perticular cell type
int MESHQACellQualityVector::getNumberCellsAboveMaxQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0;

    return nCellsAboveMaxQuality[qindex];
}

//! Returns the Number of Cells below Minimum Value of the Quality associated with perticular cell type
int MESHQACellQualityVector::getNumberCellsBelowMinQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0;

    return nCellsBelowMinQuality[qindex];
}

//! Returns the Number of Cells failed for the Quality associated with perticular cell type
int MESHQACellQualityVector::getNumberCellsFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return 0;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return 0;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return 0;

    return (nCellsAboveMaxQuality[qindex]+nCellsBelowMinQuality[qindex]);
}

//! Returns the Cell ID with Maximum Value of the Quality associated with perticular cell type
int MESHQACellQualityVector::getMaxQualityCellId(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return -1;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return -1;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return -1;

    return MaxQualityCellId[qindex];
}

//! Returns the Cell ID with Minimum Value of the Quality associated with perticular cell type
int MESHQACellQualityVector::getMinQualityCellId(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const
{
    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return -1;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return -1;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return -1;

    return MinQualityCellId[qindex];
}

//! Returns the vector of Mesh Cell Objects Above Maximum value of Quality of Perticular Cell Type
vector<MESHQACellObject> MESHQACellQualityVector::getCellsAboveMaxFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Create the Empty Vector
    vector<MESHQACellObject> vec;

    // Nothing to Analyze
    if (!state)
        return vec;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return vec;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return vec;

    // Check if Cell Type is set for Cell Objects
    MESHQAEnums::CellType cType;
    cType = CellsVector[0].getCellType();
    if (cType == MESHQAEnums::CT_Undefined)
        return vec;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return vec;
    
    // Check the status of the given quality type
    if (QualityState[qindex] == 0)
        return vec;
    
    double MaxValue = MaxPassQuality[qindex];
    MESHQACellObject CellObject;
    int qvalue;
    for (int i = 0; i < nCells; i++) {
        qvalue = CellsVector[i].get_Quality(qualityType);
        if (qvalue > MaxValue) {
            CellObject = CellsVector[i];
            vec.push_back(CellObject);
        }
    }

    return vec;
}

//! Returns the vector of Mesh Cell Objects Below Minimum value of Quality of Perticular Cell Type
vector<MESHQACellObject> MESHQACellQualityVector::getCellsBelowMinFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Create the Empty Vector
    vector<MESHQACellObject> vec;

    // Nothing to Analyze
    if (!state)
        return vec;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return vec;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return vec;
    
    // Check if Cell Type is set for Cell Objects
    MESHQAEnums::CellType cType;
    cType = CellsVector[0].getCellType();
    if (cType == MESHQAEnums::CT_Undefined)
        return vec;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return vec;

    // Check the status of the given quality type
    if (QualityState[qindex] == 0)
        return vec;

    double MinValue = MinPassQuality[qindex];
    MESHQACellObject CellObject;
    int qvalue;
    for (int i = 0; i < nCells; i++) {
        qvalue = CellsVector[i].get_Quality(qualityType);
        if (qvalue < MinValue) {
            CellObject = CellsVector[i];
            vec.push_back(CellObject);
        }
    }

    return vec;
}

//! Returns the vector of Mesh Cell Objects with Failed Quality of Perticular Cell Type
//! First Failed Cells of Above Maximum Pass Value are Stored
//! Next Failed Cells Below Minimum Pass value are Stored
vector<MESHQACellObject> MESHQACellQualityVector::getCellsFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Create the Empty Vector
    vector<MESHQACellObject> vec;

    // Nothing to Analyze
    if (!state)
        return vec;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return vec;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return vec;
    
    // Check if Cell Type is set for Cell Objects
    MESHQAEnums::CellType cType;
    cType = CellsVector[0].getCellType();
    if (cType == MESHQAEnums::CT_Undefined)
        return vec;

    int qindex;
    // Get the Quality Index according to Cell Type
    qindex = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    // No Match for Quality for Cell Type
    if (qindex == -1)
        return vec;

    // Check the status of the given quality type
    if (QualityState[qindex] == 0)
        return vec;

    // Max Value Cells First
    double Value = MaxPassQuality[qindex];
    MESHQACellObject CellObject;
    int qvalue;
    for (int i = 0; i < nCells; i++) {
        qvalue = CellsVector[i].get_Quality(qualityType);
        if (qvalue > Value) {
            CellObject = CellsVector[i];
            vec.push_back(CellObject);
        }
    }

    // Min Value Cells Next
    Value = MinPassQuality[qindex];
    for (int i = 0; i < nCells; i++) {
        qvalue = CellsVector[i].get_Quality(qualityType);
        if (qvalue < Value) {
            CellObject = CellsVector[i];
            vec.push_back(CellObject);
        }
    }

    return vec;
}

//! Updates the Statistical Data Associated with the CellsVector
void MESHQACellQualityVector::update_statistics(MESHQAEnums::CellType cellType)
{
    // Nothing to Update
    if (!state)
        return;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    int nquality = MESHQAEnums::GetNumberQualityType(cellType);
    if (nquality == 0) return;

    double qvalue = 0.0;
    MESHQAEnums::QualityType qualityType;
    
    // Check if Statistical Data Update is actually required
    for (int i = 0; i < nquality; i++) {
        if (QualityState[i] == 1) continue;

        // Get the quality Type from integer index
        qualityType = MESHQAEnums::Int2QualityType(cellType,i);

        // Get the Initial value for Min Max
        MaxQuality[i] = -MESHQA_DBL_MAX;
        MinQuality[i] = MESHQA_DBL_MAX;
        MaxQualityCellId[i] = -1;
        MinQualityCellId[i] = -1;
        
        // Set Avg and Std Dev to Zero
        AvgQuality[i] = 0.0;
        StdDevQuality[i] = 0.0;
        // Set Cells Above Below Quality to Zero
        nCellsAboveMaxQuality[i] = 0;
        nCellsBelowMinQuality[i] = 0;
        
        // Compute the Avg and Simultaniously acquire the Min-Max
        for (int j = 0; j < nCells; j++) {
            qvalue = CellsVector[j].get_Quality(qualityType);
            if (MaxQuality[i] < qvalue) {
                MaxQuality[i] = qvalue;
                MaxQualityCellId[i] = CellsVector[j].id;
            }
            if (MinQuality[i] > qvalue) {
                MinQuality[i] = qvalue;
                MinQualityCellId[i] = CellsVector[j].id;
            }
            if (qvalue > MaxPassQuality[i])
                nCellsAboveMaxQuality[i]++;
            if (qvalue < MinPassQuality[i])
                nCellsBelowMinQuality[i]++;

            AvgQuality[i] = AvgQuality[i] + qvalue;
        }

        AvgQuality[i] = AvgQuality[i]/((double) nCells);

        // Now Compute the Standard Deviation
        for (int j = 0; j < nCells; j++) {
            qvalue = CellsVector[j].get_Quality(qualityType);
            StdDevQuality[i] = StdDevQuality[i] + ((qvalue - AvgQuality[i])*(qvalue - AvgQuality[i]));
        }
        StdDevQuality[i] = StdDevQuality[i]/((double) (nCells-1));
        StdDevQuality[i] = sqrt(StdDevQuality[i]);
    }
    return;
}

//! Updates the Statistical Data Associated with the CellsVector of Perticular Quality
void MESHQACellQualityVector::update_statistics_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Update
    if (!state)
        return;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    // Verify the Quality Type
    if (qualityType == MESHQAEnums::QT_Undefined)
        return;

    int qualId = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    if (qualId == -1) return;
    
    // Check if Statistical Data Update is actually required
    if (QualityState[qualId] == 1) return;
    
    // Get the Initial value for Min Max
    MaxQuality[qualId] = -MESHQA_DBL_MAX;
    MinQuality[qualId] = MESHQA_DBL_MAX;
    MaxQualityCellId[qualId] = -1;
    MinQualityCellId[qualId] = -1;

    // Set Avg and Std Dev to Zero
    AvgQuality[qualId] = 0.0;
    StdDevQuality[qualId] = 0.0;
    // Set Cells Above Below Quality to Zero
    nCellsAboveMaxQuality[qualId] = 0;
    nCellsBelowMinQuality[qualId] = 0;
    
    double qvalue;
    // Compute the Avg and Simultaniously acquire the Min-Max
    for (int j = 0; j < nCells; j++) {
        qvalue = CellsVector[j].get_Quality(qualityType);
        if (MaxQuality[qualId] < qvalue) {
            MaxQuality[qualId] = qvalue;
            MaxQualityCellId[qualId] = CellsVector[j].id;
        }
        if (MinQuality[qualId] > qvalue) {
            MinQuality[qualId] = qvalue;
            MinQualityCellId[qualId] = CellsVector[j].id;
        }
        if (qvalue > MaxPassQuality[qualId])
            nCellsAboveMaxQuality[qualId]++;
        if (qvalue < MinPassQuality[qualId])
            nCellsBelowMinQuality[qualId]++;

        AvgQuality[qualId] = AvgQuality[qualId] + qvalue;
    }
    
    AvgQuality[qualId] = AvgQuality[qualId]/((float) nCells);

    // Now Compute the Standard Deviation
    for (int j = 0; j < nCells; j++) {
        qvalue = CellsVector[j].get_Quality(qualityType);
        StdDevQuality[qualId] = StdDevQuality[qualId] + ((qvalue - AvgQuality[qualId])*(qvalue - AvgQuality[qualId]));
    }
    StdDevQuality[qualId] = StdDevQuality[qualId]/((float) (nCells-1));
    StdDevQuality[qualId] = sqrt(StdDevQuality[qualId]);
    
    return;
}

//! Sets the Default value of Qualities Associated with a Cell Type
void MESHQACellQualityVector::setDefaultPassValue(MESHQAEnums::CellType cellType)
{
    // Nothing to Update
    if (!state)
        return;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    // Get the number of Qualities of a perticular cell type
    int nqual = MESHQAEnums::GetNumberQualityType(cellType);
    if (nqual == 0) return;

    MESHQAEnums::QualityType qtype;
    for (int nq = 0; nq < nqual; nq++) {
        qtype = MESHQAEnums::Int2QualityType(cellType, nq);
        setDefaultPassValue_quality(cellType, qtype);
    }
    
    return;
}

//! Sets the Default Value of a Quality of Cell
void MESHQACellQualityVector::setDefaultPassValue_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Update
    if (!state)
        return;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    int qualId = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    if (qualId == -1) return;
    
    switch (cellType) {
        case MESHQAEnums::CT_Edge2:
            switch (qualityType) {
                case MESHQAEnums::QT_Length:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                default:
                    MaxPassQuality[qualId] = 0.0;
                    MinPassQuality[qualId] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Tri3:
            switch (qualityType) {
                case MESHQAEnums::QT_Aspect:
                    MaxPassQuality[qualId] = 1.3;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_Area:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_SmallestAngle:
                    MaxPassQuality[qualId] = 60.0;
                    MinPassQuality[qualId] = 30.0;
                    break;
                case MESHQAEnums::QT_LargestAngle:
                    MaxPassQuality[qualId] = 90.0;
                    MinPassQuality[qualId] = 30.0;
                    break;
                case MESHQAEnums::QT_Condition:
                    MaxPassQuality[qualId] = 1.3;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_NormalizedJacobian:
                    MaxPassQuality[qualId] = 1.155;
                    MinPassQuality[qualId] = 0.5;
                    break;
                case MESHQAEnums::QT_Shear:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.75;
                    break;
                case MESHQAEnums::QT_Shape:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.75;
                    break;
                case MESHQAEnums::QT_RelativeSize:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.75;
                    break;
                case MESHQAEnums::QT_ShapeSize:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.75;
                    break;
                default:
                    MaxPassQuality[qualId] = 0.0;
                    MinPassQuality[qualId] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Quad4:
            switch (qualityType) {
                case MESHQAEnums::QT_Aspect:
                    MaxPassQuality[qualId] = 4.0;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_Area:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_SmallestAngle:
                    MaxPassQuality[qualId] = 90.0;
                    MinPassQuality[qualId] = 45.0;
                    break;
                case MESHQAEnums::QT_LargestAngle:
                    MaxPassQuality[qualId] = 135.0;
                    MinPassQuality[qualId] = 90.0;
                    break;
                case MESHQAEnums::QT_Condition:
                    MaxPassQuality[qualId] = 4.0;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_Jacobian:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_NormalizedJacobian:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.5;
                    break;
                case MESHQAEnums::QT_Shear:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.3;
                    break;
                case MESHQAEnums::QT_Shape:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.3;
                    break;
                case MESHQAEnums::QT_RelativeSize:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.5;
                    break;
                case MESHQAEnums::QT_ShapeSize:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.2;
                    break;
                case MESHQAEnums::QT_Skew:
                    MaxPassQuality[qualId] = 0.5;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_Taper:
                    MaxPassQuality[qualId] = 0.7;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_Warpage:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.9;
                    break;
                case MESHQAEnums::QT_Stretch:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.25;
                    break;
                case MESHQAEnums::QT_Oddy:
                    MaxPassQuality[qualId] = 16.0;
                    MinPassQuality[qualId] = 0.0;
                    break;
                default:
                    MaxPassQuality[qualId] = 0.0;
                    MinPassQuality[qualId] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Tetra4:
            switch (qualityType) {
                case MESHQAEnums::QT_Aspect:
                    MaxPassQuality[qualId] = 3.0;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_AspectGamma:
                    MaxPassQuality[qualId] = 3.0;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_Condition:
                    MaxPassQuality[qualId] = 3.0;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_Jacobian:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_NormalizedJacobian:
                    MaxPassQuality[qualId] = 1.414;
                    MinPassQuality[qualId] = 0.5;
                    break;
                case MESHQAEnums::QT_Shear:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.3;
                    break;
                case MESHQAEnums::QT_Shape:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.3;
                    break;
                case MESHQAEnums::QT_RelativeSize:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.5;
                    break;
                case MESHQAEnums::QT_ShapeSize:
                    MaxPassQuality[qualId] = 1.01;
                    MinPassQuality[qualId] = 0.2;
                    break;
                case MESHQAEnums::QT_Volume:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                default:
                    MaxPassQuality[qualId] = 0.0;
                    MinPassQuality[qualId] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Pyra5:
            switch (qualityType) {
                case MESHQAEnums::QT_Volume:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                default:
                    MaxPassQuality[qualId] = 0.0;
                    MinPassQuality[qualId] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Prism6:
            switch (qualityType) {
                case MESHQAEnums::QT_Volume:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                default:
                    MaxPassQuality[qualId] = 0.0;
                    MinPassQuality[qualId] = 0.0;
            }
            break;
        case MESHQAEnums::CT_Hexa8:
            switch (qualityType) {
                case MESHQAEnums::QT_Aspect:
                    MaxPassQuality[qualId] = 4.0;
                    MinPassQuality[qualId] = 1.0;
                    break;
                case MESHQAEnums::QT_Condition:
                    MaxPassQuality[qualId] = 8.0;
                    MinPassQuality[qualId] = 0.999999;
                    break;
                case MESHQAEnums::QT_Jacobian:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_NormalizedJacobian:
                    MaxPassQuality[qualId] = 1.01;
                    MinPassQuality[qualId] = 0.5;
                    break;
                case MESHQAEnums::QT_Shear:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.3;
                    break;
                case MESHQAEnums::QT_Shape:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.3;
                    break;
                case MESHQAEnums::QT_RelativeSize:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.5;
                    break;
                case MESHQAEnums::QT_ShapeSize:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.2;
                    break;
                case MESHQAEnums::QT_Skew:
                    MaxPassQuality[qualId] = 0.5;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_Taper:
                    MaxPassQuality[qualId] = 0.4;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_Stretch:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.25;
                    break;
                case MESHQAEnums::QT_Oddy:
                    MaxPassQuality[qualId] = 20.0;
                    MinPassQuality[qualId] = -1e-06;
                    break;
                case MESHQAEnums::QT_Volume:
                    MaxPassQuality[qualId] = 1.0;
                    MinPassQuality[qualId] = 0.65;
                    break;
                case MESHQAEnums::QT_Diagonal:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                case MESHQAEnums::QT_Dimension:
                    MaxPassQuality[qualId] = 1e+30;
                    MinPassQuality[qualId] = 0.0;
                    break;
                default:
                    MaxPassQuality[qualId] = 0.0;
                    MinPassQuality[qualId] = 0.0;
            }
            break;
        default:
            MaxPassQuality[qualId] = 0.0;
            MinPassQuality[qualId] = 0.0;
            break;
    }
    return;
}

//! Prints the Analyzed Quality of Perticular Cell Type
void MESHQACellQualityVector::print_statistics(MESHQAEnums::CellType cellType)
{
    // Nothing to Update
    if (!state)
        return;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    if (nCells <= 0) return;

    // Get the number of Qualities of a perticular cell type
    int nqual = MESHQAEnums::GetNumberQualityType(cellType);
    if (nqual == 0) return;

    // This Method of printing should be replaced
    cout.width(110);
    cout.fill('=');
    cout << "\n";
    cout << "Summary of " << MESHQAEnums::CellTypeToString(cellType);
    cout << " Quality Metrics: Total = " << nCells << " " << MESHQAEnums::CellTypeToString(cellType);
    cout << "\n";
    cout.width(110);
    cout.fill('=');
    cout << "\n";
    cout.setf(ios::right);
    cout.width(20);
    cout.fill(' ');
    cout << "Function Name";
    cout.width(15);
    cout << "Average";
    cout.width(15);
    cout << "Std Dev";
    cout.width(15);
    cout << "Minimum";
    cout.width(15);
    cout << "(Cell ID)";
    cout.width(15);
    cout << "Maximum";
    cout.width(15);
    cout << "(Cell ID)";
    cout << endl;
    cout.width(110);
    cout.fill('-');
    cout << '\n';
    MESHQAEnums::QualityType qtype;
    for (int nq = 0; nq < nqual; nq++) {
        qtype = MESHQAEnums::Int2QualityType(cellType, nq);
        print_statistics_quality(cellType, qtype);
    }
    cout.width(110);
    cout.fill('-');
    cout << '\n';
    return;
}

//! Prints the Analyzed Perticular Quality of Perticular Cell Type
void MESHQACellQualityVector::print_statistics_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType)
{
    // Nothing to Update
    if (!state)
        return;

    // Verify the Cell Type
    if (cellType == MESHQAEnums::CT_Undefined)
        return;

    int qualId = MESHQAEnums::GetQualityIndex(cellType, qualityType);
    if (qualId == -1) return;

    cout.width(20);
    cout.fill(' ');
    cout.setf(ios::right);
    cout << MESHQAEnums::QualityTypeToString(qualityType);
    cout.precision(6);
    cout.width(15);
    cout << AvgQuality[qualId];
    cout.width(15);
    cout << StdDevQuality[qualId];
    cout.width(15);
    cout << MinQuality[qualId];
    cout.width(15);
    cout << MinQualityCellId[qualId];
    cout.width(15);
    cout << MaxQuality[qualId];
    cout.width(15);
    cout << MaxQualityCellId[qualId];
    cout << endl;

    return;
}

