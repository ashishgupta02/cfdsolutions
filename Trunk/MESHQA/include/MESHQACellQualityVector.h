/*******************************************************************************
 * File:        MESHQACellQualityVector.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _MESHQACELLQUALITYVECTOR_H
#define	_MESHQACELLQUALITYVECTOR_H

#include "MESHQACommon.h"
#include "MESHQAEnums.h"

using namespace std;

// Forward Declaration
class MESHQACellObject;

class MESHQACellQualityVector {
// Private Data - Attributes
private:
    int state;
// Protected Data - Attributes
protected:
    int nCells;
    vector<MESHQACellObject> CellsVector;
    int    nCellsAboveMaxQuality[MESHQA_NQUALITYTYPES];
    int    nCellsBelowMinQuality[MESHQA_NQUALITYTYPES];
    int    MaxQualityCellId[MESHQA_NQUALITYTYPES];
    int    MinQualityCellId[MESHQA_NQUALITYTYPES];
    double MaxQuality[MESHQA_NQUALITYTYPES];
    double MinQuality[MESHQA_NQUALITYTYPES];
    double StdDevQuality[MESHQA_NQUALITYTYPES];
    double AvgQuality[MESHQA_NQUALITYTYPES];
    double MaxPassQuality[MESHQA_NQUALITYTYPES];
    double MinPassQuality[MESHQA_NQUALITYTYPES];
    int    QualityState[MESHQA_NQUALITYTYPES];
public:
    // Constructors and Distructor
    MESHQACellQualityVector();
    MESHQACellQualityVector(const MESHQACellQualityVector& orig);
    virtual ~MESHQACellQualityVector();
    
    // Basic Vector Operations
    void   clear();
    bool   empty() const;
    size_t size() const;
    void   add(const MESHQACellObject &cellObject);

    // Global Initialize
    void   initialize();
    // Intialize Perticular Quality Type
    void   initialize_quality(MESHQAEnums::QualityType qualityType);
    
    // Global Analyze
    void   analyze();
    // Analyze Perticular Quality Type
    void   analyze_quality(MESHQAEnums::QualityType qualityType);

    // Quality Assesment
    void      setMaxPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, double maxValue);
    void      setMinPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType, double minValue);
    double    getMaxPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    double    getMinPassQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    double    getMaxQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    double    getMinQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    double    getAvgQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    double    getStdDevQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    int       getNumberCellsAboveMaxQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    int       getNumberCellsBelowMinQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    int       getNumberCellsFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    int       getMaxQualityCellId(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    int       getMinQualityCellId(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType) const;
    vector<MESHQACellObject> getCellsFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    vector<MESHQACellObject> getCellsAboveMaxFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    vector<MESHQACellObject> getCellsBelowMinFailQuality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
protected:
    // Set Default the Max - Min Pass Value Range for Quality
    void   setDefaultPassValue(MESHQAEnums::CellType cellType);
    void   setDefaultPassValue_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    // Update the Quality Statistics of a Cell
    void   update_statistics(MESHQAEnums::CellType cellType);
    void   update_statistics_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
    // Print Quality Statistics
    void   print_statistics(MESHQAEnums::CellType cellType);
    void   print_statistics_quality(MESHQAEnums::CellType cellType, MESHQAEnums::QualityType qualityType);
};

#endif	/* _MESHQACELLQUALITYVECTOR_H */

