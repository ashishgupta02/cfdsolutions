/*******************************************************************************
 * File:        MESHQACellObject.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _MESHQACELLOBJECT_H
#define	_MESHQACELLOBJECT_H

#include "MESHQACommon.h"
#include "MESHQAEnums.h"
#include <corestruct.h>
#include <vector>
#include <map>

using namespace std;

class MESHQACellObject {
// Private Data - Attributes
private:
    int state;
    int normalstate;
// Public Data - Attributes
public:
    int id;
// Protected Data - Attributes
protected:
    MESHQAEnums::CellType cellType;
    vector<VERTEX*>  ver;
    map<int, double> quality;
    //vector<double>  quality;
    vector<double>   normal;
// Constructors, Distructor and Operators
public:
    MESHQACellObject();
    MESHQACellObject(MESHQAEnums::CellType ctype);
//    MESHQACellObject(const MESHQACellObject& orig);
    virtual ~MESHQACellObject();
// Basic Quality Operations
public:
    void reset();
    void addVertex(VERTEX* pver);
    void setCellType(MESHQAEnums::CellType ctype);
    MESHQAEnums::CellType getCellType() const;

    // Global Initialize
    void intialize();
    // Intialize Perticular Quality Type
    void intialize_quality(MESHQAEnums::QualityType qualityType);
    // Intialize Normals
    void intialize_normal();
    
    // Global Analyze
    void analyze();
    // Analyze Perticular Quality Type
    void analyze_quality(MESHQAEnums::QualityType qualityType);
    // Analyze Normals
    void analyze_normal();
    // Analyze Orientation
    bool analyze_orientation();

    // Return the Requested Quality Value
    double get_Quality(MESHQAEnums::QualityType qualityType);
    // Return the Computed Normal Vector
    vector<double> get_normal() const;
};

#endif	/* _MESHQACELLOBJECT_H */

