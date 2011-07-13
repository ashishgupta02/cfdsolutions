/*******************************************************************************
 * File:        MESHQAQualityMetrics.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _MESHQAQUALITYMETRICS_H
#define	_MESHQAQUALITYMETRICS_H

#include "MESHQACommon.h"
#include "MESHQAEnums.h"
#include "verdict.h"

class MESHQAQualityMetrics {
public:
    // Egde Quality
    static double analyze_Edge2(MESHQAEnums::QualityType qt, double coordinates[][3]);
    // Triangle Quality
    static double analyze_Tri3(MESHQAEnums::QualityType qt, double coordinates[][3]);
    // Qualdrilateral Quality
    static double analyze_Quad4(MESHQAEnums::QualityType qt, double coordinates[][3]);
    // Tetra Quality
    static double analyze_Tetra4(MESHQAEnums::QualityType qt, double coordinates[][3]);
    // Pyramid Quality
    static double analyze_Pyra5(MESHQAEnums::QualityType qt, double coordinates[][3]);
    // Prism Quality
    static double analyze_Prism6(MESHQAEnums::QualityType qt, double coordinates[][3]);
    // Hexa Quality
    static double analyze_Hexa8(MESHQAEnums::QualityType qt, double coordinates[][3]);
};

#endif	/* _MESHQAQUALITYMETRICS_H */

