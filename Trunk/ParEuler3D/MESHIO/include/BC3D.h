/*******************************************************************************
 * File:        BC3D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef BC3D_H
#define	BC3D_H

#include <vector>
#include <set>

using namespace std;

class BC3D {
public:
    vector<int>         boundaryFaceCount;
    vector<int>         boundaryNodeCount;
    vector<int>         globalBoundaryFaceCount;
    vector<int>         globalBoundaryNodeCount;
    vector< set<int> >  bocoNodes;
public:
    BC3D();
    ~BC3D();
};

#endif	/* BC3D_H */

