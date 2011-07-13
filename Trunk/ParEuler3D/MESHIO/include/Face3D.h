/*******************************************************************************
 * File:        Face3D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _FACE3D_H
#define	_FACE3D_H

#include <vector>
#include <map>

#include "Point3D.h"
#include "Vector3D.h"

class Node3D;

using namespace std;

class Face3D {
public:
    int              id;
    int              bc;
    int              parentIndex;
    int              parent;
    int              neighbor;
    int              nodeCount;
    Point3D          centroid;
    Vector3D         normal;
    vector<int>      nodes;
    map<int, double> average;
    double           area;
    double           mdot, weightL;
    
public:
    Face3D();
    ~Face3D();
};

#endif	/* _FACE3D_H */

