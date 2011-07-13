/*******************************************************************************
 * File:        Node3D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _NODE3D_H
#define	_NODE3D_H

#include <vector>
#include <map>
#include <set>

#include "Point3D.h"

using namespace std;

class Node3D : public Point3D {
public:
    int              id;
    int              globalId;
    vector<int>      cells;
    vector<int>      ghosts;
    map<int, double> average;
    set<int>         bcs;
public:
    Node3D();
    //Node3D(double x = 0., double y = 0., double z = 0.);
    ~Node3D();
};

#endif	/* _NODE3D_H */

