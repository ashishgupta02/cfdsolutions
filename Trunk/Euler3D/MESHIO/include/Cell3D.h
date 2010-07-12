/*******************************************************************************
 * File:        Cell3D.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _CELL3D_H
#define	_CELL3D_H

#include <vector>
#include <map>
#include <set>

#include "Point3D.h"
#include "Vector3D.h"

class Node3D;
class Face3D;

using namespace std;

class Cell3D {
public:
    int         globalId;
    int         nodeCount;
    int         faceCount;
    int         neighborCellCount;
    int         ghostCount;
    int         globalCellCount;
    double      volume;
    double      lengthScale;
    double      closest_wall_distance;
    Point3D     centroid;
    vector<int> nodes;
    vector<int> faces;
    vector<int> ghosts;
    vector<int> neighborCells;
    double      p, T, rho, dt;
    Vector3D    v, grad[5];
    // Gradients are stored as p,u,v,w,T in order
    map<int, Vector3D> gradMap;
    double      update[5];
public:
    Cell3D();
    ~Cell3D();
    bool HaveNodes(int &nodelistsize, int nodelist[]);
};

#endif	/* _CELL3D_H */

