/*******************************************************************************
 * File:        Ghost3D.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _GHOST3D_H
#define	_GHOST3D_H

#include <vector>

#include "Point3D.h"
#include "Vector3D.h"

using namespace std;

class Ghost3D {
public:
    int         globalId;
    int         partition;
    int         matrix_id;
    int         id_in_owner;
    Point3D     centroid;
    vector<int> cells;
    double      p, T, rho, closest_wall_distance;
    // Gradients are stored as p,u,v,w,T in order
    Vector3D    v, grad[5];
    // TODO do we need these updates for ghosts?
    double      update[5];
public:
    Ghost3D();
    ~Ghost3D();
};

#endif	/* _GHOST3D_H */

