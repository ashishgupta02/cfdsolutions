/*******************************************************************************
 * File:        Face3D.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include "Utils.h"
#include "Node3D.h"
#include "Face3D.h"
#include "Commons.h"

//------------------------------------------------------------------------------
//! Default Constructor
//------------------------------------------------------------------------------
Face3D::Face3D() {
    id              = -1;
    bc              = INTERNAL;
    parentIndex     = -1;
    parent          = -1;
    neighbor        = -1;
    nodeCount       = -1;
    centroid.pos[0] = 0.0;
    centroid.pos[1] = 0.0;
    centroid.pos[2] = 0.0;
    normal.vec[0]   = 0.0;
    normal.vec[1]   = 0.0;
    normal.vec[2]   = 0.0;
    area            = 0.0;
    mdot            = 0.0;
    weightL         = 0.0;
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
Face3D::~Face3D() {

}

