/*******************************************************************************
 * File:        Ghost3D.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include "Utils.h"
#include "Ghost3D.h"

//------------------------------------------------------------------------------
//! Default Constructor
//------------------------------------------------------------------------------
Ghost3D::Ghost3D() {
    globalId        = -1;
    partition       = -1;
    matrix_id       = -1;
    id_in_owner     = -1;
    centroid.pos[0] = 0.0;
    centroid.pos[1] = 0.0;
    centroid.pos[2] = 0.0;
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
Ghost3D::~Ghost3D() {
    
}
