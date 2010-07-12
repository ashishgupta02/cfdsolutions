/*******************************************************************************
 * File:        Cell3D.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include "Utils.h"
#include "Node3D.h"
#include "Face3D.h"
#include "Cell3D.h"

//------------------------------------------------------------------------------
//! Default Constructor
//------------------------------------------------------------------------------
Cell3D::Cell3D() {
    globalId              = -1;
    nodeCount             = 0;
    faceCount             = 0;
    neighborCellCount     = 0;
    ghostCount            = 0;
    globalCellCount       = 0;
    volume                = 0.0;
    lengthScale           = 0.0;
    closest_wall_distance = 0.0;
    centroid.pos[0]       = 0.0;
    centroid.pos[1]       = 0.0;
    centroid.pos[2]       = 0.0;
}

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
Cell3D::~Cell3D() {

}

//------------------------------------------------------------------------------
//! Check if Cell has nodes in the list
//------------------------------------------------------------------------------
bool Cell3D::HaveNodes(int &nodelistsize, int nodelist []) {
    bool match;
    for (int i = 0; i < nodelistsize; ++i) {
        match = false;
        for (int j = 0; j < nodeCount; ++j) {
            if (nodelist[i] == nodes[j]) {
                match = true;
                break;
            }
        }
        if (!match)
            return false;
    }
    return true;
}

