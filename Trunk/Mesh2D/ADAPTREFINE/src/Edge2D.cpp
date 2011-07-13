/*******************************************************************************
 * File:        Edge2D.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include <stdlib.h>
#include "Edge2D.h"

// *****************************************************************************
// *****************************************************************************
Edge2D::Edge2D() {
    node1 = -1;
    node2 = -1;
    cell1 = -1;
    cell2 = -1;
    Flag  = 0;
    nData = 0;
    Data  = NULL;
}

// *****************************************************************************
// *****************************************************************************
Edge2D::~Edge2D() {
    if (Data != NULL)
        delete[] Data;
}

