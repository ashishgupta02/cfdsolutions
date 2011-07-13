/*******************************************************************************
 * File:        Edge2D.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#ifndef _EDGE2D_H
#define	_EDGE2D_H

class Edge2D {
public:
    int node1;
    int node2;
    int cell1;
    int cell2;
    int Flag;
    int nData;
    double *Data;
public:
    Edge2D();
    ~Edge2D();
};

#endif	/* _EDGE2D_H */

