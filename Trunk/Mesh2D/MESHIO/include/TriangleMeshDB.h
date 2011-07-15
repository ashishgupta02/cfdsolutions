/*******************************************************************************
 * File:        TriangleMeshDB.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _TRIANGLEMESHDB_H
#define	_TRIANGLEMESHDB_H

#include "List.h"

class TriangleMeshDB {
protected:
    int NNode;
    int NBlock;
    int NTri;
    int NConstant;
    int NVariable;
    int NBoundary;
    int *NBoundarySegments;
    int ***BoundarySegments;
    // Coordinates
    double *X;
    double *Y;
    double **Variable;
    double *Constant;
    char **VariableName;
    // Cell2Node Connectivity
    int (*Tri)[3];
    // Connectivity Data Structure
    int  **Cell2Cell;
    List **Node2Cell;
    List **Node2Node;
    // Tag Interior and Boundary Nodes
    int *IBTag;
public:
    TriangleMeshDB();
    virtual ~TriangleMeshDB();
    // File IO for Steve L. Karman
    void SLK_MeshReader(const char* FileName);
    void SLK_MeshWriter(const char* FileName);
    void SLK_GnuplotWriter(const char* FileName);
    // File IO for W. Kyle Anderson
    void WKA_MeshReader(const char* FileName);
    void WKA_MeshWriter(const char* FileName);
    void WKA_GnuplotWriter(const char* FileName);
    // Topology Operations
    void Create_Interior_Boundary_Tag();
    void Delete_Interior_Boundary_Tag();
    void Create_Connectivity(int Sort);
protected:
    void Create_Node2Cell_Connectivity();
    void Create_Cell2Cell_Connectivity();
    void Create_Node2Node_Connectivity(int Sort);
    void Delete_Connectivity();
private:
    void Init();
};


#endif	/* _TRIANGLEMESHDB_H */

