/* 
 * File:   Euler2D_Mesh.h
 * Author: Ashish Gupta
 *
 * Created on March 21, 2010, 7:50 PM
 */

#ifndef _EULER2D_MESH_H
#define	_EULER2D_MESH_H

#include "Mesh.h"

class Euler2D_Mesh {
protected:
    int restart;
    int changeIndex;
    MESH mesh;
    CELL *cell;
    EDGE *edge;
    NODE *node;
    BOUNDARYEDGE *boundaryEdge;
    BOUNDARYNODE *boundaryNode;
    int *BNTag;
public:
    Euler2D_Mesh();
    virtual ~Euler2D_Mesh();
    // Read Edge Based Mesh File
    void WKA_MeshReader(const char* FileName);
    // Read Restart File
    void Read_RestartFile(const char* FileName);
    // Write Restart File
    void Write_RestartFile(const char* FileName);
    // Export to SLK File
    void SLK_MeshWriter(const char* FileName);
    // Export Gnuplot Solution Data
    void Write_Solution_GnuplotFile(const char* FileName);
    // Export Tecplot File
    void Write_TecplotFile(const char* FileName);
    // Export VTK Unstructured .vtu file
    void Write_VTK_Unstructured_File(const char* FileName);
protected:
    void WKA_ExtractCells(void);
    void Tag_Boundary_Nodes(void);
    void Compute_Geometric_Properties();
private:
    void Init();
};

#endif	/* _EULER2D_MESH_H */

