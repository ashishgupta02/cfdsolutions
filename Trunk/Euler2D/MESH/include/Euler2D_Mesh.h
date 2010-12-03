/*******************************************************************************
 * File:        Euler2D_Mesh.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EULER2D_MESH_H
#define	_EULER2D_MESH_H

#include "Mesh.h"
#include "Euler2D_Mesh_LinearElasticSmoother.h"

class Euler2D_Mesh {
protected:
    int changeIndex;
    MESH mesh;
    CELL *cell;
    EDGE *edge;
    NODE *node;
    BOUNDARYEDGE *boundaryEdge;
    BOUNDARYNODE *boundaryNode;
    int *BNTag;
    // Linear Elastic Smoother
    Euler2D_Mesh_LinearElasticSmoother LESmooth;
public:
    Euler2D_Mesh();
    virtual ~Euler2D_Mesh();
    // Read Edge Based Mesh File
    void WKA_MeshReader(const char* FileName);
    // Write Edge Based Mesh File
    void WKA_MeshWriter(const char* FileName);
    // Read Restart File
    void Read_RestartFile(const char* FileName);
    // Write Restart File
    void Write_RestartFile(const char* FileName);
    // Read Cell Based Mesh File
    void SLK_MeshReader(const char* FileName);
    // Export to SLK File
    void SLK_MeshWriter(const char* FileName);
    // Export Gnuplot Solution Data
    void Write_Solution_GnuplotFile(const char* FileName);
    // Export Tecplot File
    void Write_TecplotFile(const char* FileName);
    // Export VTK Unstructured .vtu file
    void Write_VTK_Unstructured_File(const char* FileName);
    // Export VTK Unstructured .vtu file for Debug Variables
    void Write_VTK_Unstructured_DebugFile(const char* FileName, double *Data);
protected:
    void WKA_ExtractCells(void);
    void Tag_Boundary_Nodes(void);
    void Compute_Geometric_Properties();
private:
    void Init();
};

#endif	/* _EULER2D_MESH_H */

