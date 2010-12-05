/*******************************************************************************
 * File:        Euler2D_Mesh_LinearElasticSmoother.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EULER2D_MESH_LINEARELASTICSMOOTHER_H
#define	_EULER2D_MESH_LINEARELASTICSMOOTHER_H

#include "MC.h"
#include "Mesh.h"

class Euler2D_Mesh_LinearElasticSmoother {
protected:
    int BlockSize;
    int VectorSize;   // No of Nodes
    int changeIndex;
    int Parent;
    MESH mesh;
    CELL *cell;
    EDGE *edge;
    NODE *node;
    BOUNDARYEDGE *boundaryEdge;
    BOUNDARYNODE *boundaryNode;
    int *BNTag;
    // Matrix Computation Compressed Row Storage
    MC_CRS LESmoothBlockMatrix;
public:
    Euler2D_Mesh_LinearElasticSmoother();
    virtual ~Euler2D_Mesh_LinearElasticSmoother();
    void Mesh_Smoother_Prepare();
    void Initialize_Mesh_Smoother(const char* FileName, MC_CRS *Object);
    void Initialize_Mesh_Smoother(MESH inputMesh, CELL *ptrCell, EDGE *ptrEdge,
            NODE *ptrNode, BOUNDARYEDGE *ptrBoundaryEdge, BOUNDARYNODE *ptrBoundaryNode,
            int *ptrBNTag, MC_CRS *Object);
    void Mesh_Smoother(int MaxIter, double Relaxation, double *U, double *V);
    void Mesh_Smoother_Finalize();
    int IsParent();
protected:
    // Create and Initialize CRS LESmoother Block Matrix
    void Create_CRS_LESmoothBlockMatrix(MC_CRS *Object);

    // Compute the Matrix for Linear Elastic Equation
    void Compute_CRS_LESmoothBlockMatrix();

    // Compute the Perturbation Vector V and U
    double Solve_CRS_LESmoothBlockMatrix(int Iteration, double Relax, double *U, double *V);

    // Geometric Properties
    double Get_Poission_Ratio(int CellID);
    double Get_Youngs_Modulus(int CellID);
    double Compute_Tri_ConditionNumber(int n0, int n1, int n2);
    double Compute_Tri_AspectRatio(int n0, int n1, int n2);
    double Compute_Tri_Area(int n0, int n1, int n2);

    // Copy From Mesh Code
    void WKA_MeshReader(const char* FileName);
    void WKA_ExtractCells(void);
    void Tag_Boundary_Nodes(void);
    void Compute_Geometric_Properties();
private:
    void Init();
    void Reset();
};

#endif	/* _EULER2D_MESH_LINEARELASTICSMOOTHER_H */

