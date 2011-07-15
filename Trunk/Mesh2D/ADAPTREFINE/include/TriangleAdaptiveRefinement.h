/*******************************************************************************
 * File:        TriangleAdaptiveRefinement.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _TRIANGLEADAPTIVEREFINEMENT_H
#define	_TRIANGLEADAPTIVEREFINEMENT_H

#include "TriangleMeshDB.h"
#include "Edge2D.h"

class TriangleAdaptiveRefinement: virtual public TriangleMeshDB {
protected:
    int    NNode_Refine;
    int    NNode_Old;
    int    NTri_Old;
    int    **Cell2Node_6;
    int    *NodeCoarsen;
    // Solution Fields
    int    NField;
    int    *FieldTag;
    double **Field;
    char   **FieldName;
    // Adaptation Function
    int    AFType;
    int    Frefine;
    int    Fcoarsen;
    int    Fcoarsen_Old;
    int    FuncType;
    int    NAdapt;
    int    NFieldAdapt;
    int    ForceTriMesh;
    int    NLimitCoarsen;
    int    NLimitRefine;
    double AFrefine;
    double AFcoarsen;
    double AFavg;
    double Sigma;
    double Cr;
    double Cc;
    double P;
public:
    TriangleAdaptiveRefinement();
    virtual ~TriangleAdaptiveRefinement();
    void   AdaptiveRefinement();
    void   Get_Input_Parameters();
protected:
    virtual double Analytic_Function(double CoordX, double CoordY);
    double Compute_Tri_AspectRatio(double x1, double y1, double x2, double y2, double x3, double y3);
    void   Compute_Gauss_Gradient(int IBndy, double *F, double *Fx, double *Fy);
    double Adaptation_Function(int Node1, int Node2, double *F, double *Fx, double *Fy);
    void   Compute_Adaptation_Threshold(double *F, double *Fx, double *Fy);
    void   Create_Tri6_Connectivity();
    void   Delete_Tri6_Connectivity();
    void   Delete_Mesh_Connectivity();
    void   Get_Field_Data();
    void   Compute_Derived_FlowField();
    void   Tag_FlowField();
    int    Get_NumberOfEdges();
    // Refinement Functions
    void   Refine_Mesh();
    void   Refine_Get_Edge_Data(int NEdge, Edge2D *RefineEdge, double *RefineTreshold);
    int    Refine_Triangle(int conn[6], int tri[][3]);
    void   Adaptation_Refine_Mark_Edges(int NEdge, Edge2D *RefineEdge, List& CandiEdge);
    void   Adaptation_Refine_Mark_Edges_Old(double *F, double *Fx, double *Fy);
    void   Get_Potential_Refine_Edges(int NEdge, Edge2D *RefineEdge, double *RefineTreshold);
    void   Generate_Refine_New_Nodes_And_Solution();
    void   Update_Refine_Boundaries();
    void   Update_Refine_Triangles();
    // Coarsening Functions
    void   Coarsen_Mesh();
    void   Adaptation_Coarsen_Mark_Nodes(double *F, double *Fx, double *Fy);
    int    Get_Potential_Coarsen_Nodes(List& CandiNodes, double *F, double *Fx, double *Fy);
    int    Finalize_Coarsen_Mark_Nodes();
    int    Compress_Coarsen_Nodes(List& FinalCandiNodes);
    void   Generate_BowyerWatson_Delaunay_TriMesh();
private:
    void Init();
};

#endif	/* _TRIANGLEADAPTIVEREFINEMENT_H */

