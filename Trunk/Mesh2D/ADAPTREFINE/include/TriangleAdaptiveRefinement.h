/* 
 * File:   TriangleAdaptiveRefinement.h
 * Author: Ashish Gupta
 *
 * Created on March 15, 2010, 8:06 PM
 */

#ifndef _TRIANGLEADAPTIVEREFINEMENT_H
#define	_TRIANGLEADAPTIVEREFINEMENT_H

#include "TriangleMeshDB.h"

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
    int    FuncType;
    int    NAdapt;
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
    void AdaptiveRefinement();
    void Get_Input_Parameters();
protected:
    virtual double Analytic_Function(double CoordX, double CoordY);
    void Create_Tri6_Connectivity();
    void Delete_Tri6_Connectivity();
    int Refine_Triangle(int conn[6], int tri[][3]);
    double Compute_Tri_AspectRatio(double x1, double y1, double x2, double y2, double x3, double y3);
    void Compute_Gauss_Gradient(int IBndy, double *F, double *Fx, double *Fy);
    double Adaptation_Function(int Node1, int Node2, double *F, double *Fx, double *Fy);
    void Compute_Adaptation_Threshold(double *F, double *Fx, double *Fy);
    void Adaptation_Refine_Mark_Edges(double *F, double *Fx, double *Fy);
    void Adaptation_Coarsen_Mark_Nodes(double *F, double *Fx, double *Fy);
    void Generate_Refine_New_Nodes_And_Solution();
    void Update_Boundaries();
    void Update_Triangles();
    void Delete_Mesh_Connectivity();
    void Get_Field_Data();
    void Compress_Coarsen_Nodes();
    void Generate_BowyerWatson_Delaunay_TriMesh();
    void Compute_Derived_FlowField();
    void Tag_FlowField();
private:
    void Init();
};

#endif	/* _TRIANGLEADAPTIVEREFINEMENT_H */

