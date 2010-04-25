/*
 * File:   TriangleAdaptiveRefinement.cpp
 * Author: Ashish Gupta
 *
 * Created on March 15, 2010, 8:06 PM
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string.h>

#include "Utils.h"
#include "MUtils.h"
#include "TriangleAdaptiveRefinement.h"
#include "TriangleMesh.h"
#include "TriangleMeshSLK.h"

// *****************************************************************************
// *****************************************************************************
TriangleAdaptiveRefinement::TriangleAdaptiveRefinement() {
    printf("\n====================================================");
    printf("\n    Mesh Module : Adaptive Refinement               ");
    printf("\n====================================================\n");

    // Initialize the Data
    Init();
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Init() {
    NNode_Refine = 0;
    NNode_Old    = 0;
    NTri_Old     = 0;
    Cell2Node_6  = NULL;
    NodeCoarsen  = NULL;
    NField       = 0;
    FieldTag     = NULL;
    Field        = NULL;
    FieldName    = NULL;
    AFType       = 1;
    Frefine      = 0;
    Fcoarsen     = 0;
    Fcoarsen_Old = 0;
    FuncType     = 0;
    NAdapt       = 0;
    NFieldAdapt  = 0;
    ForceTriMesh = 0;
    NLimitCoarsen= -1;
    NLimitRefine = -1;
    AFrefine     = 0.0;
    AFcoarsen    = 0.0;
    AFavg        = 0.0;
    Sigma        = 0.0;
    Cr           = 1.0;
    Cc           = 1.0;
    P            = 1.0;
}

// *****************************************************************************
// *****************************************************************************
TriangleAdaptiveRefinement::~TriangleAdaptiveRefinement() {
    int i;

    if (Cell2Node_6 != NULL) {
        for (i = 0; i < NTri_Old; i++) {
            if (Cell2Node_6[i] != NULL)
                delete[] Cell2Node_6[i];
        }
        delete[] Cell2Node_6;
    }

    if (NodeCoarsen != NULL)
        delete[] NodeCoarsen;

    if (FieldTag != NULL)
        delete[] FieldTag;
    
    if (NField > 0) {
        if (Field != NULL) {
            for (i = 0; i < NField; i++) {
                if (Field[i] != NULL)
                    delete[] Field[i];
            }
            delete[] Field;
        }
    }

    if (FieldName != NULL) {
        for (i = 0; i < NField; i++) {
            if (FieldName[i] != NULL)
                delete FieldName[i];
        }
        delete FieldName;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Get_Input_Parameters() {
    
    printf("---------------------------------------------------\n");
    std::cout << "Perform Refinement (0: No, 1:Yes)              : ";
    std::cin  >> Frefine;

    std::cout << "Perform Coarsening (0: No, 1:Yes)              : ";
    std::cin  >> Fcoarsen;

    if (Frefine == 1) {
        std::cout << "Input Refinement Coefficient (Cr)              : ";
        std::cin  >> Cr;
        std::cout << "Input Max Refinement Node Creation Limit       : ";
        std::cin  >> NLimitRefine;
    }

    if (Fcoarsen == 1) {
        std::cout << "Input Coarsening Coefficient (Cc)              : ";
        std::cin  >> Cc;
        std::cout << "Input Max Coarsening Node Deletion Limit       : ";
        std::cin  >> NLimitCoarsen;
    }

    if ((Frefine != 1) && (Fcoarsen !=1))
        return;

    std::cout << "Input Adaptation Exponent (P)                  : ";
    std::cin  >> P;

    NAdapt = 0;
    std::cout << "Input Adaptation Cycle (>= 0)                  : ";
    std::cin  >> NAdapt;
    
    std::cout << "Input Adaptation Method (0:Linear, 1:Gradient) : ";
    std::cin  >> AFType;
    if ((AFType < 0) || (AFType > 1))
        error("AdaptiveRefinement: %s\n", "Invalid Adaptation Method");

    FuncType = -1;
    std::cout << "Input Function To Adapt Options                ? " << std::endl;
    std::cout << "Flow Field = 0 " << std::endl;
    std::cout << "Circle     = 1 " << std::endl;
    std::cout << "Cosine     = 2 " << std::endl;
    std::cout << "CosineSine = 3 " << std::endl;
    std::cout << "Function Type                                  : ";
    std::cin  >> FuncType;
    if ((FuncType < 0) || (FuncType > 3))
        error("AdaptiveRefinement: %s\n", "Invalid Analytic Function");

    Tag_FlowField();

    if (Fcoarsen != 1) {
        ForceTriMesh = 0;
        std::cout << "Force Triangular Mesher (0: No, 1: Yes)        : ";
        std::cin  >> ForceTriMesh;
        if ((ForceTriMesh < 0) || (ForceTriMesh > 1))
            error("AdaptiveRefinement: %s\n", "Force Triangular Mesher - Invalid Valid");
    } else
        ForceTriMesh = 1;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Tag_FlowField() {
    int i;

    if (FuncType != 0)
        return;

    // Allocate Memory
    FieldTag = new int[NVariable+3];
#ifdef DEBUG
    if (FieldTag == NULL)
        error("Tag_FlowField: %s\n", "Error Allocating Memory 1");
#endif

    // Select the Adaptaion Function
    std::cout << "Select Flow Field Variables for Adaptation     ? " << std::endl;
    for (i = 0; i < NVariable; i++) {
        FieldTag[i] = 0;
        std::cout << "Field: " << i+1 << " = " << VariableName[i] << " : (0/1) : ";
        std::cin  >> FieldTag[i];
    }

    // Derived Variables
    FieldTag[NVariable+0] = 0;
    std::cout << "Field: " << NVariable+1 << " = Velocity_Magnitude " << ": (0/1) : ";
    std::cin  >> FieldTag[NVariable+0];
    FieldTag[NVariable+1] = 0;
    std::cout << "Field: " << NVariable+2 << " = Pressure " << ": (0/1) : ";
    std::cin  >> FieldTag[NVariable+1];
    FieldTag[NVariable+2] = 0;
    std::cout << "Field: " << NVariable+3 << " = Mach_Number " << ": (0/1) : ";
    std::cin  >> FieldTag[NVariable+2];
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::AdaptiveRefinement() {
    int iAdapt;

    printf("---------------------------------------------------\n");
    if ((Frefine != 1) && (Fcoarsen != 1))
        return;
    
    info("Adaptive Refine: Ini: Nodes %7.0d - Triangles %7.0d", NNode, NTri);
    
    // Save Coarsening Flag: This is done because Coarsening may be deactivated
    Fcoarsen_Old = Fcoarsen;

    for (iAdapt = 0; iAdapt < NAdapt; iAdapt++) {
        Fcoarsen  = Fcoarsen_Old;

        // Coarsening Starts
        if (Fcoarsen == 1)
            Coarsen_Mesh();
        
        // Refinement Starts
        if (Frefine == 1)
            Refine_Mesh();
        
        info("Adaptive Refine: %2.0d - Nodes %7.0d - Triangles %7.0d", iAdapt+1, NNode, NTri);
        printf("---------------------------------------------------\n");
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Refine_Mesh() {
    int iField, NEdge;
    double *RefineTreshold;
    Edge2D *RefineEdge;

    NNode_Old = NNode;
    NTri_Old  = NTri;
    
    // Create Interior and Boundary Node Tagging
    Create_Interior_Boundary_Tag();

    // Create All Possible Connectivities and Sorted
    // Sorting is done to facilate Local Triangulation while Coarsening
    Create_Connectivity(1);

    // Expand Connectivity to Higher Order
    Create_Tri6_Connectivity();

    // Free the Triangle Memory
    if (Tri != NULL) {
        delete[] Tri;
        Tri = NULL;
    }

    // Get the Solution Field
    Get_Field_Data();

    // Get the Number of Data Field
    NFieldAdapt = 0;
    for (iField = 0; iField < NField; iField++) {
        if (FuncType == 0) {
            if (FieldTag[iField] == 0)
                continue;
        }
        NFieldAdapt++;
    }

    if (NFieldAdapt == 0)
        error("Refine_Mesh: %s", "No Data Field Found to Adapt !");

    // Allocate memory to store Refinement Treshold for Each Field
    RefineTreshold = NULL;
    RefineTreshold = new double[NFieldAdapt];
    for (iField = 0; iField < NFieldAdapt; iField++)
        RefineTreshold[iField] = 0.0;
    
    // Get the Number of Edges in Mesh
    NEdge = Get_NumberOfEdges();
    if (NEdge <= 0)
        error("Refine_Mesh: %s", "Anamoly Occured !");

    RefineEdge = NULL;
    RefineEdge = new Edge2D[NEdge];
    
    // Get the Refinement Edge Data
    Refine_Get_Edge_Data(NEdge, RefineEdge, RefineTreshold);

    // Initialize
    NNode_Refine = NNode;
    // Get the Most Potential Edges
    Get_Potential_Refine_Edges(NEdge, RefineEdge, RefineTreshold);
    
    // Generate New Nodes and Interpolate Solution Field
    Generate_Refine_New_Nodes_And_Solution();

    // Update Boundaries with New Segements
    Update_Refine_Boundaries();

    // Create New Triangles if Trimesher conditions is not forced
    if (ForceTriMesh == 0) {
        // Update and Expand the Triangle Connectivities
        Update_Refine_Triangles();
    } else {
        // Call Mesher
        Generate_BowyerWatson_Delaunay_TriMesh();
    }

    // Delete Expanded Connectivity
    Delete_Tri6_Connectivity();

    // Free the Mesh Connectivities as they are changed
    Delete_Mesh_Connectivity();

    Delete_Interior_Boundary_Tag();
    
    // Free the Used Memory
    if (RefineEdge != NULL)
        delete[] RefineEdge;
    if (RefineTreshold != NULL)
        delete[] RefineTreshold;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Coarsen_Mesh() {
    int i, iField, iNode, iTri, NLC, iCNodes, nodeid, itmp;
    double *Fx, *Fy;
    List CandiNodes;
    List FinalCandiNodes;
    
    // Coarsening Starts
    if (Fcoarsen != 1)
        return;
    
    // Get the Solution Field
    Get_Field_Data();
    
    // Allocate Memory to Store Gradients
    Fx = Fy = NULL;
    Fx = new double[NNode];
    Fy = new double[NNode];
#ifdef DEBUG
    if ((Fx == NULL) || (Fy == NULL))
        error("Coaren_Mesh: %s\n", "Error Allocating Memory 1");
#endif
    
    // Create Coarsen Tag Array
    NodeCoarsen = new int [NNode];
#ifdef DEBUG
    if (NodeCoarsen == NULL)
        error("Coaren_Mesh: %s\n", "Error Allocating Memory 2");
#endif

    NLC = NLimitCoarsen;
    iCNodes = 0;
    while (NLC > 0) {
        NNode_Old = NNode;
        NTri_Old  = NTri;

        CandiNodes.Reset(5);
        FinalCandiNodes.Reset(5);
        
        // Create Interior and Boundary Node Tagging
        Create_Interior_Boundary_Tag();

        // Create All Possible Connectivities and Sorted
        // Sorting is done to facilate Local Triangulation while Coarsening
        Create_Connectivity(1);

        // Expand Connectivity to Higher Order
        Create_Tri6_Connectivity();

        // Free the Triangle Memory
        if (Tri != NULL) {
            delete[] Tri;
            Tri = NULL;
        }
        
        // Initialize the Coarsen Node Tag Array
        for (iNode = 0; iNode < NNode; iNode++)
            NodeCoarsen[iNode] = 0;
        
        NFieldAdapt = 0;
        // Function Loop
        for (iField = 0; iField < NField; iField++) {
            if (FuncType == 0) {
                if (FieldTag[iField] == 0)
                    continue;
            }

            NFieldAdapt++;
            // Compute Gradient also on Boundary
            Compute_Gauss_Gradient(1, Field[iField], Fx, Fy);

            // Compute Refinement Threshold
            Compute_Adaptation_Threshold(Field[iField], Fx, Fy);

            // Mark Nodes for Coarsening
            Adaptation_Coarsen_Mark_Nodes(Field[iField], Fx, Fy);
        }
        
        // Check if Coarsening is required
        iCNodes = Finalize_Coarsen_Mark_Nodes();
        // Exit While Loop if No Nodes are Found
        if (iCNodes == 0) {
            // Get the Triangle 3 Connectivity Back
            Tri = new int[NTri][3];
#ifdef DEBUG
            if (Tri == NULL)
                error("Coaren_Mesh: %s\n", "Error Allocating Memory 3");
#endif
            for (iTri = 0; iTri < NTri; iTri++) {
                // Copy Values
                for (i = 0; i < 3; i++)
                    Tri[iTri][i] = Cell2Node_6[iTri][i];
            }

            // Finally Delete Tri6 Connectivity
            Delete_Tri6_Connectivity();
            break;
        }

        // Get the List of Nodes to be Deleted
        for (iNode = 0; iNode < NNode; iNode++) {
            if (NodeCoarsen[iNode] < 0)
                CandiNodes.Add_To_List(iNode);
        }
        
        if (CandiNodes.max != iCNodes)
            error("Coaren_Mesh: %s", "Anomaly Occured!");
        
        // Find the Most Potential Nodes for each Field
        for (iField = 0; iField < NField; iField++) {
            if (FuncType == 0) {
                if (FieldTag[iField] == 0)
                    continue;
            }

            // Compute Gradient also on Boundary
            Compute_Gauss_Gradient(1, Field[iField], Fx, Fy);
            nodeid = -1;
            nodeid = Get_Potential_Coarsen_Nodes(CandiNodes, Field[iField], Fx, Fy);
            if (nodeid != -1)
                FinalCandiNodes.Check_List(nodeid);
        }
        
        // Now Delete the Nodes Potential of Potential Nodes
        iCNodes = 0;
        iCNodes = Compress_Coarsen_Nodes(FinalCandiNodes);
        
        // Get the Triangle 3 Connectivity Back
        Tri =  new int[NTri][3];
#ifdef DEBUG
        if (Tri == NULL)
            error("Coaren_Mesh: %s\n", "Error Allocating Memory 4");
#endif
        for (iTri = 0; iTri < NTri; iTri++) {
            // Copy Values
            for (i = 0; i < 3; i++)
                Tri[iTri][i] = Cell2Node_6[iTri][i];
        }
        
        // Finally Delete Tri6 Connectivity
        Delete_Tri6_Connectivity();

        // Now Delete Connectivity As they are no longer Valid
        itmp = NNode;
        NNode = NNode_Old;
        Delete_Connectivity();
        NNode = itmp;
        
        // Delete Interior and Boundary Node Tagging
        Delete_Interior_Boundary_Tag();

        info("Coaren: Deleted %4d Nodes - %6d NODES %6d TRIANGLE", (NLimitCoarsen-NLC-iCNodes+2), NNode, NTri);
        // Finally Decrement
        NLC = NLC - iCNodes;
    }

    // Delete the Coarsen Node Tag
    if (NodeCoarsen != NULL)
        delete[] NodeCoarsen;
    NodeCoarsen = NULL;

    // Free the Used Memory
    if (Fx != NULL)
        delete[] Fx;
    if (Fy != NULL)
        delete[] Fy;
    
    Delete_Mesh_Connectivity();
    info("Coarsen_Mesh: Deleted Nodes     %d", NLimitCoarsen-NLC);
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Create_Tri6_Connectivity() {
    int i, iTri;

    Cell2Node_6 = new int*[NTri];
#ifdef DEBUG
    if (Cell2Node_6 == NULL)
        error("Create_Tri6_Connectivity: %s\n", "Error Allocating Memory 1");
#endif
    for (iTri = 0; iTri < NTri; iTri++) {
        Cell2Node_6[iTri] = NULL;
        Cell2Node_6[iTri] = new int[6];
#ifdef DEBUG
        if (Cell2Node_6[iTri] == NULL)
            error("Create_Tri6_Connectivity: %s\n", "Error Allocating Memory 2");
#endif
        // Copy Values and Initialize
        for (i = 0; i < 3; i++)
            Cell2Node_6[iTri][i] = Tri[iTri][i];
        for (i = 3; i < 6; i++)
            Cell2Node_6[iTri][i] = -1;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Delete_Tri6_Connectivity() {
    int i;

    if (Cell2Node_6 != NULL) {
        for (i = 0; i < NTri_Old; i++) {
            if (Cell2Node_6[i] != NULL)
                delete[] Cell2Node_6[i];
        }
        delete[] Cell2Node_6;
    }
    
    Cell2Node_6 = NULL;
}

// *****************************************************************************
// *****************************************************************************
double TriangleAdaptiveRefinement::Analytic_Function(double CoordX, double CoordY) {
    double a, b, value = 0.0;
    switch (FuncType) {
        case 1: // Circle Function
            CoordX -= 0.5;
            CoordY -= 0.5;
            a = sqrt(CoordX * CoordX + CoordY * CoordY);
            if (a <= 0.25 && a >= 0.125)
                value = 1.0;
            break;
        case 2: // Cos Function
            a = CoordX * M_PI * 2.0;
            b = cos(a)*0.25 + 0.5;
            if (CoordY > b)
                value = 1.0;
            break;
        case 3: // Cos and Sin Function
            value = (cos(2.0*M_PI*CoordX) + sin(2.0*M_PI*CoordY));
            break;
    }
    return value;
}

// *****************************************************************************
// *****************************************************************************
int TriangleAdaptiveRefinement::Refine_Triangle(int conn[6], int tri[][3]) {
    double asp1, asp2;

    // This routine returns the number of sub-triangles
    // The sub-triangle connectivity will be stored in the "tri" array

    int nt = 0;

    int mask;
    mask = 0;
    if (conn[3] >= 0) mask += 1;
    if (conn[4] >= 0) mask += 2;
    if (conn[5] >= 0) mask += 4;

    switch (mask) {
        case 0: // 000
            tri[nt][0] = conn[0];
            tri[nt][1] = conn[1];
            tri[nt][2] = conn[2];
            nt++;
            break;
        case 1: // 001
            tri[nt][0] = conn[0];
            tri[nt][1] = conn[3];
            tri[nt][2] = conn[2];
            nt++;
            tri[nt][0] = conn[3];
            tri[nt][1] = conn[1];
            tri[nt][2] = conn[2];
            nt++;
            break;
        case 2: // 010
            tri[nt][0] = conn[0];
            tri[nt][1] = conn[1];
            tri[nt][2] = conn[4];
            nt++;
            tri[nt][0] = conn[0];
            tri[nt][1] = conn[4];
            tri[nt][2] = conn[2];
            nt++;
            break;
        case 3: // 011
            tri[nt][0] = conn[3];
            tri[nt][1] = conn[1];
            tri[nt][2] = conn[4];
            nt++;
            asp1 = Compute_Tri_AspectRatio(X[conn[0]], Y[conn[0]], X[conn[3]], Y[conn[3]], X[conn[4]], Y[conn[4]]);
            asp1 = MAX(asp1, Compute_Tri_AspectRatio(X[conn[0]], Y[conn[0]], X[conn[4]], Y[conn[4]], X[conn[2]], Y[conn[2]]));
            asp2 = Compute_Tri_AspectRatio(X[conn[0]], Y[conn[0]], X[conn[3]], Y[conn[3]], X[conn[2]], Y[conn[2]]);
            asp2 = MAX(asp2, Compute_Tri_AspectRatio(X[conn[3]], Y[conn[3]], X[conn[4]], Y[conn[4]], X[conn[2]], Y[conn[2]]));
            if (asp1 < asp2) {
                tri[nt][0] = conn[0];
                tri[nt][1] = conn[3];
                tri[nt][2] = conn[4];
                nt++;
                tri[nt][0] = conn[0];
                tri[nt][1] = conn[4];
                tri[nt][2] = conn[2];
                nt++;
            } else {
                tri[nt][0] = conn[0];
                tri[nt][1] = conn[3];
                tri[nt][2] = conn[2];
                nt++;
                tri[nt][0] = conn[3];
                tri[nt][1] = conn[4];
                tri[nt][2] = conn[2];
                nt++;
            }
            break;
        case 4: // 100
            tri[nt][0] = conn[0];
            tri[nt][1] = conn[1];
            tri[nt][2] = conn[5];
            nt++;
            tri[nt][0] = conn[5];
            tri[nt][1] = conn[1];
            tri[nt][2] = conn[2];
            nt++;
            break;
        case 5: // 101
            tri[nt][0] = conn[0];
            tri[nt][1] = conn[3];
            tri[nt][2] = conn[5];
            nt++;
            asp1 = Compute_Tri_AspectRatio(X[conn[3]], Y[conn[3]], X[conn[1]], Y[conn[1]], X[conn[2]], Y[conn[2]]);
            asp1 = MAX(asp1, Compute_Tri_AspectRatio(X[conn[3]], Y[conn[3]], X[conn[2]], Y[conn[2]], X[conn[5]], Y[conn[5]]));
            asp2 = Compute_Tri_AspectRatio(X[conn[3]], Y[conn[3]], X[conn[1]], Y[conn[1]], X[conn[5]], Y[conn[5]]);
            asp2 = MAX(asp2, Compute_Tri_AspectRatio(X[conn[1]], Y[conn[1]], X[conn[2]], Y[conn[2]], X[conn[5]], Y[conn[5]]));
            if (asp1 < asp2) {
                tri[nt][0] = conn[3];
                tri[nt][1] = conn[1];
                tri[nt][2] = conn[2];
                nt++;
                tri[nt][0] = conn[3];
                tri[nt][1] = conn[2];
                tri[nt][2] = conn[5];
                nt++;
            } else {
                tri[nt][0] = conn[3];
                tri[nt][1] = conn[1];
                tri[nt][2] = conn[5];
                nt++;
                tri[nt][0] = conn[1];
                tri[nt][1] = conn[2];
                tri[nt][2] = conn[5];
                nt++;
            }
            break;
        case 6: // 110
            tri[nt][0] = conn[5];
            tri[nt][1] = conn[4];
            tri[nt][2] = conn[2];
            nt++;
            asp1 = Compute_Tri_AspectRatio(X[conn[0]], Y[conn[0]], X[conn[1]], Y[conn[1]], X[conn[4]], Y[conn[4]]);
            asp1 = MAX(asp1, Compute_Tri_AspectRatio(X[conn[0]], Y[conn[0]], X[conn[4]], Y[conn[4]], X[conn[5]], Y[conn[5]]));
            asp2 = Compute_Tri_AspectRatio(X[conn[0]], Y[conn[0]], X[conn[1]], Y[conn[1]], X[conn[5]], Y[conn[5]]);
            asp2 = MAX(asp2, Compute_Tri_AspectRatio(X[conn[1]], Y[conn[1]], X[conn[4]], Y[conn[4]], X[conn[5]], Y[conn[5]]));
            if (asp1 < asp2) {
                tri[nt][0] = conn[0];
                tri[nt][1] = conn[1];
                tri[nt][2] = conn[4];
                nt++;
                tri[nt][0] = conn[0];
                tri[nt][1] = conn[4];
                tri[nt][2] = conn[5];
                nt++;
            } else {
                tri[nt][0] = conn[0];
                tri[nt][1] = conn[1];
                tri[nt][2] = conn[5];
                nt++;
                tri[nt][0] = conn[1];
                tri[nt][1] = conn[4];
                tri[nt][2] = conn[5];
                nt++;
            }
            break;
        case 7: // 111
            tri[nt][0] = conn[0];
            tri[nt][1] = conn[3];
            tri[nt][2] = conn[5];
            nt++;
            tri[nt][0] = conn[3];
            tri[nt][1] = conn[1];
            tri[nt][2] = conn[4];
            nt++;
            tri[nt][0] = conn[5];
            tri[nt][1] = conn[4];
            tri[nt][2] = conn[2];
            nt++;
            tri[nt][0] = conn[3];
            tri[nt][1] = conn[4];
            tri[nt][2] = conn[5];
            nt++;
            break;
        default:
            error("Refine_Triangle: %s\n", "Found Anomaly 1");
            break;
    }
    return (nt);
}

// *****************************************************************************
// *****************************************************************************
double TriangleAdaptiveRefinement::Compute_Tri_AspectRatio(double x1, double y1, double x2, double y2, double x3, double y3) {
    double a, b, c, s, asp;

    a = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    b = sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2));
    c = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
    s = 0.5 * (a + b + c);
    asp = (a * b * c) / (8.0 * (s - a)*(s - b)*(s - c));
    
    return (asp);
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Compute_Gauss_Gradient(int IBndy, double *F, double *Fx, double *Fy) {
    int iNode, iTri, TriID, n0, n1, n2;
    double nx0, ny0, nx1, ny1, nx2, ny2, area;

    if ((F == NULL) || (Fx == NULL) || (Fy == NULL) || (NNode <= 0))
        return;

    for (iNode = 0; iNode < NNode; iNode++) {
        if (IBndy == 0) {
            // Ignore Boundary Nodes
            if (IBTag[iNode] > -1)
                continue;
        }
        
        area  = 0.0;
        // Initialize Gradients
        Fx[iNode] = Fy[iNode] = 0.0;
        // Loop over all Neigbouring Triangles
        for (iTri = 0; iTri < Node2Cell[iNode]->max; iTri++) {
            TriID = Node2Cell[iNode]->list[iTri];

            n0 = Cell2Node_6[TriID][0];
            n1 = Cell2Node_6[TriID][1];
            n2 = Cell2Node_6[TriID][2];

            // Compute the Normals
            nx0 = (Y[n2] - Y[n1]);
            ny0 = -(X[n2] - X[n1]);
            nx1 = (Y[n0] - Y[n2]);
            ny1 = -(X[n0] - X[n2]);
            nx2 = (Y[n1] - Y[n0]);
            ny2 = -(X[n1] - X[n0]);

            // Compute Area of the Triangle
            area += fabs(0.5*(nx1*ny2 - nx2*ny1));

            // Compute the Gradient
            Fx[iNode] -= 0.5*(F[n0]*nx0 + F[n1]*nx1 + F[n2]*nx2);
            Fy[iNode] -= 0.5*(F[n0]*ny0 + F[n1]*ny1 + F[n2]*ny2);
        }
        Fx[iNode] /= area;
        Fy[iNode] /= area;
    }
}

// *****************************************************************************
// *****************************************************************************
double TriangleAdaptiveRefinement::Adaptation_Function(int Node1, int Node2,
        double *F, double *Fx, double *Fy) {
    double AF = 0.0;
    double length;
    double Ex, Ey;

    if ((F == NULL) || (Fx == NULL) || (Fy == NULL) || (NNode <= 0))
        return 0.0;

    // Edge and Length
    Ex = X[Node2]-X[Node1];
    Ey = Y[Node2]-Y[Node1];
    length = sqrt(Ex*Ex + Ey*Ey);
    Ex /= length;
    Ey /= length;
    
    switch (AFType) {
        case 0:
            AF = fabs(F[Node2] - F[Node1])*pow(length, P-1.0);
            break;
        case 1:
            AF  = fabs(Fx[Node1]*Ex + Fy[Node1]*Ey);
            AF += fabs(Fx[Node2]*Ex + Fy[Node2]*Ey);
            AF = 0.5*AF*pow(length, P);
            break;
    }
    
    return AF;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Compute_Adaptation_Threshold(double *F, double *Fx, double *Fy) {
    int i, iNode, nodeid, nEdge;
    double var;

    nEdge = 0;
    AFrefine  = 0.0;
    AFcoarsen = 0.0;
    AFavg     = 0.0;
    Sigma     = 0.0;

    // Compute the No of Edges and Get the Average Adaptive Function
    for (iNode = 0; iNode < NNode; iNode++) {
        for (i = 0; i < Node2Node[iNode]->max; i++) {
            nodeid = Node2Node[iNode]->list[i];
            if (iNode < nodeid) {
                nEdge++;
                // Compute the Adaptation Function
                AFavg += Adaptation_Function(iNode, nodeid, F, Fx, Fy);
            }
        }
    }
    AFavg /= ((double)nEdge);

    // Compute the Sigma or Standard Deviation
    for (iNode = 0; iNode < NNode; iNode++) {
        for (i = 0; i < Node2Node[iNode]->max; i++) {
            nodeid = Node2Node[iNode]->list[i];
            if (iNode < nodeid) {
                var = Adaptation_Function(iNode, nodeid, F, Fx, Fy) - AFavg;
                Sigma += var*var;
            }
        }
    }
    Sigma = sqrt(Sigma/((double)nEdge));
    AFrefine = AFavg  + Cr*Sigma;
    AFcoarsen = AFavg - Cc*Sigma;
}

// *****************************************************************************
// *****************************************************************************
int TriangleAdaptiveRefinement::Get_NumberOfEdges() {
    int i, iNode, nodeid, count;

    // Edge Loop
    count = 0;
    for (iNode = 0; iNode < NNode; iNode++) {
        for (i = 0; i < Node2Node[iNode]->max; i++) {
            nodeid = Node2Node[iNode]->list[i];
            if (iNode < nodeid)
                count++;
        }
    }

    return count;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Refine_Get_Edge_Data(int NEdge, Edge2D *RefineEdge, double *RefineTreshold) {
    int i, j, k, iNode, iField, iEdge, nodeid, cellid, icell, iFA;
    double *Fx, *Fy;
    
    if ((RefineEdge == NULL) || (NFieldAdapt <= 0))
        return;

    // Edge Loop to Get the Connectivity Data
    iEdge = 0;
    for (iNode = 0; iNode < NNode; iNode++) {
        for (i = 0; i < Node2Node[iNode]->max; i++) {
            nodeid = Node2Node[iNode]->list[i];
            if (iNode < nodeid) {
                // Get the Edge Nodes
                RefineEdge[iEdge].node1 = iNode;
                RefineEdge[iEdge].node2 = nodeid;

                RefineEdge[iEdge].Flag = 0;

                if (RefineEdge[iEdge].nData == 0) {
                    RefineEdge[iEdge].nData = NFieldAdapt;
                    RefineEdge[iEdge].Data = new double[NFieldAdapt];
                    for (j = 0; j < RefineEdge[iEdge].nData; j++)
                        RefineEdge[iEdge].Data[j] = 0.0;
                }

                // Get the Cells IDs
                // Loop Over Triangles Around iNode
                icell = 0;
                for (j = 0; j < Node2Cell[iNode]->max; j++) {
                    cellid = Node2Cell[iNode]->list[j];
                    // Check if cell contains the other node of edge
                    for (k = 0; k < 3; k++) {
                        if (nodeid == Cell2Node_6[cellid][k]) {
                            if (icell == 0) {
                                RefineEdge[iEdge].cell1 = cellid;
                                icell++;
                            } else
                                RefineEdge[iEdge].cell2 = cellid;
                            // exit loop
                            k = 4;
                        }
                    } // k - loop
                } // iNode Surround Triangle Loop

                if (icell > 1)
                   error("Refine_Get_Edge_Data: %s", "Anomaly Occured 1!");

                // Finally Increase the Edge Count
                iEdge++;
            }
        }
    }

    // Allocate Memory to Store Gradients
    Fx = Fy = NULL;
    Fx = new double[NNode];
    Fy = new double[NNode];
#ifdef DEBUG
    if ((Fx == NULL) || (Fy == NULL))
        error("Refine_Get_Edge_Data: %s\n", "Error Allocating Memory 1");
#endif

    // Function Loop for Refinement
    iFA = 0;
    for (iField = 0; iField < NField; iField++) {
        if (FuncType == 0) {
            if (FieldTag[iField] == 0)
                continue;
        }

        // Compute Gradient also on Boundary
        Compute_Gauss_Gradient(1, Field[iField], Fx, Fy);

        // Compute Refinement Threshold
        Compute_Adaptation_Threshold(Field[iField], Fx, Fy);
        RefineTreshold[iFA] = AFrefine;
        
        // Edge Loop
        for (iEdge = 0; iEdge < NEdge; iEdge++) {
            // Get the Adaptation Function Value
            iNode  = RefineEdge[iEdge].node1;
            nodeid = RefineEdge[iEdge].node2;
            RefineEdge[iEdge].Data[iFA] = Adaptation_Function(iNode, nodeid, Field[iField], Fx, Fy);
            // Tag Edge if Only if Refine Condition Exceeds
            if (RefineEdge[iEdge].Data[iFA] > AFrefine)
                RefineEdge[iEdge].Flag++;
        }

        // Increment the Field Count
        iFA++;
    }

    // Free the Used Memory
    if (Fx != NULL)
        delete[] Fx;
    if (Fy != NULL)
        delete[] Fy;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Adaptation_Refine_Mark_Edges(int NEdge, Edge2D *RefineEdge, List& CandiEdge) {
    int i, j, k, cellid, edgeid, redge, iEdge, n1, n2;
    int NNode_New;

    // Check if Refinement Conditions
    if ((NFieldAdapt <= 0) && (RefineEdge != NULL) && (NEdge <= 0))
        return;
    
    NNode_New = NNode_Refine;
    // Edge Loop: Test, Mark Each Edge and Notify Neighbour
    for (iEdge = 0; iEdge < CandiEdge.max; iEdge++) {
        redge = CandiEdge.list[iEdge];

        // Get the node and cell ids
        n1 = RefineEdge[redge].node1;
        n2 = RefineEdge[redge].node2;

        for (i = 0; i < 2; i++) {
            if (i == 0)
                cellid = RefineEdge[redge].cell1;
            else
                cellid = RefineEdge[redge].cell2;

            // For Boundary Edge
            if (cellid == -1)
                continue;
            
            // Insert the new Node in cell1 edge
            // Find the Edge ID of the Triangle with n1 - n2
            for (j = 0; j < 3; j++) {
                if (n2 == Cell2Node_6[cellid][j]) {
                    edgeid = -1;
                    for (k = 0; k < 3; k++) {
                        if (n1 == Cell2Node_6[cellid][k]) {
                            // Mark the Triangle Edge For Refinement
                            switch (k + j) {
                                case 1: // Edge1
                                    edgeid = 3;
                                    break;
                                case 3: // Edge2
                                    edgeid = 4;
                                    break;
                                case 2: // Edge3
                                    edgeid = 5;
                                    break;
                                default:
                                    error("Adaptation_Refine_Mark_Edges: %s\n", "Found Anomaly 1");
                            }
                            // Assign the new node number
                            Cell2Node_6[cellid][edgeid] = NNode_New;
                            k = 4; // exit loop
                        }
                    } // k - loop
                    j = 4; // exit loop
                }
            } // j - loop
        }

        // Increment the Node
        NNode_New++;
    }
    
    NNode_Refine = NNode_New;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Get_Potential_Refine_Edges(int NEdge, Edge2D *RefineEdge, double *RefineTreshold) {
    int iEdge, iREdge, NLR, edgeid, range;
    int *LGEdgeID;
    double *AdaptVar;
    List PotEdge, FinalEdge;

    // Check if Refinement Conditions
    if ((NFieldAdapt <= 0) && (RefineEdge != NULL))
        return;

    NLR    = NLimitRefine;
    iREdge = 0;
    while (NLR > 0) {
        // Get the Most Potential Edges
        for (iEdge = 0; iEdge < NEdge; iEdge++) {
            if (RefineEdge[iEdge].Flag > 0)
                PotEdge.Add_To_List(iEdge);
            if (RefineEdge[iEdge].Flag == NFieldAdapt)
                FinalEdge.Add_To_List(iEdge);
        }
        
        if (NFieldAdapt == 1) {
            LGEdgeID = NULL;
            LGEdgeID = new int[FinalEdge.max];
            AdaptVar = NULL;
            AdaptVar = new double[FinalEdge.max];
#ifdef DEBUG
            if ((AdaptVar == NULL) || (LGEdgeID == NULL))
                error("Get_Potential_Refine_Edges: %s\n", "Error Allocating Memory 1");
#endif
            for (iEdge = 0; iEdge < FinalEdge.max; iEdge++) {
                LGEdgeID[iEdge] = iEdge;
                edgeid = FinalEdge.list[iEdge];
                AdaptVar[iEdge] = RefineEdge[edgeid].Data[0];
            }

            // Now Sort the List of Edges to Get the Maximum to Minimum Value Distribution of Adaptation Function
            // All Values of Adaptation Function function will be greater then Adaptation Threshold
            sort_high2low(FinalEdge.max, LGEdgeID, AdaptVar);
            
            // Resize to get the Final Edges
            range = MIN(NLR, FinalEdge.max);
            PotEdge.Reset(range);
            for (iEdge = 0; iEdge < range; iEdge++) {
                edgeid = LGEdgeID[iEdge];
                PotEdge.Add_To_List(FinalEdge.list[edgeid]);
            }
            
            // Free Memory
            if (LGEdgeID != NULL)
                delete[] LGEdgeID;
            if (AdaptVar != NULL)
                delete[] AdaptVar;
        } else {
            error("Get_Potential_Refine_Edges: %s", "Multiple Field Not implemented - TODO");
        }

        // Edge Loop: Test, Mark Each Edge and Notify Neighbour
        Adaptation_Refine_Mark_Edges(NEdge, RefineEdge, PotEdge);
        
        NLR = NLR - PotEdge.max;
        if (NFieldAdapt == 1)
            break;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Adaptation_Refine_Mark_Edges_Old(double *F, double *Fx, double *Fy) {
    int i, j, k, l, iNode, nodeid, cellid, edgeid;
    int NNode_New, rflag;
    double var;

    // Check if Refinement is Asked
    if (Frefine == 1) {
        NNode_New = NNode_Refine;
        // Edge Loop: Test, Mark Each Edge and Notify Neighbour
        for (iNode = 0; iNode < NNode; iNode++) {
            for (i = 0; i < Node2Node[iNode]->max; i++) {
                nodeid = Node2Node[iNode]->list[i];
                if (iNode < nodeid) {
                    var = Adaptation_Function(iNode, nodeid, F, Fx, Fy);
                    // Refine Only Condition Exceeds
                    if (var > AFrefine) {
                        rflag = 0; // 1: Not Refined, 0: Already Refined
                        // Loop Over Triangles Around iNode
                        for (j = 0; j < Node2Cell[iNode]->max; j++) {
                            cellid = Node2Cell[iNode]->list[j];
                            // Check if cell contains the other node of edge
                            for (k = 0; k < 3; k++) {
                                if (nodeid == Cell2Node_6[cellid][k]) {
                                    // Find the Edge ID of the Triangle
                                    edgeid = -1;
                                    for (l = 0; l < 3; l++) {
                                        if (iNode == Cell2Node_6[cellid][l]) {
                                            // Mark the Triangle Edge For Refinement
                                            switch (l+k) {
                                                case 1: // Edge1
                                                    edgeid = 3;
                                                    break;
                                                case 3: // Edge2
                                                    edgeid = 4;
                                                    break;
                                                case 2: // Edge3
                                                    edgeid = 5;
                                                    break;
                                                default:
                                                    error("Adaptation_Mark_Edges: %s\n", "Found Anomaly 1");
                                            }
                                            // This condition is set to avoid if multiple time setting of edge to get refined
                                            // and evently counting wrong no of node
                                            // This will happen if for multi function senario: example: pressure + temperature
                                            if (Cell2Node_6[cellid][edgeid] == -1) {
                                                Cell2Node_6[cellid][edgeid] = NNode_New;
                                                rflag = 1;
                                            }
                                            l = 4; // exit loop
                                        }
                                    } // l - loop
                                    k = 4; // exit loop
                                }
                            } // k - loop
                        } // iNode Surround Triangle Loop
                        // Update New NNodes and Triangles
                        if (rflag == 1) {
                            NNode_New++;
                        }
                    }
                }
            }
        }
        NNode_Refine = NNode_New;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Generate_Refine_New_Nodes_And_Solution() {
    int i, iNode, iTri, iField, node1, node2;
    double *X_new, *Y_new, *Field_new;

    if ((NNode_Refine <= NNode) && (NNode_Refine != 0))
        error("Generate_New_Nodes_And_Solution: %s\n", "Found Anomaly 1");

    // Initialize
    node1 = -1;
    node2 = -1;
    
    // Allocate Memory for coordinates
    X_new = NULL;
    Y_new = NULL;
    X_new = new double[NNode_Refine];
    Y_new = new double[NNode_Refine];
#ifdef DEBUG
    if ((X_new == NULL) || (Y_new == NULL))
        error("Generate_Refine_New_Nodes_And_Solution: %s\n", "Error Allocating Memory 1");
#endif
    
    // Copy the Old Coordinates
    for (iNode = 0; iNode < NNode; iNode++) {
        X_new[iNode] = X[iNode];
        Y_new[iNode] = Y[iNode];
    }

    // Now Get the Value for Newly Created Nodes Coordinates
    for (iTri = 0; iTri < NTri; iTri++) {
        for (i = 3; i < 6; i++) {
            // Check if Edge is set for refinement
            if (Cell2Node_6[iTri][i] < 0)
                continue;
            switch (i) {
                case 3:
                    node1 = Cell2Node_6[iTri][0];
                    node2 = Cell2Node_6[iTri][1];
                    break;
                case 4:
                    node1 = Cell2Node_6[iTri][1];
                    node2 = Cell2Node_6[iTri][2];
                    break;
                case 5:
                    node1 = Cell2Node_6[iTri][2];
                    node2 = Cell2Node_6[iTri][0];
                    break;
            }
            X_new[Cell2Node_6[iTri][i]] = 0.5*(X[node1] + X[node2]);
            Y_new[Cell2Node_6[iTri][i]] = 0.5*(Y[node1] + Y[node2]);
        }
    }

    // Free the memory for Old Coordinates
    if (X != NULL)
        delete[] X;
    if (Y != NULL)
        delete[] Y;
    
    // Replace the Pointers for with New Nodes
    X = X_new;
    Y = Y_new;
    X_new = NULL;
    Y_new = NULL;

    // Now Interpolate the New Solution Field
    Field_new = NULL;
    for (iField = 0; iField < NField; iField++) {
        Field_new = new double[NNode_Refine];
#ifdef DEBUG
        if (Field_new == NULL)
            error("Generate_Refine_New_Nodes_And_Solution: %s\n", "Error Allocating Memory 2");
#endif
        // Copy the Field Values
        for (iNode = 0; iNode < NNode; iNode++)
            Field_new[iNode] = Field[iField][iNode];
        
        for (iTri = 0; iTri < NTri; iTri++) {
            for (i = 3; i < 6; i++) {
                // Check if Edge is set for refinement
                if (Cell2Node_6[iTri][i] < 0)
                    continue;
                switch (i) {
                    case 3:
                        node1 = Cell2Node_6[iTri][0];
                        node2 = Cell2Node_6[iTri][1];
                        break;
                    case 4:
                        node1 = Cell2Node_6[iTri][1];
                        node2 = Cell2Node_6[iTri][2];
                        break;
                    case 5:
                        node1 = Cell2Node_6[iTri][2];
                        node2 = Cell2Node_6[iTri][0];
                        break;
                }
                Field_new[Cell2Node_6[iTri][i]] = 0.5*(Field[iField][node1] + Field[iField][node2]);
            }
        }
        
        // Free the memory Field Memory
        if (Field[iField] != NULL)
            delete[] Field[iField];

        Field[iField] = Field_new;
        Field_new = NULL;
    }

    // Update the Number of Nodes
    NNode = NNode_Refine;
    info("Refinement:   Created Nodes:  %5d Total %5d", NNode-NNode_Old, NNode);
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Update_Refine_Boundaries() {
    int i, j, k, ibnd, ibs, node1, node2, cellid, edgeid;
    int nbs_new, nodeid;

    // Initialize
    edgeid = -1;
    
    for (ibnd = 0; ibnd < NBoundary; ibnd++) {
        nbs_new = NBoundarySegments[ibnd];
        for (ibs = 0; ibs < NBoundarySegments[ibnd]; ibs++) {
            node1 = BoundarySegments[ibnd][ibs][0];
            node2 = BoundarySegments[ibnd][ibs][1];
            // Get the Boundary Triangle
            // Loop Over Triangles Around node1
            for (i = 0; i < Node2Cell[node1]->max; i++) {
                cellid = Node2Cell[node1]->list[i];
                // Check if cell contains the other node of edge
                for (j = 0; j < 3; j++) {
                    if (node2 == Cell2Node_6[cellid][j]) {
                        // Find the Edge ID of the Triangle
                        for (k = 0; k < 3; k++) {
                            if (node1 == Cell2Node_6[cellid][k]) {
                                // Mark the Triangle Edge For Refinement
                                switch (j+k) {
                                    case 1: // Edge1
                                        edgeid = 3;
                                        break;
                                    case 3: // Edge2
                                        edgeid = 4;
                                        break;
                                    case 2: // Edge3
                                        edgeid = 5;
                                        break;
                                    default:
                                        error("Update_Boundaries: %s\n", "Found Anomaly 1");
                                }
                                // Check if Boundary Segment is Refined
                                if (Cell2Node_6[cellid][edgeid] != -1) {
                                    nodeid = Cell2Node_6[cellid][edgeid];
                                    // Expand Boundary Segment array for new edge
                                    nbs_new++;
                                    BoundarySegments[ibnd] = (int **) realloc((void *)BoundarySegments[ibnd], nbs_new*sizeof(int*));
                                    BoundarySegments[ibnd][nbs_new - 1] = (int *)malloc(2*sizeof(int));
                                    // Reconstruct Boundary Segment Edge with first half, and add 2nd half of new edge to end
                                    BoundarySegments[ibnd][ibs][0] = node1;
                                    BoundarySegments[ibnd][ibs][1] = nodeid;
                                    BoundarySegments[ibnd][nbs_new - 1][0] = nodeid;
                                    BoundarySegments[ibnd][nbs_new - 1][1] = node2;
                                }
                                k = 4; // exit loop
                            }
                        }
                        j = 4; // exit loop
                    }
                }
            }
        } // BoundarySegment Loop
        NBoundarySegments[ibnd] = nbs_new;
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Update_Refine_Triangles() {
    int i, j, iTri, nsubtri, count, NTri_New;
    int (*subtri)[3];
    
    // Free the Old Triangle Memory if exits
    if (Tri != NULL) {
        delete[] Tri;
        Tri = NULL;
    }

    // Count the number of new triangles
    count = 0;
    for (iTri = 0; iTri < NTri; iTri++) {
        for (i = 3; i < 6; i++)
            if (Cell2Node_6[iTri][i] > -1)
                count++;
    }

    NTri_New = NTri+count;
    
    // Allocate Memory to store expanded triangles
    Tri = new int[NTri_New][3];
#ifdef DEBUG
    if ((Tri == NULL))
        error("Update_Triangles: %s\n", "Error Allocating Memory 1");
#endif

    // Allocate Memory to hold subtriangles : 4 for Maximum division allowed
    subtri = new int[4][3];
    count = NTri;
    for (iTri = 0; iTri < NTri; iTri++) {
        nsubtri = 0;
        for (i = 3; i < 6; i++) {
            if (Cell2Node_6[iTri][i] != -1)
                nsubtri++;
        }
        // Check if Triangle is subdivided
        if (nsubtri == 0) {
            // Copy the Old Connectivity
            for (i = 0; i < 3; i++)
                Tri[iTri][i] = Cell2Node_6[iTri][i];
            continue;
        }
        
        // Get the SubTriangle connectivities
        if (Refine_Triangle(Cell2Node_6[iTri], subtri) != (nsubtri+1))
           error("Update_Triangles: %s\n", "Found Anomaly 1");
        
        // Copy the New Subdivided Triangle Connectivity
        for (i = 0; i < 3; i++)
            Tri[iTri][i] = subtri[0][i];
        for (i = 1; i < (nsubtri+1); i++) {
            for (j = 0; j < 3; j++)
                Tri[count][j] = subtri[i][j];
            count++;
        }
    }
    
    // Check if Count Matched Total New Triangles
    if (count == NTri_New)
        NTri = NTri_New;
    else
        error("Update_Triangles: %s %d\n", "Found Anomaly 2", count);

    // Free the resources Used
    if (subtri != NULL)
        delete[] subtri;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Delete_Mesh_Connectivity() {
    int i;
    
    if (FuncType != 0) {
        if (NField > 0) {
            if (Field != NULL) {
                for (i = 0; i < NField; i++) {
                    if (Field[i] != NULL)
                        delete[] Field[i];
                }
                delete[] Field;
            }
        }
    }

    if (FuncType == 0) {
        if (NField > NVariable) {
            for (i = 0; i < NVariable; i++) {
                Variable[i] = Field[i];
                Field[i] = NULL;
            }
            for (i = NVariable; i < NField; i++) {
                if (Field[i] != NULL)
                    delete[] Field[i];
            }
            delete[] Field;
        }
        
        if (FieldName != NULL) {
            for (i = 0; i < NField; i++) {
                if (FieldName[i] != NULL)
                    delete FieldName[i];
            }
            delete FieldName;
        }
    }
    
    if (Node2Cell != NULL) {
        for (i = 0; i < NNode_Old; i++) {
            if (Node2Cell[i] != NULL)
                delete Node2Cell[i];
        }
        delete Node2Cell;
    }

    if (Cell2Cell != NULL) {
        for (i = 0; i < NTri_Old; i++) {
            if (Cell2Cell[i] != NULL)
                delete Cell2Cell[i];
        }
        delete Cell2Cell;
    }

    if (Node2Node != NULL) {
        for (i = 0; i < NNode_Old; i++) {
            if (Node2Node[i] != NULL)
                delete Node2Node[i];
        }
        delete Node2Node;
    }

    Field = NULL;
    FieldName = NULL;
    Node2Cell = NULL;
    Cell2Cell = NULL;
    Node2Node = NULL;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Get_Field_Data() {
    int iField, iNode;

    // Check if Solution Field is Asked for Adaptation
    if (FuncType == 0) {
        if (NVariable <= 0)
            error("AdaptiveRefinement: %s\n", "No Field Variables Found");
        Compute_Derived_FlowField();
    } else {
        // Allocate the Memory for Functions
        NField = 1;
        Field = new double*[NField];
#ifdef DEBUG
        if (Field == NULL)
            error("Get_Field_Data: %s\n", "Error Allocating Memory 1");
#endif
        for (iField = 0; iField < NField; iField++) {
            Field[iField] = NULL;
            Field[iField] = new double[NNode];
#ifdef DEBUG
            if (Field[iField] == NULL)
                error("Get_Field_Data: %s\n", "Error Allocating Memory 2");
#endif
        }

        // Compute the values of Field Function
        for (iField = 0; iField < NField; iField++) {
            for (iNode = 0; iNode < NNode; iNode++) {
                Field[iField][iNode] = Analytic_Function(X[iNode], Y[iNode]);
            }
        }
    }
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Adaptation_Coarsen_Mark_Nodes(double* F, double* Fx, double* Fy) {
    int i, iNode, nodeid, cflag;
    double var;

    // Check if Coarsen is Asked
    if (Fcoarsen == 1) {
        // Edge Loop: Test, Mark Each Edge and Notify Neighbour
        for (iNode = 0; iNode < NNode; iNode++) {
            // Check if Node is Boundary
            if (IBTag[iNode] != -1)
                continue;

            cflag = 0;
            for (i = 0; i < Node2Node[iNode]->max; i++) {
                nodeid = Node2Node[iNode]->list[i];
                var = Adaptation_Function(iNode, nodeid, F, Fx, Fy);
                // Coarsen if Only Condition Below Limits
                if (var < AFcoarsen)
                   cflag++;
            }
            // Mark Node for Deletion only if all connected edge fails
            if (cflag == Node2Node[iNode]->max)
                NodeCoarsen[iNode] += -1;
        }
    }
}

// *****************************************************************************
// *****************************************************************************
int TriangleAdaptiveRefinement::Finalize_Coarsen_Mark_Nodes() {
    int iNode, count;

    count = 0;
    // Check if Coarsen is Asked
    if (Fcoarsen == 1) {
        // Edge Loop: Test, Mark Each Edge and Notify Neighbour
        for (iNode = 0; iNode < NNode; iNode++) {
            // Check if Node is Boundary
            if (IBTag[iNode] != -1)
                continue;

            // Mark Node for Deletion only if all Adaptive Field Fails
            if (NodeCoarsen[iNode] == -NFieldAdapt) {
                NodeCoarsen[iNode] = -1;
                count++;
            } else
                NodeCoarsen[iNode] = 0;
        }
        // Check if Coarsening is Required
        // This will avoid call to Tri Mesher
        if (count == 0) {
            Fcoarsen = 0;
            // Delete the Coarsen Node Tag
            if (NodeCoarsen != NULL)
                delete[] NodeCoarsen;
            NodeCoarsen = NULL;
        }
    }

    // Return Final Number of Failed Nodes
    return count;
}

// *****************************************************************************
// *****************************************************************************
int TriangleAdaptiveRefinement::Get_Potential_Coarsen_Nodes(List& CandiNodes, double *F, double *Fx, double *Fy) {
    int i, iNode, n1, n2;
    int *LGNodeID;
    double *AdaptVar;
    double var;

    if (CandiNodes.max <= 0)
        return -1;
    
    // Allocate Memory to Store Min Value
    LGNodeID = NULL;
    LGNodeID = new int[CandiNodes.max];
    AdaptVar = NULL;
    AdaptVar = new double[CandiNodes.max];
#ifdef DEBUG
    if ((AdaptVar == NULL) || (LGNodeID == NULL))
        error("Get_Potential_Coarsen_Nodes: %s\n", "Error Allocating Memory 1");
#endif

    // Loop Over the flagged nodes for deletion and Get the Max Adaptation Function
    // Value out of Edges surrounds the Node
    for (iNode = 0; iNode < CandiNodes.max; iNode++) {
        LGNodeID[iNode] = iNode;
        AdaptVar[iNode] = 0.0;
        n1 = CandiNodes.list[iNode];
        for (i = 0; i < Node2Node[n1]->max; i++) {
            n2 = Node2Node[n1]->list[i];
            var = Adaptation_Function(n1, n2, F, Fx, Fy);
            AdaptVar[iNode] = MAX(AdaptVar[iNode], var);
        }
    }
    
    // Now Sort the List of Nodes to Get the Minimum of Maximum Value of Adaptation Function
    // All Values of Adaptation Function function will be less then Adaptation Threshold
    sort_low2high(CandiNodes.max, LGNodeID, AdaptVar);

    // Get the Most Potential Node ID
    n1 = CandiNodes.list[LGNodeID[0]];
    
    // Free Memory
    if (LGNodeID != NULL)
        delete[] LGNodeID;
    if (AdaptVar != NULL)
        delete[] AdaptVar;

    // Get the NodeID which will be deleted
    return n1;
}

// *****************************************************************************
// *****************************************************************************
int TriangleAdaptiveRefinement::Compress_Coarsen_Nodes(List& FinalCandiNodes) {
    int i, j, ibnd, iNode, iField, count, n1, n2, nodeid;
    double *X_New, *Y_New, *Field_New;
    double *Lx, *Ly;
    int Ltdim, Lnt, Lnbnd, Lnnode;
    int ibs, itri, triid;
    int *Lnbs, ***Lbs, *LGNodeID;
    int (*Ltri)[3];
    TriangleMeshSLK Mesher;
    
    if (FinalCandiNodes.max <= 0)
        return 0;

    // Find the List of Nodes which are not Connected to each Other
    // This Operation will result minimum of one node to be deleted
    if (FinalCandiNodes.max > 1) {
        for (iNode = 0; iNode < FinalCandiNodes.max; iNode++) {
            n1 = FinalCandiNodes.list[iNode];
            for (i = 0; i < Node2Node[n1]->max; i++) {
                n2 = Node2Node[n1]->list[i];
                if (FinalCandiNodes.Is_In_List(n2))
                    FinalCandiNodes.Delete_From_List(n2);
            }
        }
    }

    // Now Start Deleting Nodes
    Lnbnd = 1;
    Lnbs  = new int[Lnbnd];
    Lbs   = (int ***) malloc(Lnbnd * sizeof (int**));
    for (i = 0; i < FinalCandiNodes.max; i++) {
        nodeid = FinalCandiNodes.list[i];
        Lnnode  = Node2Node[nodeid]->max;

        // Get the Coordinates and Local to Global Node Index
        LGNodeID = NULL;
        LGNodeID = new int[Lnnode];
        Lx = NULL;
        Ly = NULL;
        Lx = new double[Lnnode];
        Ly = new double[Lnnode];
#ifdef DEBUG
        if ((LGNodeID == NULL) || (Lx == NULL) || (Ly == NULL))
            error("Compress_Coarsen_Nodes: %s\n", "Error Allocating Memory 1");
#endif
        for (j = 0; j < Lnnode; j++) {
            LGNodeID[j] = Node2Node[nodeid]->list[j];
            Lx[j] = X[Node2Node[nodeid]->list[j]];
            Ly[j] = Y[Node2Node[nodeid]->list[j]];
        }

        // Get the Boundary Segments with local node numbering
        Lnbs[0] = Lnnode;
        Lbs[0]  = (int**) malloc(Lnbs[0] * sizeof (int*));
        for (ibs = 0; ibs < Lnbs[0]-1; ibs++) {
            Lbs[0][ibs] = (int*) malloc(2 * sizeof (int));
            Lbs[0][ibs][0] = ibs;
            Lbs[0][ibs][1] = ibs+1;
        }
        Lbs[0][Lnbs[0]-1] = (int*) malloc(2 * sizeof (int));
        Lbs[0][Lnbs[0]-1][0] = Lnbs[0]-1;
        Lbs[0][Lnbs[0]-1][1] = 0;
        
        // initial guess for number of triangles
        Ltdim = Lnnode * 3;
        Ltri  = new int[Ltdim][3];
#ifdef DEBUG
        if (Ltri == NULL)
            error("Compress_Coarsen_Nodes: %s\n", "Error Allocating Memory 2");
#endif

        // iterative loop to allocate space for triangle indices
        Lnt = -1;
        while (Lnt < 0) {
            //Lnt  = Mesher.trimesh(Lnnode, Ltdim, Lnbnd, Lnbs, Lbs, Lx, Ly, Ltri);
            Lnt = trimesh(Lnnode, Ltdim, Lnbnd, Lnbs, Lbs, Lx, Ly, Ltri);

            if (Lnt < 0) {
                info("Expanding New Local Triangle Dimension from %d to %d", Ltdim, Ltdim + Lnnode);

                // Delete current memory for triangle
                if (Ltri != NULL)
                    delete[] Ltri;

                // Increment triangle dimension by number of nodes
                Ltdim += Lnnode;
                // Allocate new space for triangles
                Ltri = new int[Ltdim][3];
#ifdef DEBUG
                if (Ltri == NULL)
                    error("Compress_Coarsen_Nodes: %s\n", "Error Allocating Memory 3");
#endif
            }
        }

        // Now Modify the Triangle to Node Connectivity
        if (Lnt >= Node2Cell[FinalCandiNodes.list[i]]->max)
            error("Compress_Coarsen_Nodes: %s", "Anamoly Occured !");
        
        for (itri = 0; itri < Node2Cell[FinalCandiNodes.list[i]]->max; itri++) {
            triid = Node2Cell[FinalCandiNodes.list[i]]->list[itri];
            if (itri < Lnt) {
                Cell2Node_6[triid][0] = LGNodeID[Ltri[itri][0]];
                Cell2Node_6[triid][1] = LGNodeID[Ltri[itri][1]];
                Cell2Node_6[triid][2] = LGNodeID[Ltri[itri][2]];
            } else {
                Cell2Node_6[triid][0] = -1;
                Cell2Node_6[triid][1] = -1;
                Cell2Node_6[triid][2] = -1;
            }
        }

        // Free the Memory
        if (Ltri != NULL)
            delete[] Ltri;
        if (Lx != NULL)
            delete Lx;
        if (Ly != NULL)
            delete Ly;
        if (LGNodeID != NULL)
            delete LGNodeID;
        if (Lbs[0] != NULL) {
            for (int ibs = 0; ibs < Lnbs[0]; ibs++)
                free(Lbs[0][ibs]);
            free(Lbs[0]);
        }
    }
    if (Lnbs != NULL)
        delete[] Lnbs;
    if (Lbs != NULL)
        free(Lbs);

    // Get the New Node IDs
    count = 0;
    for (iNode = 0; iNode < NNode; iNode++) {
        if (FinalCandiNodes.Is_In_List(iNode))
            NodeCoarsen[iNode] = -1;
        else {
            NodeCoarsen[iNode] = count;
            count++;
        }
    }

    // Now Compress the Coordinates X, Y
    X_New = NULL;
    Y_New = NULL;
    X_New = new double[count];
    Y_New = new double[count];
#ifdef DEBUG
    if ((X_New == NULL) || (Y_New == NULL))
        error("Compress_Coarsen_Nodes: %s\n", "Error Allocating Memory 4");
#endif
    count = 0;
    for (iNode = 0; iNode < NNode; iNode++) {
        if (NodeCoarsen[iNode] < 0)
            continue;
        X_New[count] = X[iNode];
        Y_New[count] = Y[iNode];
        count++;
    }
    // Now Swap the memory
    if (X != NULL)
        delete[] X;
    if (Y != NULL)
        delete[] Y;
    X = X_New;
    Y = Y_New;
    X_New = NULL;
    Y_New = NULL;
    
    // Now Compress the Fields
    for (iField = 0; iField < NField; iField++) {
        Field_New = new double[count];
#ifdef DEBUG
        if (Field_New == NULL)
            error("Compress_Coarsen_Nodes: %s\n", "Error Allocating Memory 5");
#endif
        count = 0;
        for (iNode = 0; iNode < NNode; iNode++) {
            if (NodeCoarsen[iNode] < 0)
                continue;
            Field_New[count] = Field[iField][iNode];
            count++;
        }

        // Now Swap the memory
        if (Field[iField] != NULL)
            delete[] Field[iField];
        Field[iField] = Field_New;
        Field_New = NULL;
    }

    // Update the Boundaries
    for (ibnd = 0; ibnd < NBoundary; ibnd++) {
        for (ibs = 0; ibs < NBoundarySegments[ibnd]; ibs++) {
            // Node 1
            iNode = BoundarySegments[ibnd][ibs][0];
            BoundarySegments[ibnd][ibs][0] = NodeCoarsen[iNode];
            // Node 2
            iNode = BoundarySegments[ibnd][ibs][1];
            BoundarySegments[ibnd][ibs][1] = NodeCoarsen[iNode];
        }
    }

    // Update the NNode
    NNode = count;

    // Compress the Triangle Connectivity
    count = 0;
    for (itri = 0; itri < NTri; itri++) {
        if (Cell2Node_6[itri][0] < 0)
            continue;
        Cell2Node_6[count][0] = NodeCoarsen[Cell2Node_6[itri][0]];
        Cell2Node_6[count][1] = NodeCoarsen[Cell2Node_6[itri][1]];
        Cell2Node_6[count][2] = NodeCoarsen[Cell2Node_6[itri][2]];
        count++;
    }
    NTri = count;

    // Return the Number of Nodes Deleted
    return FinalCandiNodes.max;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Generate_BowyerWatson_Delaunay_TriMesh() {
    int tdim, nt;
    int (*Tri_New)[3];
    TriangleMeshSLK Mesher;

    // initial guess for number of triangles
    tdim = NNode * 3;
    Tri_New = new int[tdim][3];
#ifdef DEBUG
    if (Tri_New == NULL)
        error("Generate_BowyerWatson_Delaunay_TriMesh: %s\n", "Error Allocating Memory 1");
#endif

    // iterative loop to allocate space for triangle indices
    nt = -1;
    while (nt < 0) {
        nt  = Mesher.trimesh(NNode, tdim, NBoundary, NBoundarySegments, BoundarySegments, X, Y, Tri_New);
        //nt = trimesh(NNode, tdim, NBoundary, NBoundarySegments, BoundarySegments, X, Y, Tri_New);

        if (nt < 0) {
            info("Expanding New Triangle Dimension from %d to %d", tdim, tdim + NNode);

            // Delete current memory for triangle
            if (Tri_New != NULL)
                delete[] Tri_New;
            
            // Increment triangle dimension by number of nodes
            tdim += NNode;
            // Allocate new space for triangles
            Tri_New = new int[tdim][3];
#ifdef DEBUG
            if (Tri_New == NULL)
                error("Generate_BowyerWatson_Delaunay_TriMesh: %s\n", "Error Allocating Memory 2");
#endif
        }
    }
    info("Number of New Triangles Generated = %d", nt);

    // Now Swap the Memory
    if (Tri != NULL)
        delete[] Tri;
    Tri = Tri_New;
    NTri = nt;
    Tri_New = NULL;
}

// *****************************************************************************
// *****************************************************************************
void TriangleAdaptiveRefinement::Compute_Derived_FlowField() {
    int  i, iNode;
    double Q1, Q2, Q3, Q4, Gamma;
    double var1;
    
    // Check if Flow Field is active
    if (FuncType != 0)
        return;

    if (NVariable == 0)
        error("Compute_Derived_FlowField: %s", "No Flow Field Variables Found");

    // Check if Flow Field has atleast 4 Variables
    if (NVariable != 4) {
        warn("Compute_Derived_FlowField: %s", "Unable to Compute Derived Quantities");
        return;
    }

    // Allocate Memory to Store Derived Variables
    NField = NVariable + 3;
    Field = NULL;
    Field = new double*[NField];
#ifdef DEBUG
    if (Field == NULL)
        error("Compute_Derived_FlowField: %s\n", "Error Allocating Memory 1");
#endif
    for (i = NVariable; i < NField; i++) {
        Field[i] = NULL;
        Field[i] = new double[NNode];
#ifdef DEBUG
        if (Field[i] == NULL)
            error("Compute_Derived_FlowField: %s\n", "Error Allocating Memory 2");
#endif
    }

    // Get the conserved variables
    for (i = 0; i < NVariable; i++)
        Field[i] = Variable[i];

    // Compute Derived Variables
    Gamma = Constant[3];
    for (iNode = 0; iNode < NNode; iNode++) {
        Q1 = Variable[0][iNode];
        Q2 = Variable[1][iNode];
        Q3 = Variable[2][iNode];
        Q4 = Variable[3][iNode];
        
        var1 = (Q2*Q2 + Q3*Q3)/(Q1*Q1);

        // Compute Velocity Magnitude
        Field[NVariable+0][iNode] = sqrt(var1);

        // Compute Pressure
        Field[NVariable+1][iNode] = (Gamma - 1.0)*Q1*((Q4/Q1) - 0.5*var1);

        // Compute Mach Number
        var1 = sqrt(Gamma*Field[NVariable+1][iNode]/Q1);
        Field[NVariable+2][iNode] = Field[NVariable+0][iNode]/var1;
    }

    // Get the Name of Variables
    FieldName = NULL;
    FieldName = new char*[NField];
#ifdef DEBUG
    if (FieldName == NULL)
        error("Compute_Derived_FlowField: %s\n", "Error Allocating Memory 3");
#endif
    for (i = 0; i < NField; i++) {
        FieldName[i] = NULL;
        FieldName[i] = new char[33];
#ifdef DEBUG
        if (FieldName[i] == NULL)
            error("Compute_Derived_FlowField: %s\n", "Error Allocating Memory 4");
#endif
    }

    // Get the conserved variables Names
    for (i = 0; i < NVariable; i++)
        strcpy(FieldName[i], VariableName[i]);
    
    // Get the Name of Derived Variables
    strcpy(FieldName[NVariable+0], "Velocity_Magnitude");
    strcpy(FieldName[NVariable+1], "Pressure");
    strcpy(FieldName[NVariable+2], "Mach_Number");
}

