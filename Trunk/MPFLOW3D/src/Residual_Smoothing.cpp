/*******************************************************************************
 * File:        Residual_Smoothing.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifdef DEBUG
#include <assert.h>
#endif

// Custom header files
#include "Trim_Utils.h"
#include "RestartIO.h"
#include "MeshIO.h"
#include "Commons.h"
#include "Solver.h"
#include "MC.h"
#include "Gradient.h"
#include "CompressibleUtils.h"
#include "Material.h"
#include "DebugSolver.h"

// Static Variable for Speed Up
static int    RM_DB             = 0;
static int    *nodeType         = NULL;
static double *RM_Res_Conv      = NULL;
static double *RM_Res_Conv_Old  = NULL;
static double *RM_Res_Diss      = NULL;
static double *RM_Res_Diss_Old  = NULL;
static double *RM_Relaxation    = NULL;
static double *cArea            = NULL;

// Shared Variables
double *RM_MaxEigenValue        = NULL;
double *RM_SumMaxEigenValue     = NULL;

//------------------------------------------------------------------------------
//! Create Residual Smoothing Data Structure
//------------------------------------------------------------------------------
void Residual_Smoothing_Init(void) {
    int node_L, node_R;
    double area;
    
    if (RM_DB == 0) {
        nodeType    = new int[nNode];
        
        // Allocate based on Residual Smoothing Scheme
        if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION)
            RM_Res_Conv = new double[NEQUATIONS*nNode];
        else 
            RM_Res_Conv = NULL;
        RM_Res_Diss = new double[NEQUATIONS*nNode];
        
        // Allocate for Implicit Method
        if ((ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT)
                || (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT_CONSTANT)) {
            // Allocate based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION)
                RM_Res_Conv_Old = new double[NEQUATIONS*nNode];
            else
                RM_Res_Conv_Old = NULL;
            RM_Res_Diss_Old = new double[NEQUATIONS*nNode];
        } else {
            RM_Res_Conv_Old = NULL;
            RM_Res_Diss_Old = NULL;
        }

        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
            RM_SumMaxEigenValue = new double[nNode];
            // Allocate Sum of Number of Nodes Connected to each nodes
            RM_MaxEigenValue    = new double[crs_IA_Node2Node[nNode] - crs_IA_Node2Node[0]];
            RM_Relaxation       = new double[crs_IA_Node2Node[nNode] - crs_IA_Node2Node[0]];
        } else {
            RM_Relaxation       = NULL;
            RM_MaxEigenValue    = NULL;
            RM_SumMaxEigenValue = NULL;
        }
        
        // Allocate memory to store control volume area
        cArea = new double[nNode];
        
        // Initialize and Classify the Node Type Based on Boundary Type
        for (int inode = 0; inode < nNode; inode++) {
            nodeType[inode] = -1;
            cArea[inode]    = 0.0;
        }
        for (int i = 0; i < nBEdge; i++)
            nodeType[bndEdge[i].node[0]] = bndEdge[i].type;

        // Compute the Control Volume Surface Area
        // Internal Edges
        for (int i = 0; i < nEdge; i++) {
            // Get two nodes of edge
            node_L = intEdge[i].node[0];
            node_R = intEdge[i].node[1];

            // Get area vector
            area   = intEdge[i].areav.magnitude();

            // Left Node
            cArea[node_L] += area;

            // Right Node
            cArea[node_R] += area;
        }
        // Boundary Edges
        for (int i = 0; i < nBEdge; i++) {
            // Get two nodes of edge
            node_L = bndEdge[i].node[0];
            node_R = bndEdge[i].node[1];

            // Get area vector
            area   = bndEdge[i].areav.magnitude();

            // Left Node
            cArea[node_L] += area;
        }
        RM_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete Residual Smoothing Data Structure
//------------------------------------------------------------------------------
void Residual_Smoothing_Finalize(void) {
    // Static Variables
    if (nodeType != NULL)
        delete[] nodeType;
    if (RM_Res_Conv != NULL)
        delete[] RM_Res_Conv;
    if (RM_Res_Conv_Old != NULL)
        delete[] RM_Res_Conv_Old;
    if (RM_Res_Diss != NULL)
        delete[] RM_Res_Diss;
    if (RM_Res_Diss_Old != NULL)
        delete[] RM_Res_Diss_Old;
    if (cArea != NULL)
        delete[] cArea;
    
    // Shared Variables
    if (RM_Relaxation != NULL)
        delete[] RM_Relaxation;
    if (RM_MaxEigenValue != NULL)
        delete[] RM_MaxEigenValue;
    if (RM_SumMaxEigenValue != NULL)
        delete[] RM_SumMaxEigenValue;
    RM_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset Residual Smoothing Data Structure
//------------------------------------------------------------------------------
void Residual_Smoothing_Reset(void) { 
    if (RM_DB == 0)
        Residual_Smoothing_Init();
    
    if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
        for (int i = 0; i < nNode; i++)
            RM_SumMaxEigenValue[i] = 0.0;
        for (int i = 0; i < (crs_IA_Node2Node[nNode] - crs_IA_Node2Node[0]); i++) {
            RM_MaxEigenValue[i]    = 0.0;
            RM_Relaxation[i]       = 0.0;
        }
    }
}

//------------------------------------------------------------------------------
//! Explicit Residual Smoothing 
//------------------------------------------------------------------------------
void Residual_Smoothing_Explicit(void) {
    double Coeff1, Coeff2;
    double Q[5], mach;
    int nid;
    
    // Compute Unweighted Residual Smooth
    if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_EXPLICIT) {
        // Compute Res_t
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Terms
                Res1[inode] = DeltaT[inode]*Res1[inode];
                Res2[inode] = DeltaT[inode]*Res2[inode];
                Res3[inode] = DeltaT[inode]*Res3[inode];
                Res4[inode] = DeltaT[inode]*Res4[inode];
                Res5[inode] = DeltaT[inode]*Res5[inode];
            }
            // Dissipation Terms
            Res1_Diss[inode] = DeltaT[inode]*Res1_Diss[inode];
            Res2_Diss[inode] = DeltaT[inode]*Res2_Diss[inode];
            Res3_Diss[inode] = DeltaT[inode]*Res3_Diss[inode];
            Res4_Diss[inode] = DeltaT[inode]*Res4_Diss[inode];
            Res5_Diss[inode] = DeltaT[inode]*Res5_Diss[inode];
        }
        
        Coeff1 = ResidualSmoothRelaxation;
        for (int iSmooth = 0; iSmooth < ResidualSmoothNIteration; iSmooth++) {
            for (int inode = 0; inode < nNode; inode++) {
                // Only Internal Nodes
                if (nodeType[inode] != -1)
                    continue;
                
                // Only Local Residual Smoothing
                if (ResidualSmoothType == RESIDUAL_SMOOTH_TYPE_LOCAL) {
                    // Get the Variables
                    Q[0] = Q1[inode];
                    Q[1] = Q2[inode];
                    Q[2] = Q3[inode];
                    Q[3] = Q4[inode];
                    Q[4] = Q5[inode];

                    // Compute Mach Number
                    mach = Material_Get_Mach(Q);
                    if ((mach > 0.01*Ref_Mach) && (mach > 1.0e-5))
                        continue;
                }
                
                Coeff2 = ((1.0 - ResidualSmoothRelaxation)/(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    RM_Res_Conv[NEQUATIONS*inode + 0] = Coeff1*Res1[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 1] = Coeff1*Res2[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 2] = Coeff1*Res3[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 3] = Coeff1*Res4[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 4] = Coeff1*Res5[inode];
                }
                // Dissipative Term
                RM_Res_Diss[NEQUATIONS*inode + 0] = Coeff1*Res1_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 1] = Coeff1*Res2_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 2] = Coeff1*Res3_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 3] = Coeff1*Res4_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 4] = Coeff1*Res5_Diss[inode];
                for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                    nid = crs_JA_Node2Node[i];
                    // Compute based on Residual Smoothing Scheme
                    if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                        // Convective Term
                        RM_Res_Conv[NEQUATIONS*inode + 0] += Coeff2*Res1[nid];
                        RM_Res_Conv[NEQUATIONS*inode + 1] += Coeff2*Res2[nid];
                        RM_Res_Conv[NEQUATIONS*inode + 2] += Coeff2*Res3[nid];
                        RM_Res_Conv[NEQUATIONS*inode + 3] += Coeff2*Res4[nid];
                        RM_Res_Conv[NEQUATIONS*inode + 4] += Coeff2*Res5[nid];
                    }
                    // Dissipative Term
                    RM_Res_Diss[NEQUATIONS*inode + 0] += Coeff2*Res1_Diss[nid];
                    RM_Res_Diss[NEQUATIONS*inode + 1] += Coeff2*Res2_Diss[nid];
                    RM_Res_Diss[NEQUATIONS*inode + 2] += Coeff2*Res3_Diss[nid];
                    RM_Res_Diss[NEQUATIONS*inode + 3] += Coeff2*Res4_Diss[nid];
                    RM_Res_Diss[NEQUATIONS*inode + 4] += Coeff2*Res5_Diss[nid];
                }
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    Res1[inode] = RM_Res_Conv[NEQUATIONS*inode + 0];
                    Res2[inode] = RM_Res_Conv[NEQUATIONS*inode + 1];
                    Res3[inode] = RM_Res_Conv[NEQUATIONS*inode + 2];
                    Res4[inode] = RM_Res_Conv[NEQUATIONS*inode + 3];
                    Res5[inode] = RM_Res_Conv[NEQUATIONS*inode + 4];
                }
                // Dissipative Term
                Res1_Diss[inode] = RM_Res_Diss[NEQUATIONS*inode + 0];
                Res2_Diss[inode] = RM_Res_Diss[NEQUATIONS*inode + 1];
                Res3_Diss[inode] = RM_Res_Diss[NEQUATIONS*inode + 2];
                Res4_Diss[inode] = RM_Res_Diss[NEQUATIONS*inode + 3];
                Res5_Diss[inode] = RM_Res_Diss[NEQUATIONS*inode + 4];
            }
        }
        
        // Compute back Res
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Term
                Res1[inode] = Res1[inode]/DeltaT[inode];
                Res2[inode] = Res2[inode]/DeltaT[inode];
                Res3[inode] = Res3[inode]/DeltaT[inode];
                Res4[inode] = Res4[inode]/DeltaT[inode];
                Res5[inode] = Res5[inode]/DeltaT[inode];
            }
            // Dissipative Term
            Res1_Diss[inode] = Res1_Diss[inode]/DeltaT[inode];
            Res2_Diss[inode] = Res2_Diss[inode]/DeltaT[inode];
            Res3_Diss[inode] = Res3_Diss[inode]/DeltaT[inode];
            Res4_Diss[inode] = Res4_Diss[inode]/DeltaT[inode];
            Res5_Diss[inode] = Res5_Diss[inode]/DeltaT[inode];
        }
    }
    
    // Weighted Residual Smoothing for Anisotropic Grid with large Face sizes
    // the surface of a control volume for a point
    if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_EXPLICIT_WEIGHTED) {
        int node_L, node_R;
        double area;
        
        // Compute Res_t
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Term
                Res1[inode] = DeltaT[inode]*Res1[inode];
                Res2[inode] = DeltaT[inode]*Res2[inode];
                Res3[inode] = DeltaT[inode]*Res3[inode];
                Res4[inode] = DeltaT[inode]*Res4[inode];
                Res5[inode] = DeltaT[inode]*Res5[inode];
            }
            // Dissipation Terms
            Res1_Diss[inode] = DeltaT[inode]*Res1_Diss[inode];
            Res2_Diss[inode] = DeltaT[inode]*Res2_Diss[inode];
            Res3_Diss[inode] = DeltaT[inode]*Res3_Diss[inode];
            Res4_Diss[inode] = DeltaT[inode]*Res4_Diss[inode];
            Res5_Diss[inode] = DeltaT[inode]*Res5_Diss[inode];
        }
        
        Coeff1 = ResidualSmoothRelaxation;
        for (int iSmooth = 0; iSmooth < ResidualSmoothNIteration; iSmooth++) {
            for (int inode = 0; inode < nNode; inode++) {
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    RM_Res_Conv[NEQUATIONS*inode + 0] = 0.0;
                    RM_Res_Conv[NEQUATIONS*inode + 1] = 0.0;
                    RM_Res_Conv[NEQUATIONS*inode + 2] = 0.0;
                    RM_Res_Conv[NEQUATIONS*inode + 3] = 0.0;
                    RM_Res_Conv[NEQUATIONS*inode + 4] = 0.0;
                }
                // Dissipation Terms
                RM_Res_Diss[NEQUATIONS*inode + 0] = 0.0;
                RM_Res_Diss[NEQUATIONS*inode + 1] = 0.0;
                RM_Res_Diss[NEQUATIONS*inode + 2] = 0.0;
                RM_Res_Diss[NEQUATIONS*inode + 3] = 0.0;
                RM_Res_Diss[NEQUATIONS*inode + 4] = 0.0;
            }

            // Internal Edges
            for (int i = 0; i < nEdge; i++) {
                // Get two nodes of edge
                node_L = intEdge[i].node[0];
                node_R = intEdge[i].node[1];

                // Get area vector
                area   = intEdge[i].areav.magnitude();
                
                // Left Node
                if (nodeType[node_L] == -1) {
                    // Compute based on Residual Smoothing Scheme
                    if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                        // Convective Term
                        RM_Res_Conv[NEQUATIONS*node_L + 0] += area*Res1[node_R];
                        RM_Res_Conv[NEQUATIONS*node_L + 1] += area*Res2[node_R];
                        RM_Res_Conv[NEQUATIONS*node_L + 2] += area*Res3[node_R];
                        RM_Res_Conv[NEQUATIONS*node_L + 3] += area*Res4[node_R];
                        RM_Res_Conv[NEQUATIONS*node_L + 4] += area*Res5[node_R];
                    }
                    // Dissipation Terms
                    RM_Res_Diss[NEQUATIONS*node_L + 0] += area*Res1_Diss[node_R];
                    RM_Res_Diss[NEQUATIONS*node_L + 1] += area*Res2_Diss[node_R];
                    RM_Res_Diss[NEQUATIONS*node_L + 2] += area*Res3_Diss[node_R];
                    RM_Res_Diss[NEQUATIONS*node_L + 3] += area*Res4_Diss[node_R];
                    RM_Res_Diss[NEQUATIONS*node_L + 4] += area*Res5_Diss[node_R];
                }
                
                // Right Node
                if (nodeType[node_R] == -1) {
                    // Compute based on Residual Smoothing Scheme
                    if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                        // Convective Term
                        RM_Res_Conv[NEQUATIONS*node_R + 0] += area*Res1[node_L];
                        RM_Res_Conv[NEQUATIONS*node_R + 1] += area*Res2[node_L];
                        RM_Res_Conv[NEQUATIONS*node_R + 2] += area*Res3[node_L];
                        RM_Res_Conv[NEQUATIONS*node_R + 3] += area*Res4[node_L];
                        RM_Res_Conv[NEQUATIONS*node_R + 4] += area*Res5[node_L];
                    }
                    // Dissipation Terms
                    RM_Res_Diss[NEQUATIONS*node_R + 0] += area*Res1_Diss[node_L];
                    RM_Res_Diss[NEQUATIONS*node_R + 1] += area*Res2_Diss[node_L];
                    RM_Res_Diss[NEQUATIONS*node_R + 2] += area*Res3_Diss[node_L];
                    RM_Res_Diss[NEQUATIONS*node_R + 3] += area*Res4_Diss[node_L];
                    RM_Res_Diss[NEQUATIONS*node_R + 4] += area*Res5_Diss[node_L];
                }
            }

            for (int inode = 0; inode < nNode; inode++) {
                // Only Internal Nodes
                if (nodeType[inode] != -1)
                    continue;
                
                // Only Local Residual Smoothing
                if (ResidualSmoothType == RESIDUAL_SMOOTH_TYPE_LOCAL) {
                    // Get the Variables
                    Q[0] = Q1[inode];
                    Q[1] = Q2[inode];
                    Q[2] = Q3[inode];
                    Q[3] = Q4[inode];
                    Q[4] = Q5[inode];

                    // Compute Mach Number
                    mach = Material_Get_Mach(Q);
                    if ((mach > 0.01*Ref_Mach) && (mach > 1.0e-5))
                        continue;
                }
                
                Coeff2 = (1.0 - ResidualSmoothRelaxation)/cArea[inode];
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    Res1[inode] = Coeff1*Res1[inode] + Coeff2*RM_Res_Conv[NEQUATIONS*inode + 0];
                    Res2[inode] = Coeff1*Res2[inode] + Coeff2*RM_Res_Conv[NEQUATIONS*inode + 1];
                    Res3[inode] = Coeff1*Res3[inode] + Coeff2*RM_Res_Conv[NEQUATIONS*inode + 2];
                    Res4[inode] = Coeff1*Res4[inode] + Coeff2*RM_Res_Conv[NEQUATIONS*inode + 3];
                    Res5[inode] = Coeff1*Res5[inode] + Coeff2*RM_Res_Conv[NEQUATIONS*inode + 4];
                }
                // Dissipation Terms
                Res1_Diss[inode] = Coeff1*Res1_Diss[inode] + Coeff2*RM_Res_Diss[NEQUATIONS*inode + 0];
                Res2_Diss[inode] = Coeff1*Res2_Diss[inode] + Coeff2*RM_Res_Diss[NEQUATIONS*inode + 1];
                Res3_Diss[inode] = Coeff1*Res3_Diss[inode] + Coeff2*RM_Res_Diss[NEQUATIONS*inode + 2];
                Res4_Diss[inode] = Coeff1*Res4_Diss[inode] + Coeff2*RM_Res_Diss[NEQUATIONS*inode + 3];
                Res5_Diss[inode] = Coeff1*Res5_Diss[inode] + Coeff2*RM_Res_Diss[NEQUATIONS*inode + 4];
            }
        }
        
        // Compute back Res
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Term
                Res1[inode] = Res1[inode]/DeltaT[inode];
                Res2[inode] = Res2[inode]/DeltaT[inode];
                Res3[inode] = Res3[inode]/DeltaT[inode];
                Res4[inode] = Res4[inode]/DeltaT[inode];
                Res5[inode] = Res5[inode]/DeltaT[inode];
            }
            // Dissipation Term
            Res1_Diss[inode] = Res1_Diss[inode]/DeltaT[inode];
            Res2_Diss[inode] = Res2_Diss[inode]/DeltaT[inode];
            Res3_Diss[inode] = Res3_Diss[inode]/DeltaT[inode];
            Res4_Diss[inode] = Res4_Diss[inode]/DeltaT[inode];
            Res5_Diss[inode] = Res5_Diss[inode]/DeltaT[inode];
        }
    }
}

//------------------------------------------------------------------------------
//! Implicit Residual Smoothing 
//------------------------------------------------------------------------------
void Residual_Smoothing_Implicit(void) {
    int nid, nneigh;
    double Coeff1, Coeff2;
    double Q[5], mach;
    
    // Compute Smooth Residual using Variable Epsilon
    if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
        // Compute Res_t
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Term
                Res1[inode] = DeltaT[inode]*Res1[inode];
                Res2[inode] = DeltaT[inode]*Res2[inode];
                Res3[inode] = DeltaT[inode]*Res3[inode];
                Res4[inode] = DeltaT[inode]*Res4[inode];
                Res5[inode] = DeltaT[inode]*Res5[inode];

                RM_Res_Conv_Old[NEQUATIONS*inode + 0] = Res1[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 1] = Res2[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 2] = Res3[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 3] = Res4[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 4] = Res5[inode];
            }
            // Dissipation Term
            Res1_Diss[inode] = DeltaT[inode]*Res1_Diss[inode];
            Res2_Diss[inode] = DeltaT[inode]*Res2_Diss[inode];
            Res3_Diss[inode] = DeltaT[inode]*Res3_Diss[inode];
            Res4_Diss[inode] = DeltaT[inode]*Res4_Diss[inode];
            Res5_Diss[inode] = DeltaT[inode]*Res5_Diss[inode];
            
            RM_Res_Diss_Old[NEQUATIONS*inode + 0] = Res1_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 1] = Res2_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 2] = Res3_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 3] = Res4_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 4] = Res5_Diss[inode];
        }
        
        // Compute the Variable Relaxation Factor
        for (int inode = 0; inode < nNode; inode++) {
            for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode + 1]; i++) {
                nneigh = crs_IA_Node2Node[inode + 1] - crs_IA_Node2Node[inode];  
                RM_Relaxation[i] = 1.0 + pow(0.5*(RM_SumMaxEigenValue[inode] - RM_MaxEigenValue[i])/(RM_MaxEigenValue[i]), (2.0/3.0));
                RM_Relaxation[i] = (CFL/CFL_MAX)*((2.0*RM_MaxEigenValue[i])/(RM_MaxEigenValue[i] + RM_SumMaxEigenValue[inode]/nneigh))*RM_Relaxation[i];
                RM_Relaxation[i] = MIN(MAX(0.0, 0.25*(RM_Relaxation[i]*RM_Relaxation[i] - 1.0)), 1.0);
            }
        }
        
        for (int iSmooth = 0; iSmooth < ResidualSmoothNIteration; iSmooth++) {
            for (int inode = 0; inode < nNode; inode++) {
                // Only Internal Nodes
                if (nodeType[inode] != -1)
                    continue;
                
                // Only Local Residual Smoothing
                if (ResidualSmoothType == RESIDUAL_SMOOTH_TYPE_LOCAL) {
                    // Get the Variables
                    Q[0] = Q1[inode];
                    Q[1] = Q2[inode];
                    Q[2] = Q3[inode];
                    Q[3] = Q4[inode];
                    Q[4] = Q5[inode];

                    // Compute Mach Number
                    mach = Material_Get_Mach(Q);
                    if ((mach > 0.01*Ref_Mach) && (mach > 1.0e-5))
                        continue;
                }
                
                Coeff2 = 1.0;
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    RM_Res_Conv[NEQUATIONS*inode + 0] = Res1[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 1] = Res2[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 2] = Res3[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 3] = Res4[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 4] = Res5[inode];
                }
                // Dissipation Term
                RM_Res_Diss[NEQUATIONS*inode + 0] = Res1_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 1] = Res2_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 2] = Res3_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 3] = Res4_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 4] = Res5_Diss[inode];
                
                for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                    nid     = crs_JA_Node2Node[i];
                    Coeff1  = RM_Relaxation[i];
                    Coeff2 += Coeff1;
                    
                    // Compute based on Residual Smoothing Scheme
                    if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                        // Convective Term
                        RM_Res_Conv[NEQUATIONS*inode + 0] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 0];
                        RM_Res_Conv[NEQUATIONS*inode + 1] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 1];
                        RM_Res_Conv[NEQUATIONS*inode + 2] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 2];
                        RM_Res_Conv[NEQUATIONS*inode + 3] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 3];
                        RM_Res_Conv[NEQUATIONS*inode + 4] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 4];
                    }
                    // Dissipation Term
                    RM_Res_Diss[NEQUATIONS*inode + 0] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 0];
                    RM_Res_Diss[NEQUATIONS*inode + 1] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 1];
                    RM_Res_Diss[NEQUATIONS*inode + 2] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 2];
                    RM_Res_Diss[NEQUATIONS*inode + 3] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 3];
                    RM_Res_Diss[NEQUATIONS*inode + 4] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 4];
                }
                // Compute the New Residual and Copy it to Old
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    RM_Res_Conv_Old[NEQUATIONS*inode + 0] = RM_Res_Conv[NEQUATIONS*inode + 0]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 1] = RM_Res_Conv[NEQUATIONS*inode + 1]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 2] = RM_Res_Conv[NEQUATIONS*inode + 2]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 3] = RM_Res_Conv[NEQUATIONS*inode + 3]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 4] = RM_Res_Conv[NEQUATIONS*inode + 4]/Coeff2;
                }
                // Dissipation Term
                RM_Res_Diss_Old[NEQUATIONS*inode + 0] = RM_Res_Diss[NEQUATIONS*inode + 0]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 1] = RM_Res_Diss[NEQUATIONS*inode + 1]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 2] = RM_Res_Diss[NEQUATIONS*inode + 2]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 3] = RM_Res_Diss[NEQUATIONS*inode + 3]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 4] = RM_Res_Diss[NEQUATIONS*inode + 4]/Coeff2;
            }
        }
        
        // Compute back Res
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Term
                Res1[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 0]/DeltaT[inode];
                Res2[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 1]/DeltaT[inode];
                Res3[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 2]/DeltaT[inode];
                Res4[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 3]/DeltaT[inode];
                Res5[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 4]/DeltaT[inode];
            }
            // Dissipation Term
            Res1_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 0]/DeltaT[inode];
            Res2_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 1]/DeltaT[inode];
            Res3_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 2]/DeltaT[inode];
            Res4_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 3]/DeltaT[inode];
            Res5_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 4]/DeltaT[inode];
        }
    }
    
    // Compute Smooth Residual using Constant Epsilon
    if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT_CONSTANT) {
        // Compute Res_t
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Term
                Res1[inode] = DeltaT[inode]*Res1[inode];
                Res2[inode] = DeltaT[inode]*Res2[inode];
                Res3[inode] = DeltaT[inode]*Res3[inode];
                Res4[inode] = DeltaT[inode]*Res4[inode];
                Res5[inode] = DeltaT[inode]*Res5[inode];

                RM_Res_Conv_Old[NEQUATIONS*inode + 0] = Res1[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 1] = Res2[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 2] = Res3[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 3] = Res4[inode];
                RM_Res_Conv_Old[NEQUATIONS*inode + 4] = Res5[inode];
            }
            // Dissipation Term
            Res1_Diss[inode] = DeltaT[inode]*Res1_Diss[inode];
            Res2_Diss[inode] = DeltaT[inode]*Res2_Diss[inode];
            Res3_Diss[inode] = DeltaT[inode]*Res3_Diss[inode];
            Res4_Diss[inode] = DeltaT[inode]*Res4_Diss[inode];
            Res5_Diss[inode] = DeltaT[inode]*Res5_Diss[inode];
            
            RM_Res_Diss_Old[NEQUATIONS*inode + 0] = Res1_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 1] = Res2_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 2] = Res3_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 3] = Res4_Diss[inode];
            RM_Res_Diss_Old[NEQUATIONS*inode + 4] = Res5_Diss[inode];
        }
        
        Coeff1 = ResidualSmoothRelaxation;
        for (int iSmooth = 0; iSmooth < ResidualSmoothNIteration; iSmooth++) {
            for (int inode = 0; inode < nNode; inode++) {
                // Only Internal Nodes
                if (nodeType[inode] != -1)
                    continue;
                
                // Only Local Residual Smoothing
                if (ResidualSmoothType == RESIDUAL_SMOOTH_TYPE_LOCAL) {
                    // Get the Variables
                    Q[0] = Q1[inode];
                    Q[1] = Q2[inode];
                    Q[2] = Q3[inode];
                    Q[3] = Q4[inode];
                    Q[4] = Q5[inode];

                    // Compute Mach Number
                    mach = Material_Get_Mach(Q);
                    if ((mach > 0.01*Ref_Mach) && (mach > 1.0e-5))
                        continue;
                }
                
                Coeff2 = (1.0 + ResidualSmoothRelaxation*(crs_IA_Node2Node[inode+1] - crs_IA_Node2Node[inode]));
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    RM_Res_Conv[NEQUATIONS*inode + 0] = Res1[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 1] = Res2[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 2] = Res3[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 3] = Res4[inode];
                    RM_Res_Conv[NEQUATIONS*inode + 4] = Res5[inode];
                }
                // Dissipation Term
                RM_Res_Diss[NEQUATIONS*inode + 0] = Res1_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 1] = Res2_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 2] = Res3_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 3] = Res4_Diss[inode];
                RM_Res_Diss[NEQUATIONS*inode + 4] = Res5_Diss[inode];
                
                for (int i = crs_IA_Node2Node[inode]; i < crs_IA_Node2Node[inode+1]; i++) {
                    nid   = crs_JA_Node2Node[i];
                    
                    // Compute based on Residual Smoothing Scheme
                    if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                        // Convective Term
                        RM_Res_Conv[NEQUATIONS*inode + 0] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 0];
                        RM_Res_Conv[NEQUATIONS*inode + 1] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 1];
                        RM_Res_Conv[NEQUATIONS*inode + 2] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 2];
                        RM_Res_Conv[NEQUATIONS*inode + 3] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 3];
                        RM_Res_Conv[NEQUATIONS*inode + 4] += Coeff1*RM_Res_Conv_Old[NEQUATIONS*nid + 4];
                    }
                    // Dissipation Term
                    RM_Res_Diss[NEQUATIONS*inode + 0] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 0];
                    RM_Res_Diss[NEQUATIONS*inode + 1] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 1];
                    RM_Res_Diss[NEQUATIONS*inode + 2] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 2];
                    RM_Res_Diss[NEQUATIONS*inode + 3] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 3];
                    RM_Res_Diss[NEQUATIONS*inode + 4] += Coeff1*RM_Res_Diss_Old[NEQUATIONS*nid + 4];
                }
                // Compute the New Residual and Copy it to Old
                // Compute based on Residual Smoothing Scheme
                if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                    // Convective Term
                    RM_Res_Conv_Old[NEQUATIONS*inode + 0] = RM_Res_Conv[NEQUATIONS*inode + 0]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 1] = RM_Res_Conv[NEQUATIONS*inode + 1]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 2] = RM_Res_Conv[NEQUATIONS*inode + 2]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 3] = RM_Res_Conv[NEQUATIONS*inode + 3]/Coeff2;
                    RM_Res_Conv_Old[NEQUATIONS*inode + 4] = RM_Res_Conv[NEQUATIONS*inode + 4]/Coeff2;
                }
                // Dissipation Term
                RM_Res_Diss_Old[NEQUATIONS*inode + 0] = RM_Res_Diss[NEQUATIONS*inode + 0]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 1] = RM_Res_Diss[NEQUATIONS*inode + 1]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 2] = RM_Res_Diss[NEQUATIONS*inode + 2]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 3] = RM_Res_Diss[NEQUATIONS*inode + 3]/Coeff2;
                RM_Res_Diss_Old[NEQUATIONS*inode + 4] = RM_Res_Diss[NEQUATIONS*inode + 4]/Coeff2;
            }
        }
        
        // Compute back Res
        for (int inode = 0; inode < nNode; inode++) {
            // Compute based on Residual Smoothing Scheme
            if (ResidualSmoothScheme == RESIDUAL_SMOOTH_SCHEME_CONVECTIVE_DISSIPATION) {
                // Convective Term
                Res1[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 0]/DeltaT[inode];
                Res2[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 1]/DeltaT[inode];
                Res3[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 2]/DeltaT[inode];
                Res4[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 3]/DeltaT[inode];
                Res5[inode] = RM_Res_Conv_Old[NEQUATIONS*inode + 4]/DeltaT[inode];
            }
            // Dissipation Term
            Res1_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 0]/DeltaT[inode];
            Res2_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 1]/DeltaT[inode];
            Res3_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 2]/DeltaT[inode];
            Res4_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 3]/DeltaT[inode];
            Res5_Diss[inode] = RM_Res_Diss_Old[NEQUATIONS*inode + 4]/DeltaT[inode];
        }
    }
}

//------------------------------------------------------------------------------
//! Perform Residual Smoothing
//------------------------------------------------------------------------------
void Residual_Smoothing(void) {
    switch (ResidualSmoothMethod) {
        case RESIDUAL_SMOOTH_METHOD_EXPLICIT:
            Residual_Smoothing_Explicit();
            break;
        case RESIDUAL_SMOOTH_METHOD_EXPLICIT_WEIGHTED:
            Residual_Smoothing_Explicit();
            break;
        case RESIDUAL_SMOOTH_METHOD_IMPLICIT:
            Residual_Smoothing_Implicit();
            break;
        case RESIDUAL_SMOOTH_METHOD_IMPLICIT_CONSTANT:
            Residual_Smoothing_Implicit();
            break;
    }
}

