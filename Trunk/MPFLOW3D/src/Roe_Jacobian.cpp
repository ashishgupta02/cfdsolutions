/*******************************************************************************
 * File:        Roe_Jacobian.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "CompressibleUtils.h"
#include "Material.h"
#include "Solver.h"

// Static Variable for Speed Up
static int      RoeJac_DB   = 0;
static double **RoeJac_dFdL = NULL;
static double **RoeJac_dFdR = NULL;
static double **RoeJac_ARoe = NULL;

//------------------------------------------------------------------------------
//! Create Roe Jacobian Data Structure
//------------------------------------------------------------------------------
void RoeJac_Init(void) {
    int i;
    double *tmp = NULL;
    
    // Check if Roe Data Structure is required
    if (RoeJac_DB == 0) {
        RoeJac_dFdL = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp         = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            RoeJac_dFdL[i] = &tmp[i*NEQUATIONS];
        tmp         = NULL;
        RoeJac_dFdR = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp         = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            RoeJac_dFdR[i] = &tmp[i*NEQUATIONS];
        tmp         = NULL;
        RoeJac_ARoe = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp         = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            RoeJac_ARoe[i] = &tmp[i*NEQUATIONS];
        tmp         = NULL;
        RoeJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete Roe Jacobian Data Structure
//------------------------------------------------------------------------------
void RoeJac_Finalize(void) {
    double *tmp = NULL;
    
    tmp = RoeJac_dFdL[0];
    free(tmp);
    tmp = RoeJac_dFdR[0];
    free(tmp);
    tmp = RoeJac_ARoe[0];
    free(tmp);
    tmp = NULL;
    free(RoeJac_dFdL);
    free(RoeJac_dFdR);
    free(RoeJac_ARoe);
    RoeJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset Roe Jacobian Data Structure
//------------------------------------------------------------------------------
void RoeJac_Reset(void) {
    int i, j;
    
    if (RoeJac_DB == 0)
        RoeJac_Init();
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        for (j = 0; j < NEQUATIONS; j++) {
            RoeJac_dFdL[i][j] = 0.0;
            RoeJac_dFdR[i][j] = 0.0;
            RoeJac_ARoe[i][j] = 0.0;
        }
    }
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_Roe_Exact(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes Roe Jacobian and Build System of Equation: Roe Approximate
//! Note: Jacobian are computed using first order Q's
//------------------------------------------------------------------------------
void Compute_Jacobian_Roe_Approximate(int AddTime, int Iteration) {
    int i, j, k, iNode, iEdge, ibEdge;
    int node_L, node_R;
    int idgn, idgnL, idgnR, ofdgnL, ofdgnR;
    double area, tmp;
    Vector3D areavec;
    
    // Initialize/Reset the Data Structure
    RoeJac_Reset();
    
    // Initialize the CRS Matrix
    for (i = 0; i < SolverBlockMatrix.DIM; i++) {
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                SolverBlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Internal Edges
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
        // Get two nodes of edge
        node_L = intEdge[iEdge].node[0];
        node_R = intEdge[iEdge].node[1];
        
        // Get area vector
        areavec = intEdge[iEdge].areav;
        area    = areavec.magnitude();
        
        // Get the diagonal Locations
        idgnL = SolverBlockMatrix.IAU[node_L];
        idgnR = SolverBlockMatrix.IAU[node_R];
        
        // Get the Off-Diagonal Locations
        // node_L: ofdgnL-> node_R;
        // node_R: ofdgnR-> node_L;
        ofdgnL = -1;
        for (i = SolverBlockMatrix.IA[node_L]; i < SolverBlockMatrix.IA[node_L+1]; i++) {
            if (SolverBlockMatrix.JA[i] == node_R) {
                ofdgnL = i;
                break;
            }
        }
        ofdgnR = -1;
        for (i = SolverBlockMatrix.IA[node_R]; i < SolverBlockMatrix.IA[node_R+1]; i++) {
            if (SolverBlockMatrix.JA[i] == node_L) {
                ofdgnR = i;
                break;
            }
        }
        
        // Initialize the Helper Matrix
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                RoeJac_dFdL[i][j] = 0.0;
                RoeJac_dFdR[i][j] = 0.0;
                RoeJac_ARoe[i][j] = 0.0;
            }
        }
        
        // Compute RoeJac_dFdL
        Compute_Flux_Jacobian_Euler_Convective(node_L, areavec, RoeJac_dFdL);
        
        // Compute RoeJac_dFdR
        Compute_Flux_Jacobian_Euler_Convective(node_R, areavec, RoeJac_dFdR);
        
        // Compute Dissipation Matrix Roe
        Compute_Dissipation_Matrix_Roe(node_L, node_R, areavec, RoeJac_ARoe);
        
        // Finally Compute the dFlux/dQ_L and dFlux/dQ_R
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                RoeJac_dFdL[i][j] = 0.5*(RoeJac_dFdL[i][j] + RoeJac_ARoe[i][j])*area;
                RoeJac_dFdR[i][j] = 0.5*(RoeJac_dFdR[i][j] - RoeJac_ARoe[i][j])*area;
            }
        }
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                // Diagonal
                SolverBlockMatrix.A[idgnL][i][j] += RoeJac_dFdL[i][j];
                SolverBlockMatrix.A[idgnR][i][j] -= RoeJac_dFdR[i][j];
                // Off-Diagonal
                SolverBlockMatrix.A[ofdgnL][i][j] += RoeJac_dFdR[i][j];
                SolverBlockMatrix.A[ofdgnR][i][j] -= RoeJac_dFdL[i][j];
            }
        }
    }
    
    // Boundary Edges
    for (ibEdge = 0; ibEdge < nBEdge; ibEdge++) {
        // Get two nodes of edge
        node_L = bndEdge[ibEdge].node[0];
        node_R = bndEdge[ibEdge].node[1];
        
        // Get area vector
        areavec = bndEdge[ibEdge].areav;
        area    = areavec.magnitude();
        
        // Get the diagonal Locations - Only Physical
        // No diagonal and off-diagonal locations exists for Ghost Nodes
        idgnL = SolverBlockMatrix.IAU[node_L];
        
        // Initialize the Helper Matrix - Only Physical
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                RoeJac_dFdL[i][j] = 0.0;
                RoeJac_ARoe[i][j] = 0.0;
            }
        }
        
        // Compute RoeJac_dFdL
        Compute_Flux_Jacobian_Euler_Convective(node_L, areavec, RoeJac_dFdL);
        
        // Compute Dissipation Matrix Roe
        Compute_Dissipation_Matrix_Roe(node_L, node_R, areavec, RoeJac_ARoe);
        
        
        // Finally Compute the dFlux/dQ_L
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                RoeJac_dFdL[i][j] = 0.5*(RoeJac_dFdL[i][j] + RoeJac_ARoe[i][j])*area;
        }
        
        // Update the Diagonal Term of Physical Node Only
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                SolverBlockMatrix.A[idgnL][i][j] += RoeJac_dFdL[i][j];
        }
    }
    
    // Copy the Residuals to Block Matrix which is B
    // And Copy Cinv/DeltaT
    for (iNode = 0; iNode < nNode; iNode++) {
        // Get the LHS
        SolverBlockMatrix.B[iNode][0] = -(Res1_Conv[iNode] + Res1_Diss[iNode]);
        SolverBlockMatrix.B[iNode][1] = -(Res2_Conv[iNode] + Res2_Diss[iNode]);
        SolverBlockMatrix.B[iNode][2] = -(Res3_Conv[iNode] + Res3_Diss[iNode]);
        SolverBlockMatrix.B[iNode][3] = -(Res4_Conv[iNode] + Res4_Diss[iNode]);
        SolverBlockMatrix.B[iNode][4] = -(Res5_Conv[iNode] + Res5_Diss[iNode]);
        
        if (AddTime == TRUE) {
            // Get the diagonal location
            idgn = SolverBlockMatrix.IAU[iNode];
            // Initialize the Helper Matrix
            for (i = 0; i < NEQUATIONS; i++)
                for (j = 0; j < NEQUATIONS; j++)
                    RoeJac_ARoe[i][j] = 0.0;

            // Compute the Transformation Matrix: Cinv = Binv.A
            Compute_Transformed_Preconditioner_Matrix_Roe(iNode, 1, RoeJac_ARoe);

            // Scale the Diagonal Matrix
            tmp = cVolume[iNode]/DeltaT[iNode];
            for (i = 0; i < NEQUATIONS; i++)
                for (j = 0; j < NEQUATIONS; j++)
                    RoeJac_ARoe[i][j] *= tmp;

            // Add to the Diagonal SOE
            for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
                for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                        SolverBlockMatrix.A[idgn][j][k] += RoeJac_ARoe[j][k];
            }
        }
    }
}
