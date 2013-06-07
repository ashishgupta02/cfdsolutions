/*******************************************************************************
 * File:        FiniteDifference_Jacobian.cpp
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
static int      FDJac_DB          = 0;
static double  *FDJac_Q_L         = NULL;
static double  *FDJac_Q_R         = NULL;
static double  *FDJac_FluxP_Conv  = NULL;
static double  *FDJac_FluxP_Diss  = NULL;
static double  *FDJac_FluxM_Conv  = NULL;
static double  *FDJac_FluxM_Diss  = NULL;
static double **FDJac_dFdL        = NULL;
static double **FDJac_dFdR        = NULL;

//------------------------------------------------------------------------------
//! Create Finite Difference Jacobian Data Structure
//------------------------------------------------------------------------------
void FDJac_Init(void) {
    int i;
    double *tmp = NULL;
    
    // Check if Roe Data Structure is required
    if (FDJac_DB == 0) {
        FDJac_Q_L        = (double *)  malloc(NEQUATIONS*sizeof(double));
        FDJac_Q_R        = (double *)  malloc(NEQUATIONS*sizeof(double));
        FDJac_FluxP_Conv = (double *)  malloc(NEQUATIONS*sizeof(double));
        FDJac_FluxP_Diss = (double *)  malloc(NEQUATIONS*sizeof(double));
        FDJac_FluxM_Conv = (double *)  malloc(NEQUATIONS*sizeof(double));
        FDJac_FluxM_Diss = (double *)  malloc(NEQUATIONS*sizeof(double));
        FDJac_dFdL       = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp            = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i] = &tmp[i*NEQUATIONS];
        tmp            = NULL;
        FDJac_dFdR       = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp            = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdR[i] = &tmp[i*NEQUATIONS];
        tmp            = NULL;
        
        FDJac_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete Finite Difference Jacobian Data Structure
//------------------------------------------------------------------------------
void FDJac_Finalize(void) {
    double *tmp = NULL;
    
    tmp = FDJac_dFdL[0];
    free(tmp);
    tmp = FDJac_dFdR[0];
    free(tmp);
    tmp = NULL;
    free(FDJac_Q_L);
    free(FDJac_Q_R);
    free(FDJac_FluxP_Conv);
    free(FDJac_FluxP_Diss);
    free(FDJac_FluxM_Conv);
    free(FDJac_FluxM_Diss);
    free(FDJac_dFdL);
    free(FDJac_dFdR);
    FDJac_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset Finite Difference Jacobian Data Structure
//------------------------------------------------------------------------------
void FDJac_Reset(void) {
    int i, j;
    
    if (FDJac_DB == 0)
        FDJac_Init();
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        FDJac_Q_L[i]        = 0.0;
        FDJac_Q_R[i]        = 0.0;
        FDJac_FluxP_Conv[i] = 0.0;
        FDJac_FluxP_Diss[i] = 0.0;
        FDJac_FluxM_Conv[i] = 0.0;
        FDJac_FluxM_Diss[i] = 0.0;
        for (j = 0; j < NEQUATIONS; j++) {
            FDJac_dFdL[i][j] = 0.0;
            FDJac_dFdR[i][j] = 0.0;
        }
    }
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//! Note Currently Internal Edges works : Forward, Backward or Central
//! Boundary Edges Only Central -- Some more research needed.
//! Note: Jacobian is computed using first order Q's
//------------------------------------------------------------------------------
void Compute_Jacobian_FiniteDifference(int AddTime, int Iteration) {
    int i, j, k, iNode, iEdge, ibEdge;
    int node_L, node_R;
    int idgn, idgnL, idgnR, ofdgnL, ofdgnR;
    Vector3D areavec;
    double eps = 1.0E-8;
    int FiniteDifference = JacobianMethod;
    double Coeff = 0.5;
    double tmp;
    
    // Initialize the CRS Matrix
    for (i = 0; i < SolverBlockMatrix.DIM; i++) {
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                SolverBlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Copy the Residuals to Block Matrix which is B
    // And Copy Cinv.Vol/DeltaT
    for (iNode = 0; iNode < nNode; iNode++) {
        // Get the diagonal location
        idgn = SolverBlockMatrix.IAU[iNode];
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
                    FDJac_dFdL[i][j] = 0.0;

            // Compute the Transformation Matrix: Cinv = Binv.A
            Compute_Transformed_Preconditioner_Matrix(iNode, 1, FDJac_dFdL);

            // Scale the Diagonal Matrix
            tmp = cVolume[iNode]/DeltaT[iNode];
            for (i = 0; i < NEQUATIONS; i++)
                for (j = 0; j < NEQUATIONS; j++)
                    FDJac_dFdL[i][j] *= tmp;

            // Add to the Diagonal SOE
            for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
                for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                        SolverBlockMatrix.A[idgn][j][k] = FDJac_dFdL[j][k];
            }
        }
    }
    
    // Set Flux Recompute Flag
    CogSolver.FluxRecomputeFlag = TRUE;
    
    // Check if Alternating
    if (FiniteDifference == JACOBIAN_METHOD_ALTERNATE) {
        if (Iteration%2 == 0)
            FiniteDifference = JACOBIAN_METHOD_FORWARD;
        else
            FiniteDifference = JACOBIAN_METHOD_BACKWARD;
    }
    
    // Set the value of Coeff based on Finite Difference Method
    if (FiniteDifference == JACOBIAN_METHOD_CENTRAL) // Central
        Coeff =  0.5;
    else // Forward or Backward
        Coeff =  1.0;
    
    // Internal Edges
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
        // Get two nodes of edge
        node_L = intEdge[iEdge].node[0];
        node_R = intEdge[iEdge].node[1];
        
        // Get area vector
        areavec = intEdge[iEdge].areav;
        
        // Backup the Copy of Q_L and Q_R
        // Left
        FDJac_Q_L[0] = Q1[node_L];
        FDJac_Q_L[1] = Q2[node_L];
        FDJac_Q_L[2] = Q3[node_L];
        FDJac_Q_L[3] = Q4[node_L];
        FDJac_Q_L[4] = Q5[node_L];
        // Right
        FDJac_Q_R[0] = Q1[node_R];
        FDJac_Q_R[1] = Q2[node_R];
        FDJac_Q_R[2] = Q3[node_R];
        FDJac_Q_R[3] = Q4[node_R];
        FDJac_Q_R[4] = Q5[node_R];
        
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
                FDJac_dFdL[i][j] = 0.0;
                FDJac_dFdR[i][j] = 0.0;
            }
        }
        
        // ---------------------------------------------------------------------
        // - Left Node ---------------------------------------------------------
        // -- Q1 ---------------------------------------------------------------
        // Now Perturb the Left Node Q1 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) // Central or Forward
            Q1[node_L] = FDJac_Q_L[0] + eps;
        else
            Q1[node_L] = FDJac_Q_L[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        
        // Now Perturb the Left Node Q1 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD))  // Central or Backward
            Q1[node_L] = FDJac_Q_L[0] - eps;
        else
            Q1[node_L] = FDJac_Q_L[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][0] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q1
        Q1[node_L] = FDJac_Q_L[0];
        
        // -- Q2 ---------------------------------------------------------------
        // Now Perturb the Left Node Q2 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q2[node_L] = FDJac_Q_L[1] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Left Node Q2 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q2[node_L] = FDJac_Q_L[1] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][1] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q2
        Q2[node_L] = FDJac_Q_L[1];
        
        // -- Q3 ---------------------------------------------------------------
        // Now Perturb the Left Node Q3 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q3[node_L] = FDJac_Q_L[2] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Left Node Q3 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q3[node_L] = FDJac_Q_L[2] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][2] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q3
        Q3[node_L] = FDJac_Q_L[2];
        
        // -- Q4 ---------------------------------------------------------------
        // Now Perturb the Left Node Q4 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q4[node_L] = FDJac_Q_L[3] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Left Node Q4 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q4[node_L] = FDJac_Q_L[3] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][3] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q4
        Q4[node_L] = FDJac_Q_L[3];
        
        // -- Q5 ---------------------------------------------------------------
        // Now Perturb the Left Node Q5 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q5[node_L] = FDJac_Q_L[4] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Left Node Q5 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q5[node_L] = FDJac_Q_L[4] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][4] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q5
        Q5[node_L] = FDJac_Q_L[4];
        
        // ---------------------------------------------------------------------
        // - Right Node --------------------------------------------------------
        // ---------------------------------------------------------------------
        // -- Q1 --
        // Now Perturb the Right Node Q1 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) // Central or Forward
            Q1[node_R] = FDJac_Q_R[0] + eps;
        else
            Q1[node_R] = FDJac_Q_R[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        
        // Now Perturb the Right Node Q1 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) // Central or Backward
            Q1[node_R] = FDJac_Q_R[0] - eps;
        else
            Q1[node_R] = FDJac_Q_R[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdR[i][0] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q1
        Q1[node_R] = FDJac_Q_R[0];
        
        // -- Q2 ---------------------------------------------------------------
        // Now Perturb the Right Node Q2 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q2[node_R] = FDJac_Q_R[1] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Right Node Q2 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q2[node_R] = FDJac_Q_R[1] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdR[i][1] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q2
        Q2[node_R] = FDJac_Q_R[1];
        
        // -- Q3 ---------------------------------------------------------------
        // Now Perturb the Right Node Q3 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q3[node_R] = FDJac_Q_R[2] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Right Node Q3 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q3[node_R] = FDJac_Q_R[2] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdR[i][2] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q3
        Q3[node_R] = FDJac_Q_R[2];
        
        // -- Q4 ---------------------------------------------------------------
        // Now Perturb the Right Node Q4 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q4[node_R] = FDJac_Q_R[3] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Right Node Q4 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q4[node_R] = FDJac_Q_R[3] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdR[i][3] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q4
        Q4[node_R] = FDJac_Q_R[3];
        
        // -- Q5 ---------------------------------------------------------------
        // Now Perturb the Right Node Q5 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q5[node_R] = FDJac_Q_R[4] + eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        }
        
        // Now Perturb the Right Node Q5 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q5[node_R] = FDJac_Q_R[4] - eps;
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdR[i][4] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q5
        Q5[node_R] = FDJac_Q_R[4];
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                // Diagonal
                SolverBlockMatrix.A[idgnL][i][j] += FDJac_dFdL[i][j];
                SolverBlockMatrix.A[idgnR][i][j] -= FDJac_dFdR[i][j];
                // Off-Diagonal
                SolverBlockMatrix.A[ofdgnL][i][j] += FDJac_dFdR[i][j];
                SolverBlockMatrix.A[ofdgnR][i][j] -= FDJac_dFdL[i][j];
            }
        }
    }
    
    // Note: Currently Forcing Boundary Contribution as Central
    // More Research has to be done for Forward and Backward to work here
    FiniteDifference = JACOBIAN_METHOD_CENTRAL;
    Coeff = 0.5;
    // Boundary Edges
    for (ibEdge = 0; ibEdge < nBEdge; ibEdge++) {
        // Get two nodes of edge
        node_L = bndEdge[ibEdge].node[0];
        node_R = bndEdge[ibEdge].node[1];
        
        // Get area vector
        areavec = bndEdge[ibEdge].areav;
        
        // Backup the Copy of Q_L and Q_R
        // Left - Physical
        FDJac_Q_L[0] = Q1[node_L];
        FDJac_Q_L[1] = Q2[node_L];
        FDJac_Q_L[2] = Q3[node_L];
        FDJac_Q_L[3] = Q4[node_L];
        FDJac_Q_L[4] = Q5[node_L];
        // Right - Ghost
        FDJac_Q_R[0] = Q1[node_R];
        FDJac_Q_R[1] = Q2[node_R];
        FDJac_Q_R[2] = Q3[node_R];
        FDJac_Q_R[3] = Q4[node_R];
        FDJac_Q_R[4] = Q5[node_R];
        
        // Get the diagonal Locations - Only Physical
        // No diagonal and off-diagonal locations exists for Ghost Nodes
        idgnL = SolverBlockMatrix.IAU[node_L];
        
        // Initialize the Helper Matrix - Only Physical
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                FDJac_dFdL[i][j] = 0.0;
        }
        
        // ---------------------------------------------------------------------
        // Left Node - Physical ------------------------------------------------
        // -- Q1 ---------------------------------------------------------------
        // Now Perturb the Left Node Q1 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q1[node_L] = FDJac_Q_L[0] + eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
        } else {
            Q1[node_L] = FDJac_Q_L[0];
        }
        // Based on Finite Difference Method, One Sided Do only once
        Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
        // Reset the Ghost Node Values
        Q1[node_R] = FDJac_Q_R[0];
        Q2[node_R] = FDJac_Q_R[1];
        Q3[node_R] = FDJac_Q_R[2];
        Q4[node_R] = FDJac_Q_R[3];
        Q5[node_R] = FDJac_Q_R[4];
        
        // Now Perturb the Left Node Q1 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q1[node_L] = FDJac_Q_L[0] - eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
        } else {
            Q1[node_L] = FDJac_Q_L[0];
        }
        // Based on Finite Difference Method, One Sided Do only once
        Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
        // Reset the Ghost Node Values
        Q1[node_R] = FDJac_Q_R[0];
        Q2[node_R] = FDJac_Q_R[1];
        Q3[node_R] = FDJac_Q_R[2];
        Q4[node_R] = FDJac_Q_R[3];
        Q5[node_R] = FDJac_Q_R[4];
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][0] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q1
        Q1[node_L] = FDJac_Q_L[0];
        
        // -- Q2 ---------------------------------------------------------------
        // Now Perturb the Left Node Q2 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q2[node_L] = FDJac_Q_L[1] + eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        // Now Perturb the Left Node Q2 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q2[node_L] = FDJac_Q_L[1] - eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][1] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q2
        Q2[node_L] = FDJac_Q_L[1];
        
        // -- Q3 ---------------------------------------------------------------
        // Now Perturb the Left Node Q3 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q3[node_L] = FDJac_Q_L[2] + eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        // Now Perturb the Left Node Q3 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q3[node_L] = FDJac_Q_L[2] - eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][2] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q3
        Q3[node_L] = FDJac_Q_L[2];
        
        // -- Q4 ---------------------------------------------------------------
        // Now Perturb the Left Node Q4 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q4[node_L] = FDJac_Q_L[3] + eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        // Now Perturb the Left Node Q4 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q4[node_L] = FDJac_Q_L[3] - eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][3] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q4
        Q4[node_L] = FDJac_Q_L[3];
        
        // -- Q5 ---------------------------------------------------------------
        // Now Perturb the Left Node Q5 +eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_FORWARD)) { // Central or Forward
            Q5[node_L] = FDJac_Q_L[4] + eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxP_Conv, FDJac_FluxP_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        // Now Perturb the Left Node Q5 -eps
        if ((FiniteDifference == JACOBIAN_METHOD_CENTRAL) || (FiniteDifference == JACOBIAN_METHOD_BACKWARD)) { // Central or Backward
            Q5[node_L] = FDJac_Q_L[4] - eps;
            // Update the Boundary Condition Ghost Node 
            Apply_Boundary_Condition(ibEdge, Iteration);
            Compute_Flux(node_L, node_R, areavec, FDJac_FluxM_Conv, FDJac_FluxM_Diss, FALSE);
            // Reset the Ghost Node Values
            Q1[node_R] = FDJac_Q_R[0];
            Q2[node_R] = FDJac_Q_R[1];
            Q3[node_R] = FDJac_Q_R[2];
            Q4[node_R] = FDJac_Q_R[3];
            Q5[node_R] = FDJac_Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            FDJac_dFdL[i][4] = Coeff*((FDJac_FluxP_Conv[i] + FDJac_FluxP_Diss[i]) - (FDJac_FluxM_Conv[i] + FDJac_FluxM_Diss[i]))/eps;
        // Reset the Q5
        Q5[node_L] = FDJac_Q_L[4];
        
        // Update the Diagonal Term of Physical Node Only
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                SolverBlockMatrix.A[idgnL][i][j] += FDJac_dFdL[i][j];
        }
    }
    
    // Reset Flux Recompute Flag
    CogSolver.FluxRecomputeFlag = FALSE;
}

