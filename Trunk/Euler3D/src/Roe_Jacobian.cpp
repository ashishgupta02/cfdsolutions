/*******************************************************************************
 * File:        Roe_Jacobian.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "CompressibleUtils.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_Exact_Roe(int AddTime, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//! Note: Jacobian is computed using first order Q's
//! Note: This Routine Does not Work
//------------------------------------------------------------------------------
void Compute_Jacobian_Approximate_Roe(int AddTime, int Iteration) {
    int i, j, k, iNode, iEdge, ibEdge;
    int node_L, node_R;
    int idgn, idgnL, idgnR, ofdgnL, ofdgnR;
    double area;
    double Q_L[NEQUATIONS];
    double Q_R[NEQUATIONS];
    double **dFdL;
    double **dFdR;
    double **ARoe;
    Vector3D areavec;
    
    // Create the Helper Matrix
    dFdL = (double **) malloc(NEQUATIONS*sizeof(double*));
    dFdR = (double **) malloc(NEQUATIONS*sizeof(double*));
    ARoe = (double **) malloc(NEQUATIONS*sizeof(double*));
    for (i = 0; i < NEQUATIONS; i++) {
        dFdL[i] = (double *) malloc(NEQUATIONS*sizeof(double));
        dFdR[i] = (double *) malloc(NEQUATIONS*sizeof(double));
        ARoe[i] = (double *) malloc(NEQUATIONS*sizeof(double));
    }
    
    // Initialize the CRS Matrix
    for (i = 0; i < SolverBlockMatrix.DIM; i++) {
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                SolverBlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Copy the Residuals to Block Matrix which is B
    // And Copy I/DeltaT
    for (iNode = 0; iNode < nNode; iNode++) {
        // Get the diagonal location
        idgn = SolverBlockMatrix.IAU[iNode];
        // Get the LHS
        SolverBlockMatrix.B[iNode][0] = -Res1[iNode];
        SolverBlockMatrix.B[iNode][1] = -Res2[iNode];
        SolverBlockMatrix.B[iNode][2] = -Res3[iNode];
        SolverBlockMatrix.B[iNode][3] = -Res4[iNode];
        SolverBlockMatrix.B[iNode][4] = -Res5[iNode];
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            if (AddTime == 1) {
                for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                    if (k == j)
                        SolverBlockMatrix.A[idgn][j][k] = cVolume[iNode]/DeltaT[iNode];
            }
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
        
        // Backup the Copy of Q_L and Q_R
        // Left
        Q_L[0] = Q1[node_L];
        Q_L[1] = Q2[node_L];
        Q_L[2] = Q3[node_L];
        Q_L[3] = Q4[node_L];
        Q_L[4] = Q5[node_L];
        // Right
        Q_R[0] = Q1[node_R];
        Q_R[1] = Q2[node_R];
        Q_R[2] = Q3[node_R];
        Q_R[3] = Q4[node_R];
        Q_R[4] = Q5[node_R];
        
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
                dFdL[i][j] = 0.0;
                dFdR[i][j] = 0.0;
                ARoe[i][j] = 0.0;
            }
        }
        
        // Compute dFdL
        ConservativeEulerFluxJacobian(Q_L, areavec, dFdL, Gamma);
        
        // Compute dFdR
        ConservativeEulerFluxJacobian(Q_R, areavec, dFdR, Gamma);
        
        // Compute Roe Jacobian
        Compute_RoeAJacobian(Q_L, Q_R, areavec, ARoe);
        
        // Finally Compute the dFlux/dQ_L and dFlux/dQ_R
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                dFdL[i][j] = 0.5*(dFdL[i][j] - ARoe[i][j])*area;
                dFdR[i][j] = 0.5*(dFdR[i][j] + ARoe[i][j])*area;
            }
        }
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                // Diagonal
                SolverBlockMatrix.A[idgnL][i][j] += dFdL[i][j];
                SolverBlockMatrix.A[idgnR][i][j] -= dFdR[i][j];
                // Off-Diagonal
                SolverBlockMatrix.A[ofdgnL][i][j] += dFdR[i][j];
                SolverBlockMatrix.A[ofdgnR][i][j] -= dFdL[i][j];
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
        
        // Backup the Copy of Q_L and Q_R
        // Left - Physical
        Q_L[0] = Q1[node_L];
        Q_L[1] = Q2[node_L];
        Q_L[2] = Q3[node_L];
        Q_L[3] = Q4[node_L];
        Q_L[4] = Q5[node_L];
        // Right - Ghost
        Q_R[0] = Q1[node_R];
        Q_R[1] = Q2[node_R];
        Q_R[2] = Q3[node_R];
        Q_R[3] = Q4[node_R];
        Q_R[4] = Q5[node_R];
        
        // Get the diagonal Locations - Only Physical
        // No diagonal and off-diagonal locations exists for Ghost Nodes
        idgnL = SolverBlockMatrix.IAU[node_L];
        
        // Initialize the Helper Matrix - Only Physical
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                dFdL[i][j] = 0.0;
                ARoe[i][j] = 0.0;
            }
        }
        
        // Compute dFdL
        ConservativeEulerFluxJacobian(Q_L, areavec, dFdL, Gamma);
        
        // Compute Roe Jacobian
        Compute_RoeAJacobian(Q_L, Q_R, areavec, ARoe);
        
        // Finally Compute the dFlux/dQ_L
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                dFdL[i][j] = 0.5*(dFdL[i][j] - ARoe[i][j])*area;
        }
        
        // Update the Diagonal Term of Physical Node Only
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                SolverBlockMatrix.A[idgnL][i][j] += dFdL[i][j];
        }
    }
    
    // Delete the Helper Matrix
    for (i = 0; i < NEQUATIONS; i++) {
        free(ARoe[i]);
        free(dFdL[i]);
        free(dFdR[i]);
    }
    free(ARoe);
    free(dFdL);
    free(dFdR);
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//! Note Currently Internal Edges works : Forward, Backward or Central
//! Boundary Edges Only Central -- Some more research needed.
//! Note: Jacobian is computed using first order Q's
//------------------------------------------------------------------------------
void Compute_Jacobian_FiniteDifference_Roe(int AddTime, int Iteration) {
    int i, j, k, iNode, iEdge, ibEdge;
    int node_L, node_R;
    int idgn, idgnL, idgnR, ofdgnL, ofdgnR;
    double Q_L[NEQUATIONS];
    double Q_R[NEQUATIONS];
    double FluxP[NEQUATIONS];
    double FluxM[NEQUATIONS];
    double dFdL[NEQUATIONS][NEQUATIONS];
    double dFdR[NEQUATIONS][NEQUATIONS];
    Vector3D areavec;
    double eps = 1.0E-8;
    int FiniteDifference = JacobianMethod;
    double Coeff = 0.5;
    
    // Initialize the CRS Matrix
    for (i = 0; i < SolverBlockMatrix.DIM; i++) {
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                SolverBlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Copy the Residuals to Block Matrix which is B
    // And Copy I/DeltaT
    for (iNode = 0; iNode < nNode; iNode++) {
        // Get the diagonal location
        idgn = SolverBlockMatrix.IAU[iNode];
        // Get the LHS
        SolverBlockMatrix.B[iNode][0] = -Res1[iNode];
        SolverBlockMatrix.B[iNode][1] = -Res2[iNode];
        SolverBlockMatrix.B[iNode][2] = -Res3[iNode];
        SolverBlockMatrix.B[iNode][3] = -Res4[iNode];
        SolverBlockMatrix.B[iNode][4] = -Res5[iNode];
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            if (AddTime == 1) {
                for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                    if (k == j)
                        SolverBlockMatrix.A[idgn][j][k] = cVolume[iNode]/DeltaT[iNode];
            }
        }
    }
    
    // Check if Alternating
    if (FiniteDifference == SOLVER_JACOBIAN_ALTERNATE) {
        if (Iteration%2 == 0)
            FiniteDifference = SOLVER_JACOBIAN_FORWARD;
        else
            FiniteDifference = SOLVER_JACOBIAN_BACKWARD;
    }
    
    // Set the value of Coeff based on Finite Difference Method
    if (FiniteDifference == SOLVER_JACOBIAN_CENTRAL) // Central
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
        Q_L[0] = Q1[node_L];
        Q_L[1] = Q2[node_L];
        Q_L[2] = Q3[node_L];
        Q_L[3] = Q4[node_L];
        Q_L[4] = Q5[node_L];
        // Right
        Q_R[0] = Q1[node_R];
        Q_R[1] = Q2[node_R];
        Q_R[2] = Q3[node_R];
        Q_R[3] = Q4[node_R];
        Q_R[4] = Q5[node_R];
        
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
                dFdL[i][j] = 0.0;
                dFdR[i][j] = 0.0;
            }
        }
        
        // ---------------------------------------------------------------------
        // - Left Node ---------------------------------------------------------
        // -- Q1 ---------------------------------------------------------------
        // Now Perturb the Left Node Q1 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) // Central or Forward
            Q1[node_L] = Q_L[0] + eps;
        else
            Q1[node_L] = Q_L[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        
        // Now Perturb the Left Node Q1 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD))  // Central or Backward
            Q1[node_L] = Q_L[0] - eps;
        else
            Q1[node_L] = Q_L[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][0] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q1
        Q1[node_L] = Q_L[0];
        
        // -- Q2 ---------------------------------------------------------------
        // Now Perturb the Left Node Q2 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q2[node_L] = Q_L[1] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Left Node Q2 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q2[node_L] = Q_L[1] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][1] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q2
        Q2[node_L] = Q_L[1];
        
        // -- Q3 ---------------------------------------------------------------
        // Now Perturb the Left Node Q3 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q3[node_L] = Q_L[2] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Left Node Q3 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q3[node_L] = Q_L[2] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][2] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q3
        Q3[node_L] = Q_L[2];
        
        // -- Q4 ---------------------------------------------------------------
        // Now Perturb the Left Node Q4 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q4[node_L] = Q_L[3] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Left Node Q4 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q4[node_L] = Q_L[3] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][3] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q4
        Q4[node_L] = Q_L[3];
        
        // -- Q5 ---------------------------------------------------------------
        // Now Perturb the Left Node Q5 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q5[node_L] = Q_L[4] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Left Node Q5 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q5[node_L] = Q_L[4] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][4] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q5
        Q5[node_L] = Q_L[4];
        
        // ---------------------------------------------------------------------
        // - Right Node --------------------------------------------------------
        // ---------------------------------------------------------------------
        // -- Q1 --
        // Now Perturb the Right Node Q1 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) // Central or Forward
            Q1[node_R] = Q_R[0] + eps;
        else
            Q1[node_R] = Q_R[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        
        // Now Perturb the Right Node Q1 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) // Central or Backward
            Q1[node_R] = Q_R[0] - eps;
        else
            Q1[node_R] = Q_R[0];
        // Based on Finite Difference Method, One Sided Do only once
        Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdR[i][0] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q1
        Q1[node_R] = Q_R[0];
        
        // -- Q2 ---------------------------------------------------------------
        // Now Perturb the Right Node Q2 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q2[node_R] = Q_R[1] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Right Node Q2 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q2[node_R] = Q_R[1] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdR[i][1] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q2
        Q2[node_R] = Q_R[1];
        
        // -- Q3 ---------------------------------------------------------------
        // Now Perturb the Right Node Q3 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q3[node_R] = Q_R[2] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Right Node Q3 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q3[node_R] = Q_R[2] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdR[i][2] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q3
        Q3[node_R] = Q_R[2];
        
        // -- Q4 ---------------------------------------------------------------
        // Now Perturb the Right Node Q4 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q4[node_R] = Q_R[3] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Right Node Q4 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q4[node_R] = Q_R[3] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdR[i][3] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q4
        Q4[node_R] = Q_R[3];
        
        // -- Q5 ---------------------------------------------------------------
        // Now Perturb the Right Node Q5 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q5[node_R] = Q_R[4] + eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        }
        
        // Now Perturb the Right Node Q5 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q5[node_R] = Q_R[4] - eps;
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdR[i][4] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q5
        Q5[node_R] = Q_R[4];
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                // Diagonal
                SolverBlockMatrix.A[idgnL][i][j] += dFdL[i][j];
                SolverBlockMatrix.A[idgnR][i][j] -= dFdR[i][j];
                // Off-Diagonal
                SolverBlockMatrix.A[ofdgnL][i][j] += dFdR[i][j];
                SolverBlockMatrix.A[ofdgnR][i][j] -= dFdL[i][j];
            }
        }
    }
    
    // Note: Currently Forcing Boundary Contribution as Central
    // More Research has to be done for Forward and Backward to work here
    FiniteDifference = SOLVER_JACOBIAN_CENTRAL;
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
        Q_L[0] = Q1[node_L];
        Q_L[1] = Q2[node_L];
        Q_L[2] = Q3[node_L];
        Q_L[3] = Q4[node_L];
        Q_L[4] = Q5[node_L];
        // Right - Ghost
        Q_R[0] = Q1[node_R];
        Q_R[1] = Q2[node_R];
        Q_R[2] = Q3[node_R];
        Q_R[3] = Q4[node_R];
        Q_R[4] = Q5[node_R];
        
        // Get the diagonal Locations - Only Physical
        // No diagonal and off-diagonal locations exists for Ghost Nodes
        idgnL = SolverBlockMatrix.IAU[node_L];
        
        // Initialize the Helper Matrix - Only Physical
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                dFdL[i][j] = 0.0;
        }
        
        // ---------------------------------------------------------------------
        // Left Node - Physical ------------------------------------------------
        // -- Q1 ---------------------------------------------------------------
        // Now Perturb the Left Node Q1 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q1[node_L] = Q_L[0] + eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
        } else {
            Q1[node_L] = Q_L[0];
        }
        // Based on Finite Difference Method, One Sided Do only once
        Compute_RoeFlux(node_L, node_R, areavec, FluxP);
        // Reset the Ghost Node Values
        Q1[node_R] = Q_R[0];
        Q2[node_R] = Q_R[1];
        Q3[node_R] = Q_R[2];
        Q4[node_R] = Q_R[3];
        Q5[node_R] = Q_R[4];
        
        // Now Perturb the Left Node Q1 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q1[node_L] = Q_L[0] - eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
        } else {
            Q1[node_L] = Q_L[0];
        }
        // Based on Finite Difference Method, One Sided Do only once
        Compute_RoeFlux(node_L, node_R, areavec, FluxM);
        // Reset the Ghost Node Values
        Q1[node_R] = Q_R[0];
        Q2[node_R] = Q_R[1];
        Q3[node_R] = Q_R[2];
        Q4[node_R] = Q_R[3];
        Q5[node_R] = Q_R[4];
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][0] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q1
        Q1[node_L] = Q_L[0];
        
        // -- Q2 ---------------------------------------------------------------
        // Now Perturb the Left Node Q2 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q2[node_L] = Q_L[1] + eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        // Now Perturb the Left Node Q2 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q2[node_L] = Q_L[1] - eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][1] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q2
        Q2[node_L] = Q_L[1];
        
        // -- Q3 ---------------------------------------------------------------
        // Now Perturb the Left Node Q3 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q3[node_L] = Q_L[2] + eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        // Now Perturb the Left Node Q3 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q3[node_L] = Q_L[2] - eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][2] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q3
        Q3[node_L] = Q_L[2];
        
        // -- Q4 ---------------------------------------------------------------
        // Now Perturb the Left Node Q4 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q4[node_L] = Q_L[3] + eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        // Now Perturb the Left Node Q4 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q4[node_L] = Q_L[3] - eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][3] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q4
        Q4[node_L] = Q_L[3];
        
        // -- Q5 ---------------------------------------------------------------
        // Now Perturb the Left Node Q5 +eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_FORWARD)) { // Central or Forward
            Q5[node_L] = Q_L[4] + eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxP);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        // Now Perturb the Left Node Q5 -eps
        if ((FiniteDifference == SOLVER_JACOBIAN_CENTRAL) || (FiniteDifference == SOLVER_JACOBIAN_BACKWARD)) { // Central or Backward
            Q5[node_L] = Q_L[4] - eps;
            // Update the Boundary Condition Ghost Node 
            Compute_Characteristic_BoundaryCondition(ibEdge, Iteration);
            Compute_RoeFlux(node_L, node_R, areavec, FluxM);
            // Reset the Ghost Node Values
            Q1[node_R] = Q_R[0];
            Q2[node_R] = Q_R[1];
            Q3[node_R] = Q_R[2];
            Q4[node_R] = Q_R[3];
            Q5[node_R] = Q_R[4];
        }
        
        for (i = 0; i < NEQUATIONS; i++)
            dFdL[i][4] = Coeff*(FluxP[i] - FluxM[i])/eps;
        // Reset the Q5
        Q5[node_L] = Q_L[4];
        
        // Update the Diagonal Term of Physical Node Only
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                SolverBlockMatrix.A[idgnL][i][j] += dFdL[i][j];
        }
    }
}

//------------------------------------------------------------------------------
//! Computes the Jacobian for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Jacobian_Roe(int AddTime, int Iteration) {
    switch (JacobianMethod) {
        case SOLVER_JACOBIAN_CENTRAL:
            Compute_Jacobian_FiniteDifference_Roe(AddTime, Iteration);
            break;
        case SOLVER_JACOBIAN_FORWARD:
            Compute_Jacobian_FiniteDifference_Roe(AddTime, Iteration);
            break;
        case SOLVER_JACOBIAN_BACKWARD:
            Compute_Jacobian_FiniteDifference_Roe(AddTime, Iteration);
            break;
        case SOLVER_JACOBIAN_ALTERNATE:
            Compute_Jacobian_FiniteDifference_Roe(AddTime, Iteration);
            break;
        case SOLVER_JACOBIAN_APPROX:
            Compute_Jacobian_Approximate_Roe(AddTime, Iteration);
            break;
        case SOLVER_JACOBIAN_EXACT:
            Compute_Jacobian_Exact_Roe(AddTime, Iteration);
            break;
        default:
            error("Compute_Jacobian_Roe: Invalid Jacobian Method - %d", JacobianMethod);
            break;
    }
}

