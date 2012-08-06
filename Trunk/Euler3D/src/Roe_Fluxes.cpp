/*******************************************************************************
 * File:        Roe_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 * Reference:   A Low-Mach Number Fix for Roe's Approximate Riemann Solver
 *              Felix Rieper, Journal of Computational Physics, 230 (2011)
 ******************************************************************************/

#include "License.h"

#ifdef DEBUG
#include <assert.h>
#endif

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Material.h"
#include "Solver.h"
#include "Residual_Smoothing.h"

// Static Variable for Speed Up
static int     Roe_DB           = 0;
static double *Roe_fluxA        = NULL;
static double *Roe_flux_L       = NULL;
static double *Roe_flux_R       = NULL;
static double *Roe_Q_L          = NULL;
static double *Roe_Q_R          = NULL;
static double *Roe_dw           = NULL;
static double *Roe_dQ           = NULL;
static double **Roe_A           = NULL;
static double **Roe_Eigen       = NULL;
static double **Roe_M           = NULL;
static double **Roe_Minv        = NULL;
static double **Roe_P           = NULL;
static double **Roe_Pinv        = NULL;
static double **Roe_T           = NULL;
static double **Roe_Tinv        = NULL;
    
//------------------------------------------------------------------------------
//! Create Roe Scheme Data Structure
//------------------------------------------------------------------------------
void Roe_Init(void) {
    int i;
    
    // Check if Roe Data Structure is required
    if (Roe_DB == 0) {
        Roe_fluxA  = (double *)  malloc(NEQUATIONS*sizeof(double));
        Roe_flux_L = (double *)  malloc(NEQUATIONS*sizeof(double));
        Roe_flux_R = (double *)  malloc(NEQUATIONS*sizeof(double));
        Roe_Q_L    = (double *)  malloc(NEQUATIONS*sizeof(double));
        Roe_Q_R    = (double *)  malloc(NEQUATIONS*sizeof(double));
        Roe_dw     = (double *)  malloc(NEQUATIONS*sizeof(double));
        Roe_dQ     = (double *)  malloc(NEQUATIONS*sizeof(double));
        Roe_A      = (double **) malloc(NEQUATIONS*sizeof(double*));
        Roe_Eigen  = (double **) malloc(NEQUATIONS*sizeof(double*));
        Roe_M      = (double **) malloc(NEQUATIONS*sizeof(double*));
        Roe_Minv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        Roe_P      = (double **) malloc(NEQUATIONS*sizeof(double*));
        Roe_Pinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        Roe_T      = (double **) malloc(NEQUATIONS*sizeof(double*));
        Roe_Tinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        for (i = 0; i < NEQUATIONS; i++) {
            Roe_A[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            Roe_Eigen[i] = (double *) malloc(NEQUATIONS*sizeof(double));
            Roe_M[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            Roe_Minv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
            Roe_P[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            Roe_Pinv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
            Roe_T[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            Roe_Tinv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
        }
        Roe_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete Roe Scheme Data Structure
//------------------------------------------------------------------------------
void Roe_Finalize(void) {
    int i;
    
    // Free the Memory
    for (i = 0; i < NEQUATIONS; i++) {
        free(Roe_A[i]);
        free(Roe_Eigen[i]);
        free(Roe_M[i]);
        free(Roe_Minv[i]);
        free(Roe_P[i]);
        free(Roe_Pinv[i]);
        free(Roe_T[i]);
        free(Roe_Tinv[i]);
    }
    free(Roe_fluxA);
    free(Roe_flux_L);
    free(Roe_flux_R);
    free(Roe_Q_L);
    free(Roe_Q_R);
    free(Roe_dw);
    free(Roe_dQ);
    free(Roe_A);
    free(Roe_Eigen);
    free(Roe_M);
    free(Roe_Minv);
    free(Roe_P);
    free(Roe_Pinv);
    free(Roe_T);
    free(Roe_Tinv);
    Roe_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset Roe Scheme Data Structure
//------------------------------------------------------------------------------
void Roe_Reset(void) {
    int i, j;
    
    if (Roe_DB == 0)
        Roe_Init();
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        Roe_fluxA[i]  = 0.0;
        Roe_flux_L[i] = 0.0;
        Roe_flux_R[i] = 0.0;
        Roe_Q_L[i]    = 0.0;
        Roe_Q_R[i]    = 0.0;
        Roe_dw[i]     = 0.0;
        Roe_dQ[i]     = 0.0;
        for (j = 0; j < NEQUATIONS; j++) {
            Roe_A[i][j]     = 0.0;
            Roe_Eigen[i][j] = 0.0;
            Roe_M[i][j]     = 0.0;
            Roe_Minv[i][j]  = 0.0;
            Roe_P[i][j]     = 0.0;
            Roe_Pinv[i][j]  = 0.0;
            Roe_T[i][j]     = 0.0;
            Roe_Tinv[i][j]  = 0.0;
        }
    }
}

//------------------------------------------------------------------------------
//! Compute Roe Flux Jacobian A
//! Note this function should not be called from Compute_RoeFlux
//------------------------------------------------------------------------------
void Compute_RoeAJacobian(double *Q_L, double *Q_R, Vector3D areavec, double **AJacobian_Roe) {
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, nx, ny, nz;
    
    // Initialization
    Roe_Reset();
    
    // Get area vector
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    // Compute Equation of State
    Compute_EOS_Variables_Face(Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
    Compute_EOS_Variables_Face(Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
    
    // ROE AVERAGE VARIABLES
    rho   = sqrt(rho_R * rho_L);
    sigma = rho/(rho_L + rho);
    u     = u_L  + sigma*(u_R  - u_L);
    v     = v_L  + sigma*(v_R  - v_L);
    w     = w_L  + sigma*(w_R  - w_L);
    ht    = ht_L + sigma*(ht_R - ht_L);
    phi   = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
    c     = 0.0;
    switch (NonDimensional_Type) {
        case NONDIMENSIONAL_GENERIC:
            c = (Gamma - 1.0)*ht - phi;
            c = sqrt(c);
            break;
        case NONDIMENSIONAL_BTW:
            c = ht/(Ref_Mach*Ref_Mach) - phi;
            c = sqrt(c);
            break;
        case NONDIMENSIONAL_LMROE:
            c = (Gamma - 1.0)*ht - phi;
            c = sqrt(c);
            break;
        default:
            error("Compute_RoeAJacobian: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
            break;
    }
    ubar  = u*nx + v*ny + w*nz;
    
    // M
    Roe_M[0][0] = 1.0;
    Roe_M[0][1] = 0.0;
    Roe_M[0][2] = 0.0;
    Roe_M[0][3] = 0.0;
    Roe_M[0][4] = 0.0;

    Roe_M[1][0] = u;
    Roe_M[1][1] = rho;
    Roe_M[1][2] = 0.0;
    Roe_M[1][3] = 0.0;
    Roe_M[1][4] = 0.0;

    Roe_M[2][0] = v;
    Roe_M[2][1] = 0.0;
    Roe_M[2][2] = rho;
    Roe_M[2][3] = 0.0;
    Roe_M[2][4] = 0.0;

    Roe_M[3][0] = w;
    Roe_M[3][1] = 0.0;
    Roe_M[3][2] = 0.0;
    Roe_M[3][3] = rho;
    Roe_M[3][4] = 0.0;

    switch (NonDimensional_Type) {
        case NONDIMENSIONAL_GENERIC:
            Roe_M[4][0] = phi/(Gamma - 1.0);
            Roe_M[4][1] = rho * u;
            Roe_M[4][2] = rho * v;
            Roe_M[4][3] = rho * w;
            Roe_M[4][4] = 1.0/(Gamma - 1.0);
            break;
        case NONDIMENSIONAL_BTW:
            Roe_M[4][0] = phi * Ref_Mach * Ref_Mach;
            Roe_M[4][1] = (Gamma - 1.0) * rho * u * Ref_Mach * Ref_Mach;
            Roe_M[4][2] = (Gamma - 1.0) * rho * v * Ref_Mach * Ref_Mach;
            Roe_M[4][3] = (Gamma - 1.0) * rho * w * Ref_Mach * Ref_Mach;
            Roe_M[4][4] = Ref_Mach * Ref_Mach;
            break;
        case NONDIMENSIONAL_LMROE:
            Roe_M[4][0] = phi/(Gamma - 1.0);
            Roe_M[4][1] = rho * u;
            Roe_M[4][2] = rho * v;
            Roe_M[4][3] = rho * w;
            Roe_M[4][4] = 1.0/(Gamma - 1.0);
            break;
        default:
            error("Compute_RoeAJacobian: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
            break;
    }
    
    // Minv
    Roe_Minv[0][0] = 1.0;
    Roe_Minv[0][1] = 0.0;
    Roe_Minv[0][2] = 0.0;
    Roe_Minv[0][3] = 0.0;
    Roe_Minv[0][4] = 0.0;

    Roe_Minv[1][0] = -u/rho;
    Roe_Minv[1][1] = 1.0/rho;
    Roe_Minv[1][2] = 0.0;
    Roe_Minv[1][3] = 0.0;
    Roe_Minv[1][4] = 0.0;

    Roe_Minv[2][0] = -v/rho;
    Roe_Minv[2][1] = 0.0;
    Roe_Minv[2][2] = 1.0/rho;
    Roe_Minv[2][3] = 0.0;
    Roe_Minv[2][4] = 0.0;

    Roe_Minv[3][0] = -w/rho;
    Roe_Minv[3][1] = 0.0;
    Roe_Minv[3][2] = 0.0;
    Roe_Minv[3][3] = 1.0/rho;
    Roe_Minv[3][4] = 0.0;

    switch (NonDimensional_Type) {
        case NONDIMENSIONAL_GENERIC:
            Roe_Minv[4][0] = phi;
            Roe_Minv[4][1] = -u * (Gamma - 1.0);
            Roe_Minv[4][2] = -v * (Gamma - 1.0);
            Roe_Minv[4][3] = -w * (Gamma - 1.0);
            Roe_Minv[4][4] = (Gamma - 1.0);
            break;
        case NONDIMENSIONAL_BTW:
            Roe_Minv[4][0] = phi;
            Roe_Minv[4][1] = -u * (Gamma - 1.0);
            Roe_Minv[4][2] = -v * (Gamma - 1.0);
            Roe_Minv[4][3] = -w * (Gamma - 1.0);
            Roe_Minv[4][4] = 1.0/(Ref_Mach * Ref_Mach);
            break;
        case NONDIMENSIONAL_LMROE:
            Roe_Minv[4][0] = phi;
            Roe_Minv[4][1] = -u * (Gamma - 1.0);
            Roe_Minv[4][2] = -v * (Gamma - 1.0);
            Roe_Minv[4][3] = -w * (Gamma - 1.0);
            Roe_Minv[4][4] = (Gamma - 1.0);
            break;
        default:
            error("Compute_RoeAJacobian: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
            break;
    }

    // P
    Roe_P[0][0] = nx;
    Roe_P[0][1] = ny;
    Roe_P[0][2] = nz;
    Roe_P[0][3] = rho/c;
    Roe_P[0][4] = rho/c;

    Roe_P[1][0] = 0.0;
    Roe_P[1][1] = -nz;
    Roe_P[1][2] = ny;
    Roe_P[1][3] = nx;
    Roe_P[1][4] = -nx;

    Roe_P[2][0] = nz;
    Roe_P[2][1] = 0.0;
    Roe_P[2][2] = -nx;
    Roe_P[2][3] = ny;
    Roe_P[2][4] = -ny;

    Roe_P[3][0] = -ny;
    Roe_P[3][1] = nx;
    Roe_P[3][2] = 0.0;
    Roe_P[3][3] = nz;
    Roe_P[3][4] = -nz;

    Roe_P[4][0] = 0.0;
    Roe_P[4][1] = 0.0;
    Roe_P[4][2] = 0.0;
    Roe_P[4][3] = rho * c;
    Roe_P[4][4] = rho * c;

    // Pinv
    Roe_Pinv[0][0] = nx;
    Roe_Pinv[0][1] = 0.0;
    Roe_Pinv[0][2] = nz;
    Roe_Pinv[0][3] = -ny;
    Roe_Pinv[0][4] = -nx/(c * c);

    Roe_Pinv[1][0] = ny;
    Roe_Pinv[1][1] = -nz;
    Roe_Pinv[1][2] = 0.0;
    Roe_Pinv[1][3] = nx;
    Roe_Pinv[1][4] = -ny/(c * c);

    Roe_Pinv[2][0] = nz;
    Roe_Pinv[2][1] = ny;
    Roe_Pinv[2][2] = -nx;
    Roe_Pinv[2][3] = 0.0;
    Roe_Pinv[2][4] = -nz/(c * c);

    Roe_Pinv[3][0] = 0.0;
    Roe_Pinv[3][1] = 0.5 * nx;
    Roe_Pinv[3][2] = 0.5 * ny;
    Roe_Pinv[3][3] = 0.5 * nz;
    Roe_Pinv[3][4] = 1.0/(2.0 * rho * c);

    Roe_Pinv[4][0] = 0.0;
    Roe_Pinv[4][1] = -0.5 * nx;
    Roe_Pinv[4][2] = -0.5 * ny;
    Roe_Pinv[4][3] = -0.5 * nz;
    Roe_Pinv[4][4] = 1.0/(2.0 * rho * c);

    // START: Computing 
    // Roe: A = M*P*|Lambda|*Pinv*Minv

    // Calculate T = M*P
    MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_P, Roe_T);

    // Calculate EigenMatrix |Lambda|
    // Apply Entropy Fix
    if (EntropyFix != 0) {
        Roe_EntropyFix(ubar_L, c_L, ubar_R, c_R, ubar, c, Roe_Eigen);
    } else {
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar + c);
        Roe_Eigen[4][4] = fabs(ubar - c);
    }
    
    // Calculate Tinv = Pinv*Minv
    MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);

    // EigenMatrix*Tinv = |Lambda|*Pinv*Minv
    MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);

    // Get Matrix A = M*P*|Lambda|*Pinv*Minv
    MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, AJacobian_Roe);
    // END: Computing A = M*P*|Lambda|*Pinv*Minv
}

//------------------------------------------------------------------------------
//! Compute Roe Flux Optimized
//------------------------------------------------------------------------------
void Compute_RoeFlux_Optimized(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar, q2, cinv, gm1oc;
    double l1, l2, l3, l4, l5;
    double r1, r2, r3, r4, r5;
    double sigma, alpha, maxlambda;
    
    double area;
    double nx, ny, nz;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        Flux_Roe_Conv[i] = 0.0;
        Flux_Roe_Diss[i] = 0.0;
    }
    Roe_Reset();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Internal Nodes
        if (node_R < nNode) {
            // Make Solution Second Order
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, Roe_Q_L, Roe_Q_R);
            } else {
                Roe_Q_L[0] = Q1[node_L];
                Roe_Q_L[1] = Q2[node_L];
                Roe_Q_L[2] = Q3[node_L];
                Roe_Q_L[3] = Q4[node_L];
                Roe_Q_L[4] = Q5[node_L];
                Roe_Q_R[0] = Q1[node_R];
                Roe_Q_R[1] = Q2[node_R];
                Roe_Q_R[2] = Q3[node_R];
                Roe_Q_R[3] = Q4[node_R];
                Roe_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            Roe_Q_L[0] = Q1[node_L];
            Roe_Q_L[1] = Q2[node_L];
            Roe_Q_L[2] = Q3[node_L];
            Roe_Q_L[3] = Q4[node_L];
            Roe_Q_L[4] = Q5[node_L];
            Roe_Q_R[0] = Q1[node_R];
            Roe_Q_R[1] = Q2[node_R];
            Roe_Q_R[2] = Q3[node_R];
            Roe_Q_R[3] = Q4[node_R];
            Roe_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        q2    = u*u + v*v + w*w;
        phi   = 0.5*(Gamma - 1.0)*q2;
        c     = (Gamma - 1.0)*ht - phi;
        c     = sqrt(c);
        ubar  = u*nx + v*ny + w*nz;
        cinv  = 1.0/c;
        gm1oc = (Gamma - 1.0)*cinv;
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Min and Max EigenValue Value
            MinEigenLamda1 = MIN(MinEigenLamda1, fabs(ubar));
            MaxEigenLamda1 = MAX(MaxEigenLamda1, fabs(ubar));
            MinEigenLamda4 = MIN(MinEigenLamda4, fabs(ubar + c));
            MaxEigenLamda4 = MAX(MaxEigenLamda4, fabs(ubar + c));
            MinEigenLamda5 = MIN(MinEigenLamda5, fabs(ubar - c));
            MaxEigenLamda5 = MAX(MaxEigenLamda5, fabs(ubar - c));

            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(fabs(ubar), MAX(fabs(ubar + c), fabs(ubar - c)));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothType == RESIDUAL_SMOOTH_IMPLICIT) {
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L + 1]; i++) {
                // Get the Node ID and Match with node_R
                if (crs_JA_Node2Node[i] == node_R) {
                    RM_MaxEigenValue[i] = maxlambda*area;
                    break;
                }
            }
            RM_SumMaxEigenValue[node_L] += maxlambda*area;
            // Only Physical Nodes
            if (node_R < nNode) {
                for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R + 1]; i++) {
                    // Get the Node ID and Match with node_L
                    if (crs_JA_Node2Node[i] == node_L) {
                        RM_MaxEigenValue[i] = maxlambda*area;
                        break;
                    }
                }
                RM_SumMaxEigenValue[node_R] += maxlambda*area;
            }
        }
        
        // Wave Traveling From Left to Right
        if (ubar > 0.0) {
            // Compute Convective Flux: flux_L 
            Flux_Roe_Conv[0] = rho_L*ubar_L;
            Flux_Roe_Conv[1] =   u_L*Flux_Roe_Conv[0] + p_L*nx;
            Flux_Roe_Conv[2] =   v_L*Flux_Roe_Conv[0] + p_L*ny;
            Flux_Roe_Conv[3] =   w_L*Flux_Roe_Conv[0] + p_L*nz;
            Flux_Roe_Conv[4] =  ht_L*Flux_Roe_Conv[0];
            
            // Compute the Dissipation Flux
            // Subsonic
            if (ubar < c) {
                // Eigenvalues 1, 2, 3, 4 are positive, 5 is negative
                // So, compute form the left with the one eigenvector that comes into play
                
                // Compute dQ
                // Conservative Variable Formulation
                if (Variable_Type == VARIABLE_CONSERVATIVE) {
                    Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
                    Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
                    Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
                    Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
                    Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
                }
                // Primitive Variable Formulation
                if ((Variable_Type == VARIABLE_PRIMITIVE_PUT) || (Variable_Type == VARIABLE_PRIMITIVE_RUP)) {
                    Roe_dQ[0] = rho_R      - rho_L;
                    Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
                    Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
                    Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
                    Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
                }
                
                // Eigenvalues
                // Apply Entropy Fix
                if (EntropyFix != 0) {
                    Roe_EntropyFix(ubar_L, c_L, ubar_R, c_R, ubar, c, Roe_Eigen);
                } else {
                    Roe_Eigen[0][0] = fabs(ubar);
                    Roe_Eigen[1][1] = Roe_Eigen[0][0];
                    Roe_Eigen[2][2] = Roe_Eigen[0][0];
                    Roe_Eigen[3][3] = fabs(ubar + c);
                    Roe_Eigen[4][4] = fabs(ubar - c);
                }
                Roe_Eigen[4][4] = ISIGN(ubar - c)*Roe_Eigen[4][4];
        
                l1 = ubar + 0.5*gm1oc*q2;
                l2 = -nx - gm1oc*u;
                l3 = -ny - gm1oc*v;
                l4 = -nz - gm1oc*w;
                l5 = gm1oc;
                
                // multiplicative factor: 0.5/rho (left eigenvectors)
                //                        rho     (right eigenvectors)
                //                        eig5    (eigenvalue scaling)
                alpha = 0.5*(Roe_dQ[0]*l1 + Roe_dQ[1]*l2 + Roe_dQ[2]*l3 + Roe_dQ[3]*l4 + Roe_dQ[4]*l5)*Roe_Eigen[4][4];
                
                r1 = cinv;
                r2 = u*cinv - nx;
                r3 = v*cinv - ny;
                r4 = w*cinv - nz;
                r5 = 0.5*q2*cinv - ubar + 1.0/gm1oc;

                Flux_Roe_Diss[0] = alpha*r1;
                Flux_Roe_Diss[1] = alpha*r2;
                Flux_Roe_Diss[2] = alpha*r3;
                Flux_Roe_Diss[3] = alpha*r4;
                Flux_Roe_Diss[4] = alpha*r5;
            }
        } else { // Wave Traveling From Right to Left
            // Compute Convective Flux: flux_R
            Flux_Roe_Conv[0] = rho_R*ubar_R;
            Flux_Roe_Conv[1] =   u_R*Flux_Roe_Conv[0] + p_R*nx;
            Flux_Roe_Conv[2] =   v_R*Flux_Roe_Conv[0] + p_R*ny;
            Flux_Roe_Conv[3] =   w_R*Flux_Roe_Conv[0] + p_R*nz;
            Flux_Roe_Conv[4] =  ht_R*Flux_Roe_Conv[0];
            
            // Compute the Dissipation Flux
            // Subsonic
            if (ubar > -c) {
                // Eigenvalues 4 is positive, 1, 2, 3, 5 are negative
                // So, compute form the right with the one eigenvector that comes into play
                
                // Compute dQ
                // Conservative Variable Formulation
                if (Variable_Type == VARIABLE_CONSERVATIVE) {
                    Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
                    Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
                    Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
                    Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
                    Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
                }
                // Primitive Variable Formulation
                if ((Variable_Type == VARIABLE_PRIMITIVE_PUT) || (Variable_Type == VARIABLE_PRIMITIVE_RUP)) {
                    Roe_dQ[0] = rho_R      - rho_L;
                    Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
                    Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
                    Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
                    Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
                }
                
                // Eigenvalues
                // Apply Entropy Fix
                if (EntropyFix != 0) {
                    Roe_EntropyFix(ubar_L, c_L, ubar_R, c_R, ubar, c, Roe_Eigen);
                } else {
                    Roe_Eigen[0][0] = fabs(ubar);
                    Roe_Eigen[1][1] = Roe_Eigen[0][0];
                    Roe_Eigen[2][2] = Roe_Eigen[0][0];
                    Roe_Eigen[3][3] = fabs(ubar + c);
                    Roe_Eigen[4][4] = fabs(ubar - c);
                }
                Roe_Eigen[3][3] = ISIGN(ubar + c)*Roe_Eigen[3][3];
                
                l1 = -ubar + 0.5*gm1oc*q2;
                l2 = nx - gm1oc*u;
                l3 = ny - gm1oc*v;
                l4 = nz - gm1oc*w;
                l5 = gm1oc;
                
                // multiplicative factor: 0.5/rho (left eigenvectors)
                //                        rho     (right eigenvectors)
                //                        eig4    (eigenvalue scaling)
                alpha = 0.5*(Roe_dQ[0]*l1 + Roe_dQ[1]*l2 + Roe_dQ[2]*l3 + Roe_dQ[3]*l4 + Roe_dQ[4]*l5)*Roe_Eigen[3][3];
                
                r1 = cinv;
                r2 = u*cinv + nx;
                r3 = v*cinv + ny;
                r4 = w*cinv + nz;
                r5 = 0.5*q2*cinv + ubar + 1.0/gm1oc;

                Flux_Roe_Diss[0] = - alpha*r1;
                Flux_Roe_Diss[1] = - alpha*r2;
                Flux_Roe_Diss[2] = - alpha*r3;
                Flux_Roe_Diss[3] = - alpha*r4;
                Flux_Roe_Diss[4] = - alpha*r5;
            }
        }
        
        // Multiply the Roe Flux with Area
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] *= area;
            Flux_Roe_Diss[i] *= area;
        }
    } else
        error("Compute_RoeFlux_Optimized: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux
//------------------------------------------------------------------------------
void Compute_RoeFlux(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j, k, maxcount;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar, ubar1, ubar2;
    double Mach, nmax, fix, fixUn, sigma, maxlambda;
    Vector3D Vn1, Vn2;
    
    double area;
    double nx, ny, nz;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Internal Nodes
        if (node_R < nNode) {
            // Make Solution Second Order
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, Roe_Q_L, Roe_Q_R);
            } else {
                Roe_Q_L[0] = Q1[node_L];
                Roe_Q_L[1] = Q2[node_L];
                Roe_Q_L[2] = Q3[node_L];
                Roe_Q_L[3] = Q4[node_L];
                Roe_Q_L[4] = Q5[node_L];
                Roe_Q_R[0] = Q1[node_R];
                Roe_Q_R[1] = Q2[node_R];
                Roe_Q_R[2] = Q3[node_R];
                Roe_Q_R[3] = Q4[node_R];
                Roe_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            Roe_Q_L[0] = Q1[node_L];
            Roe_Q_L[1] = Q2[node_L];
            Roe_Q_L[2] = Q3[node_L];
            Roe_Q_L[3] = Q4[node_L];
            Roe_Q_L[4] = Q5[node_L];
            Roe_Q_R[0] = Q1[node_R];
            Roe_Q_R[1] = Q2[node_R];
            Roe_Q_R[2] = Q3[node_R];
            Roe_Q_R[3] = Q4[node_R];
            Roe_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // Compute flux_L
        Roe_flux_L[0] = rho_L*ubar_L;
        Roe_flux_L[1] =   u_L*Roe_flux_L[0] + p_L*nx;
        Roe_flux_L[2] =   v_L*Roe_flux_L[0] + p_L*ny;
        Roe_flux_L[3] =   w_L*Roe_flux_L[0] + p_L*nz;
        Roe_flux_L[4] =  ht_L*Roe_flux_L[0];
        
        // Compute flux_R
        Roe_flux_R[0] = rho_R*ubar_R;
        Roe_flux_R[1] =   u_R*Roe_flux_R[0] + p_R*nx;
        Roe_flux_R[2] =   v_R*Roe_flux_R[0] + p_R*ny;
        Roe_flux_R[3] =   w_R*Roe_flux_R[0] + p_R*nz;
        Roe_flux_R[4] =  ht_R*Roe_flux_R[0];
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        phi   = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c     = 0.0;
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                c = (Gamma - 1.0)*ht - phi;
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_BTW:
                c = ht/(Ref_Mach*Ref_Mach) - phi;
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_LMROE:
                c = (Gamma - 1.0)*ht - phi;
                c = sqrt(c);
                break;
            default:
                error("Compute_RoeFlux: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }
        ubar  = u*nx + v*ny + w*nz;

        // Do if Low Mach Number Fix is Requested
        if (PrecondMethod == SOLVER_PRECOND_ROE_LMFIX) {
            // Compute the point on the plane using normal
            // Get the max of nx, ny, nz
            nmax = fabs(nx);
            maxcount = 1;
            if (nmax < fabs(ny)) {
                nmax = fabs(ny);
                maxcount = 2;
            }
            if (nmax < fabs(nz)) {
                nmax = fabs(nz);
                maxcount = 3;
            }

            // Compute the vector in the plane
            if ((maxcount == 1) && (nmax > DBL_TOLERANCE)) {
                Vn1.vec[0] = -(ny + nz)/nx;
                Vn1.vec[1] = 1.0;
                Vn1.vec[2] = 1.0;
            } else if ((maxcount == 2) && (nmax > DBL_TOLERANCE)) {
                Vn1.vec[0] = 1.0;
                Vn1.vec[1] = -(nx + nz)/ny;
                Vn1.vec[2] = 1.0;
            } else if ((maxcount == 3) && (nmax > DBL_TOLERANCE)) {
                Vn1.vec[0] = 1.0;
                Vn1.vec[1] = 1.0;
                Vn1.vec[2] = -(nx + ny)/nz;
            } else {
                error("Compute_RoeFlux: Unable to compute point on the plane (nx, ny, nz): (%lf, %lf, %lf)", nx, ny, nz);
            }
            Vn1.normalize();
            ubar1 = u*Vn1.vec[0] +  v*Vn1.vec[1] +  w*Vn1.vec[2];

            // Compute the Second vector in the plain
            Vn2 = Vn1%areavec;
            Vn2.normalize();
            ubar2 = u*Vn2.vec[0] +  v*Vn2.vec[1] +  w*Vn2.vec[2];
            
            // Compute Local Mach Scaling Fix
            Mach = fabs(ubar) + fabs(ubar1) + fabs(ubar2);
            Mach = Mach/c;
            fix = MIN(Mach, 1.0);
        } else // ROE
            fix = 1.0;
              
        // M
        Roe_M[0][0] = 1.0;
        Roe_M[0][1] = 0.0;
        Roe_M[0][2] = 0.0;
        Roe_M[0][3] = 0.0;
        Roe_M[0][4] = 0.0;

        Roe_M[1][0] = u;
        Roe_M[1][1] = rho;
        Roe_M[1][2] = 0.0;
        Roe_M[1][3] = 0.0;
        Roe_M[1][4] = 0.0;

        Roe_M[2][0] = v;
        Roe_M[2][1] = 0.0;
        Roe_M[2][2] = rho;
        Roe_M[2][3] = 0.0;
        Roe_M[2][4] = 0.0;

        Roe_M[3][0] = w;
        Roe_M[3][1] = 0.0;
        Roe_M[3][2] = 0.0;
        Roe_M[3][3] = rho;
        Roe_M[3][4] = 0.0;

        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                Roe_M[4][0] = phi/(Gamma - 1.0);
                Roe_M[4][1] = rho * u;
                Roe_M[4][2] = rho * v;
                Roe_M[4][3] = rho * w;
                Roe_M[4][4] = 1.0/(Gamma - 1.0);
                break;
            case NONDIMENSIONAL_BTW:
                Roe_M[4][0] = phi * Ref_Mach * Ref_Mach;
                Roe_M[4][1] = (Gamma - 1.0) * rho * u * Ref_Mach * Ref_Mach;
                Roe_M[4][2] = (Gamma - 1.0) * rho * v * Ref_Mach * Ref_Mach;
                Roe_M[4][3] = (Gamma - 1.0) * rho * w * Ref_Mach * Ref_Mach;
                Roe_M[4][4] = Ref_Mach * Ref_Mach;
                break;
            case NONDIMENSIONAL_LMROE:
                Roe_M[4][0] = phi/(Gamma - 1.0);
                Roe_M[4][1] = rho * u;
                Roe_M[4][2] = rho * v;
                Roe_M[4][3] = rho * w;
                Roe_M[4][4] = 1.0/(Gamma - 1.0);
                break;
            default:
                error("Compute_RoeFlux: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }

        // Minv
        Roe_Minv[0][0] = 1.0;
        Roe_Minv[0][1] = 0.0;
        Roe_Minv[0][2] = 0.0;
        Roe_Minv[0][3] = 0.0;
        Roe_Minv[0][4] = 0.0;

        Roe_Minv[1][0] = -u/rho;
        Roe_Minv[1][1] = 1.0/rho;
        Roe_Minv[1][2] = 0.0;
        Roe_Minv[1][3] = 0.0;
        Roe_Minv[1][4] = 0.0;

        Roe_Minv[2][0] = -v/rho;
        Roe_Minv[2][1] = 0.0;
        Roe_Minv[2][2] = 1.0/rho;
        Roe_Minv[2][3] = 0.0;
        Roe_Minv[2][4] = 0.0;

        Roe_Minv[3][0] = -w/rho;
        Roe_Minv[3][1] = 0.0;
        Roe_Minv[3][2] = 0.0;
        Roe_Minv[3][3] = 1.0/rho;
        Roe_Minv[3][4] = 0.0;

        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                Roe_Minv[4][0] = phi;
                Roe_Minv[4][1] = -u * (Gamma - 1.0);
                Roe_Minv[4][2] = -v * (Gamma - 1.0);
                Roe_Minv[4][3] = -w * (Gamma - 1.0);
                Roe_Minv[4][4] = (Gamma - 1.0);
                break;
            case NONDIMENSIONAL_BTW:
                Roe_Minv[4][0] = phi;
                Roe_Minv[4][1] = -u * (Gamma - 1.0);
                Roe_Minv[4][2] = -v * (Gamma - 1.0);
                Roe_Minv[4][3] = -w * (Gamma - 1.0);
                Roe_Minv[4][4] = 1.0/(Ref_Mach * Ref_Mach);
                break;
            case NONDIMENSIONAL_LMROE:
                Roe_Minv[4][0] = phi;
                Roe_Minv[4][1] = -u * (Gamma - 1.0);
                Roe_Minv[4][2] = -v * (Gamma - 1.0);
                Roe_Minv[4][3] = -w * (Gamma - 1.0);
                Roe_Minv[4][4] = (Gamma - 1.0);
                break;
            default:
                error("Compute_RoeFlux: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }
        
        // P
        Roe_P[0][0] = nx;
        Roe_P[0][1] = ny;
        Roe_P[0][2] = nz;
        Roe_P[0][3] = rho/c;
        Roe_P[0][4] = rho/c;

        Roe_P[1][0] = 0.0;
        Roe_P[1][1] = -nz;
        Roe_P[1][2] = ny;
        Roe_P[1][3] = nx;
        Roe_P[1][4] = -nx;

        Roe_P[2][0] = nz;
        Roe_P[2][1] = 0.0;
        Roe_P[2][2] = -nx;
        Roe_P[2][3] = ny;
        Roe_P[2][4] = -ny;

        Roe_P[3][0] = -ny;
        Roe_P[3][1] = nx;
        Roe_P[3][2] = 0.0;
        Roe_P[3][3] = nz;
        Roe_P[3][4] = -nz;

        Roe_P[4][0] = 0.0;
        Roe_P[4][1] = 0.0;
        Roe_P[4][2] = 0.0;
        Roe_P[4][3] = rho * c;
        Roe_P[4][4] = rho * c;

        // Pinv
        Roe_Pinv[0][0] = nx;
        Roe_Pinv[0][1] = 0.0;
        Roe_Pinv[0][2] = nz;
        Roe_Pinv[0][3] = -ny;
        Roe_Pinv[0][4] = -nx/(c * c);

        Roe_Pinv[1][0] = ny;
        Roe_Pinv[1][1] = -nz;
        Roe_Pinv[1][2] = 0.0;
        Roe_Pinv[1][3] = nx;
        Roe_Pinv[1][4] = -ny/(c * c);

        Roe_Pinv[2][0] = nz;
        Roe_Pinv[2][1] = ny;
        Roe_Pinv[2][2] = -nx;
        Roe_Pinv[2][3] = 0.0;
        Roe_Pinv[2][4] = -nz/(c * c);

        Roe_Pinv[3][0] = 0.0;
        Roe_Pinv[3][1] = 0.5 * nx;
        Roe_Pinv[3][2] = 0.5 * ny;
        Roe_Pinv[3][3] = 0.5 * nz;
        Roe_Pinv[3][4] = 1.0/(2.0 * rho * c);

        Roe_Pinv[4][0] = 0.0;
        Roe_Pinv[4][1] = -0.5 * nx;
        Roe_Pinv[4][2] = -0.5 * ny;
        Roe_Pinv[4][3] = -0.5 * nz;
        Roe_Pinv[4][4] = 1.0/(2.0 * rho * c);

        // START: Computing 
        // if Roe:      A*dQ = M*P*|Lambda|*Pinv*Minv*dQ
        // if LMRoeFix: A*dw = M*P*|Lambda|*Pinv*dw
        // Note: w = non-conservative {rho, v, u, w, p}
        // Where: dQ = Minv*dw
        
        // Calculate T = M*P
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_P, Roe_T);
        
        // Calculate EigenMatrix |Lambda|
        for (j = 0; j < 5; j++)
            for (k = 0; k < 5; k++)
                Roe_Eigen[j][k] = 0.0;

        // Apply Entropy Fix
        if (EntropyFix != 0) {
            Roe_EntropyFix(ubar_L, c_L, ubar_R, c_R, ubar, c, Roe_Eigen);
        } else {
            Roe_Eigen[0][0] = fabs(ubar);
            Roe_Eigen[1][1] = Roe_Eigen[0][0];
            Roe_Eigen[2][2] = Roe_Eigen[0][0];
            Roe_Eigen[3][3] = fabs(ubar + c);
            Roe_Eigen[4][4] = fabs(ubar - c);
        }
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Min and Max EigenValue Value
            MinEigenLamda1 = MIN(MinEigenLamda1, Roe_Eigen[0][0]);
            MaxEigenLamda1 = MAX(MaxEigenLamda1, Roe_Eigen[0][0]);
            MinEigenLamda4 = MIN(MinEigenLamda4, Roe_Eigen[3][3]);
            MaxEigenLamda4 = MAX(MaxEigenLamda4, Roe_Eigen[3][3]);
            MinEigenLamda5 = MIN(MinEigenLamda5, Roe_Eigen[4][4]);
            MaxEigenLamda5 = MAX(MaxEigenLamda5, Roe_Eigen[4][4]);

            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothType == RESIDUAL_SMOOTH_IMPLICIT) {
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L + 1]; i++) {
                // Get the Node ID and Match with node_R
                if (crs_JA_Node2Node[i] == node_R) {
                    RM_MaxEigenValue[i] = maxlambda*area;
                    break;
                }
            }
            RM_SumMaxEigenValue[node_L] += maxlambda*area;
            // Only Physical Nodes
            if (node_R < nNode) {
                for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R + 1]; i++) {
                    // Get the Node ID and Match with node_L
                    if (crs_JA_Node2Node[i] == node_L) {
                        RM_MaxEigenValue[i] = maxlambda*area;
                        break;
                    }
                }
                RM_SumMaxEigenValue[node_R] += maxlambda*area;
            }
        }
        
        // Apply Low Mach Fix if Requested
        if (PrecondMethod == SOLVER_PRECOND_ROE_LMFIX) {
            // Temporary = Tinv = EigenMatrix*Pinv
            MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Pinv, Roe_Tinv);

            // Get Matrix A = M*P*|Lambda|*Pinv
            MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_Tinv, Roe_A);

            // Compute dQ1
            Roe_dQ[0] = rho_R - rho_L;
            Roe_dQ[1] = 0.0;
            Roe_dQ[2] = 0.0;
            Roe_dQ[3] = 0.0;
            Roe_dQ[4] = p_R - p_L;
            // Compute |A|*dQ1
            MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
            
            fixUn = fix*((u_R - u_L)*nx + (v_R - v_L)*ny + (w_R - w_L)*nz);
            
            // Compute dQ3
            Roe_dw[0] = (0.5*rho*fixUn*(fabs(ubar + c) - fabs(ubar - c)))/c;
            Roe_dw[1] = fabs(ubar)*((u_R - u_L) - fixUn*nx) + 0.5*nx*fixUn*(fabs(ubar + c) + fabs(ubar - c));
            Roe_dw[2] = fabs(ubar)*((v_R - v_L) - fixUn*ny) + 0.5*ny*fixUn*(fabs(ubar + c) + fabs(ubar - c));
            Roe_dw[3] = fabs(ubar)*((w_R - w_L) - fixUn*nz) + 0.5*nz*fixUn*(fabs(ubar + c) + fabs(ubar - c));
            Roe_dw[4] = 0.5*rho*c*fixUn*(fabs(ubar + c) - fabs(ubar - c));

            // Compute M*dQ3
            MC_Matrix_Mul_Vector(5, 5, Roe_M, Roe_dw, Roe_dQ);
            
            // Finally Compute A*dw = M*dQ3 + M*P*|Lambda|*Pinv*dQ1 =  M*P*|Lambda|*Pinv*dw
            for (i = 0; i < NEQUATIONS; i++)
                Roe_fluxA[i] += Roe_dQ[i];
            // END: Computing A*dw = M*P*|Lambda|*Pinv*dw
        } else { // ROE
            // Calculate Tinv = Pinv*Minv
            MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
            
            // EigenMatrix*Tinv = |Lambda|*Pinv*Minv
            MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);

            // Tinv = EigenMatrix*Tinv
            for (j = 0; j < 5; j++) {
                for (k = 0; k < 5; k++) {
                    Roe_Tinv[j][k] = Roe_A[j][k];
                    Roe_A[j][k] = 0.0;
                }
            }

            // Get Matrix A = M*P*|Lambda|*Pinv*Minv
            MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_Tinv, Roe_A);

            // Set up matrix vector multiplication
            // Compute dQ
            // Conservative Variable Formulation
            if (Variable_Type == VARIABLE_CONSERVATIVE) {
                Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
                Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
                Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
                Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
                Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
            }
            // Primitive Variable Formulation
            if ((Variable_Type == VARIABLE_PRIMITIVE_PUT) || (Variable_Type == VARIABLE_PRIMITIVE_RUP)) {
                Roe_dQ[0] = rho_R      - rho_L;
                Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
                Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
                Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
                Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
            }

            // Compute |A|*dQ
            MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
            // END: Computing A*dQ = M*P*|Lambda|*Pinv*Minv*dQ
        }
        
        // Compute the Roe Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Weiss and Smith Pre-Conditioner 
//------------------------------------------------------------------------------
void Compute_RoeFlux_WeissSmith_Old(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, area, nx, ny, nz;
    
    double eps, alpha, beta, pubar, pc, Ur, Mstar, Cstar, deltaUbar;
    double deltaU, deltaP, mach;
    double dtmp;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Internal Nodes
        if (node_R < nNode) {
            // Make Solution Second Order
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, Roe_Q_L, Roe_Q_R);
            } else {
                Roe_Q_L[0] = Q1[node_L];
                Roe_Q_L[1] = Q2[node_L];
                Roe_Q_L[2] = Q3[node_L];
                Roe_Q_L[3] = Q4[node_L];
                Roe_Q_L[4] = Q5[node_L];
                Roe_Q_R[0] = Q1[node_R];
                Roe_Q_R[1] = Q2[node_R];
                Roe_Q_R[2] = Q3[node_R];
                Roe_Q_R[3] = Q4[node_R];
                Roe_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            Roe_Q_L[0] = Q1[node_L];
            Roe_Q_L[1] = Q2[node_L];
            Roe_Q_L[2] = Q3[node_L];
            Roe_Q_L[3] = Q4[node_L];
            Roe_Q_L[4] = Q5[node_L];
            Roe_Q_R[0] = Q1[node_R];
            Roe_Q_R[1] = Q2[node_R];
            Roe_Q_R[2] = Q3[node_R];
            Roe_Q_R[3] = Q4[node_R];
            Roe_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // Compute flux_L
        Roe_flux_L[0] = rho_L*ubar_L;
        Roe_flux_L[1] =   u_L*Roe_flux_L[0] + p_L*nx;
        Roe_flux_L[2] =   v_L*Roe_flux_L[0] + p_L*ny;
        Roe_flux_L[3] =   w_L*Roe_flux_L[0] + p_L*nz;
        Roe_flux_L[4] =  ht_L*Roe_flux_L[0];
        
        // Compute flux_R
        Roe_flux_R[0] = rho_R*ubar_R;
        Roe_flux_R[1] =   u_R*Roe_flux_R[0] + p_R*nx;
        Roe_flux_R[2] =   v_R*Roe_flux_R[0] + p_R*ny;
        Roe_flux_R[3] =   w_R*Roe_flux_R[0] + p_R*nz;
        Roe_flux_R[4] =  ht_R*Roe_flux_R[0];
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        // Numerical Average Variables
//        rho = 0.5*(rho_R + rho_L);
//        u   = 0.5*(u_R   + u_L);
//        v   = 0.5*(v_R   + v_L);
//        w   = 0.5*(w_R   + w_L);
//        ht  = 0.5*(ht_R  + ht_L);
        // Other Quantities
        phi  = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c    = (Gamma - 1.0)*ht - phi;
        c    = sqrt(c);
        ubar = u*nx + v*ny + w*nz;
        
        // Compute Precondition Variables
        mach = sqrt(u*u + v*v + w*w)/c;
        eps  = MIN(1.0, MAX(1e-10, MAX(mach*mach, 1.4*fabs(p_R - p_L)/(rho*c*c))));
        
        // Compute beta for Ideal Gas
        beta = 1.0/(c*c);
        
        // Compute Ur for Ideal Gas
        dtmp = sqrt(u*u + v*v + w*w);
        if (dtmp < eps*c)
            Ur = eps*c;
        else if ((eps*c < dtmp) && (dtmp < c))
            Ur = dtmp;
        else
            Ur = c;
        
        // Compute alpha
        alpha  = 0.5*(1.0 - beta*Ur*Ur);
        
        pubar  = ubar*(1.0 - alpha);
        pc     = sqrt(alpha*alpha*ubar*ubar + Ur*Ur);
        Mstar  = 0.5*(fabs(pubar + pc) + fabs(pubar - pc));
        Cstar  = 0.5*(fabs(pubar + pc) - fabs(pubar - pc))/c;
        deltaUbar = (u_R - u_L)*nx + (v_R - v_L)*ny + (w_R - w_L)*nz;
        deltaU    = Mstar*deltaUbar + (Cstar - (1.0 - 2.0*alpha)*fabs(ubar) - alpha*ubar*Mstar)*(p_R - p_L)/(rho*Ur*Ur);
        deltaP    = Mstar*(p_R - p_L) + (Cstar - fabs(ubar) + alpha*ubar*Mstar)*rho*deltaUbar;
        
        // Compute |Ubar|*(Q_R - Q_L)*(n.n)
        dtmp = fabs(ubar)*(nx*nx + ny*ny + nz*nz);
        // Conservative Variable Formulation
        if (Variable_Type == VARIABLE_CONSERVATIVE) {
            for (i = 0; i < NEQUATIONS; i++)
                Roe_fluxA[i] = dtmp*(Roe_Q_R[i] - Roe_Q_L[i]);
        }
        // Primitive Variable Formulation
        if ((Variable_Type == VARIABLE_PRIMITIVE_PUT) || (Variable_Type == VARIABLE_PRIMITIVE_RUP)) {
            Roe_fluxA[0] = dtmp*(rho_R      - rho_L);
            Roe_fluxA[1] = dtmp*(rho_R*u_R  - rho_L*u_L);
            Roe_fluxA[2] = dtmp*(rho_R*v_R  - rho_L*v_L);
            Roe_fluxA[3] = dtmp*(rho_R*w_R  - rho_L*w_L);
            Roe_fluxA[4] = dtmp*(rho_R*et_R - rho_L*et_L);
        }
        
        // Compute deltaU*Q_Roe
        dtmp = deltaU*rho*(nx*nx + ny*ny + nz*nz);
        Roe_fluxA[0] += dtmp;
        Roe_fluxA[1] += dtmp*u;
        Roe_fluxA[2] += dtmp*v;
        Roe_fluxA[3] += dtmp*w;
        Roe_fluxA[4] += dtmp*ht;
        
        // Compute deltaP*n
        Roe_fluxA[0] += 0.0;
        Roe_fluxA[1] += deltaP*nx;
        Roe_fluxA[2] += deltaP*ny;
        Roe_fluxA[3] += deltaP*nz;
        Roe_fluxA[4] += deltaP*ubar;
        
        // Compute the Roe Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Weiss and Smith Pre-Conditioner 
//------------------------------------------------------------------------------
void Compute_RoeFlux_WeissSmith(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, area, nx, ny, nz;
    double lambda1, lambda2, lambda3, lambda4, lambda5, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Internal Nodes
        if (node_R < nNode) {
            // Make Solution Second Order
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, Roe_Q_L, Roe_Q_R);
            } else {
                Roe_Q_L[0] = Q1[node_L];
                Roe_Q_L[1] = Q2[node_L];
                Roe_Q_L[2] = Q3[node_L];
                Roe_Q_L[3] = Q4[node_L];
                Roe_Q_L[4] = Q5[node_L];
                Roe_Q_R[0] = Q1[node_R];
                Roe_Q_R[1] = Q2[node_R];
                Roe_Q_R[2] = Q3[node_R];
                Roe_Q_R[3] = Q4[node_R];
                Roe_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            Roe_Q_L[0] = Q1[node_L];
            Roe_Q_L[1] = Q2[node_L];
            Roe_Q_L[2] = Q3[node_L];
            Roe_Q_L[3] = Q4[node_L];
            Roe_Q_L[4] = Q5[node_L];
            Roe_Q_R[0] = Q1[node_R];
            Roe_Q_R[1] = Q2[node_R];
            Roe_Q_R[2] = Q3[node_R];
            Roe_Q_R[3] = Q4[node_R];
            Roe_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // Compute flux_L
        Roe_flux_L[0] = rho_L*ubar_L;
        Roe_flux_L[1] =   u_L*Roe_flux_L[0] + p_L*nx;
        Roe_flux_L[2] =   v_L*Roe_flux_L[0] + p_L*ny;
        Roe_flux_L[3] =   w_L*Roe_flux_L[0] + p_L*nz;
        Roe_flux_L[4] =  ht_L*Roe_flux_L[0];
        
        // Compute flux_R
        Roe_flux_R[0] = rho_R*ubar_R;
        Roe_flux_R[1] =   u_R*Roe_flux_R[0] + p_R*nx;
        Roe_flux_R[2] =   v_R*Roe_flux_R[0] + p_R*ny;
        Roe_flux_R[3] =   w_R*Roe_flux_R[0] + p_R*nz;
        Roe_flux_R[4] =  ht_R*Roe_flux_R[0];
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        phi  = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c    = (Gamma - 1.0)*ht - phi;
        c    = sqrt(c);
        ubar = u*nx + v*ny + w*nz;
        
        // Compute dQ
        // Conservative Variable Formulation
        if (Variable_Type == VARIABLE_CONSERVATIVE) {
            Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
            Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
            Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
            Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
            Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        }
        // Primitive Variable Formulation
        if ((Variable_Type == VARIABLE_PRIMITIVE_PUT) || (Variable_Type == VARIABLE_PRIMITIVE_RUP)) {
            Roe_dQ[0] = rho_R      - rho_L;
            Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
            Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
            Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
            Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
        }

        double beta, alpha, alpha1, alpha2, alpha3, alpha4, alpha5, A, B, r, s, t;
        double dtmp;
        
        // Compute the Eigenvalues of the Precondition Dissipation Term
        beta    = MIN(1.0, sqrt(u*u + v*v + w*w)/c);
        beta    = 1.0;
        lambda1 = ubar;
        lambda2 = lambda1;
        lambda3 = lambda1;
        dtmp = (1.0 - beta*beta)*(1.0 - beta*beta)*ubar*ubar + 4.0*beta*beta*c*c;
        lambda4 = 0.5*((1.0 + beta*beta)*ubar + sqrt(dtmp));
        lambda5 = 0.5*((1.0 + beta*beta)*ubar - sqrt(dtmp));
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(lambda1, MAX(lambda4, lambda5));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Compute Precondition variables
        A      = p_R - p_L;
        B      = -rho*((u_R - u_L)*nx + (v_R - v_L)*ny + (w_R - w_L)*nz);
        r      = lambda3 - lambda1*beta*beta;
        s      = lambda4 - lambda1*beta*beta;
        t      = 0.5*(lambda5 - lambda4);
        alpha1 = Roe_dQ[0] - A/c*c;
        alpha2 = (ubar*nx - u)*Roe_dQ[0] + (ny*ny + nz*nz)*Roe_dQ[1] - nx*ny*Roe_dQ[2] - nx*nz*Roe_dQ[3];
        alpha3 = (ubar*nz - w)*Roe_dQ[0] - nx*nz*Roe_dQ[1] - nz*ny*Roe_dQ[2] + (ny*ny + nz*nz)*Roe_dQ[3];
        alpha4 =  (s*A + beta*beta*c*c*B)/(2.0*beta*beta*c*c*t);
        alpha5 = -(r*A + beta*beta*c*c*B)/(2.0*beta*beta*c*c*t);
        alpha  = -(alpha2*nx + alpha3*nz)/ny;
        
        // Compute Dissipation Term
        Roe_fluxA[0] = fabs(lambda1)*alpha1 + fabs(lambda4)*alpha4 + fabs(lambda5)*alpha5;
        Roe_fluxA[1] = fabs(lambda1)*(alpha1*u + alpha2) + fabs(lambda4)*alpha4*(u + r*nx) + fabs(lambda5)*alpha5*(u + s*nx);
        Roe_fluxA[2] = fabs(lambda1)*(alpha1*v + alpha)  + fabs(lambda4)*alpha4*(v + r*ny) + fabs(lambda5)*alpha5*(v + s*ny);
        Roe_fluxA[3] = fabs(lambda1)*(alpha1*w + alpha3) + fabs(lambda4)*alpha4*(w + r*nz) + fabs(lambda5)*alpha5*(w + s*nz);
        Roe_fluxA[4] = fabs(lambda1)*(0.5*alpha1*(u*u + v*v + w*w) + alpha2*u + alpha*v + alpha3*w)
                        + fabs(lambda4)*alpha4*(rho*ht + r*lambda1) + fabs(lambda5)*alpha5*(rho*ht + s*lambda1);
        
        // Compute the Roe Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Preconditioned Roe Residual with Weiss Smith Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Precondition_Residual_Roe_WeissSmith(void) {
    double rho, u, v, w, ht, et, p, T, c, q2, mach;
    double Ur, eps, Cp, Theta, drhodt;
    double Q[5];
    double flux_roe[5];
    double flux_roe_diss[5];
    double dtmp;
    
    // Multiply by Precondition Matrix and Construct the Flux
    for (int i = 0; i < nNode; i++) {
        // Get the Variables
        Q[0] = Q1[i];
        Q[1] = Q2[i];
        Q[2] = Q3[i];
        Q[3] = Q4[i];
        Q[4] = Q5[i];
        
        // Compute Equation of State
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
        // Compute the EPS
        eps  = MIN(1.0, MAX(1e-10, mach*mach));

        // Compute Ur for Ideal Gas
        dtmp = sqrt(q2);
        if (dtmp < eps*c)
            Ur = eps*c;
        else if ((eps*c < dtmp) && (dtmp < c))
            Ur = dtmp;
        else
            Ur = c;

        Theta  = (1.0/(Ur*Ur) + (Gamma - 1.0)/(c*c));
        Cp     = 1/(Gamma - 1.0);
        drhodt = -(rho*rho)/(Gamma*(p + Gauge_Pressure));

        // Get the Convective Flux
        flux_roe_diss[0] = Res1[i];
        flux_roe_diss[1] = Res2[i];
        flux_roe_diss[2] = Res3[i];
        flux_roe_diss[3] = Res4[i];
        flux_roe_diss[4] = Res5[i];
        
        // Finally Compute precondition residual
        flux_roe[0] = 0.0;
        flux_roe[1] = 0.0;
        flux_roe[2] = 0.0;
        flux_roe[3] = 0.0;
        flux_roe[4] = 0.0;
        // Manual Matrix Vector Multiplication
        flux_roe[0] = ht*flux_roe_diss[0] - flux_roe_diss[4] + u*flux_roe_diss[1] - u*u*flux_roe_diss[0] 
                        + v*flux_roe_diss[2] - v*v*flux_roe_diss[0] + w*flux_roe_diss[3] - w*w*flux_roe_diss[0];
        flux_roe[0] = (drhodt*flux_roe[0] + Cp*rho*flux_roe_diss[0])/(drhodt + Cp*Theta*rho);
        flux_roe[1] = (flux_roe_diss[1] - u*flux_roe_diss[0])/rho;
        flux_roe[2] = (flux_roe_diss[2] - v*flux_roe_diss[0])/rho;
        flux_roe[3] = (flux_roe_diss[3] - w*flux_roe_diss[0])/rho;
        flux_roe[4] = Theta*(flux_roe_diss[4] - u*flux_roe_diss[1] - v*flux_roe_diss[2] - w*flux_roe_diss[3]) 
                        + flux_roe_diss[0]*(1.0 + Theta*(u*u + v*v + w*w - ht));
        flux_roe[4] = flux_roe[4]/(drhodt + Cp*Theta*rho);

        // Compute the Preconditioned Roe Flux Residual
        Res1[i] = flux_roe[0] + Res1_Diss[i];
        Res2[i] = flux_roe[1] + Res2_Diss[i];
        Res3[i] = flux_roe[2] + Res3_Diss[i];
        Res4[i] = flux_roe[3] + Res4_Diss[i];
        Res5[i] = flux_roe[4] + Res5_Diss[i];
        
        // Transform the Residual in Conservative form
        // Primitive Variable Formulation
        if (Variable_Type == VARIABLE_CONSERVATIVE) {
            q2 = 0.5*(Gamma - 1.0)*q2;
            
            // Compute the Transformation Matrix
            Roe_P[0][0] = Gamma/T;
            Roe_P[0][1] = 0.0;
            Roe_P[0][2] = 0.0;
            Roe_P[0][3] = 0.0;
            Roe_P[0][4] = -rho/T;

            Roe_P[1][0] = Gamma*u/T;
            Roe_P[1][1] = rho;
            Roe_P[1][2] = 0.0;
            Roe_P[1][3] = 0.0;
            Roe_P[1][4] = -u*rho/T;

            Roe_P[2][0] = Gamma*v/T;
            Roe_P[2][1] = 0.0;
            Roe_P[2][2] = rho;
            Roe_P[2][3] = 0.0;
            Roe_P[2][4] = -v*rho/T;

            Roe_P[3][0] = Gamma*w/T;
            Roe_P[3][1] = 0.0;
            Roe_P[3][2] = 0.0;
            Roe_P[3][3] = rho;
            Roe_P[3][4] = -w*rho/T;

            Roe_P[4][0] = (1.0 + Gamma*q2/T)/(Gamma - 1.0);
            Roe_P[4][1] = rho*u;
            Roe_P[4][2] = rho*v;
            Roe_P[4][3] = rho*v;
            Roe_P[4][4] = -q2*rho/((Gamma - 1.0)*T);
            
            // Get the Conservative Residual
            flux_roe[0] = Res1[i];
            flux_roe[1] = Res2[i];
            flux_roe[2] = Res3[i];
            flux_roe[3] = Res4[i];
            flux_roe[4] = Res5[i];
            MC_Matrix_Mul_Vector(5, 5, Roe_P, flux_roe, flux_roe_diss);
            
            // Update the Residual in Primitive Form
            Res1[i] = flux_roe_diss[0];
            Res2[i] = flux_roe_diss[1];
            Res3[i] = flux_roe_diss[2];
            Res4[i] = flux_roe_diss[3];
            Res5[i] = flux_roe_diss[4];
        }
    }
}

//------------------------------------------------------------------------------
//! Compute Preconditioned Roe Residual with Cecile Voizat Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Precondition_Residual_Roe_CecileVoizat(void) {
    int i, j, nid;
    double Q[5], lQ[5];
    double flux_roe[5];
    double flux_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, phi, sigma, reta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Multiply by Precondition Matrix and Construct the Flux
    for (i = 0; i < nNode; i++) {
        // Get the Variables
        Q[0] = Q1[i];
        Q[1] = Q2[i];
        Q[2] = Q3[i];
        Q[3] = Q4[i];
        Q[4] = Q5[i];
        
        // Compute Equation of State
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        phi = 0.5*(Gamma - 1.0)*q2;
        
        //======================================================================
        // Compute Precondition of Convective Flux and Assemble Total Flux
        //======================================================================
        // STEP - 1:
        // Compute the Local Max Mach and Max Change in Pressure around a node
        dp_max   = 0.0;
        mach_max = 0.0;
        // Check if Local Precondition is Required
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            // Check if Precondition Smoother is Requested
            if (PrecondSmooth != 0) {
                for (j = crs_IA_Node2Node[i]; j < crs_IA_Node2Node[i+1]; j++) {
                    nid   = crs_JA_Node2Node[j];
                    // Get the local Q's
                    lQ[0] = Q1[nid];
                    lQ[1] = Q2[nid];
                    lQ[2] = Q3[nid];
                    lQ[3] = Q4[nid];
                    lQ[4] = Q5[nid];
                    // Compute Local Equation of State
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    mach_max = MAX(mach_max, lmach);
                    dp_max   = MAX(dp_max, fabs(p - lp));
                }
            }
        }
        
        // STEP 2:
        // Compute the required parameters
        sigma = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            sigma = MAX(mach, mach_max);
            if (sigma < Ref_Mach)
                sigma = sqrt(Ref_Mach*sigma);
            sigma = MIN(1.0, sigma);
//            //if (Inf_Mach < 0.5) {
//                sigma = MIN(sqrt(1.0e-8/sqrt(q2)), Inf_Mach);
//                // Check if Precondition Smoother is Requested
//                if (PrecondSmooth != 0)
//                    sigma = MIN(1.0, MAX(MAX(MAX(sigma, mach), mach_max), 2.0*dp_max/(rho*c*c)));
//                else
//                    sigma = MIN(1.0, MAX(sigma, mach));
//            //sigma = MAX(mach, dp_max/(rho*c*c));
//            //sigma = MAX(1.0e-6, MIN(mach_max, dp_max/(Inf_Mach*Inf_Mach*rho*c*c)));
//                if (sigma > Inf_Mach)
//                    sigma  = Inf_Mach;
//            //}
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            sigma = MIN(1.0, Ref_Mach);
        
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sigma);
        MaxPrecondSigma = MAX(MaxPrecondSigma, sigma);
        
        // Compute the Precondition in Conservative form
        reta = (Gamma - 1.0)*(sigma*sigma - 1.0)/(c*c);

        // STEP 3:
        Roe_P[0][0] = 1.0 + 0.5*q2*reta;
        Roe_P[0][1] = -u*reta;
        Roe_P[0][2] = -v*reta;
        Roe_P[0][3] = -w*reta;
        Roe_P[0][4] = reta;

        Roe_P[1][0] = 0.5*u*q2*reta;
        Roe_P[1][1] = 1.0 - u*u*reta;
        Roe_P[1][2] = -u*v*reta;
        Roe_P[1][3] = -u*w*reta;
        Roe_P[1][4] = u*reta;

        Roe_P[2][0] = 0.5*v*q2*reta;
        Roe_P[2][1] = -u*v*reta;
        Roe_P[2][2] = 1.0 - v*v*reta;
        Roe_P[2][3] = -v*w*reta;
        Roe_P[2][4] = v*reta;

        Roe_P[3][0] = 0.5*w*q2*reta;
        Roe_P[3][1] = -u*w*reta;
        Roe_P[3][2] = -v*w*reta;
        Roe_P[3][3] = 1.0 - w*w*reta;
        Roe_P[3][4] = w*reta;

        Roe_P[4][0] = 0.5*q2*ht*reta;
        Roe_P[4][1] = -u*ht*reta;
        Roe_P[4][2] = -v*ht*reta;
        Roe_P[4][3] = -w*ht*reta;
        Roe_P[4][4] = 1.0 + ht*reta;
        
        // STEP 4:
        // Compute Preconditioned Convective Flux
        // Get the Convective Flux
        flux_roe_diss[0] = Res1[i];
        flux_roe_diss[1] = Res2[i];
        flux_roe_diss[2] = Res3[i];
        flux_roe_diss[3] = Res4[i];
        flux_roe_diss[4] = Res5[i];
        if (PrecondMethod == SOLVER_PRECOND_ROE_CV_ORIGINAL) {
            flux_roe_diss[0] += Res1_Diss[i];
            flux_roe_diss[1] += Res2_Diss[i];
            flux_roe_diss[2] += Res3_Diss[i];
            flux_roe_diss[3] += Res4_Diss[i];
            flux_roe_diss[4] += Res5_Diss[i];
        }

        // Finally Compute precondition residual
        flux_roe[0] = 0.0;
        flux_roe[1] = 0.0;
        flux_roe[2] = 0.0;
        flux_roe[3] = 0.0;
        flux_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_P, flux_roe_diss, flux_roe);

        // STEP 5:
        // Compute the Preconditioned Roe Flux Residual
        Res1[i] = flux_roe[0];
        Res2[i] = flux_roe[1];
        Res3[i] = flux_roe[2];
        Res4[i] = flux_roe[3];
        Res5[i] = flux_roe[4];
        // Marco Modification
        if (PrecondMethod == SOLVER_PRECOND_ROE_CV) {
            Res1[i] += Res1_Diss[i];
            Res2[i] += Res2_Diss[i];
            Res3[i] += Res3_Diss[i];
            Res4[i] += Res4_Diss[i];
            Res5[i] += Res5_Diss[i];
        }
        
        // Transform the Residual in Primitive form
        // Primitive Variable Formulation
        if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
            q2 = 0.5*(Gamma - 1.0)*q2;
            
            // Compute the Transformation Matrix
            Roe_P[0][0] = q2;
            Roe_P[0][1] = (1.0 - Gamma)*u;
            Roe_P[0][2] = (1.0 - Gamma)*v;
            Roe_P[0][3] = (1.0 - Gamma)*w;
            Roe_P[0][4] = (Gamma - 1.0);

            Roe_P[1][0] = -u/rho;
            Roe_P[1][1] = 1.0/rho;
            Roe_P[1][2] = 0.0;
            Roe_P[1][3] = 0.0;
            Roe_P[1][4] = 0.0;

            Roe_P[2][0] = -v/rho;
            Roe_P[2][1] = 0.0;
            Roe_P[2][2] = 1.0/rho;
            Roe_P[2][3] = 0.0;
            Roe_P[2][4] = 0.0;

            Roe_P[3][0] = -w/rho;
            Roe_P[3][1] = 0.0;
            Roe_P[3][2] = 0.0;
            Roe_P[3][3] = 1.0/rho;
            Roe_P[3][4] = 0.0;

            Roe_P[4][0] = (q2*Gamma - T)/rho;
            Roe_P[4][1] = (1.0 - Gamma)*Gamma*u/rho;
            Roe_P[4][2] = (1.0 - Gamma)*Gamma*v/rho;
            Roe_P[4][3] = (1.0 - Gamma)*Gamma*w/rho;
            Roe_P[4][4] = (1.0 - Gamma)*Gamma/rho;
            
            // Get the Conservative Residual
            flux_roe[0] = Res1[i];
            flux_roe[1] = Res2[i];
            flux_roe[2] = Res3[i];
            flux_roe[3] = Res4[i];
            flux_roe[4] = Res5[i];
            MC_Matrix_Mul_Vector(5, 5, Roe_P, flux_roe, flux_roe_diss);
            
            // Update the Residual in Primitive Form
            Res1[i] = flux_roe_diss[0];
            Res2[i] = flux_roe_diss[1];
            Res3[i] = flux_roe_diss[2];
            Res4[i] = flux_roe_diss[3];
            Res5[i] = flux_roe_diss[4];
        }
    }
}

//------------------------------------------------------------------------------
//! Compute Preconditioned Roe Residual with Briley Taylor Whitfield Pre-Conditioner
//! Note Variables are in Primitive Form: Density Velocity Pressure
//------------------------------------------------------------------------------
void Compute_Precondition_Residual_Roe_BTW(void) {
    int i, j, nid;
    double Q[5], lQ[5];
    double flux_roe[5];
    double flux_roe_conv[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, phi, beta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Multiply by Precondition Matrix and Construct the Flux
    for (i = 0; i < nNode; i++) {
        // Get the Variables
        Q[0] = Q1[i];
        Q[1] = Q2[i];
        Q[2] = Q3[i];
        Q[3] = Q4[i];
        Q[4] = Q5[i];
        
        // Compute Equation of State
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        phi = 0.5*(Gamma - 1.0)*q2;
        
        //======================================================================
        // Compute Precondition of Convective Flux and Assemble Total Flux
        //======================================================================
        // STEP - 1:
        // Compute the Local Max Mach and Max Change in Pressure around a node
        dp_max   = 0.0;
        mach_max = 0.0;
        // Check if Local Precondition is Required
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            // Check if Precondition Smoother is Requested
            if (PrecondSmooth != 0) {
                for (j = crs_IA_Node2Node[i]; j < crs_IA_Node2Node[i+1]; j++) {
                    nid   = crs_JA_Node2Node[j];
                    // Get the local Q's
                    lQ[0] = Q1[nid];
                    lQ[1] = Q2[nid];
                    lQ[2] = Q3[nid];
                    lQ[3] = Q4[nid];
                    lQ[4] = Q5[nid];
                    // Compute Local Equation of State
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    mach_max = MAX(mach_max, lmach);
                    dp_max   = MAX(dp_max, fabs(p - lp));
                }
            }
        }

        // STEP 2:
        // Compute the required parameters
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            beta = MAX(mach, mach_max);
            if (beta < Ref_Mach) {
                beta = MAX(5.0e-5, MAX(sqrt(1.0e-11/sqrt(q2)), sqrt(Ref_Mach*beta)));
                beta = MIN(beta, Ref_Mach);
            }
            beta = MIN(1.0, beta);
//            beta = MAX(beta, sqrt(1.0e-12/sqrt(q2)));
//            beta = MAX(beta, 2.0*dp_max/(rho*c*c));
//            beta = MAX(MAX(beta, 1.0e-5), 0.01*Ref_Mach);
//            beta = MIN(1.0, MIN(beta, Ref_Mach));
//            if (beta > 0.7) {
//                beta = 1.0;
//            } else {
//                beta = 2.0*beta*beta;
//                beta = beta/(1.0 - beta);
//                beta = sqrt(beta/(1 + (Gamma - 1.0)*beta));
//            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        beta = beta*beta;
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
        MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));
        
        
        // STEP 3:
        // Compute P*Minv
        Roe_Minv[0][0] = 1.0;
        Roe_Minv[0][1] = 0.0;
        Roe_Minv[0][2] = 0.0;
        Roe_Minv[0][3] = 0.0;
        Roe_Minv[0][4] = 0.0;

        Roe_Minv[1][0] = -u/rho;
        Roe_Minv[1][1] = 1.0/rho;
        Roe_Minv[1][2] = 0.0;
        Roe_Minv[1][3] = 0.0;
        Roe_Minv[1][4] = 0.0;

        Roe_Minv[2][0] = -v/rho;
        Roe_Minv[2][1] = 0.0;
        Roe_Minv[2][2] = 1.0/rho;
        Roe_Minv[2][3] = 0.0;
        Roe_Minv[2][4] = 0.0;

        Roe_Minv[3][0] = -w/rho;
        Roe_Minv[3][1] = 0.0;
        Roe_Minv[3][2] = 0.0;
        Roe_Minv[3][3] = 1.0/rho;
        Roe_Minv[3][4] = 0.0;

        Roe_Minv[4][0] = beta*phi;
        Roe_Minv[4][1] = -beta*u*(Gamma - 1.0);
        Roe_Minv[4][2] = -beta*v*(Gamma - 1.0);
        Roe_Minv[4][3] = -beta*w*(Gamma - 1.0);
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                Roe_Minv[4][4] = beta*(Gamma - 1.0);
                break;
            case NONDIMENSIONAL_BTW:
                Roe_Minv[4][4] = beta/(Ref_Mach*Ref_Mach);
                break;
            case NONDIMENSIONAL_LMROE:
                Roe_Minv[4][4] = beta*(Gamma - 1.0);
                break;
            default:
                error("Compute_Precondition_Residual_Roe_BTW: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }
        
        // STEP 4:
        // Compute Preconditioned Convective Flux
        // Get the Convective Flux
        flux_roe_conv[0] = Res1[i];
        flux_roe_conv[1] = Res2[i];
        flux_roe_conv[2] = Res3[i];
        flux_roe_conv[3] = Res4[i];
        flux_roe_conv[4] = Res5[i];
        if (PrecondMethod == SOLVER_PRECOND_ROE_BTW_ORIGINAL) {
            flux_roe_conv[0] += Res1_Diss[i];
            flux_roe_conv[1] += Res2_Diss[i];
            flux_roe_conv[2] += Res3_Diss[i];
            flux_roe_conv[3] += Res4_Diss[i];
            flux_roe_conv[4] += Res5_Diss[i];
        }

        // Finally Compute precondition residual
        flux_roe[0] = 0.0;
        flux_roe[1] = 0.0;
        flux_roe[2] = 0.0;
        flux_roe[3] = 0.0;
        flux_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_Minv, flux_roe_conv, flux_roe);
        
        // STEP 5:
        // Compute the Preconditioned Roe Flux Residual
        Res1[i] = flux_roe[0];
        Res2[i] = flux_roe[1];
        Res3[i] = flux_roe[2];
        Res4[i] = flux_roe[3];
        Res5[i] = flux_roe[4];
        // Marco Modification
        if (PrecondMethod == SOLVER_PRECOND_ROE_BTW) {
            Res1[i] += Res1_Diss[i];
            Res2[i] += Res2_Diss[i];
            Res3[i] += Res3_Diss[i];
            Res4[i] += Res4_Diss[i];
            Res5[i] += Res5_Diss[i];
        }
    }
}

//------------------------------------------------------------------------------
//! Compute Preconditioned Roe Residual with Eriksson Pre-Conditioner
//! Note Variables are in Primitive Form: Density Velocity Pressure
//------------------------------------------------------------------------------
void Compute_Precondition_Residual_Roe_ERIKSSON(void) {
    int i, j, nid;
    double Q[5], lQ[5];
    double flux_roe[5];
    double flux_roe_conv[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, phi, beta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Multiply by Precondition Matrix and Construct the Flux
    for (i = 0; i < nNode; i++) {
        // Get the Variables
        Q[0] = Q1[i];
        Q[1] = Q2[i];
        Q[2] = Q3[i];
        Q[3] = Q4[i];
        Q[4] = Q5[i];
        
        // Compute Equation of State
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        phi = 0.5*(Gamma - 1.0)*q2;
        
        //======================================================================
        // Compute Precondition of Convective Flux and Assemble Total Flux
        //======================================================================
        // STEP - 1:
        // Compute the Local Max Mach and Max Change in Pressure around a node
        dp_max   = 0.0;
        mach_max = 0.0;
        // Check if Local Precondition is Required
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            // Check if Precondition Smoother is Requested
            if (PrecondSmooth != 0) {
                for (j = crs_IA_Node2Node[i]; j < crs_IA_Node2Node[i+1]; j++) {
                    nid   = crs_JA_Node2Node[j];
                    // Get the local Q's
                    lQ[0] = Q1[nid];
                    lQ[1] = Q2[nid];
                    lQ[2] = Q3[nid];
                    lQ[3] = Q4[nid];
                    lQ[4] = Q5[nid];
                    // Compute Local Equation of State
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    mach_max = MAX(mach_max, lmach);
                    dp_max   = MAX(dp_max, fabs(p - lp));
                }
            }
        }
        
        // STEP 2:
        // Compute the required parameters
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            beta = MAX(mach, mach_max);
            if (beta < Ref_Mach)
                beta = sqrt(Ref_Mach*beta);
            beta = MIN(1.0, beta);
//            beta = MAX(beta, sqrt(1.0e-12/sqrt(q2)));
//            beta = MAX(beta, 2.0*dp_max/(rho*c*c));
//            beta = MAX(MAX(beta, 1.0e-5), 0.01*Ref_Mach);
//            beta = MIN(1.0, MIN(beta, Ref_Mach));
//            if (beta > 0.7) {
//                beta = 1.0;
//            } else {
//                beta = 2.0*beta*beta;
//                beta = beta/(1.0 - beta);
//                beta = sqrt(beta/(1 + (Gamma - 1.0)*beta));
//            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        beta = beta*beta;
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
        MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));
        
        
        // STEP 3:
        // Compute P*Minv
        Roe_Minv[0][0] =  1.0 + (beta - 1.0)*(Gamma - 1.0)*q2/(2.0*c*c);
        Roe_Minv[0][1] = -(beta - 1.0)*(Gamma - 1.0)*u/(c*c);
        Roe_Minv[0][2] = -(beta - 1.0)*(Gamma - 1.0)*v/(c*c);
        Roe_Minv[0][3] = -(beta - 1.0)*(Gamma - 1.0)*w/(c*c);
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                Roe_Minv[0][4] = (beta - 1.0)*(Gamma - 1.0)/(c*c);
                break;
            case NONDIMENSIONAL_BTW:
                Roe_Minv[0][4] = (beta - 1.0)/(c*c*Ref_Mach*Ref_Mach);
                break;
            case NONDIMENSIONAL_LMROE:
                Roe_Minv[0][4] = (beta - 1.0)*(Gamma - 1.0)/(c*c);
                break;
            default:
                error("Compute_Precondition_Residual_Roe_ERIKSSON: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }

        Roe_Minv[1][0] = -u/rho;
        Roe_Minv[1][1] = 1.0/rho;
        Roe_Minv[1][2] = 0.0;
        Roe_Minv[1][3] = 0.0;
        Roe_Minv[1][4] = 0.0;

        Roe_Minv[2][0] = -v/rho;
        Roe_Minv[2][1] = 0.0;
        Roe_Minv[2][2] = 1.0/rho;
        Roe_Minv[2][3] = 0.0;
        Roe_Minv[2][4] = 0.0;

        Roe_Minv[3][0] = -w/rho;
        Roe_Minv[3][1] = 0.0;
        Roe_Minv[3][2] = 0.0;
        Roe_Minv[3][3] = 1.0/rho;
        Roe_Minv[3][4] = 0.0;

        Roe_Minv[4][0] =  beta*phi;
        Roe_Minv[4][1] = -beta*u*(Gamma - 1.0);
        Roe_Minv[4][2] = -beta*v*(Gamma - 1.0);
        Roe_Minv[4][3] = -beta*w*(Gamma - 1.0);
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                Roe_Minv[4][4] =  beta*(Gamma - 1.0);
                break;
            case NONDIMENSIONAL_BTW:
                Roe_Minv[4][4] =  beta/(Ref_Mach*Ref_Mach);
                break;
            case NONDIMENSIONAL_LMROE:
                Roe_Minv[4][4] =  beta*(Gamma - 1.0);
                break;
            default:
                error("Compute_Precondition_Residual_Roe_ERIKSSON: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }
        
        // STEP 4:
        // Compute Preconditioned Convective Flux
        // Get the Convective Flux
        flux_roe_conv[0] = Res1[i];
        flux_roe_conv[1] = Res2[i];
        flux_roe_conv[2] = Res3[i];
        flux_roe_conv[3] = Res4[i];
        flux_roe_conv[4] = Res5[i];
        if (PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL) {
            flux_roe_conv[0] += Res1_Diss[i];
            flux_roe_conv[1] += Res2_Diss[i];
            flux_roe_conv[2] += Res3_Diss[i];
            flux_roe_conv[3] += Res4_Diss[i];
            flux_roe_conv[4] += Res5_Diss[i];
        }

        // Finally Compute precondition residual
        flux_roe[0] = 0.0;
        flux_roe[1] = 0.0;
        flux_roe[2] = 0.0;
        flux_roe[3] = 0.0;
        flux_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_Minv, flux_roe_conv, flux_roe);
        
        // STEP 5:
        // Compute the Preconditioned Roe Flux Residual
        Res1[i] = flux_roe[0];
        Res2[i] = flux_roe[1];
        Res3[i] = flux_roe[2];
        Res4[i] = flux_roe[3];
        Res5[i] = flux_roe[4];
        // Marco Modification
        if (PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON) {
            Res1[i] += Res1_Diss[i];
            Res2[i] += Res2_Diss[i];
            Res3[i] += Res3_Diss[i];
            Res4[i] += Res4_Diss[i];
            Res5[i] += Res5_Diss[i];
        }
    }
}

//------------------------------------------------------------------------------
//! Compute Transformed Roe Residual in Conservative to Primitive Form
//------------------------------------------------------------------------------
void Compute_Transform_Residual_Roe_ConservativeToPrimitive(void) {
    int i;
    double Q[5];
    double flux_roe[5];
    double flux_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach;
    double dtmp;
    
    // Transform the Residual in Conservative Form to Primitive Form
    // Primitive Variable Formulation
    if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
        for (i = 0; i < nNode; i++) {
            // Get the Variables
            Q[0] = Q1[i];
            Q[1] = Q2[i];
            Q[2] = Q3[i];
            Q[3] = Q4[i];
            Q[4] = Q5[i];

            // Compute Equation of State
            Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
            
            // Compute the Transformation Matrix
            dtmp = 0.5*(Gamma - 1.0)*q2;
            Roe_P[0][0] = dtmp;
            Roe_P[0][1] = (1.0 - Gamma)*u;
            Roe_P[0][2] = (1.0 - Gamma)*v;
            Roe_P[0][3] = (1.0 - Gamma)*w;
            Roe_P[0][4] = (Gamma - 1.0);

            Roe_P[1][0] = -u/rho;
            Roe_P[1][1] = 1.0/rho;
            Roe_P[1][2] = 0.0;
            Roe_P[1][3] = 0.0;
            Roe_P[1][4] = 0.0;

            Roe_P[2][0] = -v/rho;
            Roe_P[2][1] = 0.0;
            Roe_P[2][2] = 1.0/rho;
            Roe_P[2][3] = 0.0;
            Roe_P[2][4] = 0.0;

            Roe_P[3][0] = -w/rho;
            Roe_P[3][1] = 0.0;
            Roe_P[3][2] = 0.0;
            Roe_P[3][3] = 1.0/rho;
            Roe_P[3][4] = 0.0;

            Roe_P[4][0] = (dtmp*Gamma - T)/rho;
            Roe_P[4][1] = (1.0 - Gamma)*Gamma*u/rho;
            Roe_P[4][2] = (1.0 - Gamma)*Gamma*v/rho;
            Roe_P[4][3] = (1.0 - Gamma)*Gamma*w/rho;
            Roe_P[4][4] = (Gamma - 1.0)*Gamma/rho;

            // Get the Conservative Residual
            flux_roe[0] = Res1[i];
            flux_roe[1] = Res2[i];
            flux_roe[2] = Res3[i];
            flux_roe[3] = Res4[i];
            flux_roe[4] = Res5[i];
            MC_Matrix_Mul_Vector(5, 5, Roe_P, flux_roe, flux_roe_diss);

            // Update the Residual in Primitive Form
            Res1[i] = flux_roe_diss[0];
            Res2[i] = flux_roe_diss[1];
            Res3[i] = flux_roe_diss[2];
            Res4[i] = flux_roe_diss[3];
            Res5[i] = flux_roe_diss[4]; 
        }
    } else
        error("Compute_Transform_Residual_Roe_ConservativeToPrimitive: Invalid Transformation - %d", Variable_Type);
}

//------------------------------------------------------------------------------
//! Compute Precondition Residual and Transform the Residual Accordingly
//------------------------------------------------------------------------------
void Compute_Precondition_Residual_Roe(void) {
    switch (PrecondMethod) {
        case SOLVER_PRECOND_NONE: // Roe
            // No Precondition but Transformation is needed for primitive variables
            // Primitive Variable Formulation
            if (Variable_Type == VARIABLE_PRIMITIVE_PUT)
                Compute_Transform_Residual_Roe_ConservativeToPrimitive();
            break;
        case SOLVER_PRECOND_ROE_LMFIX: // LMRoe
            // No Precondition but Transformation is needed for primitive variables
            // Primitive Variable Formulation
            if (Variable_Type == VARIABLE_PRIMITIVE_PUT)
                Compute_Transform_Residual_Roe_ConservativeToPrimitive();
            break;
        case SOLVER_PRECOND_ROE_WS: // Roe Weiss Smith Pre-Conditioner and Transformation
            Compute_Precondition_Residual_Roe_WeissSmith();
            break;
        case SOLVER_PRECOND_ROE_CV: // Roe Cecile Voizat Pre-Conditioner and Transformation
            Compute_Precondition_Residual_Roe_CecileVoizat();
            break;
        case SOLVER_PRECOND_ROE_CV_ORIGINAL: // Roe Cecile Voizat Pre-Conditioner and Transformation
            Compute_Precondition_Residual_Roe_CecileVoizat();
            break;
        case SOLVER_PRECOND_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner and Transformation
            Compute_Precondition_Residual_Roe_BTW();
            break;
        case SOLVER_PRECOND_ROE_BTW_ORIGINAL: // Roe Briley Taylor Whitfield Pre-Conditioner and Transformation
            Compute_Precondition_Residual_Roe_BTW();
            break;
        case SOLVER_PRECOND_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner and Transformation
            Compute_Precondition_Residual_Roe_ERIKSSON();
            break;
        case SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL: // Roe Eriksson Pre-Conditioner and Transformation
            Compute_Precondition_Residual_Roe_ERIKSSON();
            break;
        default:
            error("Compute_Precondition_Residual_Roe: Invalid Solver Precondition Scheme - %d", PrecondMethod);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Cecile Voizat Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_RoeFlux_CecileVoizat(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Internal Nodes
        if (node_R < nNode) {
            // Make Solution Second Order
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, Roe_Q_L, Roe_Q_R);
            } else {
                Roe_Q_L[0] = Q1[node_L];
                Roe_Q_L[1] = Q2[node_L];
                Roe_Q_L[2] = Q3[node_L];
                Roe_Q_L[3] = Q4[node_L];
                Roe_Q_L[4] = Q5[node_L];
                Roe_Q_R[0] = Q1[node_R];
                Roe_Q_R[1] = Q2[node_R];
                Roe_Q_R[2] = Q3[node_R];
                Roe_Q_R[3] = Q4[node_R];
                Roe_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            Roe_Q_L[0] = Q1[node_L];
            Roe_Q_L[1] = Q2[node_L];
            Roe_Q_L[2] = Q3[node_L];
            Roe_Q_L[3] = Q4[node_L];
            Roe_Q_L[4] = Q5[node_L];
            Roe_Q_R[0] = Q1[node_R];
            Roe_Q_R[1] = Q2[node_R];
            Roe_Q_R[2] = Q3[node_R];
            Roe_Q_R[3] = Q4[node_R];
            Roe_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // Compute flux_L
        Roe_flux_L[0] = rho_L*ubar_L;
        Roe_flux_L[1] =   u_L*Roe_flux_L[0] + p_L*nx;
        Roe_flux_L[2] =   v_L*Roe_flux_L[0] + p_L*ny;
        Roe_flux_L[3] =   w_L*Roe_flux_L[0] + p_L*nz;
        Roe_flux_L[4] =  ht_L*Roe_flux_L[0];
        
        // Compute flux_R
        Roe_flux_R[0] = rho_R*ubar_R;
        Roe_flux_R[1] =   u_R*Roe_flux_R[0] + p_R*nx;
        Roe_flux_R[2] =   v_R*Roe_flux_R[0] + p_R*ny;
        Roe_flux_R[3] =   w_R*Roe_flux_R[0] + p_R*nz;
        Roe_flux_R[4] =  ht_R*Roe_flux_R[0];
        
        AvgType = 1;
        if (AvgType == 1) {
            // ROE AVERAGE VARIABLES
            rho   = sqrt(rho_R * rho_L);
            sigma = rho/(rho_L + rho);
            u     = u_L  + sigma*(u_R  - u_L);
            v     = v_L  + sigma*(v_R  - v_L);
            w     = w_L  + sigma*(w_R  - w_L);
            ht    = ht_L + sigma*(ht_R - ht_L);
        } else {
            // SIMPLE AVERAGE VARIABLES
            rho   = 0.5*(rho_R + rho_L);
            u     = 0.5*(u_R  + u_L);
            v     = 0.5*(v_R  + v_L);
            w     = 0.5*(w_R  + w_L);
            ht    = 0.5*(ht_R + ht_L);
        }
        q2   = u*u + v*v + w*w;
        phi  = 0.5*(Gamma - 1.0)*q2;
        ubar = u*nx + v*ny + w*nz;
        c    = 0.0;
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                c = (Gamma - 1.0)*ht - phi;
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_BTW:
                c = ht/(Ref_Mach*Ref_Mach) - phi;
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_LMROE:
                c = (Gamma - 1.0)*ht - phi;
                c = sqrt(c);
                break;
            default:
                error("Compute_RoeFlux_CecileVoizat: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }
        mach = sqrt(q2)/c;
        
        if (isnan(c))
            printf("CV Hell %10.5e\n", ubar);
        
        // Compute dQ
        // Conservative Variable Formulation
        if (Variable_Type == VARIABLE_CONSERVATIVE) {
            Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
            Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
            Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
            Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
            Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        }
        // Primitive Variable Formulation
        if ((Variable_Type == VARIABLE_PRIMITIVE_PUT) || (Variable_Type == VARIABLE_PRIMITIVE_RUP)) {
            Roe_dQ[0] = rho_R      - rho_L;
            Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
            Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
            Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
            Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
        }
        
        int nid;
        double lQ[NEQUATIONS];
        double alpha, beta, zeta;
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        
        //======================================================================
        // Compute Precondition Dissipation Flux
        //======================================================================
        // STEP - 1:
        // Compute the Local Max Mach and Max Change in Pressure around a node L and R
        // If node R is not ghost compute average
        dp_L     = dp_R    = 0.0;
        lmach_L  = lmach_R = 0.0;
        dp_max   = 0.0;
        mach_max = 0.0;
        // Check if Local Precondition is Required
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            // Check of Precondition Smooth is Required
            if (PrecondSmooth != 0) {
                for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L+1]; i++) {
                    nid   = crs_JA_Node2Node[i];
                    // Get the local Q's
                    lQ[0] = Q1[nid];
                    lQ[1] = Q2[nid];
                    lQ[2] = Q3[nid];
                    lQ[3] = Q4[nid];
                    lQ[4] = Q5[nid];
                    // Compute Local Equation of State
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    lmach_L = MAX(lmach_L, lmach);
                    dp_L    = MAX(dp_L, fabs(p_L - lp));
                }
                dp_max   = dp_L;
                mach_max = lmach_L;
                // Check if Right Node is not a Ghost Boundary Node
                if (node_R < nNode) {
                    for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R+1]; i++) {
                        nid  = crs_JA_Node2Node[i];
                        // Get the local Q's
                        lQ[0] = Q1[nid];
                        lQ[1] = Q2[nid];
                        lQ[2] = Q3[nid];
                        lQ[3] = Q4[nid];
                        lQ[4] = Q5[nid];
                        // Compute Local Equation of State
                        Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                        // Compute Smoothing Parameters
                        lmach_R = MAX(lmach_R, lmach);
                        dp_R    = MAX(dp_R, fabs(p_R - lp));
                    }
                    dp_max   = 0.5*(dp_L + dp_R);
                    mach_max = 0.5*(lmach_L + lmach_R);
                }
            } else {
                mach_max = MAX(mach_L, mach_R);
                dp_max   = fabs(p_L - p_R);
            }
        }
        
        // STEP 2:
        // Compute the required parameters
        sigma = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            sigma = MAX(mach, mach_max);
            if (sigma < Ref_Mach)
                sigma = sqrt(Ref_Mach*sigma);
            sigma = MIN(1.0, sigma);
//            //if (Ref_Mach < 0.5) {
//                sigma = MIN(sqrt(1.0e-8/fabs(ubar)), Ref_Mach);
//                // Check if Precondition Smoother is Requested
//                //if (PrecondSmooth != 0)
//                    sigma = MIN(1.0, MAX(MAX(MAX(sigma, sqrt(q2)/c), mach_max), 2.0*dp_max/(rho*c*c)));
//                //else
//                //    sigma = MIN(1.0, MAX(sigma, sqrt(q2)/c));
//            //sigma = MAX(sqrt(q2)/c, dp_max/(rho*c*c));
//            //sigma = MAX(1.0e-6, MIN(mach_max, dp_max/(Inf_Mach*Inf_Mach*rho*c*c)));
//                if (sigma > Ref_Mach)
//                    sigma = Ref_Mach;
//            //}
            if (node_L < nNode && node_R < nNode) {
                PrecondSigma[node_L] += sigma;
                PrecondSigma[node_R] += sigma;
            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            sigma = MIN(1.0, Ref_Mach);
        
        alpha = ubar*ubar*(sigma*sigma*sigma*sigma - 2.0*sigma*sigma + 1.0) + 4.0*c*c*sigma*sigma;
        alpha = sqrt(alpha);
        beta  = (sigma*sigma - 1.0)*ubar + alpha;
        zeta  = (1.0 - sigma*sigma)*ubar + alpha;
        
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(0.5*(ubar + sigma*sigma*ubar + alpha));
        Roe_Eigen[4][4] = fabs(0.5*(ubar + sigma*sigma*ubar - alpha));
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Min and Max EigenValue Value
            MinEigenLamda1 = MIN(MinEigenLamda1, Roe_Eigen[0][0]);
            MaxEigenLamda1 = MAX(MaxEigenLamda1, Roe_Eigen[0][0]);
            MinEigenLamda4 = MIN(MinEigenLamda4, Roe_Eigen[3][3]);
            MaxEigenLamda4 = MAX(MaxEigenLamda4, Roe_Eigen[3][3]);
            MinEigenLamda5 = MIN(MinEigenLamda5, Roe_Eigen[4][4]);
            MaxEigenLamda5 = MAX(MaxEigenLamda5, Roe_Eigen[4][4]);

            // Min and Max Precondition Variable
            MinPrecondSigma = MIN(MinPrecondSigma, sigma);
            MaxPrecondSigma = MAX(MaxPrecondSigma, sigma);

            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothType == RESIDUAL_SMOOTH_IMPLICIT) {
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L + 1]; i++) {
                // Get the Node ID and Match with node_R
                if (crs_JA_Node2Node[i] == node_R) {
                    RM_MaxEigenValue[i] = maxlambda*area;
                    break;
                }
            }
            RM_SumMaxEigenValue[node_L] += maxlambda*area;
            // Only Physical Nodes
            if (node_R < nNode) {
                for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R + 1]; i++) {
                    // Get the Node ID and Match with node_L
                    if (crs_JA_Node2Node[i] == node_L) {
                        RM_MaxEigenValue[i] = maxlambda*area;
                        break;
                    }
                }
                RM_SumMaxEigenValue[node_R] += maxlambda*area;
            }
        }
        
        // STEP 4:
        // Marco Modification
        if (PrecondMethod == SOLVER_PRECOND_ROE_CV) {
            // Compute Tmod_g = M3.P
            Roe_T[0][0] = nx;
            Roe_T[0][1] = ny;
            Roe_T[0][2] = nz;
            Roe_T[0][3] = 2.0*rho*sigma*sigma/zeta;
            Roe_T[0][4] = 2.0*rho*sigma*sigma/beta;

            Roe_T[1][0] = u*nx;
            Roe_T[1][1] = u*ny - rho*nz;
            Roe_T[1][2] = rho*ny + u*nz;
            Roe_T[1][3] = rho*(nx + 2.0*u*sigma*sigma/zeta);
            Roe_T[1][4] = rho*(-nx + 2.0*u*sigma*sigma/beta);

            Roe_T[2][0] = v*nx + rho*nz;
            Roe_T[2][1] = v*ny;
            Roe_T[2][2] = -rho*nx + v*nz;
            Roe_T[2][3] = rho*(ny + 2.0*v*sigma*sigma/zeta);
            Roe_T[2][4] = rho*(-ny + 2.0*v*sigma*sigma/beta);

            Roe_T[3][0] = w*nx - rho*ny;
            Roe_T[3][1] = rho*nx + w*ny;
            Roe_T[3][2] = w*nz;
            Roe_T[3][3] = rho*(nz+ 2.0*w*sigma*sigma/zeta);
            Roe_T[3][4] = rho*(-nz + 2.0*w*sigma*sigma/beta);

            Roe_T[4][0] = -w*rho*ny + v*rho*nz + 0.5*nx*q2;
            Roe_T[4][1] =  w*rho*nx - u*rho*nz + 0.5*ny*q2;
            Roe_T[4][2] = -v*rho*nx + u*rho*ny + 0.5*nz*q2;
            Roe_T[4][3] = rho*(ubar + sigma*sigma*(2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*zeta));
            Roe_T[4][4] = rho*(-ubar + sigma*sigma*(2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*beta));
        }
        
        // Original 
        if (PrecondMethod == SOLVER_PRECOND_ROE_CV_ORIGINAL) {
            // Compute T_g = M3.Sinv.Pinv
            Roe_T[0][0] = nx;
            Roe_T[0][1] = ny;
            Roe_T[0][2] = nz;
            Roe_T[0][3] = 2.0*rho/zeta;
            Roe_T[0][4] = 2.0*rho/beta;

            Roe_T[1][0] = u*nx;
            Roe_T[1][1] = u*ny - rho*nz;
            Roe_T[1][2] = rho*ny + u*nz;
            Roe_T[1][3] = rho*(nx + 2.0*u/zeta);
            Roe_T[1][4] = rho*(-nx + 2.0*u/beta);

            Roe_T[2][0] = v*nx + rho*nz;
            Roe_T[2][1] = v*ny;
            Roe_T[2][2] = -rho*nx + v*nz;
            Roe_T[2][3] = rho*(ny + 2.0*v/zeta);
            Roe_T[2][4] = rho*(-ny + 2.0*v/beta);

            Roe_T[3][0] = w*nx - rho*ny;
            Roe_T[3][1] = rho*nx + w*ny;
            Roe_T[3][2] = w*nz;
            Roe_T[3][3] = rho*(nz+ 2.0*w/zeta);
            Roe_T[3][4] = rho*(-nz + 2.0*w/beta);

            Roe_T[4][0] = -w*rho*ny + v*rho*nz + 0.5*nx*q2;
            Roe_T[4][1] =  w*rho*nx - u*rho*nz + 0.5*ny*q2;
            Roe_T[4][2] = -v*rho*nx + u*rho*ny + 0.5*nz*q2;
            Roe_T[4][3] = rho*(ubar + (2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*zeta));
            Roe_T[4][4] = rho*(-ubar + (2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*beta));
        }
        
        // STEP 5:
        // Compute T_d = Pinv.M3inv
        Roe_Tinv[0][0] = (w*ny - v*nz)/rho + (nx*(2.0*c*c - (Gamma - 1.0)*q2))/(2.0*c*c);
        Roe_Tinv[0][1] = (u*(Gamma - 1.0)*nx)/(c*c);
        Roe_Tinv[0][2] = (v*(Gamma - 1.0)*nx)/(c*c) + nz/rho;
        Roe_Tinv[0][3] = (w*(Gamma - 1.0)*nx)/(c*c) - ny/rho;
        Roe_Tinv[0][4] = -(Gamma - 1.0)*nx/(c*c);
        
        Roe_Tinv[1][0] = (u*nz - w*nx)/rho + (ny*(2.0*c*c - (Gamma - 1.0)*q2))/(2.0*c*c);
        Roe_Tinv[1][1] = (u*(Gamma - 1.0)*ny)/(c*c) - nz/rho;
        Roe_Tinv[1][2] = (v*(Gamma - 1.0)*ny)/(c*c);
        Roe_Tinv[1][3] = (w*(Gamma - 1.0)*ny)/(c*c) + nx/rho;
        Roe_Tinv[1][4] = -(Gamma - 1.0)*ny/(c*c);
        
        Roe_Tinv[2][0] = (v*nx - u*ny)/rho + (nz*(2.0*c*c - (Gamma - 1.0)*q2))/(2.0*c*c);
        Roe_Tinv[2][1] = (u*(Gamma - 1.0)*nz)/(c*c) + ny/rho;
        Roe_Tinv[2][2] = (v*(Gamma - 1.0)*nz)/(c*c) - nx/rho;
        Roe_Tinv[2][3] = (w*(Gamma - 1.0)*nz)/(c*c);
        Roe_Tinv[2][4] = -(Gamma - 1.0)*nz/(c*c);
        
        Roe_Tinv[3][0] = ((Gamma - 1.0)*q2 - ubar*zeta)/(2.0*rho*alpha);
        Roe_Tinv[3][1] = (-2.0*u*(Gamma - 1.0) + nx*zeta)/(2.0*rho*alpha);
        Roe_Tinv[3][2] = (-2.0*v*(Gamma - 1.0) + ny*zeta)/(2.0*rho*alpha);
        Roe_Tinv[3][3] = (-2.0*w*(Gamma - 1.0) + nz*zeta)/(2.0*rho*alpha);
        Roe_Tinv[3][4] = (Gamma - 1.0)/(rho*alpha);
        
        Roe_Tinv[4][0] = ((Gamma - 1.0)*q2 + ubar*beta)/(2.0*rho*alpha);
        Roe_Tinv[4][1] = (-2.0*u*(Gamma - 1.0) - nx*beta)/(2.0*rho*alpha);
        Roe_Tinv[4][2] = (-2.0*v*(Gamma - 1.0) - ny*beta)/(2.0*rho*alpha);
        Roe_Tinv[4][3] = (-2.0*w*(Gamma - 1.0) - nz*beta)/(2.0*rho*alpha);
        Roe_Tinv[4][4] = (Gamma - 1.0)/(rho*alpha);
        
        // STEP 6:
        // Compute the Dissipation Flux
        // Tmp => P = |Lambda|*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);

        // Get Matrix |A| = T*|Lambda|*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);

        // Compute |A|*dQ
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 7:
        // Compute the Flux_Roe_Conv and Flux_Roe_Diss
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux_CecileVoizat: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Briley Taylor Whitfield Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_RoeFlux_Precondition_BTW(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Internal Nodes
        if (node_R < nNode) {
            // Make Solution Second Order
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, Roe_Q_L, Roe_Q_R);
            } else {
                Roe_Q_L[0] = Q1[node_L];
                Roe_Q_L[1] = Q2[node_L];
                Roe_Q_L[2] = Q3[node_L];
                Roe_Q_L[3] = Q4[node_L];
                Roe_Q_L[4] = Q5[node_L];
                Roe_Q_R[0] = Q1[node_R];
                Roe_Q_R[1] = Q2[node_R];
                Roe_Q_R[2] = Q3[node_R];
                Roe_Q_R[3] = Q4[node_R];
                Roe_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            Roe_Q_L[0] = Q1[node_L];
            Roe_Q_L[1] = Q2[node_L];
            Roe_Q_L[2] = Q3[node_L];
            Roe_Q_L[3] = Q4[node_L];
            Roe_Q_L[4] = Q5[node_L];
            Roe_Q_R[0] = Q1[node_R];
            Roe_Q_R[1] = Q2[node_R];
            Roe_Q_R[2] = Q3[node_R];
            Roe_Q_R[3] = Q4[node_R];
            Roe_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // Compute flux_L
        Roe_flux_L[0] = rho_L*ubar_L;
        Roe_flux_L[1] =   u_L*Roe_flux_L[0] + p_L*nx;
        Roe_flux_L[2] =   v_L*Roe_flux_L[0] + p_L*ny;
        Roe_flux_L[3] =   w_L*Roe_flux_L[0] + p_L*nz;
        Roe_flux_L[4] =  ht_L*Roe_flux_L[0];
        
        // Compute flux_R
        Roe_flux_R[0] = rho_R*ubar_R;
        Roe_flux_R[1] =   u_R*Roe_flux_R[0] + p_R*nx;
        Roe_flux_R[2] =   v_R*Roe_flux_R[0] + p_R*ny;
        Roe_flux_R[3] =   w_R*Roe_flux_R[0] + p_R*nz;
        Roe_flux_R[4] =  ht_R*Roe_flux_R[0];
        
        AvgType = 0;
        if (AvgType == 1) {
            // ROE AVERAGE VARIABLES
            rho   = sqrt(rho_R * rho_L);
            sigma = rho/(rho_L + rho);
            u     = u_L  + sigma*(u_R  - u_L);
            v     = v_L  + sigma*(v_R  - v_L);
            w     = w_L  + sigma*(w_R  - w_L);
            ht    = ht_L + sigma*(ht_R - ht_L);
        } else {
            // SIMPLE AVERAGE VARIABLES
            rho   = 0.5*(rho_R + rho_L);
            u     = 0.5*(u_R  + u_L);
            v     = 0.5*(v_R  + v_L);
            w     = 0.5*(w_R  + w_L);
            ht    = 0.5*(ht_R + ht_L);
        }
        q2   = u*u + v*v + w*w;
        phi  = 0.5*(Gamma - 1.0)*q2;
        ubar = u*nx + v*ny + w*nz;
        c    = 0.0;
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                c = (Gamma - 1.0)*ht - phi;
                if (c < 0.0 && SolverIteration < 2000)
                    c = fabs(c);
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_BTW:
                c = ht/(Ref_Mach*Ref_Mach) - phi;
                if (c < 0.0 && SolverIteration < 2000)
                    c = fabs(c);
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_LMROE:
                c = (Gamma - 1.0)*ht - phi;
                c = sqrt(c);
                break;
            default:
                error("Compute_RoeFlux_Precondition_BTW: Undefined Non-Dimensional Type - %d -1", NonDimensional_Type);
                break;
        }
        mach = sqrt(q2)/c;
        
        if (isnan(c))
            printf("Hell %10.5e\n", ubar);
        
        // Compute dQ
        // Conservative Variable Formulation
        if (Variable_Type == VARIABLE_CONSERVATIVE) {
            Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
            Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
            Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
            Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
            Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        }
        // Primitive Variable Formulation
        if (Variable_Type == VARIABLE_PRIMITIVE_PUT || Variable_Type == VARIABLE_PRIMITIVE_RUP) {
            Roe_dQ[0] = rho_R      - rho_L;
            Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
            Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
            Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
            Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
        }
        
        // Compute dw
        // Primitive Variable Formulation Pressure Velocity Temperature
        if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
            Roe_dw[0] = p_R - p_L;
            Roe_dw[1] = u_R - u_L;
            Roe_dw[2] = v_R - v_L;
            Roe_dw[3] = w_R - w_L;
            Roe_dw[4] = T_R - T_L;
        }
        // Primitive Variable Formulation Density Velocity Pressure
        if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
            Roe_dw[0] = rho_R - rho_L;
            Roe_dw[1] = u_R   - u_L;
            Roe_dw[2] = v_R   - v_L;
            Roe_dw[3] = w_R   - w_L;
            Roe_dw[4] = p_R   - p_L;
        }
        
        int nid;
        double lQ[NEQUATIONS];
        double tau, omega, beta, beta_m, beta_p;
        double chi_m, chi_p, psi_m, psi_p;
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        
        //======================================================================
        // Compute Precondition Dissipation Flux
        //======================================================================
        // STEP - 1:
        // Compute the Local Max Mach and Max Change in Pressure around a node L and R
        // If node R is not ghost compute average
        dp_L     = dp_R    = 0.0;
        lmach_L  = lmach_R = 0.0;
        dp_max   = 0.0;
        mach_max = 0.0;
        mach_max = 0.0;
        // Check if Local Precondition is Required
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            // Check of Precondition Smooth is Required
            if (PrecondSmooth != 0) {
                for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L+1]; i++) {
                    nid   = crs_JA_Node2Node[i];
                    // Get the local Q's
                    lQ[0] = Q1[nid];
                    lQ[1] = Q2[nid];
                    lQ[2] = Q3[nid];
                    lQ[3] = Q4[nid];
                    lQ[4] = Q5[nid];
                    // Compute Local Equation of State
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    lmach_L = MAX(lmach_L, lmach);
                    dp_L    = MAX(dp_L, fabs(p_L - lp));
                }
                dp_max   = dp_L;
                mach_max = lmach_L;
                // Check if Right Node is not a Ghost Boundary Node
                if (node_R < nNode) {
                    for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R+1]; i++) {
                        nid  = crs_JA_Node2Node[i];
                        // Get the local Q's
                        lQ[0] = Q1[nid];
                        lQ[1] = Q2[nid];
                        lQ[2] = Q3[nid];
                        lQ[3] = Q4[nid];
                        lQ[4] = Q5[nid];
                        // Compute Local Equation of State
                        Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                        // Compute Smoothing Parameters
                        lmach_R = MAX(lmach_R, lmach);
                        dp_R    = MAX(dp_R, fabs(p_R - lp));
                    }
                    dp_max   = 0.5*(dp_L + dp_R);
                    mach_max = 0.5*(lmach_L + lmach_R);
                }
            } else {
                mach_max = MAX(mach_L, mach_R);
                dp_max   = fabs(p_L - p_R);
            }
        }
        
        // STEP 2:
        // Compute the required parameters
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            beta = MAX(mach, mach_max);
            if (beta < Ref_Mach) {
                beta = MAX(5.0e-5, MAX(sqrt(1.0e-11/sqrt(q2)), sqrt(Ref_Mach*beta)));
                beta = MIN(beta, Ref_Mach);
            }
            beta = MIN(1.0, beta);
//            beta = MAX(beta, sqrt(1.0e-12/sqrt(q2)));
//            beta = MAX(beta, 2.0*dp_max/(rho*c*c));
//            beta = MAX(MAX(beta, 1.0e-5), 0.01*Ref_Mach);
//            beta = MIN(1.0, MIN(beta, Ref_Mach));
//            if (beta > 0.7) {
//                beta = 1.0;
//            } else {
//                beta = 2.0*beta*beta;
//                beta = beta/(1.0 - beta);
//                beta = sqrt(beta/(1 + (Gamma - 1.0)*beta));
//            }
            if (node_L < nNode && node_R < nNode) {
                PrecondSigma[node_L] += beta;
                PrecondSigma[node_R] += beta;
            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        beta   = beta*beta;
        beta_m = 0.5*(1.0 - beta);
        beta_p = 0.5*(1.0 + beta);
        tau    = beta*c*c;
        sigma  = sqrt(ubar*ubar*beta_m*beta_m + tau);
        chi_m  = ubar*beta_m - sigma;
        chi_p  = ubar*beta_m + sigma;
        omega  = 2.0*ubar*rho*beta_m/tau;
        psi_m  = 0.5*chi_m/sigma;
        psi_p  = 0.5*chi_p/sigma;
        
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar*beta_p + sigma);
        Roe_Eigen[4][4] = fabs(ubar*beta_p - sigma);
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Min and Max EigenValue Value
            MinEigenLamda1 = MIN(MinEigenLamda1, Roe_Eigen[0][0]);
            MaxEigenLamda1 = MAX(MaxEigenLamda1, Roe_Eigen[0][0]);
            MinEigenLamda4 = MIN(MinEigenLamda4, Roe_Eigen[3][3]);
            MaxEigenLamda4 = MAX(MaxEigenLamda4, Roe_Eigen[3][3]);
            MinEigenLamda5 = MIN(MinEigenLamda5, Roe_Eigen[4][4]);
            MaxEigenLamda5 = MAX(MaxEigenLamda5, Roe_Eigen[4][4]);

            // Min and Max Precondition Variable
            MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
            MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));

            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothType == RESIDUAL_SMOOTH_IMPLICIT) {
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L + 1]; i++) {
                // Get the Node ID and Match with node_R
                if (crs_JA_Node2Node[i] == node_R) {
                    RM_MaxEigenValue[i] = maxlambda*area;
                    break;
                }
            }
            RM_SumMaxEigenValue[node_L] += maxlambda*area;
            // Only Physical Nodes
            if (node_R < nNode) {
                for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R + 1]; i++) {
                    // Get the Node ID and Match with node_L
                    if (crs_JA_Node2Node[i] == node_L) {
                        RM_MaxEigenValue[i] = maxlambda*area;
                        break;
                    }
                }
                RM_SumMaxEigenValue[node_R] += maxlambda*area;
            }
        }
        
        // STEP 4:
        // Marco Modification
        if (PrecondMethod == SOLVER_PRECOND_ROE_BTW) {
            // Compute Tmod_g = R
            Roe_T[0][0] =  nx;
            Roe_T[0][1] =  ny;
            Roe_T[0][2] =  nz;
            Roe_T[0][3] =  (rho*chi_p)/tau;
            Roe_T[0][4] = -(rho*chi_m)/tau;

            Roe_T[1][0] =  0.0;
            Roe_T[1][1] = -nz;
            Roe_T[1][2] =  ny;
            Roe_T[1][3] =  nx;
            Roe_T[1][4] = -nx;

            Roe_T[2][0] =  nz;
            Roe_T[2][1] =  0.0;
            Roe_T[2][2] = -nx;
            Roe_T[2][3] =  ny;
            Roe_T[2][4] = -ny;

            Roe_T[3][0] = -ny;
            Roe_T[3][1] =  nx;
            Roe_T[3][2] =  0.0;
            Roe_T[3][3] =  nz;
            Roe_T[3][4] = -nz;

            Roe_T[4][0] =  0.0;
            Roe_T[4][1] =  0.0;
            Roe_T[4][2] =  0.0;
            Roe_T[4][3] = -rho*chi_m;
            Roe_T[4][4] =  rho*chi_p;
        }
        
        // Original
        if (PrecondMethod == SOLVER_PRECOND_ROE_BTW_ORIGINAL) {
            // Compute T_g = M1.Sinv.R
            Roe_T[0][0] =  nx;
            Roe_T[0][1] =  ny;
            Roe_T[0][2] =  nz;
            Roe_T[0][3] =  (rho*chi_p)/tau;
            Roe_T[0][4] = -(rho*chi_m)/tau;

            Roe_T[1][0] =  u*nx;
            Roe_T[1][1] =  u*ny - rho*nz;
            Roe_T[1][2] =  rho*ny + u*nz;
            Roe_T[1][3] =  rho*((u*chi_p)/tau + nx);
            Roe_T[1][4] = -rho*((u*chi_m)/tau + nx);

            Roe_T[2][0] =  v*nx + rho*nz;
            Roe_T[2][1] =  v*ny;
            Roe_T[2][2] = -rho*nx + v*nz;
            Roe_T[2][3] =  rho*((v*chi_p)/tau + ny);
            Roe_T[2][4] = -rho*((v*chi_m)/tau + ny);

            Roe_T[3][0] =  w*nx - rho*ny;
            Roe_T[3][1] =  rho*nx + w*ny;
            Roe_T[3][2] =  w*nz;
            Roe_T[3][3] =  rho*((w*chi_p)/tau + nz);
            Roe_T[3][4] = -rho*((w*chi_m)/tau + nz);

            switch (NonDimensional_Type) {
                case NONDIMENSIONAL_GENERIC:
                    Roe_T[4][0] = -w*rho*ny + v*rho*nz + 0.5*nx*q2;
                    Roe_T[4][1] =  w*rho*nx - u*rho*nz + 0.5*ny*q2;
                    Roe_T[4][2] = -v*rho*nx + u*rho*ny + 0.5*nz*q2;
                    Roe_T[4][3] =  rho*(ubar + (0.5*q2*chi_p)/tau - chi_m/(beta*(Gamma - 1.0)));
                    Roe_T[4][4] =  rho*(chi_p/(beta*(Gamma - 1.0)) - ubar - (0.5*q2*chi_m)/tau);
                    break;
                case NONDIMENSIONAL_BTW:
                    Roe_T[4][0] =  0.5*(Gamma - 1.0)*Ref_Mach*Ref_Mach*(-2.0*w*rho*ny + 2.0*v*rho*nz + nx*q2);
                    Roe_T[4][1] =  0.5*(Gamma - 1.0)*Ref_Mach*Ref_Mach*( 2.0*w*rho*nx - 2.0*u*rho*nz + ny*q2);
                    Roe_T[4][2] =  0.5*(Gamma - 1.0)*Ref_Mach*Ref_Mach*(-2.0*v*rho*nx + 2.0*u*rho*ny + nz*q2);
                    Roe_T[4][3] =  0.5*rho*Ref_Mach*Ref_Mach*((Gamma - 1.0)*(2.0*ubar + (q2*chi_p)/tau) - (2.0*chi_m)/beta);
                    Roe_T[4][4] =  0.5*rho*Ref_Mach*Ref_Mach*((2.0*chi_p)/beta - (Gamma - 1.0)*(2.0*ubar + (q2*chi_m)/tau));
                    break;
                case NONDIMENSIONAL_LMROE:
                    Roe_T[4][0] = -w*rho*ny + v*rho*nz + 0.5*nx*q2;
                    Roe_T[4][1] =  w*rho*nx - u*rho*nz + 0.5*ny*q2;
                    Roe_T[4][2] = -v*rho*nx + u*rho*ny + 0.5*nz*q2;
                    Roe_T[4][3] =  rho*(ubar + (0.5*q2*chi_p)/tau - chi_m/(beta*(Gamma - 1.0)));
                    Roe_T[4][4] =  rho*(chi_p/(beta*(Gamma - 1.0)) - ubar - (0.5*q2*chi_m)/tau);
                    break;
                default:
                    error("Compute_RoeFlux_Precondition_BTW: Undefined Non-Dimensional Type - %d -2", NonDimensional_Type);
                    break;
            }
        }
        
        // STEP 5:
        // Marco Modification or Original
        if (PrecondMethod == SOLVER_PRECOND_ROE_BTW || PrecondMethod == SOLVER_PRECOND_ROE_BTW_ORIGINAL) {
            // Compute T_d = Rinv
            Roe_Tinv[0][0] =  nx;
            Roe_Tinv[0][1] = -omega*nx*nx;
            Roe_Tinv[0][2] = -omega*nx*ny + nz;
            Roe_Tinv[0][3] = -omega*nx*nz - ny;
            Roe_Tinv[0][4] = -nx/tau;

            Roe_Tinv[1][0] =  ny;
            Roe_Tinv[1][1] = -omega*ny*nx - nz;
            Roe_Tinv[1][2] = -omega*ny*ny;
            Roe_Tinv[1][3] = -omega*ny*nz + nx;
            Roe_Tinv[1][4] = -ny/tau;

            Roe_Tinv[2][0] =  nz;
            Roe_Tinv[2][1] = -omega*nz*nx + ny;
            Roe_Tinv[2][2] = -omega*nz*ny - nx;
            Roe_Tinv[2][3] = -omega*nz*nz;
            Roe_Tinv[2][4] = -nz/tau;

            Roe_Tinv[3][0] = 0.0;
            Roe_Tinv[3][1] = nx*psi_p;
            Roe_Tinv[3][2] = ny*psi_p;
            Roe_Tinv[3][3] = nz*psi_p;
            Roe_Tinv[3][4] = 1.0/(2.0*rho*sigma);

            Roe_Tinv[4][0] = 0.0;
            Roe_Tinv[4][1] = nx*psi_m;
            Roe_Tinv[4][2] = ny*psi_m;
            Roe_Tinv[4][3] = nz*psi_m;
            Roe_Tinv[4][4] = 1.0/(2.0*rho*sigma);
        }
        
        // STEP 6:
        // Compute the Dissipation Flux
        // Tmp => P = |Lambda|*Rinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);

        // Compute |A| = R*|Lambda|*Rinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);

        // Compute |A|*dw
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dw, Roe_fluxA);
        
        // STEP 7:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux_Precondition_BTW: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Eriksson Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_RoeFlux_Precondition_ERIKSSON(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Internal Nodes
        if (node_R < nNode) {
            // Make Solution Second Order
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, Roe_Q_L, Roe_Q_R);
            } else {
                Roe_Q_L[0] = Q1[node_L];
                Roe_Q_L[1] = Q2[node_L];
                Roe_Q_L[2] = Q3[node_L];
                Roe_Q_L[3] = Q4[node_L];
                Roe_Q_L[4] = Q5[node_L];
                Roe_Q_R[0] = Q1[node_R];
                Roe_Q_R[1] = Q2[node_R];
                Roe_Q_R[2] = Q3[node_R];
                Roe_Q_R[3] = Q4[node_R];
                Roe_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            Roe_Q_L[0] = Q1[node_L];
            Roe_Q_L[1] = Q2[node_L];
            Roe_Q_L[2] = Q3[node_L];
            Roe_Q_L[3] = Q4[node_L];
            Roe_Q_L[4] = Q5[node_L];
            Roe_Q_R[0] = Q1[node_R];
            Roe_Q_R[1] = Q2[node_R];
            Roe_Q_R[2] = Q3[node_R];
            Roe_Q_R[3] = Q4[node_R];
            Roe_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // Compute flux_L
        Roe_flux_L[0] = rho_L*ubar_L;
        Roe_flux_L[1] =   u_L*Roe_flux_L[0] + p_L*nx;
        Roe_flux_L[2] =   v_L*Roe_flux_L[0] + p_L*ny;
        Roe_flux_L[3] =   w_L*Roe_flux_L[0] + p_L*nz;
        Roe_flux_L[4] =  ht_L*Roe_flux_L[0];
        
        // Compute flux_R
        Roe_flux_R[0] = rho_R*ubar_R;
        Roe_flux_R[1] =   u_R*Roe_flux_R[0] + p_R*nx;
        Roe_flux_R[2] =   v_R*Roe_flux_R[0] + p_R*ny;
        Roe_flux_R[3] =   w_R*Roe_flux_R[0] + p_R*nz;
        Roe_flux_R[4] =  ht_R*Roe_flux_R[0];
        
        AvgType = 1;
        if (AvgType == 1) {
            // ROE AVERAGE VARIABLES
            rho   = sqrt(rho_R * rho_L);
            sigma = rho/(rho_L + rho);
            u     = u_L  + sigma*(u_R  - u_L);
            v     = v_L  + sigma*(v_R  - v_L);
            w     = w_L  + sigma*(w_R  - w_L);
            ht    = ht_L + sigma*(ht_R - ht_L);
        } else {
            // SIMPLE AVERAGE VARIABLES
            rho   = 0.5*(rho_R + rho_L);
            u     = 0.5*(u_R  + u_L);
            v     = 0.5*(v_R  + v_L);
            w     = 0.5*(w_R  + w_L);
            ht    = 0.5*(ht_R + ht_L);
        }
        q2   = u*u + v*v + w*w;
        phi  = 0.5*(Gamma - 1.0)*q2;
        ubar = u*nx + v*ny + w*nz;
        c    = 0.0;
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                c = (Gamma - 1.0)*ht - phi;
                if (c < 0.0 && SolverIteration < 2000)
                    c = fabs(c);
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_BTW:
                c = ht/(Ref_Mach*Ref_Mach) - phi;
                if (c < 0.0 && SolverIteration < 2000)
                    c = fabs(c);
                c = sqrt(c);
                break;
            case NONDIMENSIONAL_LMROE:
                c = (Gamma - 1.0)*ht - phi;
                c = sqrt(c);
                break;
            default:
                error("Compute_RoeFlux_Precondition_ERIKSSON: Undefined Non-Dimensional Type - %d -1", NonDimensional_Type);
                break;
        }
        mach = sqrt(q2)/c;
        
        if (isnan(c))
            printf("Hell %10.5e\n", ubar);
        
        // Compute dQ
        // Conservative Variable Formulation
        if (Variable_Type == VARIABLE_CONSERVATIVE) {
            Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
            Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
            Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
            Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
            Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        }
        // Primitive Variable Formulation
        if (Variable_Type == VARIABLE_PRIMITIVE_PUT || Variable_Type == VARIABLE_PRIMITIVE_RUP) {
            Roe_dQ[0] = rho_R      - rho_L;
            Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
            Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
            Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
            Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
        }
        
        // Compute dw
        // Primitive Variable Formulation Pressure Velocity Temperature
        if (Variable_Type == VARIABLE_PRIMITIVE_PUT) {
            Roe_dw[0] = p_R - p_L;
            Roe_dw[1] = u_R - u_L;
            Roe_dw[2] = v_R - v_L;
            Roe_dw[3] = w_R - w_L;
            Roe_dw[4] = T_R - T_L;
        }
        // Primitive Variable Formulation Density Velocity Pressure
        if (Variable_Type == VARIABLE_PRIMITIVE_RUP) {
            Roe_dw[0] = rho_R - rho_L;
            Roe_dw[1] = u_R   - u_L;
            Roe_dw[2] = v_R   - v_L;
            Roe_dw[3] = w_R   - w_L;
            Roe_dw[4] = p_R   - p_L;
        }
        
        int nid;
        double lQ[NEQUATIONS];
        double tau, alpha, beta, beta_m, beta_p;
        double chi_m, chi_p, psi_m, psi_p, omega_m, omega_p;
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        
        //======================================================================
        // Compute Precondition Dissipation Flux
        //======================================================================
        // STEP - 1:
        // Compute the Local Max Mach and Max Change in Pressure around a node L and R
        // If node R is not ghost compute average
        dp_L     = dp_R    = 0.0;
        lmach_L  = lmach_R = 0.0;
        dp_max   = 0.0;
        mach_max = 0.0;
        // Check if Local Precondition is Required
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            // Check of Precondition Smooth is Required
            if (PrecondSmooth != 0) {
                for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L+1]; i++) {
                    nid   = crs_JA_Node2Node[i];
                    // Get the local Q's
                    lQ[0] = Q1[nid];
                    lQ[1] = Q2[nid];
                    lQ[2] = Q3[nid];
                    lQ[3] = Q4[nid];
                    lQ[4] = Q5[nid];
                    // Compute Local Equation of State
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    lmach_L = MAX(lmach_L, lmach);
                    dp_L    = MAX(dp_L, fabs(p_L - lp));
                }
                dp_max   = dp_L;
                mach_max = lmach_L;
                // Check if Right Node is not a Ghost Boundary Node
                if (node_R < nNode) {
                    for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R+1]; i++) {
                        nid  = crs_JA_Node2Node[i];
                        // Get the local Q's
                        lQ[0] = Q1[nid];
                        lQ[1] = Q2[nid];
                        lQ[2] = Q3[nid];
                        lQ[3] = Q4[nid];
                        lQ[4] = Q5[nid];
                        // Compute Local Equation of State
                        Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                        // Compute Smoothing Parameters
                        lmach_R = MAX(lmach_R, lmach);
                        dp_R    = MAX(dp_R, fabs(p_R - lp));
                    }
                    dp_max   = 0.5*(dp_L + dp_R);
                    mach_max = 0.5*(lmach_L + lmach_R);
                }
            } else {
                mach_max = MAX(mach_L, mach_R);
                dp_max   = fabs(p_L - p_R);
            }
        }

        
        // STEP 2:
        // Compute the required parameters
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            beta = MAX(mach, mach_max);
            if (beta < Ref_Mach)
                beta = sqrt(Ref_Mach*beta);
            beta = MIN(1.0, beta);
//            beta = MAX(beta, sqrt(1.0e-12/sqrt(q2)));
//            beta = MAX(beta, 2.0*dp_max/(rho*c*c));
//            beta = MAX(MAX(beta, 1.0e-5), 0.01*Ref_Mach);
//            beta = MIN(1.0, MIN(beta, Ref_Mach));
//            if (beta > 0.7) {
//                beta = 1.0;
//            } else {
//                beta = 2.0*beta*beta;
//                beta = beta/(1.0 - beta);
//                beta = sqrt(beta/(1 + (Gamma - 1.0)*beta));
//            }
            if (node_L < nNode && node_R < nNode) {
                PrecondSigma[node_L] += beta;
                PrecondSigma[node_R] += beta;
            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        beta    = beta*beta;
        beta_m  = 0.5*(1.0 - beta);
        beta_p  = 0.5*(1.0 + beta);
        tau     = 2.0*(c*c + 2.0*beta_m*phi)/(c*c*beta);
        alpha   = 2.0*sqrt(ubar*ubar*beta_m*beta_m + c*c*beta);
        chi_m   = ubar*beta_m - 0.5*alpha;
        chi_p   = ubar*beta_m + 0.5*alpha;
        psi_m   = (c*c*beta + 2.0*beta_m*ubar*chi_m)/(c*c*chi_m);
        psi_p   = (c*c*beta + 2.0*beta_m*ubar*chi_p)/(c*c*chi_p);
        omega_m = 2.0*beta_m*chi_m/(c*c*beta);
        omega_p = 2.0*beta_m*chi_p/(c*c*beta);
        
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar*beta_p + 0.5*alpha);
        Roe_Eigen[4][4] = fabs(ubar*beta_p - 0.5*alpha);
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Min and Max EigenValue Value
            MinEigenLamda1 = MIN(MinEigenLamda1, Roe_Eigen[0][0]);
            MaxEigenLamda1 = MAX(MaxEigenLamda1, Roe_Eigen[0][0]);
            MinEigenLamda4 = MIN(MinEigenLamda4, Roe_Eigen[3][3]);
            MaxEigenLamda4 = MAX(MaxEigenLamda4, Roe_Eigen[3][3]);
            MinEigenLamda5 = MIN(MinEigenLamda5, Roe_Eigen[4][4]);
            MaxEigenLamda5 = MAX(MaxEigenLamda5, Roe_Eigen[4][4]);

            // Min and Max Precondition Variable
            MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
            MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));

            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothType == RESIDUAL_SMOOTH_IMPLICIT) {
            maxlambda = MAX(Roe_Eigen[0][0], MAX(Roe_Eigen[3][3], Roe_Eigen[4][4]));
            for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L + 1]; i++) {
                // Get the Node ID and Match with node_R
                if (crs_JA_Node2Node[i] == node_R) {
                    RM_MaxEigenValue[i] = maxlambda*area;
                    break;
                }
            }
            RM_SumMaxEigenValue[node_L] += maxlambda*area;
            // Only Physical Nodes
            if (node_R < nNode) {
                for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R + 1]; i++) {
                    // Get the Node ID and Match with node_L
                    if (crs_JA_Node2Node[i] == node_L) {
                        RM_MaxEigenValue[i] = maxlambda*area;
                        break;
                    }
                }
                RM_SumMaxEigenValue[node_R] += maxlambda*area;
            }
        }
        
        // STEP 4:
        // Marco Modification (VariableType = RUP)
        if (PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON) {
            // Compute Tmod_g = R
            Roe_T[0][0] =  nx;
            Roe_T[0][1] =  ny;
            Roe_T[0][2] =  nz;
            Roe_T[0][3] = -rho*psi_m;
            Roe_T[0][4] =  rho*psi_p;

            Roe_T[1][0] =  0.0;
            Roe_T[1][1] = -nz;
            Roe_T[1][2] =  ny;
            Roe_T[1][3] =  nx;
            Roe_T[1][4] = -nx;

            Roe_T[2][0] =  nz;
            Roe_T[2][1] =  0.0;
            Roe_T[2][2] = -nx;
            Roe_T[2][3] =  ny;
            Roe_T[2][4] = -ny;

            Roe_T[3][0] = -ny;
            Roe_T[3][1] =  nx;
            Roe_T[3][2] =  0.0;
            Roe_T[3][3] =  nz;
            Roe_T[3][4] = -nz;

            Roe_T[4][0] =  0.0;
            Roe_T[4][1] =  0.0;
            Roe_T[4][2] =  0.0;
            Roe_T[4][3] = -rho*chi_m;
            Roe_T[4][4] =  rho*chi_p;
        }
        
        // Original (VariableType = RUP)
        if (PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL) {
            // Compute T_g = M1.Sinv.R
            Roe_T[0][0] =  nx;
            Roe_T[0][1] =  ny;
            Roe_T[0][2] =  nz;
            Roe_T[0][3] = -rho*(psi_m + omega_m);
            Roe_T[0][4] =  rho*(psi_p + omega_p);

            Roe_T[1][0] =  u*nx;
            Roe_T[1][1] =  u*ny - rho*nz;
            Roe_T[1][2] =  rho*ny + u*nz;
            Roe_T[1][3] =  rho*( nx - u*(psi_m + omega_m));
            Roe_T[1][4] =  rho*(-nx + u*(psi_p + omega_p));

            Roe_T[2][0] =  v*nx + rho*nz;
            Roe_T[2][1] =  v*ny;
            Roe_T[2][2] = -rho*nx + v*nz;
            Roe_T[2][3] =  rho*( ny - v*(psi_m + omega_m));
            Roe_T[2][4] =  rho*(-ny + v*(psi_p + omega_p));

            Roe_T[3][0] =  w*nx - rho*ny;
            Roe_T[3][1] =  rho*nx + w*ny;
            Roe_T[3][2] =  w*nz;
            Roe_T[3][3] =  rho*( nz - w*(psi_m + omega_m));
            Roe_T[3][4] =  rho*(-nz + w*(psi_p + omega_p));

            switch (NonDimensional_Type) {
                case NONDIMENSIONAL_GENERIC:
                    Roe_T[4][0] = -w*rho*ny + v*rho*nz + 0.5*nx*q2;
                    Roe_T[4][1] =  w*rho*nx - u*rho*nz + 0.5*ny*q2;
                    Roe_T[4][2] = -v*rho*nx + u*rho*ny + 0.5*nz*q2;
                    Roe_T[4][3] =  0.5*rho*( 2.0*ubar - (tau*chi_m)/(Gamma - 1.0) - q2*psi_m);
                    Roe_T[4][4] =  0.5*rho*(-2.0*ubar + (tau*chi_p)/(Gamma - 1.0) + q2*psi_p);
                    break;
                case NONDIMENSIONAL_BTW:
                    Roe_T[4][0] =  0.5*(Gamma - 1.0)*Ref_Mach*Ref_Mach*(-2.0*w*rho*ny + 2.0*v*rho*nz + nx*q2);
                    Roe_T[4][1] =  0.5*(Gamma - 1.0)*Ref_Mach*Ref_Mach*( 2.0*w*rho*nx - 2.0*u*rho*nz + ny*q2);
                    Roe_T[4][2] =  0.5*(Gamma - 1.0)*Ref_Mach*Ref_Mach*(-2.0*v*rho*nx + 2.0*u*rho*ny + nz*q2);
                    Roe_T[4][3] =  0.5*rho*Ref_Mach*Ref_Mach*((Gamma - 1.0)*(2.0*ubar - q2*psi_m) - tau*chi_m);
                    Roe_T[4][4] =  0.5*rho*Ref_Mach*Ref_Mach*(tau*chi_p - (Gamma - 1.0)*(2.0*ubar - q2*psi_p));
                    break;
                case NONDIMENSIONAL_LMROE:
                    Roe_T[4][0] = -w*rho*ny + v*rho*nz + 0.5*nx*q2;
                    Roe_T[4][1] =  w*rho*nx - u*rho*nz + 0.5*ny*q2;
                    Roe_T[4][2] = -v*rho*nx + u*rho*ny + 0.5*nz*q2;
                    Roe_T[4][3] =  0.5*rho*( 2.0*ubar - (tau*chi_m)/(Gamma - 1.0) - q2*psi_m);
                    Roe_T[4][4] =  0.5*rho*(-2.0*ubar + (tau*chi_p)/(Gamma - 1.0) + q2*psi_p);
                    break;
                default:
                    error("Compute_RoeFlux_Precondition_ERIKSSON: Undefined Non-Dimensional Type - %d -2", NonDimensional_Type);
                    break;
            }
        }
        
        // STEP 5:
        // Marco Modification or Original (VariableType = RUP)
        if (PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON || PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL) {
            // Compute T_d = Rinv
            Roe_Tinv[0][0] =  nx;
            Roe_Tinv[0][1] =  0.0;
            Roe_Tinv[0][2] =  nz;
            Roe_Tinv[0][3] = -ny;
            Roe_Tinv[0][4] = -nx/(c*c);

            Roe_Tinv[1][0] =  ny;
            Roe_Tinv[1][1] = -nz;
            Roe_Tinv[1][2] =  0.0;
            Roe_Tinv[1][3] =  nx;
            Roe_Tinv[1][4] = -ny/(c*c);

            Roe_Tinv[2][0] =  nz;
            Roe_Tinv[2][1] =  ny;
            Roe_Tinv[2][2] = -nx;
            Roe_Tinv[2][3] =  0.0;
            Roe_Tinv[2][4] = -nz/(c*c);

            Roe_Tinv[3][0] = 0.0;
            Roe_Tinv[3][1] = nx*chi_p/alpha;
            Roe_Tinv[3][2] = ny*chi_p/alpha;
            Roe_Tinv[3][3] = nz*chi_p/alpha;
            Roe_Tinv[3][4] = 1.0/(rho*alpha);

            Roe_Tinv[4][0] = 0.0;
            Roe_Tinv[4][1] = nx*chi_m/alpha;
            Roe_Tinv[4][2] = ny*chi_m/alpha;
            Roe_Tinv[4][3] = nz*chi_m/alpha;
            Roe_Tinv[4][4] = 1.0/(rho*alpha);
        }
        
        // STEP 6:
        // Compute the Dissipation Flux
        // Tmp => P = |Lambda|*Rinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);

        // Compute |A| = R*|Lambda|*Rinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);

        // Compute |A|*dw
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dw, Roe_fluxA);
        
        // STEP 7:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux_Precondition_ERIKSSON: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Averaged Variables
//------------------------------------------------------------------------------
void Compute_RoeVariables(double *Q_L, double *Q_R, double *Q_Roe) {
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht;
    double sigma;

    // Compute Equation of State
    Compute_EOS_Variables_ControlVolume(Q_L, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, et_L, ht_L);
    Compute_EOS_Variables_ControlVolume(Q_R, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, et_R, ht_R);
    
    // ROE AVERAGE VARIABLES
    rho   = sqrt(rho_R * rho_L);
    sigma = rho/(rho_L + rho);
    u     = u_L  + sigma*(u_R  - u_L);
    v     = v_L  + sigma*(v_R  - v_L);
    w     = w_L  + sigma*(w_R  - w_L);
    ht    = ht_L + sigma*(ht_R - ht_L);

    Q_Roe[0] = rho;
    Q_Roe[1] = rho*u;
    Q_Roe[2] = rho*v;
    Q_Roe[3] = rho*w;
    Q_Roe[4] = (rho/Gamma)*(ht + 0.5*(Gamma - 1.0)*(u*u + v*v + w*w));
    // TODO: Get_TotalEnergy(rho, pressure, u, v, w);
}

//------------------------------------------------------------------------------
//! Computes the Roe Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_Roe(int AddTime) {
    int i;
    int node_L, node_R;
    Vector3D areavec;
    double flux_roe[5];
    double flux_roe_diss[5];
    
    // Internal Edges
    for (i = 0; i < nEdge; i++) {
        // Get two nodes of edge
        node_L = intEdge[i].node[0];
        node_R = intEdge[i].node[1];

#ifdef DEBUG
        assert(node_R > node_L);
#endif
        
        // Get area vector
        areavec = intEdge[i].areav;
        
        // Compute the Roe Flux for this edge
        switch (PrecondMethod) {
            case SOLVER_PRECOND_NONE: // Roe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_LMFIX: // LMRoe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_WS: // Roe Weiss Smith Pre-Conditioner
                Compute_RoeFlux_WeissSmith(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_CV: // Roe Cecile Voizat Pre-Conditioner
                Compute_RoeFlux_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_CV_ORIGINAL: // Roe Cecile Voizat Pre-Conditioner
                Compute_RoeFlux_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
                Compute_RoeFlux_Precondition_BTW(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_BTW_ORIGINAL: // Roe Briley Taylor Whitfield Pre-Conditioner
                Compute_RoeFlux_Precondition_BTW(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner
                Compute_RoeFlux_Precondition_ERIKSSON(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL: // Roe Eriksson Pre-Conditioner
                Compute_RoeFlux_Precondition_ERIKSSON(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            default:
                error("Compute_Residual_Roe: Invalid Solver Precondition Scheme - %d -1", PrecondMethod);
                break;
        }
        
        // Compute for LHS
        Res1[node_L] += flux_roe[0];
        Res2[node_L] += flux_roe[1];
        Res3[node_L] += flux_roe[2];
        Res4[node_L] += flux_roe[3];
        Res5[node_L] += flux_roe[4];

        Res1[node_R] -= flux_roe[0];
        Res2[node_R] -= flux_roe[1];
        Res3[node_R] -= flux_roe[2];
        Res4[node_R] -= flux_roe[3];
        Res5[node_R] -= flux_roe[4];
        
        // Dissipation Term
        if (PrecondMethod == SOLVER_PRECOND_NONE || PrecondMethod == SOLVER_PRECOND_ROE_LMFIX
                || PrecondMethod == SOLVER_PRECOND_ROE_CV || PrecondMethod == SOLVER_PRECOND_ROE_CV_ORIGINAL 
                || PrecondMethod == SOLVER_PRECOND_ROE_BTW || PrecondMethod == SOLVER_PRECOND_ROE_BTW_ORIGINAL 
                || PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON || PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL) {
            // Compute for LHS
            Res1_Diss[node_L] += flux_roe_diss[0];
            Res2_Diss[node_L] += flux_roe_diss[1];
            Res3_Diss[node_L] += flux_roe_diss[2];
            Res4_Diss[node_L] += flux_roe_diss[3];
            Res5_Diss[node_L] += flux_roe_diss[4];

            Res1_Diss[node_R] -= flux_roe_diss[0];
            Res2_Diss[node_R] -= flux_roe_diss[1];
            Res3_Diss[node_R] -= flux_roe_diss[2];
            Res4_Diss[node_R] -= flux_roe_diss[3];
            Res5_Diss[node_R] -= flux_roe_diss[4];
        }
    }

    // Boundary Edges
    for (i = 0; i < nBEdge; i++) {
        // Get two nodes of edge
        node_L = bndEdge[i].node[0];
        node_R = bndEdge[i].node[1];

#ifdef DEBUG
        assert(node_R > node_L);
#endif
        
        // Get area vector
        areavec = bndEdge[i].areav;
        
        // Compute the Roe Flux for this edge
        switch (PrecondMethod) {
            case SOLVER_PRECOND_NONE: // Roe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_LMFIX: // LMRoe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_WS: // Roe Weiss Smith Pre-Conditioner
                Compute_RoeFlux_WeissSmith(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_CV: // Roe Cecile Voizat Pre-Conditioner
                Compute_RoeFlux_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_CV_ORIGINAL: // Roe Cecile Voizat Pre-Conditioner
                Compute_RoeFlux_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
                Compute_RoeFlux_Precondition_BTW(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_BTW_ORIGINAL: // Roe Briley Taylor Whitfield Pre-Conditioner
                Compute_RoeFlux_Precondition_BTW(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner
                Compute_RoeFlux_Precondition_ERIKSSON(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL: // Roe Eriksson Pre-Conditioner
                Compute_RoeFlux_Precondition_ERIKSSON(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            default:
                error("Compute_Residual_Roe: Invalid Solver Precondition Scheme - %d -2", PrecondMethod);
                break;
        }
        
        // Compute for LHS
        Res1[node_L] += flux_roe[0];
        Res2[node_L] += flux_roe[1];
        Res3[node_L] += flux_roe[2];
        Res4[node_L] += flux_roe[3];
        Res5[node_L] += flux_roe[4];
        
        // Dissipation Term
        if (PrecondMethod == SOLVER_PRECOND_NONE || PrecondMethod == SOLVER_PRECOND_ROE_LMFIX
                || PrecondMethod == SOLVER_PRECOND_ROE_CV || PrecondMethod == SOLVER_PRECOND_ROE_CV_ORIGINAL 
                || PrecondMethod == SOLVER_PRECOND_ROE_BTW || PrecondMethod == SOLVER_PRECOND_ROE_BTW_ORIGINAL
                || PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON || PrecondMethod == SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL) {
            // Compute for LHS
            Res1_Diss[node_L] += flux_roe_diss[0];
            Res2_Diss[node_L] += flux_roe_diss[1];
            Res3_Diss[node_L] += flux_roe_diss[2];
            Res4_Diss[node_L] += flux_roe_diss[3];
            Res5_Diss[node_L] += flux_roe_diss[4];
        }
    }
    
    // Precondition and Transform the Residual in Conservative or Primitive Form
    Compute_Precondition_Residual_Roe();
    if (PrecondMethod == SOLVER_PRECOND_NONE || PrecondMethod == SOLVER_PRECOND_ROE_LMFIX) {
        for (i = 0; i < nNode; i++) {
            Res1[i] += Res1_Diss[i];
            Res2[i] += Res2_Diss[i];
            Res3[i] += Res3_Diss[i];
            Res4[i] += Res4_Diss[i];
            Res5[i] += Res5_Diss[i];
        }
    }
}

