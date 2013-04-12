/*******************************************************************************
 * File:        HLLC_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 * Reference:   
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Material.h"
#include "Solver.h"
#include "Residual_Smoothing.h"

// Static Variable for Speed Up (if any)
static int HLLC_DB               = 0;
static double *HLLC_fluxA        = NULL;
static double *HLLC_flux_L       = NULL;
static double *HLLC_flux_R       = NULL;
static double *HLLC_Q_L          = NULL;
static double *HLLC_Q_R          = NULL;
static double *HLLC_dw           = NULL;
static double *HLLC_dQ           = NULL;
static double *HLLC_Eigen        = NULL;
static double **HLLC_A           = NULL;
static double **HLLC_M           = NULL;
static double **HLLC_Minv        = NULL;
static double **HLLC_P           = NULL;
static double **HLLC_Pinv        = NULL;
static double **HLLC_T           = NULL;
static double **HLLC_Tinv        = NULL;
static double **HLLC_K           = NULL;
static double **HLLC_Kinv        = NULL;

//------------------------------------------------------------------------------
//! Create HLLC Scheme Data Structure
//------------------------------------------------------------------------------
void HLLC_Init(void) {
    int i;
    double *tmp = NULL;
    
    // Check if Roe Data Structure is required
    if (HLLC_DB == 0) {
        HLLC_fluxA  = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_flux_L = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_flux_R = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_Q_L    = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_Q_R    = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_dw     = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_dQ     = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_Eigen  = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_A      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_A[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_M      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_M[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_Minv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_Minv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_P      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_P[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_Pinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_Pinv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_T      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_T[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_Tinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_Tinv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_K      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_K[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        HLLC_Kinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            HLLC_Kinv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        
        HLLC_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete HLLC Scheme Data Structure
//------------------------------------------------------------------------------
void HLLC_Finalize(void) {
     double *tmp = NULL;
    
    tmp = HLLC_A[0];
    free(tmp);
    tmp = HLLC_M[0];
    free(tmp);
    tmp = HLLC_Minv[0];
    free(tmp);
    tmp = HLLC_P[0];
    free(tmp);
    tmp = HLLC_Pinv[0];
    free(tmp);
    tmp = HLLC_T[0];
    free(tmp);
    tmp = HLLC_Tinv[0];
    free(tmp);
    tmp = HLLC_K[0];
    free(tmp);
    tmp = HLLC_Kinv[0];
    free(tmp);
    tmp = NULL;
    free(HLLC_fluxA);
    free(HLLC_flux_L);
    free(HLLC_flux_R);
    free(HLLC_Q_L);
    free(HLLC_Q_R);
    free(HLLC_dw);
    free(HLLC_dQ);
    free(HLLC_Eigen);
    free(HLLC_A);
    free(HLLC_M);
    free(HLLC_Minv);
    free(HLLC_P);
    free(HLLC_Pinv);
    free(HLLC_T);
    free(HLLC_Tinv);
    free(HLLC_K);
    free(HLLC_Kinv);
    HLLC_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset HLLC Scheme Data Structure
//------------------------------------------------------------------------------
void HLLC_Reset(void) {
    int i, j;
    
    if (HLLC_DB == 0)
        HLLC_Init();
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        HLLC_fluxA[i]  = 0.0;
        HLLC_flux_L[i] = 0.0;
        HLLC_flux_R[i] = 0.0;
        HLLC_Q_L[i]    = 0.0;
        HLLC_Q_R[i]    = 0.0;
        HLLC_dw[i]     = 0.0;
        HLLC_dQ[i]     = 0.0;
        HLLC_Eigen[i]  = 0.0;
        for (j = 0; j < NEQUATIONS; j++) {
            HLLC_A[i][j]     = 0.0;
            HLLC_M[i][j]     = 0.0;
            HLLC_Minv[i][j]  = 0.0;
            HLLC_P[i][j]     = 0.0;
            HLLC_Pinv[i][j]  = 0.0;
            HLLC_T[i][j]     = 0.0;
            HLLC_Tinv[i][j]  = 0.0;
            HLLC_K[i][j]     = 0.0;
            HLLC_Kinv[i][j]  = 0.0;
        }
    }
}

//------------------------------------------------------------------------------
//! Compute Transformed Residual Based on Preconditioned Variable Type and 
//  Desired variable type
//------------------------------------------------------------------------------
void Compute_Transformed_Residual_HLLC(void) {
    int inode;
    double Q[5];
    double res_hllc[5];
    double res_hllc_old[5];
    
    // Multiply by Precondition Matrix and Construct the Flux
    for (inode = 0; inode < nNode; inode++) {
        // Get the Variables
        Q[0] = Q1[inode];
        Q[1] = Q2[inode];
        Q[2] = Q3[inode];
        Q[3] = Q4[inode];
        Q[4] = Q5[inode];
        
        // Compute the Transformation Matrix
        // Mro
        Material_Get_Transformation_Matrix(Q, VariableType, VARIABLE_CON, HLLC_M);
        
        // Compute Transform Residual
        // Save the Residual
        res_hllc_old[0] = Res1[inode];
        res_hllc_old[1] = Res2[inode];
        res_hllc_old[2] = Res3[inode];
        res_hllc_old[3] = Res4[inode];
        res_hllc_old[4] = Res5[inode];
        
        // Finally Compute Transformed Residual
        res_hllc[0] = 0.0;
        res_hllc[1] = 0.0;
        res_hllc[2] = 0.0;
        res_hllc[3] = 0.0;
        res_hllc[4] = 0.0;
        MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, HLLC_M, res_hllc_old, res_hllc);
        Res1[inode] = res_hllc[0];
        Res2[inode] = res_hllc[1];
        Res3[inode] = res_hllc[2];
        Res4[inode] = res_hllc[3];
        Res5[inode] = res_hllc[4];
    }
}

//------------------------------------------------------------------------------
//! Compute Steady State HLLC Residual with Merkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_HLLC_Precondition_Merkel(void) {
    error("Compute_Steady_Residual_HLLC_Precondition_Merkel: Merkel Preconditioner Not Implemented - %d", PrecondMethod);
}

//------------------------------------------------------------------------------
//! Compute Steady State HLLC Residual with Turkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_HLLC_Precondition_Turkel(void) {
    int inode, j, nid;
    double Q[5], lQ[5];
    double res_hllc[5];
    double res_hllc_old[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Multiply by Precondition Matrix and Construct the Flux
    for (inode = 0; inode < nNode; inode++) {
        // Get the Variables
        Q[0] = Q1[inode];
        Q[1] = Q2[inode];
        Q[2] = Q3[inode];
        Q[3] = Q4[inode];
        Q[4] = Q5[inode];
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure can be Perturbation of them
        Material_Get_ControlVolume_Properties(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
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
                for (j = crs_IA_Node2Node[inode]; j < crs_IA_Node2Node[inode+1]; j++) {
                    nid   = crs_JA_Node2Node[j];
                    // Get the local Q's
                    lQ[0] = Q1[nid];
                    lQ[1] = Q2[nid];
                    lQ[2] = Q3[nid];
                    lQ[3] = Q4[nid];
                    lQ[4] = Q5[nid];
                    // Compute Local Equation of State
                    // Note: Based on VariableType Pressure can be Perturbation of them
                    Material_Get_ControlVolume_Properties(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    mach_max = MAX(mach_max, lmach);
                    dp_max   = MAX(dp_max, fabs(p - lp));
                }
            }
        }
        
        double delta, alpha, beta;
        
        // STEP 2:
        // Compute the required parameters
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            mach = MAX(mach, mach_max);
            if (mach < Ref_Mach)
                beta = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
            else
                beta = 0.5;
            beta = beta*(sqrt(mach*Ref_Mach)+ mach);
            beta = MIN(1.0, beta);
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, beta);
        MaxPrecondSigma = MAX(MaxPrecondSigma, beta);
        
        // Beta is function of M*M
        beta = beta*beta;
        
        // Set Preconditioner Parameters
        alpha = 0.0;
        delta = 0.0;
        switch (PrecondMethod) {
            case PRECOND_METHOD_BTW: // HLLC Briley Taylor Whitfield Pre-Conditioner
                alpha = 0.0;
                delta = 1.0;
                break;
            case PRECOND_METHOD_ERIKSSON: // HLLC Eriksson Pre-Conditioner
                alpha = 0.0;
                delta = beta;
                break;
            case PRECOND_METHOD_TURKEL: // HLLC Turkel Pre-Conditioner
                alpha = 0.4;
                delta = beta;
                break;
            default:
                error("Compute_Steady_Residual_HLLC_Precondition_Turkel: Invalid Solver Precondition Scheme - %d", PrecondMethod);
                break;
        }
        
        // STEP 3:
        // Compute the Precondition Matrix
        HLLC_K[0][0] = 1.0;
        HLLC_K[0][1] = 0.0;
        HLLC_K[0][2] = 0.0;
        HLLC_K[0][3] = 0.0;
        HLLC_K[0][4] = (delta - 1.0)/(c*c);

        HLLC_K[1][0] = 0.0;
        HLLC_K[1][1] = 1.0;
        HLLC_K[1][2] = 0.0;
        HLLC_K[1][3] = 0.0;
        HLLC_K[1][4] = -(alpha*u)/(rho*c*c);

        HLLC_K[2][0] = 0.0;
        HLLC_K[2][1] = 0.0;
        HLLC_K[2][2] = 1.0;
        HLLC_K[2][3] = 0.0;
        HLLC_K[2][4] = -(alpha*v)/(rho*c*c);

        HLLC_K[3][0] = 0.0;
        HLLC_K[3][1] = 0.0;
        HLLC_K[3][2] = 0.0;
        HLLC_K[3][3] = 1.0;
        HLLC_K[3][4] = -(alpha*w)/(rho*c*c);

        HLLC_K[4][0] = 0.0;
        HLLC_K[4][1] = 0.0;
        HLLC_K[4][2] = 0.0;
        HLLC_K[4][3] = 0.0;
        HLLC_K[4][4] = beta;
        
        // STEP 4:
        // Compute Transformation
        // Mrp
        Material_Get_Transformation_Matrix(Q, VariableType, PrecondVariableType, HLLC_M);
        // Mpr
        Material_Get_Transformation_Matrix(Q, PrecondVariableType, VariableType, HLLC_Minv);
        
        // STEP 5:
        // Compute Transformed Precondition Matrix: Mrp.P.Mpr
        MC_Matrix_Mul_Matrix(5, 5, HLLC_M, HLLC_K, HLLC_A);
        MC_Matrix_Mul_Matrix(5, 5, HLLC_A, HLLC_Minv, HLLC_T);
        
        // STEP 6:
        // Compute Preconditioned Residual
        // Save the Residual
        res_hllc_old[0] = Res1[inode];
        res_hllc_old[1] = Res2[inode];
        res_hllc_old[2] = Res3[inode];
        res_hllc_old[3] = Res4[inode];
        res_hllc_old[4] = Res5[inode];
        // Finally Compute precondition residual
        res_hllc[0] = 0.0;
        res_hllc[1] = 0.0;
        res_hllc[2] = 0.0;
        res_hllc[3] = 0.0;
        res_hllc[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, HLLC_T, res_hllc_old, res_hllc);
        Res1[inode] = res_hllc[0];
        Res2[inode] = res_hllc[1];
        Res3[inode] = res_hllc[2];
        Res4[inode] = res_hllc[3];
        Res5[inode] = res_hllc[4];
    }
}

//------------------------------------------------------------------------------
//! Compute HLLC Flux Original: Unmodified
//------------------------------------------------------------------------------
void Compute_Flux_HLLC_Original(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, h, ht, c, ubar;
    double sigma, area, nx, ny, nz, maxlambda;
    double omega_L, rhostar_L, rhoustar_L, rhovstar_L, rhowstar_L, rhoetstar_L;
    double omega_R, rhostar_R, rhoustar_R, rhovstar_R, rhowstar_R, rhoetstar_R;
    double pstar, cnst1, cnst2;
    double SL, SR, SM;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_HLLC[i] = 0.0;
    HLLC_Reset();
    
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
            if (SolverOrder == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, HLLC_Q_L, HLLC_Q_R);
            } else {
                HLLC_Q_L[0] = Q1[node_L];
                HLLC_Q_L[1] = Q2[node_L];
                HLLC_Q_L[2] = Q3[node_L];
                HLLC_Q_L[3] = Q4[node_L];
                HLLC_Q_L[4] = Q5[node_L];
                HLLC_Q_R[0] = Q1[node_R];
                HLLC_Q_R[1] = Q2[node_R];
                HLLC_Q_R[2] = Q3[node_R];
                HLLC_Q_R[3] = Q4[node_R];
                HLLC_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            HLLC_Q_L[0] = Q1[node_L];
            HLLC_Q_L[1] = Q2[node_L];
            HLLC_Q_L[2] = Q3[node_L];
            HLLC_Q_L[3] = Q4[node_L];
            HLLC_Q_L[4] = Q5[node_L];
            HLLC_Q_R[0] = Q1[node_R];
            HLLC_Q_R[1] = Q2[node_R];
            HLLC_Q_R[2] = Q3[node_R];
            HLLC_Q_R[3] = Q4[node_R];
            HLLC_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure can be Perturbation of them
        Material_Get_Face_Properties(HLLC_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Material_Get_Face_Properties(HLLC_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        h     = ht - 0.5*(u*u + v*v + w*w);
        c     = Material_Get_DH_SpeedSound(rho, h);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute HLLC Flux
        //======================================================================
        // STEP 1:
        // Compute the Roe Eigenvalues: Lambda
        // Note: Not the Absolute values
        HLLC_Eigen[0] = ubar;
        HLLC_Eigen[1] = HLLC_Eigen[0];
        HLLC_Eigen[2] = HLLC_Eigen[0];
        HLLC_Eigen[3] = ubar + c;
        HLLC_Eigen[4] = ubar - c;
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Min and Max EigenValue Value
            MinEigenLamda1 = MIN(MinEigenLamda1, fabs(HLLC_Eigen[0]));
            MaxEigenLamda1 = MAX(MaxEigenLamda1, fabs(HLLC_Eigen[0]));
            MinEigenLamda4 = MIN(MinEigenLamda4, fabs(HLLC_Eigen[3]));
            MaxEigenLamda4 = MAX(MaxEigenLamda4, fabs(HLLC_Eigen[3]));
            MinEigenLamda5 = MIN(MinEigenLamda5, fabs(HLLC_Eigen[4]));
            MaxEigenLamda5 = MAX(MaxEigenLamda5, fabs(HLLC_Eigen[4]));
            
            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(fabs(HLLC_Eigen[0]), MAX(fabs(HLLC_Eigen[3]), fabs(HLLC_Eigen[4])));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
            maxlambda = MAX(fabs(HLLC_Eigen[0]), MAX(fabs(HLLC_Eigen[3]), fabs(HLLC_Eigen[4])));
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
        
        // STEP 2:
        // Estimate the Wave Speeds SL, SR and SM
        SL = MIN(ubar_L - c_L, HLLC_Eigen[4]);
        SR = MAX(ubar_R + c_R, HLLC_Eigen[3]);
        SM = (p_R - p_L + rho_L*ubar_L*(SL - ubar_L) - rho_R*ubar_R*(SR - ubar_R))/
                                (rho_L*(SL - ubar_L) - rho_R*(SR - ubar_R));
        
        // STEP 3:
        // Compute Flux based on wave direction
        if (SL > 0.0) {
            Flux_HLLC[0] = rho_L*ubar_L;
            Flux_HLLC[1] =   u_L*Flux_HLLC[0] + p_L*nx;
            Flux_HLLC[2] =   v_L*Flux_HLLC[0] + p_L*ny;
            Flux_HLLC[3] =   w_L*Flux_HLLC[0] + p_L*nz;
            Flux_HLLC[4] =  ht_L*Flux_HLLC[0];
        } else if (SR < 0.0) {
            Flux_HLLC[0] = rho_R*ubar_R;
            Flux_HLLC[1] =   u_R*Flux_HLLC[0] + p_R*nx;
            Flux_HLLC[2] =   v_R*Flux_HLLC[0] + p_R*ny;
            Flux_HLLC[3] =   w_R*Flux_HLLC[0] + p_R*nz;
            Flux_HLLC[4] =  ht_R*Flux_HLLC[0];
        } else if ((SL <= 0.0) && (SM > 0.0)) {
            pstar   = p_L + rho_L*(SL - ubar_L)*(SM - ubar_L);
            omega_L = 1.0/(SL - SM);
            cnst1   = SL - ubar_L;
            cnst2   = pstar - p_L;
            
            rhostar_L   = omega_L*(cnst1*rho_L);
            rhoustar_L  = omega_L*(cnst1*rho_L*u_L + cnst2*nx);
            rhovstar_L  = omega_L*(cnst1*rho_L*v_L + cnst2*ny);
            rhowstar_L  = omega_L*(cnst1*rho_L*w_L + cnst2*nz);
            rhoetstar_L = omega_L*(cnst1*rho_L*et_L - (p_L + Gauge_Pressure)*ubar_L + (pstar + Gauge_Pressure)*SM);
            
            Flux_HLLC[0] = rhostar_L*SM;
            Flux_HLLC[1] = rhoustar_L*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_L*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_L*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_L + (pstar + Gauge_Pressure))*SM;        
        } else if ((SM <= 0.0) && (SR > 0.0)) {
            pstar   = p_R + rho_R*(SR - ubar_R)*(SM - ubar_R);
            omega_R = 1.0/(SR - SM);
            cnst1   = SR - ubar_R;
            cnst2   = pstar - p_R;
            
            rhostar_R   = omega_R*(cnst1*rho_R);
            rhoustar_R  = omega_R*(cnst1*rho_R*u_R + cnst2*nx);
            rhovstar_R  = omega_R*(cnst1*rho_R*v_R + cnst2*ny);
            rhowstar_R  = omega_R*(cnst1*rho_R*w_R + cnst2*nz);
            rhoetstar_R = omega_R*(cnst1*rho_R*et_R - (p_R + Gauge_Pressure)*ubar_R + (pstar + Gauge_Pressure)*SM);
            
            Flux_HLLC[0] = rhostar_R*SM;
            Flux_HLLC[1] = rhoustar_R*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_R*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_R*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_R + (pstar + Gauge_Pressure))*SM;
        } else
            error("Compute_Flux_HLLC_Original: Exception Anomaly");
        
        // Compute the HLLC Flux
        Flux_HLLC[0] *= area;
        Flux_HLLC[1] *= area;
        Flux_HLLC[2] *= area;
        Flux_HLLC[3] *= area;
        Flux_HLLC[4] *= area;
    } else
        error("Compute_Flux_HLLC_Original: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute HLLC Flux with Merkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_HLLC_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime) {
    error("Compute_Flux_HLLC_Precondition_Merkel: Merkel Preconditioner Not Implemented - %d", PrecondMethod);
}

//------------------------------------------------------------------------------
//! Compute HLLC Flux with Turkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_HLLC_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime) {
    int i, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, h, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    double omega_L, rhostar_L, rhoustar_L, rhovstar_L, rhowstar_L, rhoetstar_L;
    double omega_R, rhostar_R, rhoustar_R, rhovstar_R, rhowstar_R, rhoetstar_R;
    double pstar, cnst1, cnst2;
    double SL, SR, SM, EigenPre5_L, EigenPre4_R;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_HLLC[i] = 0.0;
    HLLC_Reset();
    
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
            if (SolverOrder == SOLVER_ORDER_SECOND) {
                Compute_SecondOrderReconstructQ(node_L, node_R, HLLC_Q_L, HLLC_Q_R);
            } else {
                HLLC_Q_L[0] = Q1[node_L];
                HLLC_Q_L[1] = Q2[node_L];
                HLLC_Q_L[2] = Q3[node_L];
                HLLC_Q_L[3] = Q4[node_L];
                HLLC_Q_L[4] = Q5[node_L];
                HLLC_Q_R[0] = Q1[node_R];
                HLLC_Q_R[1] = Q2[node_R];
                HLLC_Q_R[2] = Q3[node_R];
                HLLC_Q_R[3] = Q4[node_R];
                HLLC_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            HLLC_Q_L[0] = Q1[node_L];
            HLLC_Q_L[1] = Q2[node_L];
            HLLC_Q_L[2] = Q3[node_L];
            HLLC_Q_L[3] = Q4[node_L];
            HLLC_Q_L[4] = Q5[node_L];
            HLLC_Q_R[0] = Q1[node_R];
            HLLC_Q_R[1] = Q2[node_R];
            HLLC_Q_R[2] = Q3[node_R];
            HLLC_Q_R[3] = Q4[node_R];
            HLLC_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure can be Perturbation of them
        Material_Get_Face_Properties(HLLC_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Material_Get_Face_Properties(HLLC_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
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
        h    = ht - 0.5*(u*u + v*v + w*w);
        c    = Material_Get_DH_SpeedSound(rho, h);
        ubar = u*nx + v*ny + w*nz;
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        
        //======================================================================
        // Compute HLLC Preconditioned Flux
        //======================================================================
        // STEP 1:
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
                    // Note: Based on VariableType Pressure can be Perturbation of them
                    Material_Get_ControlVolume_Properties(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

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
                        // Note: Based on VariableType Pressure can be Perturbation of them
                        Material_Get_ControlVolume_Properties(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

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
        
        // Preconditioner Parameters
        double delta, alpha, beta;
        double zp;
        double cp, ubarp;
        
        // STEP 2:
        // Compute the required parameters
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            mach = MAX(mach, mach_max);
            if (mach < Ref_Mach)
                beta = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
            else
                beta = 0.5;
            beta = beta*(sqrt(mach*Ref_Mach)+ mach);
            beta = MIN(1.0, beta);
            
            if (node_L < nNode && node_R < nNode) {
                PrecondSigma[node_L] += beta;
                PrecondSigma[node_R] += beta;
            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        // Beta is function of M*M
        beta  = beta*beta;
        
        // Set the Preconditioner Parameters
        alpha = 0.0;
        delta = 0.0;
        switch (PrecondMethod) {
            case PRECOND_METHOD_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
                alpha = 0.0;
                delta = 1.0;
                break;
            case PRECOND_METHOD_ERIKSSON: // Roe Eriksson Pre-Conditioner
                alpha = 0.0;
                delta = beta;
                break;
            case PRECOND_METHOD_TURKEL: // Roe Turkel Pre-Conditioner
                alpha = 0.4;
                delta = beta;
                break;
            default:
                error("Compute_Flux_HLLC_Precondition_Turkel: Invalid Solver Precondition Scheme - %d", PrecondMethod);
                break;
        }
        
        // Compute the Modified Convective and Acoustic Wave speed
        zp    = 0.5*(1.0 + beta - alpha);
        ubarp = zp*ubar;
        cp    = sqrt(beta*(c*c - ubar*ubar) + zp*zp*ubar*ubar);
        
        // STEP 3:
        // Compute the Preconditioned Roe Eigenvalues: Lambda
        // Note: Not the Absolute values
        HLLC_Eigen[0] = ubar;
        HLLC_Eigen[1] = HLLC_Eigen[0];
        HLLC_Eigen[2] = HLLC_Eigen[0];
        HLLC_Eigen[3] = ubarp + cp;
        HLLC_Eigen[4] = ubarp - cp;
        
        // Compute the Preconditioned L and R Eigenvalues
        EigenPre5_L = zp*ubar_L - sqrt(beta*(c_L*c_L - ubar_L*ubar_L) + zp*zp*ubar_L*ubar_L);
        EigenPre4_R = zp*ubar_R + sqrt(beta*(c_R*c_R - ubar_R*ubar_R) + zp*zp*ubar_R*ubar_R);
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Min and Max EigenValue Value
            MinEigenLamda1 = MIN(MinEigenLamda1, fabs(HLLC_Eigen[0]));
            MaxEigenLamda1 = MAX(MaxEigenLamda1, fabs(HLLC_Eigen[0]));
            MinEigenLamda4 = MIN(MinEigenLamda4, fabs(HLLC_Eigen[3]));
            MaxEigenLamda4 = MAX(MaxEigenLamda4, fabs(HLLC_Eigen[3]));
            MinEigenLamda5 = MIN(MinEigenLamda5, fabs(HLLC_Eigen[4]));
            MaxEigenLamda5 = MAX(MaxEigenLamda5, fabs(HLLC_Eigen[4]));
            
            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(fabs(HLLC_Eigen[0]), MAX(fabs(HLLC_Eigen[3]), fabs(HLLC_Eigen[4])));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
            maxlambda = MAX(fabs(HLLC_Eigen[0]), MAX(fabs(HLLC_Eigen[3]), fabs(HLLC_Eigen[4])));
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
        // Estimate the Wave Speeds SL, SR and SM
        SL = MIN(EigenPre5_L, HLLC_Eigen[4]);
        SR = MAX(EigenPre4_R, HLLC_Eigen[3]);
        SM = (p_R - p_L + rho_L*ubar_L*(SL - ubar_L) - rho_R*ubar_R*(SR - ubar_R))/
                                (rho_L*(SL - ubar_L) - rho_R*(SR - ubar_R));
        
        // STEP 5:
        // Compute Preconditioned Flux based on wave direction
        if (SL >= 0.0) {
            Flux_HLLC[0] = rho_L*ubar_L;
            Flux_HLLC[1] =   u_L*Flux_HLLC[0] + p_L*nx;
            Flux_HLLC[2] =   v_L*Flux_HLLC[0] + p_L*ny;
            Flux_HLLC[3] =   w_L*Flux_HLLC[0] + p_L*nz;
            Flux_HLLC[4] =  ht_L*Flux_HLLC[0];
        } else if (SR < 0.0) {
            Flux_HLLC[0] = rho_R*ubar_R;
            Flux_HLLC[1] =   u_R*Flux_HLLC[0] + p_R*nx;
            Flux_HLLC[2] =   v_R*Flux_HLLC[0] + p_R*ny;
            Flux_HLLC[3] =   w_R*Flux_HLLC[0] + p_R*nz;
            Flux_HLLC[4] =  ht_R*Flux_HLLC[0];
        } else if ((SL <= 0.0) && (SM > 0.0)) {
            pstar   = p_L + rho_L*(SL - ubar_L)*(SM - ubar_L);
            omega_L = 1.0/(SL - SM);
            cnst1   = SL - ubar_L;
            cnst2   = pstar - p_L;
            
            rhostar_L   = omega_L*(cnst1*rho_L);
            rhoustar_L  = omega_L*(cnst1*rho_L*u_L + cnst2*nx);
            rhovstar_L  = omega_L*(cnst1*rho_L*v_L + cnst2*ny);
            rhowstar_L  = omega_L*(cnst1*rho_L*w_L + cnst2*nz);
            rhoetstar_L = omega_L*(cnst1*rho_L*et_L - (p_L + Gauge_Pressure)*ubar_L + (pstar + Gauge_Pressure)*SM);
            
            Flux_HLLC[0] = rhostar_L*SM;
            Flux_HLLC[1] = rhoustar_L*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_L*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_L*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_L + (pstar + Gauge_Pressure))*SM;        
        } else if ((SM <= 0.0) && (SR > 0.0)) {
            pstar   = p_R + rho_R*(SR - ubar_R)*(SM - ubar_R);
            omega_R = 1.0/(SR - SM);
            cnst1   = SR - ubar_R;
            cnst2   = pstar - p_R;
            
            rhostar_R   = omega_R*(cnst1*rho_R);
            rhoustar_R  = omega_R*(cnst1*rho_R*u_R + cnst2*nx);
            rhovstar_R  = omega_R*(cnst1*rho_R*v_R + cnst2*ny);
            rhowstar_R  = omega_R*(cnst1*rho_R*w_R + cnst2*nz);
            rhoetstar_R = omega_R*(cnst1*rho_R*et_R - (p_R + Gauge_Pressure)*ubar_R + (pstar + Gauge_Pressure)*SM);
            
            Flux_HLLC[0] = rhostar_R*SM;
            Flux_HLLC[1] = rhoustar_R*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_R*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_R*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_R + (pstar + Gauge_Pressure))*SM;
        } else
            error("Compute_Flux_HLLC_Precondition_Turkel: Exception Anomaly");
        
        // Compute the HLLC Flux
        Flux_HLLC[0] *= area;
        Flux_HLLC[1] *= area;
        Flux_HLLC[2] *= area;
        Flux_HLLC[3] *= area;
        Flux_HLLC[4] *= area;
    } else
        error("Compute_Flux_HLLC_Precondition_Turkel: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute HLLC Flux
//------------------------------------------------------------------------------
void Compute_Flux_HLLC(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime) {
    
    // Compute the Roe Flux for this edge
    switch (PrecondMethod) {
        case PRECOND_METHOD_NONE: // HLLC
            Compute_Flux_HLLC_Original(node_L, node_R, areavec, Flux_HLLC, AddTime);
            break;
        case PRECOND_METHOD_BTW: // HLLC Briley Taylor Whitfield Pre-Conditioner
            Compute_Flux_HLLC_Precondition_Turkel(node_L, node_R, areavec, Flux_HLLC, AddTime);
            break;
        case PRECOND_METHOD_ERIKSSON: // HLLC Eriksson Pre-Conditioner
            Compute_Flux_HLLC_Precondition_Turkel(node_L, node_R, areavec, Flux_HLLC, AddTime);
            break;
        case PRECOND_METHOD_MERKEL: // HLLC Merkel Pre-Conditioner
            Compute_Flux_HLLC_Precondition_Merkel(node_L, node_R, areavec, Flux_HLLC, AddTime);
            break;
        case PRECOND_METHOD_TURKEL: // HLLC Turkel Pre-Conditioner
            Compute_Flux_HLLC_Precondition_Turkel(node_L, node_R, areavec, Flux_HLLC, AddTime);
            break;
        default:
            error("Compute_Flux_HLLC: Invalid Solver Precondition Scheme - %d", PrecondMethod);
            break;
    }
}

//------------------------------------------------------------------------------
//! Computes the HLLC Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_HLLC(int AddTime) {
    int i;
    int node_L, node_R;
    Vector3D areavec;
    double flux_hllc[5];
    
    // Internal Edges
    for (i = 0; i < nEdge; i++) {
        // Get two nodes of edge
        node_L = intEdge[i].node[0];
        node_R = intEdge[i].node[1];
        
        // Get area vector
        areavec = intEdge[i].areav;
        
        // Compute the HLLC Flux for this edge
        Compute_Flux_HLLC(node_L, node_R, areavec, flux_hllc, AddTime);
        
        // L-Node
        Res1[node_L] += flux_hllc[0];
        Res2[node_L] += flux_hllc[1];
        Res3[node_L] += flux_hllc[2];
        Res4[node_L] += flux_hllc[3];
        Res5[node_L] += flux_hllc[4];
        
        // R-Node
        Res1[node_R] -= flux_hllc[0];
        Res2[node_R] -= flux_hllc[1];
        Res3[node_R] -= flux_hllc[2];
        Res4[node_R] -= flux_hllc[3];
        Res5[node_R] -= flux_hllc[4];
    }

    // Boundary Edges
    for (i = 0; i < nBEdge; i++) {
        // Get two nodes of edge
        node_L = bndEdge[i].node[0];
        node_R = bndEdge[i].node[1];
        
        // Get area vector
        areavec = bndEdge[i].areav;
        
        // Compute the HLLC Flux for this edge
        Compute_Flux_HLLC(node_L, node_R, areavec, flux_hllc, AddTime);
        
        // L-Node
        Res1[node_L] += flux_hllc[0];
        Res2[node_L] += flux_hllc[1];
        Res3[node_L] += flux_hllc[2];
        Res4[node_L] += flux_hllc[3];
        Res5[node_L] += flux_hllc[4];
    }
    
    // Transform the Residuals into Desired variable Form
    Compute_Transformed_Residual_HLLC();
    
    // Precondition the Residuals
    if ((SolverMethod == SOLVER_METHOD_STEADY) && (SolverScheme == SOLVER_SCHEME_EXPLICIT)) {
        switch (PrecondMethod) {
            case PRECOND_METHOD_NONE: // HLLC
                break;
            case PRECOND_METHOD_BTW: // HLLC Briley Taylor Whitfield Pre-Conditioner
                Compute_Steady_Residual_HLLC_Precondition_Turkel();
                break;
            case PRECOND_METHOD_ERIKSSON: // HLLC Eriksson Pre-Conditioner
                Compute_Steady_Residual_HLLC_Precondition_Turkel();
                break;
            case PRECOND_METHOD_MERKEL: // HLLC Merkel Pre-Conditioner
                Compute_Steady_Residual_HLLC_Precondition_Merkel();
                break;
            case PRECOND_METHOD_TURKEL: // HLLC Turkel Pre-Conditioner
                Compute_Steady_Residual_HLLC_Precondition_Turkel();
                break;
            default:
                error("Compute_Residual_HLLC: Invalid Solver Precondition Scheme - %d", PrecondMethod);
                break;
        }
    }
}

