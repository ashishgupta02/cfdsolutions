/*******************************************************************************
 * File:        Roe_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    4
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
static double **Roe_K           = NULL;
static double **Roe_Kinv        = NULL;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//------------------------------------------------------------------------------
//! Create Roe Scheme Data Structure
//------------------------------------------------------------------------------
void Roe_Init_New(void) {
    int i;
    double *tmp = NULL;
    
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
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_A[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_Eigen  = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_Eigen[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_M      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_M[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_Minv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_Minv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_P      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_P[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_Pinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_Pinv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_T      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_T[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_Tinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_Tinv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_K      = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_K[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        Roe_Kinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        tmp        = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            Roe_Kinv[i] = &tmp[i*NEQUATIONS];
        tmp        = NULL;
        
        Roe_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete Roe Scheme Data Structure
//------------------------------------------------------------------------------
void Roe_Finalize_New(void) {
    double *tmp = NULL;
    
    tmp = Roe_A[0];
    free(tmp);
    tmp = Roe_Eigen[0];
    free(tmp);
    tmp = Roe_M[0];
    free(tmp);
    tmp = Roe_Minv[0];
    free(tmp);
    tmp = Roe_P[0];
    free(tmp);
    tmp = Roe_Pinv[0];
    free(tmp);
    tmp = Roe_T[0];
    free(tmp);
    tmp = Roe_Tinv[0];
    free(tmp);
    tmp = Roe_K[0];
    free(tmp);
    tmp = Roe_Kinv[0];
    free(tmp);
    tmp = NULL;
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
    free(Roe_K);
    free(Roe_Kinv);
    Roe_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset Roe Scheme Data Structure
//------------------------------------------------------------------------------
void Roe_Reset_New(void) {
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
            Roe_K[i][j]     = 0.0;
            Roe_Kinv[i][j]  = 0.0;
        }
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//------------------------------------------------------------------------------
//! Compute Roe Precondition Transformed Matrix:
//  Note: For reqInv == 0: C = Mro
//                   else: Cinv = Mor
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix_None(int nodeID, int reqInv, double **PrecondMatrix) {
    double Q[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach;
    
    // Get the Variables
    Q[0] = Q1[nodeID];
    Q[1] = Q2[nodeID];
    Q[2] = Q3[nodeID];
    Q[3] = Q4[nodeID];
    Q[4] = Q5[nodeID];

    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
    
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformation
        // Mro
        Compute_Transformation_Matrix(VariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, PrecondMatrix);
    } else {
        // Compute Transformation
        // Mor
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VariableType, rho, u, v, w, c, PrecondMatrix);
    }
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix: Cecile Voizat
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix_CecileVoizat(int nodeID, int reqInv, double **PrecondMatrix) {
    int j, nid;
    double Q[5], lQ[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, sigma;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Initialization
    Roe_Reset_New();
    
    // Get the Variables
    Q[0] = Q1[nodeID];
    Q[1] = Q2[nodeID];
    Q[2] = Q3[nodeID];
    Q[3] = Q4[nodeID];
    Q[4] = Q5[nodeID];

    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);

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
            for (j = crs_IA_Node2Node[nodeID]; j < crs_IA_Node2Node[nodeID+1]; j++) {
                nid   = crs_JA_Node2Node[j];
                // Get the local Q's
                lQ[0] = Q1[nid];
                lQ[1] = Q2[nid];
                lQ[2] = Q3[nid];
                lQ[3] = Q4[nid];
                lQ[4] = Q5[nid];
                // Compute Local Equation of State
                // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        mach  = MAX(mach, mach_max);
        if (mach < Ref_Mach)
            sigma = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
        else
            sigma = 0.5;
        sigma = sigma*(sqrt(mach*Ref_Mach)+ mach);
        sigma = MIN(1.0, sigma);
    }
    if (PrecondType == PRECOND_TYPE_GLOBAL)
        sigma = MIN(1.0, Ref_Mach);
    
    
    // STEP 3:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute the Precondition Matrix
        Roe_K[0][0] = sigma*sigma;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
    } else {
        // Compute the Inverse Precondition Matrix
        Roe_Kinv[0][0] = 1.0/(sigma*sigma);
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
    }
    
    // STEP 4:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpo
        Compute_Transformation_Matrix(PrecondVariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, Roe_Minv);
    } else {
        // Compute Transformation
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
    }
    
    // STEP 5:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformed Precondition Matrix: Mrp.K.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
    } else {
        // Compute Transformed Inverse Precondition Matrix: Mrp.Kinv.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_Kinv, Roe_A);
    }
    MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, PrecondMatrix);
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix: Briley Taylor Whitfield
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix_BTW(int nodeID, int reqInv, double **PrecondMatrix) {
    int j, nid;
    double Q[5], lQ[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Initialization
    Roe_Reset_New();
    
    // Get the Variables
    Q[0] = Q1[nodeID];
    Q[1] = Q2[nodeID];
    Q[2] = Q3[nodeID];
    Q[3] = Q4[nodeID];
    Q[4] = Q5[nodeID];

    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);

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
            for (j = crs_IA_Node2Node[nodeID]; j < crs_IA_Node2Node[nodeID+1]; j++) {
                nid   = crs_JA_Node2Node[j];
                // Get the local Q's
                lQ[0] = Q1[nid];
                lQ[1] = Q2[nid];
                lQ[2] = Q3[nid];
                lQ[3] = Q4[nid];
                lQ[4] = Q5[nid];
                // Compute Local Equation of State
                // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
    }
    if (PrecondType == PRECOND_TYPE_GLOBAL)
        beta = MIN(1.0, Ref_Mach);

    beta = beta*beta;
    
    // STEP 3:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = beta;
    } else {
        // Compute the Inverse Precondition Matrix
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0/beta;
    }
    
    // STEP 4:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpo
        Compute_Transformation_Matrix(PrecondVariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, Roe_Minv);
    } else {
        // Compute Transformation
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
    }
    
    // STEP 5:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformed Precondition Matrix: Mrp.K.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
    } else {
        // Compute Transformed Inverse Precondition Matrix: Mrp.Kinv.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_Kinv, Roe_A);
    }
    MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, PrecondMatrix);
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix: Eriksson
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix_Eriksson(int nodeID, int reqInv, double **PrecondMatrix) {
    int j, nid;
    double Q[5], lQ[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Initialization
    Roe_Reset_New();
    
    // Get the Variables
    Q[0] = Q1[nodeID];
    Q[1] = Q2[nodeID];
    Q[2] = Q3[nodeID];
    Q[3] = Q4[nodeID];
    Q[4] = Q5[nodeID];

    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);

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
            for (j = crs_IA_Node2Node[nodeID]; j < crs_IA_Node2Node[nodeID+1]; j++) {
                nid   = crs_JA_Node2Node[j];
                // Get the local Q's
                lQ[0] = Q1[nid];
                lQ[1] = Q2[nid];
                lQ[2] = Q3[nid];
                lQ[3] = Q4[nid];
                lQ[4] = Q5[nid];
                // Compute Local Equation of State
                // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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

    beta = beta*beta;
    
    // STEP 3:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = (beta - 1.0)/(c*c);

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = beta;
    } else {
        // Compute the Inverse Precondition Matrix
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = (1.0 - beta)/(c*c*beta);

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0/beta;
    }
    
    // STEP 4:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpo
        Compute_Transformation_Matrix(PrecondVariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, Roe_Minv);
    } else {
        // Compute Transformation
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
    }
    
    // STEP 5:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformed Precondition Matrix: Mrp.K.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
    } else {
        // Compute Transformed Inverse Precondition Matrix: Mrp.Kinv.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_Kinv, Roe_A);
    }
    MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, PrecondMatrix);
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix: Merkel
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix_Merkel(int nodeID, int reqInv, double **PrecondMatrix) {
    int j, nid;
    double Q[5], lQ[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Initialization
    Roe_Reset_New();
    
    // Get the Variables
    Q[0] = Q1[nodeID];
    Q[1] = Q2[nodeID];
    Q[2] = Q3[nodeID];
    Q[3] = Q4[nodeID];
    Q[4] = Q5[nodeID];

    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);

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
            for (j = crs_IA_Node2Node[nodeID]; j < crs_IA_Node2Node[nodeID+1]; j++) {
                nid   = crs_JA_Node2Node[j];
                // Get the local Q's
                lQ[0] = Q1[nid];
                lQ[1] = Q2[nid];
                lQ[2] = Q3[nid];
                lQ[3] = Q4[nid];
                lQ[4] = Q5[nid];
                // Compute Local Equation of State
                // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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

    beta = beta*beta;
    
    // STEP 3:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
    } else {
        // Compute the Inverse Precondition Matrix
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
    }
    
    // STEP 4:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpo
        Compute_Transformation_Matrix(PrecondVariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, Roe_Minv);
    } else {
        // Compute Transformation
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
    }
    
    // STEP 5:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformed Precondition Matrix: Mrp.K.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
    } else {
        // Compute Transformed Inverse Precondition Matrix: Mrp.Kinv.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_Kinv, Roe_A);
    }
    MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, PrecondMatrix);
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix: Turkel
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix_Turkel(int nodeID, int reqInv, double **PrecondMatrix) {
    int j, nid;
    double Q[5], lQ[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Initialization
    Roe_Reset_New();
    
    // Get the Variables
    Q[0] = Q1[nodeID];
    Q[1] = Q2[nodeID];
    Q[2] = Q3[nodeID];
    Q[3] = Q4[nodeID];
    Q[4] = Q5[nodeID];

    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);

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
            for (j = crs_IA_Node2Node[nodeID]; j < crs_IA_Node2Node[nodeID+1]; j++) {
                nid   = crs_JA_Node2Node[j];
                // Get the local Q's
                lQ[0] = Q1[nid];
                lQ[1] = Q2[nid];
                lQ[2] = Q3[nid];
                lQ[3] = Q4[nid];
                lQ[4] = Q5[nid];
                // Compute Local Equation of State
                // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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

    beta = beta*beta;
    
    // STEP 3:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
    } else {
        // Compute the Inverse Precondition Matrix
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
    }
    
    // STEP 4:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpo
        Compute_Transformation_Matrix(PrecondVariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, Roe_Minv);
    } else {
        // Compute Transformation
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
    }
    
    // STEP 5:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformed Precondition Matrix: Mrp.K.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
    } else {
        // Compute Transformed Inverse Precondition Matrix: Mrp.Kinv.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_Kinv, Roe_A);
    }
    MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, PrecondMatrix);
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix: Weiss Smith
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix_WeissSmith(int nodeID, int reqInv, double **PrecondMatrix) {
    int j, nid;
    double Q[5], lQ[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;
    
    // Initialization
    Roe_Reset_New();
    
    // Get the Variables
    Q[0] = Q1[nodeID];
    Q[1] = Q2[nodeID];
    Q[2] = Q3[nodeID];
    Q[3] = Q4[nodeID];
    Q[4] = Q5[nodeID];

    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);

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
            for (j = crs_IA_Node2Node[nodeID]; j < crs_IA_Node2Node[nodeID+1]; j++) {
                nid   = crs_JA_Node2Node[j];
                // Get the local Q's
                lQ[0] = Q1[nid];
                lQ[1] = Q2[nid];
                lQ[2] = Q3[nid];
                lQ[3] = Q4[nid];
                lQ[4] = Q5[nid];
                // Compute Local Equation of State
                // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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

    beta = beta*beta;
    
    // STEP 3:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
    } else {
        // Compute the Inverse Precondition Matrix
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
    }
    
    // STEP 4:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpo
        Compute_Transformation_Matrix(PrecondVariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, Roe_Minv);
    } else {
        // Compute Transformation
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
    }
    
    // STEP 5:
    // Check if Inverse Preconditioner is Requested
    if (reqInv == 0) {
        // Compute Transformed Precondition Matrix: Mrp.K.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
    } else {
        // Compute Transformed Inverse Precondition Matrix: Mrp.Kinv.Mpo
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_Kinv, Roe_A);
    }
    MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, PrecondMatrix);
}

//------------------------------------------------------------------------------
//! Compute Roe Transformed Preconditioner Matrix
//  Note: For reqInv == 0: C = Mrp.K.Mpo
//                   else: Cinv = Mop.Kinv.Mpr
//------------------------------------------------------------------------------
void Compute_Roe_Transformed_Precondition_Matrix(int nodeID, int reqInv, double **PrecondMatrix) {
    // Precondition the Residuals
    switch (PrecondMethod) {
        case PRECOND_METHOD_NONE: // Roe
            Compute_Roe_Transformed_Precondition_Matrix_None(nodeID, reqInv, PrecondMatrix);
            break;
        case PRECOND_METHOD_ROE_LMFIX: // LMRoe
            Compute_Roe_Transformed_Precondition_Matrix_None(nodeID, reqInv, PrecondMatrix);
            break;
        case PRECOND_METHOD_ROE_WS: // Roe Weiss Smith Pre-Conditioner
            Compute_Roe_Transformed_Precondition_Matrix_WeissSmith(nodeID, reqInv, PrecondMatrix);
            break;
        case PRECOND_METHOD_ROE_CV: // Roe Cecile Voizat Pre-Conditioner
            Compute_Roe_Transformed_Precondition_Matrix_CecileVoizat(nodeID, reqInv, PrecondMatrix);
            break;
        case PRECOND_METHOD_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
            Compute_Roe_Transformed_Precondition_Matrix_BTW(nodeID, reqInv, PrecondMatrix);
            break;
        case PRECOND_METHOD_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner
            Compute_Roe_Transformed_Precondition_Matrix_Eriksson(nodeID, reqInv, PrecondMatrix);
            break;
        case PRECOND_METHOD_ROE_MERKEL: // Roe Merkel Pre-Conditioner
            Compute_Roe_Transformed_Precondition_Matrix_Merkel(nodeID, reqInv, PrecondMatrix);
            break;
        case PRECOND_METHOD_ROE_TURKEL: // Roe Turkel Pre-Conditioner
            Compute_Roe_Transformed_Precondition_Matrix_Turkel(nodeID, reqInv, PrecondMatrix);
            break;
        default:
            error("Compute_Roe_Transformed_Precondition_Matrix: Invalid Solver Precondition Scheme - %d", PrecondMethod);
            break;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//------------------------------------------------------------------------------
//! Compute Transformed Residual Based on Preconditioned Variable Type and 
//  Desired variable type
//------------------------------------------------------------------------------
void Compute_Transformed_Residual_Roe(void) {
    int inode;
    double Q[5];
    double res_roe[5];
    double res_roe_conv[5];
    double res_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach;
    
    // Multiply by Precondition Matrix and Construct the Flux
    for (inode = 0; inode < nNode; inode++) {
        // Get the Variables
        Q[0] = Q1[inode];
        Q[1] = Q2[inode];
        Q[2] = Q3[inode];
        Q[3] = Q4[inode];
        Q[4] = Q5[inode];
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
        // Compute Transformation
        // Mro
        Compute_Transformation_Matrix(VariableType, VARIABLE_CONSERVATIVE, rho, u, v, w, c, Roe_M);
        
        // Compute Transform Convective Residual
        // Get the Convective Residual
        res_roe_conv[0] = Res1[inode];
        res_roe_conv[1] = Res2[inode];
        res_roe_conv[2] = Res3[inode];
        res_roe_conv[3] = Res4[inode];
        res_roe_conv[4] = Res5[inode];
        
        // Finally Compute Transformed Convective Residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_M, res_roe_conv, res_roe);
        Res1[inode] = res_roe[0];
        Res2[inode] = res_roe[1];
        Res3[inode] = res_roe[2];
        Res4[inode] = res_roe[3];
        Res5[inode] = res_roe[4];
        
        // Compute Transform Dissipative Residual
        // Get the Dissipative Residual
        res_roe_diss[0] = Res1_Diss[inode];
        res_roe_diss[1] = Res2_Diss[inode];
        res_roe_diss[2] = Res3_Diss[inode];
        res_roe_diss[3] = Res4_Diss[inode];
        res_roe_diss[4] = Res5_Diss[inode];
        
        // Finally Compute Transformed Dissipative Residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_M, res_roe_diss, res_roe);
        Res1_Diss[inode] = res_roe[0];
        Res2_Diss[inode] = res_roe[1];
        Res3_Diss[inode] = res_roe[2];
        Res4_Diss[inode] = res_roe[3];
        Res5_Diss[inode] = res_roe[4];
    }
}

//------------------------------------------------------------------------------
//! Compute Steady State Roe Residual with Cecile Voizat Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_Roe_Precondition_CecileVoizat(void) {
    int inode, j, nid;
    double Q[5], lQ[5];
    double res_roe[5];
    double res_roe_conv[5];
    double res_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, sigma;
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
            mach  = MAX(mach, mach_max);
            if (mach < Ref_Mach)
                sigma = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
            else
                sigma = 0.5;
            sigma = sigma*(sqrt(mach*Ref_Mach)+ mach);
            sigma = MIN(1.0, sigma);
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            sigma = MIN(1.0, Ref_Mach);
        
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sigma);
        MaxPrecondSigma = MAX(MaxPrecondSigma, sigma);

        // STEP 3:
        // Compute the Precondition Matrix
        Roe_K[0][0] = sigma*sigma;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
        
        // STEP 4:
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // STEP 5:
        // Compute Transformed Precondition Matrix: Mrp.P.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, Roe_T);
        
        // STEP 6:
        // Compute Preconditioned Convective Residual
        // Get the Convective Residual
        res_roe_conv[0] = Res1[inode];
        res_roe_conv[1] = Res2[inode];
        res_roe_conv[2] = Res3[inode];
        res_roe_conv[3] = Res4[inode];
        res_roe_conv[4] = Res5[inode];
        // Finally Compute precondition convective residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_conv, res_roe);
        Res1[inode] = res_roe[0];
        Res2[inode] = res_roe[1];
        Res3[inode] = res_roe[2];
        Res4[inode] = res_roe[3];
        Res5[inode] = res_roe[4];
        
        // STEP 7:
        // Compute the Preconditioned Dissipative Residual
        // Get the Dissipative Residual
        res_roe_diss[0] = Res1_Diss[inode];
        res_roe_diss[1] = Res2_Diss[inode];
        res_roe_diss[2] = Res3_Diss[inode];
        res_roe_diss[3] = Res4_Diss[inode];
        res_roe_diss[4] = Res5_Diss[inode];
        // Finally Compute precondition dissipative residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_diss, res_roe);
        Res1_Diss[inode] = res_roe[0];
        Res2_Diss[inode] = res_roe[1];
        Res3_Diss[inode] = res_roe[2];
        Res4_Diss[inode] = res_roe[3];
        Res5_Diss[inode] = res_roe[4];
    }
}

//------------------------------------------------------------------------------
//! Compute Steady State Roe Residual with Briley Taylor Whitfield Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_Roe_Precondition_BTW(void) {
    int inode, j, nid;
    double Q[5], lQ[5];
    double res_roe[5];
    double res_roe_conv[5];
    double res_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        beta = beta*beta;
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
        MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));
        
        // STEP 3:
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = beta;
        
        // STEP 4:
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // STEP 5:
        // Compute Transformed Precondition Matrix: Mrp.P.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, Roe_T);
        
        // STEP 6:
        // Compute Preconditioned Convective Residual
        // Get the Convective Residual
        res_roe_conv[0] = Res1[inode];
        res_roe_conv[1] = Res2[inode];
        res_roe_conv[2] = Res3[inode];
        res_roe_conv[3] = Res4[inode];
        res_roe_conv[4] = Res5[inode];
        // Finally Compute precondition convective residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_conv, res_roe);
        Res1[inode] = res_roe[0];
        Res2[inode] = res_roe[1];
        Res3[inode] = res_roe[2];
        Res4[inode] = res_roe[3];
        Res5[inode] = res_roe[4];
        
        // STEP 7:
        // Compute the Preconditioned Dissipative Residual
        // Get the Dissipative Residual
        res_roe_diss[0] = Res1_Diss[inode];
        res_roe_diss[1] = Res2_Diss[inode];
        res_roe_diss[2] = Res3_Diss[inode];
        res_roe_diss[3] = Res4_Diss[inode];
        res_roe_diss[4] = Res5_Diss[inode];
        // Finally Compute precondition dissipative residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_diss, res_roe);
        Res1_Diss[inode] = res_roe[0];
        Res2_Diss[inode] = res_roe[1];
        Res3_Diss[inode] = res_roe[2];
        Res4_Diss[inode] = res_roe[3];
        Res5_Diss[inode] = res_roe[4];
    }
}

//------------------------------------------------------------------------------
//! Compute Steady State Roe Residual with Eriksson Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_Roe_Precondition_Eriksson(void) {
    int inode, j, nid;
    double Q[5], lQ[5];
    double res_roe[5];
    double res_roe_conv[5];
    double res_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta = beta*beta;
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
        MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));
        
        // STEP 3:
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = (beta - 1.0)/(c*c);

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = beta;
        
        // STEP 4:
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // STEP 5:
        // Compute Transformed Precondition Matrix: Mrp.P.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, Roe_T);
        
        // STEP 6:
        // Compute Preconditioned Convective Residual
        // Get the Convective Residual
        res_roe_conv[0] = Res1[inode];
        res_roe_conv[1] = Res2[inode];
        res_roe_conv[2] = Res3[inode];
        res_roe_conv[3] = Res4[inode];
        res_roe_conv[4] = Res5[inode];
        // Finally Compute precondition convective residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_conv, res_roe);
        Res1[inode] = res_roe[0];
        Res2[inode] = res_roe[1];
        Res3[inode] = res_roe[2];
        Res4[inode] = res_roe[3];
        Res5[inode] = res_roe[4];
        
        // STEP 7:
        // Compute the Preconditioned Dissipative Residual
        // Get the Dissipative Residual
        res_roe_diss[0] = Res1_Diss[inode];
        res_roe_diss[1] = Res2_Diss[inode];
        res_roe_diss[2] = Res3_Diss[inode];
        res_roe_diss[3] = Res4_Diss[inode];
        res_roe_diss[4] = Res5_Diss[inode];
        // Finally Compute precondition dissipative residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_diss, res_roe);
        Res1_Diss[inode] = res_roe[0];
        Res2_Diss[inode] = res_roe[1];
        Res3_Diss[inode] = res_roe[2];
        Res4_Diss[inode] = res_roe[3];
        Res5_Diss[inode] = res_roe[4];
    }
}

//------------------------------------------------------------------------------
//! Compute Steady State Roe Residual with Merkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_Roe_Precondition_Merkel(void) {
    int inode, j, nid;
    double Q[5], lQ[5];
    double res_roe[5];
    double res_roe_conv[5];
    double res_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta = beta*beta;
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
        MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));
        
        // STEP 3:
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
        
        // STEP 4:
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // STEP 5:
        // Compute Transformed Precondition Matrix: Mrp.P.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, Roe_T);
        
        // STEP 6:
        // Compute Preconditioned Convective Residual
        // Get the Convective Residual
        res_roe_conv[0] = Res1[inode];
        res_roe_conv[1] = Res2[inode];
        res_roe_conv[2] = Res3[inode];
        res_roe_conv[3] = Res4[inode];
        res_roe_conv[4] = Res5[inode];
        // Finally Compute precondition convective residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_conv, res_roe);
        Res1[inode] = res_roe[0];
        Res2[inode] = res_roe[1];
        Res3[inode] = res_roe[2];
        Res4[inode] = res_roe[3];
        Res5[inode] = res_roe[4];
        
        // STEP 7:
        // Compute the Preconditioned Dissipative Residual
        // Get the Dissipative Residual
        res_roe_diss[0] = Res1_Diss[inode];
        res_roe_diss[1] = Res2_Diss[inode];
        res_roe_diss[2] = Res3_Diss[inode];
        res_roe_diss[3] = Res4_Diss[inode];
        res_roe_diss[4] = Res5_Diss[inode];
        // Finally Compute precondition dissipative residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_diss, res_roe);
        Res1_Diss[inode] = res_roe[0];
        Res2_Diss[inode] = res_roe[1];
        Res3_Diss[inode] = res_roe[2];
        Res4_Diss[inode] = res_roe[3];
        Res5_Diss[inode] = res_roe[4];
    }
}

//------------------------------------------------------------------------------
//! Compute Steady State Roe Residual with Turkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_Roe_Precondition_Turkel(void) {
    int inode, j, nid;
    double Q[5], lQ[5];
    double res_roe[5];
    double res_roe_conv[5];
    double res_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta = beta*beta;
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
        MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));
        
        // STEP 3:
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
        
        // STEP 4:
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // STEP 5:
        // Compute Transformed Precondition Matrix: Mrp.P.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, Roe_T);
        
        // STEP 6:
        // Compute Preconditioned Convective Residual
        // Get the Convective Residual
        res_roe_conv[0] = Res1[inode];
        res_roe_conv[1] = Res2[inode];
        res_roe_conv[2] = Res3[inode];
        res_roe_conv[3] = Res4[inode];
        res_roe_conv[4] = Res5[inode];
        // Finally Compute precondition convective residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_conv, res_roe);
        Res1[inode] = res_roe[0];
        Res2[inode] = res_roe[1];
        Res3[inode] = res_roe[2];
        Res4[inode] = res_roe[3];
        Res5[inode] = res_roe[4];
        
        // STEP 7:
        // Compute the Preconditioned Dissipative Residual
        // Get the Dissipative Residual
        res_roe_diss[0] = Res1_Diss[inode];
        res_roe_diss[1] = Res2_Diss[inode];
        res_roe_diss[2] = Res3_Diss[inode];
        res_roe_diss[3] = Res4_Diss[inode];
        res_roe_diss[4] = Res5_Diss[inode];
        // Finally Compute precondition dissipative residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_diss, res_roe);
        Res1_Diss[inode] = res_roe[0];
        Res2_Diss[inode] = res_roe[1];
        Res3_Diss[inode] = res_roe[2];
        Res4_Diss[inode] = res_roe[3];
        Res5_Diss[inode] = res_roe[4];
    }
}

//------------------------------------------------------------------------------
//! Compute Steady State Roe Residual with Weiss Smith Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Steady_Residual_Roe_Precondition_WeissSmith(void) {
    int inode, j, nid;
    double Q[5], lQ[5];
    double res_roe[5];
    double res_roe_conv[5];
    double res_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach, beta;
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta = beta*beta;
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sqrt(beta));
        MaxPrecondSigma = MAX(MaxPrecondSigma, sqrt(beta));
        
        // STEP 3:
        // Compute the Precondition Matrix
        Roe_K[0][0] = 1.0;
        Roe_K[0][1] = 0.0;
        Roe_K[0][2] = 0.0;
        Roe_K[0][3] = 0.0;
        Roe_K[0][4] = 0.0;

        Roe_K[1][0] = 0.0;
        Roe_K[1][1] = 1.0;
        Roe_K[1][2] = 0.0;
        Roe_K[1][3] = 0.0;
        Roe_K[1][4] = 0.0;

        Roe_K[2][0] = 0.0;
        Roe_K[2][1] = 0.0;
        Roe_K[2][2] = 1.0;
        Roe_K[2][3] = 0.0;
        Roe_K[2][4] = 0.0;

        Roe_K[3][0] = 0.0;
        Roe_K[3][1] = 0.0;
        Roe_K[3][2] = 0.0;
        Roe_K[3][3] = 1.0;
        Roe_K[3][4] = 0.0;

        Roe_K[4][0] = 0.0;
        Roe_K[4][1] = 0.0;
        Roe_K[4][2] = 0.0;
        Roe_K[4][3] = 0.0;
        Roe_K[4][4] = 1.0;
        
        // STEP 4:
        // Compute Transformation
        // Mrp
        Compute_Transformation_Matrix(VariableType, PrecondVariableType, rho, u, v, w, c, Roe_M);
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // STEP 5:
        // Compute Transformed Precondition Matrix: Mrp.P.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_K, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_A, Roe_Minv, Roe_T);
        
        // STEP 6:
        // Compute Preconditioned Convective Residual
        // Get the Convective Residual
        res_roe_conv[0] = Res1[inode];
        res_roe_conv[1] = Res2[inode];
        res_roe_conv[2] = Res3[inode];
        res_roe_conv[3] = Res4[inode];
        res_roe_conv[4] = Res5[inode];
        // Finally Compute precondition convective residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_conv, res_roe);
        Res1[inode] = res_roe[0];
        Res2[inode] = res_roe[1];
        Res3[inode] = res_roe[2];
        Res4[inode] = res_roe[3];
        Res5[inode] = res_roe[4];
        
        // STEP 7:
        // Compute the Preconditioned Dissipative Residual
        // Get the Dissipative Residual
        res_roe_diss[0] = Res1_Diss[inode];
        res_roe_diss[1] = Res2_Diss[inode];
        res_roe_diss[2] = Res3_Diss[inode];
        res_roe_diss[3] = Res4_Diss[inode];
        res_roe_diss[4] = Res5_Diss[inode];
        // Finally Compute precondition dissipative residual
        res_roe[0] = 0.0;
        res_roe[1] = 0.0;
        res_roe[2] = 0.0;
        res_roe[3] = 0.0;
        res_roe[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, Roe_T, res_roe_diss, res_roe);
        Res1_Diss[inode] = res_roe[0];
        Res2_Diss[inode] = res_roe[1];
        Res3_Diss[inode] = res_roe[2];
        Res4_Diss[inode] = res_roe[3];
        Res5_Diss[inode] = res_roe[4];
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//------------------------------------------------------------------------------
//! Compute Roe Flux Original: One Sided
// Note: Will only work with conservation variables and generic non-dimensionalization 
//------------------------------------------------------------------------------
void Compute_Flux_Roe_OneSided(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, cinv, gm1oc;
    double l1, l2, l3, l4, l5;
    double r1, r2, r3, r4, r5;
    double sigma, area, nx, ny, nz, maxlambda, alpha;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        Flux_Roe_Conv[i] = 0.0;
        Flux_Roe_Diss[i] = 0.0;
    }
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        q2    = u*u + v*v + w*w;
        cinv  = 1.0/c;
        gm1oc = (Gamma - 1.0)*cinv;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
                Roe_dQ[0] = rho_R      - rho_L;
                Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
                Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
                Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
                Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
                
                // Get the signed Eigenvalues
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
                Roe_dQ[0] = rho_R      - rho_L;
                Roe_dQ[1] = rho_R*u_R  - rho_L*u_L;
                Roe_dQ[2] = rho_R*v_R  - rho_L*v_L;
                Roe_dQ[3] = rho_R*w_R  - rho_L*w_L;
                Roe_dQ[4] = rho_R*et_R - rho_L*et_L;
                
                // Get the signed Eigenvalues
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
        error("Compute_Flux_Roe_OneSided: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux: Low Dissipation
//------------------------------------------------------------------------------
void Compute_Flux_Roe_LD(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, maxcount;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, ubar1, ubar2;
    double sigma, area, nx, ny, nz, maxlambda;
    double Mach, nmax, fix;
    Vector3D Vn1, Vn2;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        
        // STEP 2:
        // Compute Low Mach Number Fix
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
            error("Compute_Flux_Roe_LD: Unable to compute point on the plane (nx, ny, nz): (%lf, %lf, %lf)", nx, ny, nz);
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
        
        // STEP 3:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VARIABLE_PRIMITIVE_RUP, rho, u, v, w, c, Roe_M);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = rho_R - rho_L;
        Roe_dQ[1] = u_R - u_L;
        Roe_dQ[2] = v_R - v_L;
        Roe_dQ[3] = w_R - w_L;
        Roe_dQ[4] = p_R - p_L;
        
        double E0, E1, E2, dU;
        E0 = Roe_Eigen[0][0];
        E1 = Roe_Eigen[3][3];
        E2 = Roe_Eigen[4][4];
        dU = Roe_dQ[1]*nx + Roe_dQ[2]*ny + Roe_dQ[3]*nz;
        Roe_dw[0] = E0*Roe_dQ[0] + (fix*(E1 + E2) - 2.0*E0)*Roe_dQ[4]/(2.0*c*c) + rho*(E1 - E2)*dU/(2.0*c);
        Roe_dw[1] = E0*Roe_dQ[1] + 0.5*nx*dU*(fix*(E1 + E2) - 2.0*E0) + nx*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[2] = E0*Roe_dQ[2] + 0.5*ny*dU*(fix*(E1 + E2) - 2.0*E0) + ny*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[3] = E0*Roe_dQ[3] + 0.5*nz*dU*(fix*(E1 + E2) - 2.0*E0) + nz*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[4] = 0.5*Roe_dQ[4]*(E1 + E2) + 0.5*rho*c*dU*(E1 - E2);
        MC_Matrix_Mul_Vector(5, 5, Roe_M, Roe_dw, Roe_fluxA);
        
        // STEP 4:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_LD: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux: LMFIX
// Reference: A Low-Mach Number Fix for Roe's Approximate Riemann Solver
//            Felix Rieper, Journal of Computational Physics, 230 (2011)
//------------------------------------------------------------------------------
void Compute_Flux_Roe_LMFIX(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, maxcount;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, ubar1, ubar2;
    double sigma, area, nx, ny, nz, maxlambda;
    double Mach, nmax, fix;
    Vector3D Vn1, Vn2;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        
        // STEP 2:
        // Compute Low Mach Number Fix
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
            error("Compute_Flux_Roe_LMFIX: Unable to compute point on the plane (nx, ny, nz): (%lf, %lf, %lf)", nx, ny, nz);
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
        
        // STEP 3:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VARIABLE_PRIMITIVE_RUP, rho, u, v, w, c, Roe_M);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = rho_R - rho_L;
        Roe_dQ[1] = u_R - u_L;
        Roe_dQ[2] = v_R - v_L;
        Roe_dQ[3] = w_R - w_L;
        Roe_dQ[4] = p_R - p_L;
        
        double E0, E1, E2, dU;
        E0 = Roe_Eigen[0][0];
        E1 = Roe_Eigen[3][3];
        E2 = Roe_Eigen[4][4];
        dU = fix*(Roe_dQ[1]*nx + Roe_dQ[2]*ny + Roe_dQ[3]*nz); // Scale the Normal Component of dU
        Roe_dw[0] = E0*Roe_dQ[0] + (E1 + E2 - 2.0*E0)*Roe_dQ[4]/(2.0*c*c) + rho*(E1 - E2)*dU/(2.0*c);
        Roe_dw[1] = E0*Roe_dQ[1] + 0.5*nx*dU*((E1 + E2) - 2.0*E0) + nx*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[2] = E0*Roe_dQ[2] + 0.5*ny*dU*((E1 + E2) - 2.0*E0) + ny*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[3] = E0*Roe_dQ[3] + 0.5*nz*dU*((E1 + E2) - 2.0*E0) + nz*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[4] = 0.5*Roe_dQ[4]*(E1 + E2) + 0.5*rho*c*dU*(E1 - E2);
        MC_Matrix_Mul_Vector(5, 5, Roe_M, Roe_dw, Roe_fluxA);
        
        // STEP 4:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_LMFIX: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux: Optimized
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Optimized(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
             
        // STEP 2:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VARIABLE_PRIMITIVE_RUP, rho, u, v, w, c, Roe_M);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = rho_R - rho_L;
        Roe_dQ[1] = u_R - u_L;
        Roe_dQ[2] = v_R - v_L;
        Roe_dQ[3] = w_R - w_L;
        Roe_dQ[4] = p_R - p_L;
        
        double E0, E1, E2, dU;
        E0 = Roe_Eigen[0][0];
        E1 = Roe_Eigen[3][3];
        E2 = Roe_Eigen[4][4];
        dU = Roe_dQ[1]*nx + Roe_dQ[2]*ny + Roe_dQ[3]*nz;
        Roe_dw[0] = E0*Roe_dQ[0] + (E1 + E2 - 2.0*E0)*Roe_dQ[4]/(2.0*c*c) + rho*(E1 - E2)*dU/(2.0*c);
        Roe_dw[1] = E0*Roe_dQ[1] + 0.5*nx*dU*((E1 + E2) - 2.0*E0) + nx*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[2] = E0*Roe_dQ[2] + 0.5*ny*dU*((E1 + E2) - 2.0*E0) + ny*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[3] = E0*Roe_dQ[3] + 0.5*nz*dU*((E1 + E2) - 2.0*E0) + nz*(E1 - E2)*Roe_dQ[4]/(2.0*rho*c);
        Roe_dw[4] = 0.5*Roe_dQ[4]*(E1 + E2) + 0.5*rho*c*dU*(E1 - E2);
        MC_Matrix_Mul_Vector(5, 5, Roe_M, Roe_dw, Roe_fluxA);
        
        // STEP 3:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Optimized: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux Original: Unmodified
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Original(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        
        // STEP 2:
        // Compute the Left Eigenvector
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
        
        // STEP 3:
        // Compute the Right Eigenvector
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
        
        // STEP 4:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VARIABLE_PRIMITIVE_RUP, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(VARIABLE_PRIMITIVE_RUP, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Calculate L = Mop*P
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_P, Roe_T);
        
        // Calculate R = Pinv*Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++)
            for (j = 0; j < 5; j++)
                Roe_P[i][j] = 0.0;
        
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];

        // Compute |A|*dqr
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 5:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Original: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Cecile Voizat Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Precondition_CecileVoizat(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
            mach  = MAX(mach, mach_max);
            if (mach < Ref_Mach)
                sigma = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
            else
                sigma = 0.5;
            sigma = sigma*(sqrt(mach*Ref_Mach)+ mach);
            sigma = MIN(1.0, sigma);
            
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0/(sigma*sigma);
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
        Roe_P[0][0] =  0.0;
        Roe_P[0][1] =  0.0;
        Roe_P[0][2] =  0.0;
        Roe_P[0][3] =  2.0*c*c*rho*sigma*sigma/zeta;
        Roe_P[0][4] =  2.0*c*c*rho*sigma*sigma/beta;

        Roe_P[1][0] =  0.0;
        Roe_P[1][1] = -nz;
        Roe_P[1][2] =  ny;
        Roe_P[1][3] =  nx;
        Roe_P[1][4] = -nx;

        Roe_P[2][0] =  nz;
        Roe_P[2][1] =  0.0;
        Roe_P[2][2] = -nx;
        Roe_P[2][3] =  ny;
        Roe_P[2][4] = -ny;

        Roe_P[3][0] = -ny;
        Roe_P[3][1] =  nx;
        Roe_P[3][2] =  0.0;
        Roe_P[3][3] =  nz;
        Roe_P[3][4] = -nz;

        Roe_P[4][0] = -Gamma*nx/rho;
        Roe_P[4][1] = -Gamma*ny/rho;
        Roe_P[4][2] = -Gamma*nz/rho;
        Roe_P[4][3] =  0.0;
        Roe_P[4][4] =  0.0;
        
        // STEP 6:
        // Compute the Right Eigenvector
        Roe_Pinv[0][0] =  0.0;
        Roe_Pinv[0][1] =  0.0;
        Roe_Pinv[0][2] =  nz;
        Roe_Pinv[0][3] = -ny;
        Roe_Pinv[0][4] = -rho*nx/Gamma;

        Roe_Pinv[1][0] =  0.0;
        Roe_Pinv[1][1] = -nz;
        Roe_Pinv[1][2] =  0.0;
        Roe_Pinv[1][3] =  nx;
        Roe_Pinv[1][4] = -rho*ny/Gamma;

        Roe_Pinv[2][0] =  0.0;
        Roe_Pinv[2][1] =  ny;
        Roe_Pinv[2][2] = -nx;
        Roe_Pinv[2][3] =  0.0;
        Roe_Pinv[2][4] = -rho*nz/Gamma;

        Roe_Pinv[3][0] =  1.0/(rho*alpha);
        Roe_Pinv[3][1] =  0.5*zeta*nx/alpha;
        Roe_Pinv[3][2] =  0.5*zeta*ny/alpha;
        Roe_Pinv[3][3] =  0.5*zeta*nz/alpha;
        Roe_Pinv[3][4] =  0.0;

        Roe_Pinv[4][0] =  1.0/(rho*alpha);
        Roe_Pinv[4][1] = -0.5*beta*nx/alpha;
        Roe_Pinv[4][2] = -0.5*beta*ny/alpha;
        Roe_Pinv[4][3] = -0.5*beta*nz/alpha;
        Roe_Pinv[4][4] =  0.0;
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
                Roe_P[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        
        // Compute |A|*dqr
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 8:
        // Compute the Flux_Roe_Conv and Flux_Roe_Diss
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Precondition_CecileVoizat: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Briley Taylor Whitfield Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Precondition_BTW(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        ubar = u*nx + v*ny + w*nz;
        // Compute Square of Speed of Sound : This is done for stability
        c    = Get_SpeedSoundSquare(u, v, w, ht);
        if (c < 0.0)
            c = fabs(c);
        c = sqrt(c);
        mach = sqrt(q2)/c;
         
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0/beta;
        
        // STEP 5:
        // Compute the Left Eigenvector
        Roe_P[0][0] =  nx;
        Roe_P[0][1] =  ny;
        Roe_P[0][2] =  nz;
        Roe_P[0][3] =  (rho*chi_p)/tau;
        Roe_P[0][4] = -(rho*chi_m)/tau;

        Roe_P[1][0] =  0.0;
        Roe_P[1][1] = -nz;
        Roe_P[1][2] =  ny;
        Roe_P[1][3] =  nx;
        Roe_P[1][4] = -nx;

        Roe_P[2][0] =  nz;
        Roe_P[2][1] =  0.0;
        Roe_P[2][2] = -nx;
        Roe_P[2][3] =  ny;
        Roe_P[2][4] = -ny;

        Roe_P[3][0] = -ny;
        Roe_P[3][1] =  nx;
        Roe_P[3][2] =  0.0;
        Roe_P[3][3] =  nz;
        Roe_P[3][4] = -nz;

        Roe_P[4][0] =  0.0;
        Roe_P[4][1] =  0.0;
        Roe_P[4][2] =  0.0;
        Roe_P[4][3] = -rho*chi_m;
        Roe_P[4][4] =  rho*chi_p;
              
        // STEP 6:
        // Compute the Right Eigenvector
        Roe_Pinv[0][0] =  nx;
        Roe_Pinv[0][1] = -omega*nx*nx;
        Roe_Pinv[0][2] = -omega*nx*ny + nz;
        Roe_Pinv[0][3] = -omega*nx*nz - ny;
        Roe_Pinv[0][4] = -nx/tau;

        Roe_Pinv[1][0] =  ny;
        Roe_Pinv[1][1] = -omega*ny*nx - nz;
        Roe_Pinv[1][2] = -omega*ny*ny;
        Roe_Pinv[1][3] = -omega*ny*nz + nx;
        Roe_Pinv[1][4] = -ny/tau;

        Roe_Pinv[2][0] =  nz;
        Roe_Pinv[2][1] = -omega*nz*nx + ny;
        Roe_Pinv[2][2] = -omega*nz*ny - nx;
        Roe_Pinv[2][3] = -omega*nz*nz;
        Roe_Pinv[2][4] = -nz/tau;

        Roe_Pinv[3][0] =  0.0;
        Roe_Pinv[3][1] =  nx*psi_p;
        Roe_Pinv[3][2] =  ny*psi_p;
        Roe_Pinv[3][3] =  nz*psi_p;
        Roe_Pinv[3][4] =  1.0/(2.0*rho*sigma);

        Roe_Pinv[4][0] =  0.0;
        Roe_Pinv[4][1] =  nx*psi_m;
        Roe_Pinv[4][2] =  ny*psi_m;
        Roe_Pinv[4][3] =  nz*psi_m;
        Roe_Pinv[4][4] =  1.0/(2.0*rho*sigma);
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
                Roe_P[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        
        // Compute |A|*dqr
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 8:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Precondition_BTW: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Eriksson Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Precondition_Eriksson(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double alpha, beta, beta_m, beta_p;
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
        beta_m  = 0.5*(1.0 - beta);
        beta_p  = 0.5*(1.0 + beta);
        alpha   = 2.0*sqrt(ubar*ubar*beta_m*beta_m + c*c*beta);
        chi_m   = ubar*beta_m - 0.5*alpha;
        chi_p   = ubar*beta_m + 0.5*alpha;
        psi_m   = (c*c*beta + 2.0*beta_m*ubar*chi_m)/(c*c*chi_m);
        psi_p   = (c*c*beta + 2.0*beta_m*ubar*chi_p)/(c*c*chi_p);
        
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = (1.0 - beta)/(c*c*beta);

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0/beta;
        
        // STEP 5:
        // Compute the Left Eigenvector
        Roe_P[0][0] =  nx;
        Roe_P[0][1] =  ny;
        Roe_P[0][2] =  nz;
        Roe_P[0][3] = -rho*psi_m;
        Roe_P[0][4] =  rho*psi_p;

        Roe_P[1][0] =  0.0;
        Roe_P[1][1] = -nz;
        Roe_P[1][2] =  ny;
        Roe_P[1][3] =  nx;
        Roe_P[1][4] = -nx;

        Roe_P[2][0] =  nz;
        Roe_P[2][1] =  0.0;
        Roe_P[2][2] = -nx;
        Roe_P[2][3] =  ny;
        Roe_P[2][4] = -ny;

        Roe_P[3][0] = -ny;
        Roe_P[3][1] =  nx;
        Roe_P[3][2] =  0.0;
        Roe_P[3][3] =  nz;
        Roe_P[3][4] = -nz;

        Roe_P[4][0] =  0.0;
        Roe_P[4][1] =  0.0;
        Roe_P[4][2] =  0.0;
        Roe_P[4][3] = -rho*chi_m;
        Roe_P[4][4] =  rho*chi_p;
        
        // STEP 6:
        // Compute the Right Eigenvector
        Roe_Pinv[0][0] =  nx;
        Roe_Pinv[0][1] =  0.0;
        Roe_Pinv[0][2] =  nz;
        Roe_Pinv[0][3] = -ny;
        Roe_Pinv[0][4] = -nx/(c*c);

        Roe_Pinv[1][0] =  ny;
        Roe_Pinv[1][1] = -nz;
        Roe_Pinv[1][2] =  0.0;
        Roe_Pinv[1][3] =  nx;
        Roe_Pinv[1][4] = -ny/(c*c);

        Roe_Pinv[2][0] =  nz;
        Roe_Pinv[2][1] =  ny;
        Roe_Pinv[2][2] = -nx;
        Roe_Pinv[2][3] =  0.0;
        Roe_Pinv[2][4] = -nz/(c*c);

        Roe_Pinv[3][0] = 0.0;
        Roe_Pinv[3][1] = nx*chi_p/alpha;
        Roe_Pinv[3][2] = ny*chi_p/alpha;
        Roe_Pinv[3][3] = nz*chi_p/alpha;
        Roe_Pinv[3][4] = 1.0/(rho*alpha);

        Roe_Pinv[4][0] = 0.0;
        Roe_Pinv[4][1] = nx*chi_m/alpha;
        Roe_Pinv[4][2] = ny*chi_m/alpha;
        Roe_Pinv[4][3] = nz*chi_m/alpha;
        Roe_Pinv[4][4] = 1.0/(rho*alpha);
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
                Roe_P[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        
        // Compute |A|*dqr
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 8:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Precondition_Eriksson: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Merkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        double beta;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
                
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar + c);
        Roe_Eigen[4][4] = fabs(ubar - c);
        
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
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
        
        // STEP 6:
        // Compute the Right Eigenvector
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
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
                Roe_P[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        
        // Compute |A|*dqr
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 8:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Precondition_Merkel: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Turkel Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        double beta;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
                
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar + c);
        Roe_Eigen[4][4] = fabs(ubar - c);
        
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
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
        
        // STEP 6:
        // Compute the Right Eigenvector
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
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
                Roe_P[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        
        // Compute |A|*dqr
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 8:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Precondition_Turkel: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Weiss Smith Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_Flux_Roe_Precondition_WeissSmith(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz, maxlambda;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe_Conv[i] = Flux_Roe_Diss[i] = 0.0;
    Roe_Reset_New();
    
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
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        double beta;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
                
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar + c);
        Roe_Eigen[4][4] = fabs(ubar - c);
        
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
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
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
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
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
        
        // STEP 6:
        // Compute the Right Eigenvector
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
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
                Roe_P[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);
        
        // Compute dqr: Based on VariableType
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        
        // Compute |A|*dqr
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // STEP 8:
        // Compute the Convective Flux and Dissipative Flux
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = -0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_Flux_Roe_Precondition_WeissSmith: Invalid Node - %d", node_L);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux: LMFIX
//! Note this function should not be called from Compute_Flux_Roe_****
// Reference: A Low-Mach Number Fix for Roe's Approximate Riemann Solver
//            Felix Rieper, Journal of Computational Physics, 230 (2011)
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_LMFIX(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    int maxcount;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, ubar1, ubar2;
    double sigma, area, nx, ny, nz;
    double Mach, nmax, fix;
    Vector3D Vn1, Vn2;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        
        // STEP 2:
        // Compute Low Mach Number Fix
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
            error("Compute_Flux_Roe_LMFIX: Unable to compute point on the plane (nx, ny, nz): (%lf, %lf, %lf)", nx, ny, nz);
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
        
        // STEP 3:
        // |A| = P.|Lambda|.Pinv
        double t, t1, t2;
        t  = Roe_Eigen[0][0];
        t1 = Roe_Eigen[3][3];
        t2 = Roe_Eigen[4][4];
        
        Roe_A[0][0] = t;
        Roe_A[0][1] = 0.5*fix*rho*nx*(t1 - t2)/c;
        Roe_A[0][2] = 0.5*fix*rho*ny*(t1 - t2)/c;
        Roe_A[0][3] = 0.5*fix*rho*nz*(t1 - t2)/c;
        Roe_A[0][4] = 0.5*(t1 + t2 - 2.0*t)/(c*c);

        Roe_A[1][0] = 0.0;
        Roe_A[1][1] = 0.5*(2.0*t + fix*nx*nx*(t1 + t2 - 2.0*t));
        Roe_A[1][2] = 0.5*fix*ny*nx*(t1 + t2 - 2.0*t);
        Roe_A[1][3] = 0.5*fix*nz*nx*(t1 + t2 - 2.0*t);
        Roe_A[1][4] = 0.5*nx*(t1 - t2)/(rho*c);

        Roe_A[2][0] = 0.0;
        Roe_A[2][1] = 0.5*fix*nx*ny*(t1 + t2 - 2.0*t);
        Roe_A[2][2] = 0.5*(2.0*t + fix*ny*ny*(t1 + t2 - 2.0*t));
        Roe_A[2][3] = 0.5*fix*nz*ny*(t1 + t2 - 2.0*t);
        Roe_A[2][4] = 0.5*ny*(t1 - t2)/(rho*c);

        Roe_A[3][0] = 0.0;
        Roe_A[3][1] = 0.5*fix*nx*nz*(t1 + t2 - 2.0*t);
        Roe_A[3][2] = 0.5*fix*ny*nz*(t1 + t2 - 2.0*t);
        Roe_A[3][3] = 0.5*(2.0*t + fix*nz*nz*(t1 + t2 - 2.0*t));
        Roe_A[3][4] = 0.5*nz*(t1 - t2)/(rho*c);

        Roe_A[4][0] = 0.0;
        Roe_A[4][1] = 0.5*fix*rho*c*nx*(t1 - t2);
        Roe_A[4][2] = 0.5*fix*rho*c*ny*(t1 - t2);
        Roe_A[4][3] = 0.5*fix*rho*c*nz*(t1 - t2);
        Roe_A[4][4] = 0.5*(t1 + t2);
        
        // STEP 4:
        // Compute the Dissipation Jacobian
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VARIABLE_PRIMITIVE_RUP, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(VARIABLE_PRIMITIVE_RUP, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Calculate Mop*|A|*Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_Minv, Dissipation_Matrix_Roe);
     } else
        error("Compute_Dissipation_Matrix_Roe_LMFIX: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux: Optimized
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Optimized(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        
        // STEP 2:
        // |A| = P.|Lambda|.Pinv
        double t, t1, t2;
        t  = Roe_Eigen[0][0];
        t1 = Roe_Eigen[3][3];
        t2 = Roe_Eigen[4][4];
        
        Roe_A[0][0] = t;
        Roe_A[0][1] = 0.5*rho*nx*(t1 - t2)/c;
        Roe_A[0][2] = 0.5*rho*ny*(t1 - t2)/c;
        Roe_A[0][3] = 0.5*rho*nz*(t1 - t2)/c;
        Roe_A[0][4] = 0.5*(t1 + t2 - 2.0*t)/(c*c);

        Roe_A[1][0] = 0.0;
        Roe_A[1][1] = 0.5*(2.0*t + nx*nx*(t1 + t2 - 2.0*t));
        Roe_A[1][2] = 0.5*ny*nx*(t1 + t2 - 2.0*t);
        Roe_A[1][3] = 0.5*nz*nx*(t1 + t2 - 2.0*t);
        Roe_A[1][4] = 0.5*nx*(t1 - t2)/(rho*c);

        Roe_A[2][0] = 0.0;
        Roe_A[2][1] = 0.5*nx*ny*(t1 + t2 - 2.0*t);
        Roe_A[2][2] = 0.5*(2.0*t + ny*ny*(t1 + t2 - 2.0*t));
        Roe_A[2][3] = 0.5*nz*ny*(t1 + t2 - 2.0*t);
        Roe_A[2][4] = 0.5*ny*(t1 - t2)/(rho*c);

        Roe_A[3][0] = 0.0;
        Roe_A[3][1] = 0.5*nx*nz*(t1 + t2 - 2.0*t);
        Roe_A[3][2] = 0.5*ny*nz*(t1 + t2 - 2.0*t);
        Roe_A[3][3] = 0.5*(2.0*t + nz*nz*(t1 + t2 - 2.0*t));
        Roe_A[3][4] = 0.5*nz*(t1 - t2)/(rho*c);

        Roe_A[4][0] = 0.0;
        Roe_A[4][1] = 0.5*rho*c*nx*(t1 - t2);
        Roe_A[4][2] = 0.5*rho*c*ny*(t1 - t2);
        Roe_A[4][3] = 0.5*rho*c*nz*(t1 - t2);
        Roe_A[4][4] = 0.5*(t1 + t2);
        
        // STEP 3:
        // Compute the Dissipation Jacobian
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VARIABLE_PRIMITIVE_RUP, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(VARIABLE_PRIMITIVE_RUP, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Calculate Mop*|A|*Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_Minv, Dissipation_Matrix_Roe);
     } else
        error("Compute_Dissipation_Matrix_Roe_Optimized: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux: Original
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Original(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        c     = Get_SpeedSound(u, v, w, ht);
        ubar  = u*nx + v*ny + w*nz;
        
        //======================================================================
        // Compute Dissipation Flux
        //======================================================================
        // STEP 1:
        // Compute the Eigenvalues of the Dissipation Term |Lambda|
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
        
        // STEP 2:
        // Compute the Left Eigenvector
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
        
        // STEP 3:
        // Compute the Right Eigenvector
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
        
        // STEP 4:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, VARIABLE_PRIMITIVE_RUP, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(VARIABLE_PRIMITIVE_RUP, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Calculate L = Mop*P
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_P, Roe_T);
        
        // Calculate R = Pinv*Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, Dissipation_Matrix_Roe);
     } else
        error("Compute_Dissipation_Matrix_Roe_Original: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux with Cecile Voizat Pre-Conditioner
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Precondition_CecileVoizat(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
                
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
             
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
            mach  = MAX(mach, mach_max);
            if (mach < Ref_Mach)
                sigma = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
            else
                sigma = 0.5;
            sigma = sigma*(sqrt(mach*Ref_Mach)+ mach);
            sigma = MIN(1.0, sigma);
            
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
        
        // STEP 4:
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0/(sigma*sigma);
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
        Roe_P[0][0] =  0.0;
        Roe_P[0][1] =  0.0;
        Roe_P[0][2] =  0.0;
        Roe_P[0][3] =  2.0*c*c*rho*sigma*sigma/zeta;
        Roe_P[0][4] =  2.0*c*c*rho*sigma*sigma/beta;

        Roe_P[1][0] =  0.0;
        Roe_P[1][1] = -nz;
        Roe_P[1][2] =  ny;
        Roe_P[1][3] =  nx;
        Roe_P[1][4] = -nx;

        Roe_P[2][0] =  nz;
        Roe_P[2][1] =  0.0;
        Roe_P[2][2] = -nx;
        Roe_P[2][3] =  ny;
        Roe_P[2][4] = -ny;

        Roe_P[3][0] = -ny;
        Roe_P[3][1] =  nx;
        Roe_P[3][2] =  0.0;
        Roe_P[3][3] =  nz;
        Roe_P[3][4] = -nz;

        Roe_P[4][0] = -Gamma*nx/rho;
        Roe_P[4][1] = -Gamma*ny/rho;
        Roe_P[4][2] = -Gamma*nz/rho;
        Roe_P[4][3] =  0.0;
        Roe_P[4][4] =  0.0;
        
        // STEP 6:
        // Compute the Right Eigenvector
        Roe_Pinv[0][0] =  0.0;
        Roe_Pinv[0][1] =  0.0;
        Roe_Pinv[0][2] =  nz;
        Roe_Pinv[0][3] = -ny;
        Roe_Pinv[0][4] = -rho*nx/Gamma;

        Roe_Pinv[1][0] =  0.0;
        Roe_Pinv[1][1] = -nz;
        Roe_Pinv[1][2] =  0.0;
        Roe_Pinv[1][3] =  nx;
        Roe_Pinv[1][4] = -rho*ny/Gamma;

        Roe_Pinv[2][0] =  0.0;
        Roe_Pinv[2][1] =  ny;
        Roe_Pinv[2][2] = -nx;
        Roe_Pinv[2][3] =  0.0;
        Roe_Pinv[2][4] = -rho*nz/Gamma;

        Roe_Pinv[3][0] =  1.0/(rho*alpha);
        Roe_Pinv[3][1] =  0.5*zeta*nx/alpha;
        Roe_Pinv[3][2] =  0.5*zeta*ny/alpha;
        Roe_Pinv[3][3] =  0.5*zeta*nz/alpha;
        Roe_Pinv[3][4] =  0.0;

        Roe_Pinv[4][0] =  1.0/(rho*alpha);
        Roe_Pinv[4][1] = -0.5*beta*nx/alpha;
        Roe_Pinv[4][2] = -0.5*beta*ny/alpha;
        Roe_Pinv[4][3] = -0.5*beta*nz/alpha;
        Roe_Pinv[4][4] =  0.0;
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, Dissipation_Matrix_Roe);
    } else
        error("Compute_Dissipation_Matrix_Roe_Precondition_CecileVoizat: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux with Briley Taylor Whitfield Pre-Conditioner
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Precondition_BTW(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
                
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
        ubar = u*nx + v*ny + w*nz;
        // Compute Square of Speed of Sound : This is done for stability
        c    = Get_SpeedSoundSquare(u, v, w, ht);
        if (c < 0.0)
            c = fabs(c);
        c = sqrt(c);
        mach = sqrt(q2)/c;
         
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        // STEP 4:
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0/beta;
        
        // STEP 5:
        // Compute the Left Eigenvector
        Roe_P[0][0] =  nx;
        Roe_P[0][1] =  ny;
        Roe_P[0][2] =  nz;
        Roe_P[0][3] =  (rho*chi_p)/tau;
        Roe_P[0][4] = -(rho*chi_m)/tau;

        Roe_P[1][0] =  0.0;
        Roe_P[1][1] = -nz;
        Roe_P[1][2] =  ny;
        Roe_P[1][3] =  nx;
        Roe_P[1][4] = -nx;

        Roe_P[2][0] =  nz;
        Roe_P[2][1] =  0.0;
        Roe_P[2][2] = -nx;
        Roe_P[2][3] =  ny;
        Roe_P[2][4] = -ny;

        Roe_P[3][0] = -ny;
        Roe_P[3][1] =  nx;
        Roe_P[3][2] =  0.0;
        Roe_P[3][3] =  nz;
        Roe_P[3][4] = -nz;

        Roe_P[4][0] =  0.0;
        Roe_P[4][1] =  0.0;
        Roe_P[4][2] =  0.0;
        Roe_P[4][3] = -rho*chi_m;
        Roe_P[4][4] =  rho*chi_p;
              
        // STEP 6:
        // Compute the Right Eigenvector
        Roe_Pinv[0][0] =  nx;
        Roe_Pinv[0][1] = -omega*nx*nx;
        Roe_Pinv[0][2] = -omega*nx*ny + nz;
        Roe_Pinv[0][3] = -omega*nx*nz - ny;
        Roe_Pinv[0][4] = -nx/tau;

        Roe_Pinv[1][0] =  ny;
        Roe_Pinv[1][1] = -omega*ny*nx - nz;
        Roe_Pinv[1][2] = -omega*ny*ny;
        Roe_Pinv[1][3] = -omega*ny*nz + nx;
        Roe_Pinv[1][4] = -ny/tau;

        Roe_Pinv[2][0] =  nz;
        Roe_Pinv[2][1] = -omega*nz*nx + ny;
        Roe_Pinv[2][2] = -omega*nz*ny - nx;
        Roe_Pinv[2][3] = -omega*nz*nz;
        Roe_Pinv[2][4] = -nz/tau;

        Roe_Pinv[3][0] =  0.0;
        Roe_Pinv[3][1] =  nx*psi_p;
        Roe_Pinv[3][2] =  ny*psi_p;
        Roe_Pinv[3][3] =  nz*psi_p;
        Roe_Pinv[3][4] =  1.0/(2.0*rho*sigma);

        Roe_Pinv[4][0] =  0.0;
        Roe_Pinv[4][1] =  nx*psi_m;
        Roe_Pinv[4][2] =  ny*psi_m;
        Roe_Pinv[4][3] =  nz*psi_m;
        Roe_Pinv[4][4] =  1.0/(2.0*rho*sigma);
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, Dissipation_Matrix_Roe);
    } else
        error("Compute_Dissipation_Matrix_Roe_Precondition_BTW: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux with Eriksson Pre-Conditioner
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Precondition_Eriksson(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double alpha, beta, beta_m, beta_p;
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
        beta_m  = 0.5*(1.0 - beta);
        beta_p  = 0.5*(1.0 + beta);
        alpha   = 2.0*sqrt(ubar*ubar*beta_m*beta_m + c*c*beta);
        chi_m   = ubar*beta_m - 0.5*alpha;
        chi_p   = ubar*beta_m + 0.5*alpha;
        psi_m   = (c*c*beta + 2.0*beta_m*ubar*chi_m)/(c*c*chi_m);
        psi_p   = (c*c*beta + 2.0*beta_m*ubar*chi_p)/(c*c*chi_p);
        
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar*beta_p + 0.5*alpha);
        Roe_Eigen[4][4] = fabs(ubar*beta_p - 0.5*alpha);
        
        // STEP 4:
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = (1.0 - beta)/(c*c*beta);

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0/beta;
        
        // STEP 5:
        // Compute the Left Eigenvector
        Roe_P[0][0] =  nx;
        Roe_P[0][1] =  ny;
        Roe_P[0][2] =  nz;
        Roe_P[0][3] = -rho*psi_m;
        Roe_P[0][4] =  rho*psi_p;

        Roe_P[1][0] =  0.0;
        Roe_P[1][1] = -nz;
        Roe_P[1][2] =  ny;
        Roe_P[1][3] =  nx;
        Roe_P[1][4] = -nx;

        Roe_P[2][0] =  nz;
        Roe_P[2][1] =  0.0;
        Roe_P[2][2] = -nx;
        Roe_P[2][3] =  ny;
        Roe_P[2][4] = -ny;

        Roe_P[3][0] = -ny;
        Roe_P[3][1] =  nx;
        Roe_P[3][2] =  0.0;
        Roe_P[3][3] =  nz;
        Roe_P[3][4] = -nz;

        Roe_P[4][0] =  0.0;
        Roe_P[4][1] =  0.0;
        Roe_P[4][2] =  0.0;
        Roe_P[4][3] = -rho*chi_m;
        Roe_P[4][4] =  rho*chi_p;
        
        // STEP 6:
        // Compute the Right Eigenvector
        Roe_Pinv[0][0] =  nx;
        Roe_Pinv[0][1] =  0.0;
        Roe_Pinv[0][2] =  nz;
        Roe_Pinv[0][3] = -ny;
        Roe_Pinv[0][4] = -nx/(c*c);

        Roe_Pinv[1][0] =  ny;
        Roe_Pinv[1][1] = -nz;
        Roe_Pinv[1][2] =  0.0;
        Roe_Pinv[1][3] =  nx;
        Roe_Pinv[1][4] = -ny/(c*c);

        Roe_Pinv[2][0] =  nz;
        Roe_Pinv[2][1] =  ny;
        Roe_Pinv[2][2] = -nx;
        Roe_Pinv[2][3] =  0.0;
        Roe_Pinv[2][4] = -nz/(c*c);

        Roe_Pinv[3][0] = 0.0;
        Roe_Pinv[3][1] = nx*chi_p/alpha;
        Roe_Pinv[3][2] = ny*chi_p/alpha;
        Roe_Pinv[3][3] = nz*chi_p/alpha;
        Roe_Pinv[3][4] = 1.0/(rho*alpha);

        Roe_Pinv[4][0] = 0.0;
        Roe_Pinv[4][1] = nx*chi_m/alpha;
        Roe_Pinv[4][2] = ny*chi_m/alpha;
        Roe_Pinv[4][3] = nz*chi_m/alpha;
        Roe_Pinv[4][4] = 1.0/(rho*alpha);
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, Dissipation_Matrix_Roe);       
    } else
        error("Compute_Dissipation_Matrix_Roe_Precondition_Eriksson: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux with Merkel Pre-Conditioner
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Precondition_Merkel(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        double beta;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
                
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar + c);
        Roe_Eigen[4][4] = fabs(ubar - c);
        
        // STEP 4:
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
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
        
        // STEP 6:
        // Compute the Right Eigenvector
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
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, Dissipation_Matrix_Roe);
    } else
        error("Compute_Dissipation_Matrix_Roe_Precondition_Merkel: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux with Turkel Pre-Conditioner
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Precondition_Turkel(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        double beta;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
                
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar + c);
        Roe_Eigen[4][4] = fabs(ubar - c);
        
        // STEP 4:
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
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
        
        // STEP 6:
        // Compute the Right Eigenvector
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
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, Dissipation_Matrix_Roe);
    } else
        error("Compute_Dissipation_Matrix_Roe_Precondition_Turkel: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Dissipation Matrix Roe Flux with Weiss Smith Pre-Conditioner
//! Note this function should not be called from Compute_Flux_Roe_****
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe_Precondition_WeissSmith(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
    int i, j, AvgType;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, ubar, q2, mach;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset_New();
    
    // Only for Physical Nodes
    if (node_L < nNode) {
        // Get area vector
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Make Solution Second Order
        // Note: For Jacobian should not be made second order for stability reasons
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
        
        // Compute Equation of State
        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
        Compute_EOS_Variables_Face(Roe_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(Roe_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        
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
        ubar = u*nx + v*ny + w*nz;
        c    = Get_SpeedSound(u, v, w, ht);
        mach = sqrt(q2)/c;
        
        int nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_L, dp_R, lmach_L, lmach_R;
        double dp_max, mach_max;
        double beta;
        
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
                    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
                        // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
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
        
        beta    = beta*beta;
                
        // STEP 3:
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(ubar + c);
        Roe_Eigen[4][4] = fabs(ubar - c);
        
        // STEP 4:
        // Compute the Precondition Matrix Inv
        Roe_Kinv[0][0] = 1.0;
        Roe_Kinv[0][1] = 0.0;
        Roe_Kinv[0][2] = 0.0;
        Roe_Kinv[0][3] = 0.0;
        Roe_Kinv[0][4] = 0.0;

        Roe_Kinv[1][0] = 0.0;
        Roe_Kinv[1][1] = 1.0;
        Roe_Kinv[1][2] = 0.0;
        Roe_Kinv[1][3] = 0.0;
        Roe_Kinv[1][4] = 0.0;

        Roe_Kinv[2][0] = 0.0;
        Roe_Kinv[2][1] = 0.0;
        Roe_Kinv[2][2] = 1.0;
        Roe_Kinv[2][3] = 0.0;
        Roe_Kinv[2][4] = 0.0;

        Roe_Kinv[3][0] = 0.0;
        Roe_Kinv[3][1] = 0.0;
        Roe_Kinv[3][2] = 0.0;
        Roe_Kinv[3][3] = 1.0;
        Roe_Kinv[3][4] = 0.0;

        Roe_Kinv[4][0] = 0.0;
        Roe_Kinv[4][1] = 0.0;
        Roe_Kinv[4][2] = 0.0;
        Roe_Kinv[4][3] = 0.0;
        Roe_Kinv[4][4] = 1.0;
        
        // STEP 5:
        // Compute the Left Eigenvector
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
        
        // STEP 6:
        // Compute the Right Eigenvector
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
        
        // STEP 7:
        // Compute the Dissipation Flux
        // Mop
        Compute_Transformation_Matrix(VARIABLE_CONSERVATIVE, PrecondVariableType, rho, u, v, w, c, Roe_M);
        
        // Mpr
        Compute_Transformation_Matrix(PrecondVariableType, VariableType, rho, u, v, w, c, Roe_Minv);
        
        // Compute the L = Mop.Kinv.P
        MC_Matrix_Mul_Matrix(5, 5, Roe_Kinv, Roe_P, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_A, Roe_T);
        
        // Compute the R = Pinv.Mpr
        MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
        // Compute |A| = L.|Eigen|.R
        for (i = 0; i < 5; i++) {
            for (j = 0; j < 5; j++) {
                Roe_A[i][j] = 0.0;
            }
        }
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_A);
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_A, Dissipation_Matrix_Roe);
    } else
        error("Compute_Dissipation_Matrix_Roe_Precondition_WeissSmith: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Computes the Dissipation Matrix Roe Flux
//------------------------------------------------------------------------------
void Compute_Dissipation_Matrix_Roe(int node_L, int node_R, Vector3D areavec, double **Dissipation_Matrix_Roe) {
        
    // Compute the Roe Flux for this edge
    switch (PrecondMethod) {
        case PRECOND_METHOD_NONE: // Roe
            Compute_Dissipation_Matrix_Roe_Original(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        case PRECOND_METHOD_ROE_LMFIX: // LMRoe
            Compute_Dissipation_Matrix_Roe_LMFIX(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        case PRECOND_METHOD_ROE_WS: // Roe Weiss Smith Pre-Conditioner
            Compute_Dissipation_Matrix_Roe_Precondition_WeissSmith(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        case PRECOND_METHOD_ROE_CV: // Roe Cecile Voizat Pre-Conditioner
            Compute_Dissipation_Matrix_Roe_Precondition_CecileVoizat(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        case PRECOND_METHOD_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
            Compute_Dissipation_Matrix_Roe_Precondition_BTW(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        case PRECOND_METHOD_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner
            Compute_Dissipation_Matrix_Roe_Precondition_Eriksson(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        case PRECOND_METHOD_ROE_MERKEL: // Roe Merkel Pre-Conditioner
            Compute_Dissipation_Matrix_Roe_Precondition_Merkel(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        case PRECOND_METHOD_ROE_TURKEL: // Roe Turkel Pre-Conditioner
            Compute_Dissipation_Matrix_Roe_Precondition_Turkel(node_L, node_R, areavec, Dissipation_Matrix_Roe);
            break;
        default:
            error("Compute_Dissipation_Matrix_Roe: Invalid Solver Precondition Scheme - %d -1", PrecondMethod);
            break;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//------------------------------------------------------------------------------
//! Computes the Roe Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_Roe_New(int AddTime) {
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
            case PRECOND_METHOD_NONE: // Roe
                Compute_Flux_Roe_Optimized(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_LMFIX: // LMRoe
                Compute_Flux_Roe_LMFIX(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_WS: // Roe Weiss Smith Pre-Conditioner
                Compute_Flux_Roe_Precondition_WeissSmith(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_CV: // Roe Cecile Voizat Pre-Conditioner
                Compute_Flux_Roe_Precondition_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
                Compute_Flux_Roe_Precondition_BTW(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner
                Compute_Flux_Roe_Precondition_Eriksson(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_MERKEL: // Roe Merkel Pre-Conditioner
                Compute_Flux_Roe_Precondition_Merkel(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_TURKEL: // Roe Turkel Pre-Conditioner
                Compute_Flux_Roe_Precondition_Turkel(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            default:
                error("Compute_Residual_Roe_New: Invalid Solver Precondition Scheme - %d -1", PrecondMethod);
                break;
        }
        
        // Convective Term
        // L-Node
        Res1[node_L] += flux_roe[0];
        Res2[node_L] += flux_roe[1];
        Res3[node_L] += flux_roe[2];
        Res4[node_L] += flux_roe[3];
        Res5[node_L] += flux_roe[4];
        // R-Node
        Res1[node_R] -= flux_roe[0];
        Res2[node_R] -= flux_roe[1];
        Res3[node_R] -= flux_roe[2];
        Res4[node_R] -= flux_roe[3];
        Res5[node_R] -= flux_roe[4];
        
        // Dissipation Term
        // L-Node
        Res1_Diss[node_L] += flux_roe_diss[0];
        Res2_Diss[node_L] += flux_roe_diss[1];
        Res3_Diss[node_L] += flux_roe_diss[2];
        Res4_Diss[node_L] += flux_roe_diss[3];
        Res5_Diss[node_L] += flux_roe_diss[4];
        // R-Node
        Res1_Diss[node_R] -= flux_roe_diss[0];
        Res2_Diss[node_R] -= flux_roe_diss[1];
        Res3_Diss[node_R] -= flux_roe_diss[2];
        Res4_Diss[node_R] -= flux_roe_diss[3];
        Res5_Diss[node_R] -= flux_roe_diss[4];
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
            case PRECOND_METHOD_NONE: // Roe
                Compute_Flux_Roe_Optimized(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_LMFIX: // LMRoe
                Compute_Flux_Roe_LMFIX(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_WS: // Roe Weiss Smith Pre-Conditioner
                Compute_Flux_Roe_Precondition_WeissSmith(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_CV: // Roe Cecile Voizat Pre-Conditioner
                Compute_Flux_Roe_Precondition_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
                Compute_Flux_Roe_Precondition_BTW(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner
                Compute_Flux_Roe_Precondition_Eriksson(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_MERKEL: // Roe Merkel Pre-Conditioner
                Compute_Flux_Roe_Precondition_Merkel(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            case PRECOND_METHOD_ROE_TURKEL: // Roe Turkel Pre-Conditioner
                Compute_Flux_Roe_Precondition_Turkel(node_L, node_R, areavec, flux_roe, flux_roe_diss, AddTime);
                break;
            default:
                error("Compute_Residual_Roe_New: Invalid Solver Precondition Scheme - %d -2", PrecondMethod);
                break;
        }
        
        // Convective Term
        // L-Node
        Res1[node_L] += flux_roe[0];
        Res2[node_L] += flux_roe[1];
        Res3[node_L] += flux_roe[2];
        Res4[node_L] += flux_roe[3];
        Res5[node_L] += flux_roe[4];
        
        // Dissipation Term
        // L-Node
        Res1_Diss[node_L] += flux_roe_diss[0];
        Res2_Diss[node_L] += flux_roe_diss[1];
        Res3_Diss[node_L] += flux_roe_diss[2];
        Res4_Diss[node_L] += flux_roe_diss[3];
        Res5_Diss[node_L] += flux_roe_diss[4];
    }
    
    // Transform the Residuals into Desired variable Form
    Compute_Transformed_Residual_Roe();
    
    // Precondition the Residuals
    if ((SolverMethod == SOLVER_METHOD_STEADY) && (SolverScheme == SOLVER_SCHEME_EXPLICIT)) {
        switch (PrecondMethod) {
            case PRECOND_METHOD_NONE: // Roe
                break;
            case PRECOND_METHOD_ROE_LMFIX: // LMRoe
                break;
            case PRECOND_METHOD_ROE_WS: // Roe Weiss Smith Pre-Conditioner
                Compute_Steady_Residual_Roe_Precondition_WeissSmith();
                break;
            case PRECOND_METHOD_ROE_CV: // Roe Cecile Voizat Pre-Conditioner
                Compute_Steady_Residual_Roe_Precondition_CecileVoizat();
                break;
            case PRECOND_METHOD_ROE_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
                Compute_Steady_Residual_Roe_Precondition_BTW();
                break;
            case PRECOND_METHOD_ROE_ERIKSSON: // Roe Eriksson Pre-Conditioner
                Compute_Steady_Residual_Roe_Precondition_Eriksson();
                break;
            case PRECOND_METHOD_ROE_MERKEL: // Roe Merkel Pre-Conditioner
                Compute_Steady_Residual_Roe_Precondition_Merkel();
                break;
            case PRECOND_METHOD_ROE_TURKEL: // Roe Turkel Pre-Conditioner
                Compute_Steady_Residual_Roe_Precondition_Turkel();
                break;
            default:
                error("Compute_Residual_Roe_New: Invalid Solver Precondition Scheme - %d -3", PrecondMethod);
                break;
        }
    }
    
    // For now just add both
//    for (i = 0; i < nNode; i++) {
//        Res1[i] += Res1_Diss[i];
//        Res2[i] += Res2_Diss[i];
//        Res3[i] += Res3_Diss[i];
//        Res4[i] += Res4_Diss[i];
//        Res5[i] += Res5_Diss[i];
//    }
}

