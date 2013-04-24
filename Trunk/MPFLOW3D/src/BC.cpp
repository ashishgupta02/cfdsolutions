/*******************************************************************************
 * File:        BC.cpp
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
#include "Vector3D.h"
#include "Commons.h"
#include "Material.h"
#include "Solver.h"
#include "LinearAlgebra.h"

// Static Variable for Speed Up
static int      BC_DB   = 0;
static int     *BC_P    = NULL;
static double  *BC_L    = NULL;
static double  *BC_X    = NULL;
static double  *BC_B    = NULL;
static double **BC_A    = NULL;

//------------------------------------------------------------------------------
//! Create Boundary Condition Data Structure
//------------------------------------------------------------------------------
void BC_Init(void) {
    int i;
    double *dtmp = NULL;
    
    if (BC_DB == 0) {
        BC_P  = (int *)     malloc(NEQUATIONS*sizeof(int));
        BC_L  = (double *)  malloc(NEQUATIONS*sizeof(double));
        BC_X  = (double *)  malloc(NEQUATIONS*sizeof(double));
        BC_B  = (double *)  malloc(NEQUATIONS*sizeof(double));
        BC_A  = (double **) malloc(NEQUATIONS*sizeof(double*));
        dtmp  = (double *)  malloc(NEQUATIONS*NEQUATIONS*sizeof(double));
        for (i = 0; i < NEQUATIONS; i++)
            BC_A[i] = &dtmp[i*NEQUATIONS];
        dtmp = NULL;
        BC_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete Boundary Condition Data Structure
//------------------------------------------------------------------------------
void BC_Finalize(void) {
    double *tmp = NULL;
    
    tmp = BC_A[0];
    free(tmp);
    tmp = NULL;
    free(BC_A);
    free(BC_P);
    free(BC_L);
    free(BC_X);
    free(BC_B);
    BC_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset Boundary Condition Data Structure
//------------------------------------------------------------------------------
void BC_Reset(void) {
    int i, j;
    
    if (BC_DB == 0)
        BC_Init();
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        BC_P[i] = i;
        BC_L[i] = 0.0;
        BC_X[i] = 0.0;
        BC_B[i] = 0.0;
        for (j = 0; j < NEQUATIONS; j++)
            BC_A[i][j] = 0.0;
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Initialize_Boundary_Condition() {
    // Get the Boundary Type
    for (int i = 0; i < nBEdge; i++)
        bndEdge[i].type = bndType[bndEdge[i].tag - 1];
}

//------------------------------------------------------------------------------
//! Preconditioned Characteristic Based Boundary Condition: Merkel Numerically Solved
//! Note: Outward Pointing Boundary Normals
//! rvalue: 0 => Success, 1: Negative Pressure
//------------------------------------------------------------------------------
int Compute_Precondition_Characteristic_BoundaryCondition_Merkel(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    error("Compute_Precondition_Characteristic_BoundaryCondition_Merkel: Merkel Preconditioner Not Implemented - %d", PrecondMethod);
    return 0;
}

//------------------------------------------------------------------------------
//! Preconditioned Characteristic Based Boundary Condition: Turkel Numerically Solved
//! Note: Outward Pointing Boundary Normals
//! rvalue: 0 => Success, 1: Negative Pressure
//------------------------------------------------------------------------------
int Compute_Precondition_Characteristic_BoundaryCondition_Turkel(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    int i;
    double nx, ny, nz;
    double   rho,   u,   v,   w,   et,   p,   c,   T,   ubar,   ht,   q2,   mach;
    double rho_i, u_i, v_i, w_i, et_i, p_i, c_i, T_i, ubar_i, ht_i, q2_i, mach_i;
    double rho_b, u_b, v_b, w_b, et_b, p_b, c_b, T_b, ubar_b, ht_b, q2_b, mach_b;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double ubar_pre, c_pre;
    int instance = BC_TYPE_NONE;
    int rvalue = 0;
    int node_L;
    
    // Initialization
    BC_Reset();
    
    // Default Initialization: Based on Variable Type Pressure will be Perturbation
    p_b_n   = Inf_Pressure;
    rho_b_n = Inf_Rho;
    u_b_n   = Inf_U;
    v_b_n   = Inf_V;
    w_b_n   = Inf_W;
    
    // Get area vector
    AreaVec.normalize();
    nx = AreaVec.vec[0];
    ny = AreaVec.vec[1];
    nz = AreaVec.vec[2];

    // Compute Equation of State
    // Physical node: Based on Variable Type Pressure will be Perturbation
    Material_Get_Face_Properties(Q_L, nx, ny, nz, rho_i, p_i, T_i, u_i, v_i, w_i, q2_i, c_i, mach_i, ubar_i, et_i, ht_i);
    
    // Ghost node: Based on Variable Type Pressure will be Perturbation
    Material_Get_Face_Properties(Q_R, nx, ny, nz, rho_b, p_b, T_b, u_b, v_b, w_b, q2_b, c_b, mach_b, ubar_b, et_b, ht_b);
    
    // Average State: Based on Variable Type Pressure will be Perturbation
    for (i = 0; i < NEQUATIONS; i++)
        Q_B[i] = 0.5*(Q_L[i] + Q_R[i]);
    Material_Get_Face_Properties(Q_B, nx, ny, nz, rho, p, T, u, v, w, q2, c, mach, ubar, et, ht);
    
    // Get two nodes of edge
    node_L = bndEdge[BEdgeID].node[0];

    int nid;
    double lQ[NEQUATIONS];
    double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
    double dp_max, mach_max;

    //==========================================================================
    // Compute Precondition Boundary Condition
    //==========================================================================
    // STEP - 1:
    // Compute the Local Max Mach and Max Change in Pressure around a node L
    // No computation for node R as it is a ghost node
    dp_max   = 0.0;
    mach_max = 0.0;
    // Check if Local Precondition is Required
    if (PrecondType == PRECOND_TYPE_LOCAL) {
        // Check if Precondition Smoother is Requested
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
                Material_Get_ControlVolume_Properties(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                // Compute Smoothing Parameters
                mach_max = MAX(mach_max, lmach);
                dp_max   = MAX(dp_max, fabs(p_i - lp));
            }
        } else {
            mach_max = MAX(mach_i, mach_b);
            dp_max   = fabs(p_i - p_b);
        }
    }
    
    double delta, alpha, beta;
    double zp, zm;
    double pv1, pv2, pv3, pv4, pv5;
    
    // STEP - 2:
    // Compute Precondition Variables
    mach = mach;
    beta = 1.0; // Corresponds to no pre-conditioning
    if (PrecondType == PRECOND_TYPE_LOCAL) {
        mach = MAX(mach, mach_max);
        if (mach < PrecondGlobalMach)
            beta = 0.5*(mach*mach)/(PrecondGlobalMach*PrecondGlobalMach) - mach/(PrecondGlobalMach) + 1.0;
        else
            beta = 0.5;
        beta = beta*(sqrt(mach*PrecondGlobalMach)+ mach);
        beta = MIN(1.0, beta);
    }
    if (PrecondType == PRECOND_TYPE_GLOBAL)
        beta = MIN(1.0, PrecondGlobalMach);

    // Min and Max Precondition Variable
    MinPrecondSigma = MIN(MinPrecondSigma, beta);
    MaxPrecondSigma = MAX(MaxPrecondSigma, beta);
    
    // Beta is function of M*M
    beta  = beta*beta;
    
    // Set Preconditioner Parameters
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
            error("Compute_Precondition_Characteristic_BoundaryCondition_Turkel: Invalid Solver Precondition Scheme - %d", PrecondMethod);
            break;
    }

    // Compute the Modified Convective and Acoustic Wave speed
    zp       = 0.5*(1.0 + beta - alpha);
    zm       = 0.5*(1.0 - beta + alpha);
    ubar_pre = zp*ubar;
    c_pre    = sqrt(beta*(c*c - ubar*ubar) + zp*zp*ubar*ubar);
    pv1      = rho*(beta - delta);
    pv2      = pv1*ubar;
    pv3      = beta*c*c - alpha*ubar*ubar;
    pv4      = (beta*c*c*(c_pre + ubar*zm))/(2.0*c_pre*pv3);
    pv5      = (beta*c*c*(c_pre - ubar*zm))/(2.0*c_pre*pv3);
    
    // STEP 3:
    // Compute the Eigenvalues of the Precondition Boundary Condition
    lamda1 = ubar;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = ubar_pre + c_pre;
    lamda5 = ubar_pre - c_pre;
    
    // STEP 4:
    // Compute the Right Eigenvector
    BC_A[0][0] =  nx;
    BC_A[0][1] =  nx*ubar*(alpha*(v*nz - w*ny) + pv1*nx)/pv3;
    BC_A[0][2] =  ( beta*c*c*nz - alpha*ubar*(w*ny*ny + u*nx*nz + w*nz*nz) + pv2*nx*ny)/pv3;
    BC_A[0][3] =  (-beta*c*c*ny + alpha*ubar*(v*ny*ny + u*nx*ny + v*nz*nz) + pv2*nx*nz)/pv3;
    BC_A[0][4] =  (-rho*delta*c*c*nx + alpha*c*c*(v*nz - w*ny) + alpha*rho*ubar*ubar*nx)/(rho*c*c*pv3);

    BC_A[1][0] = ny;
    BC_A[1][1] = (-beta*c*c*nz + alpha*ubar*(w*nx*nx + v*ny*nz + w*nz*nz) + pv2*nx*ny)/pv3;
    BC_A[1][2] = ny*ubar*(alpha*(w*nx - u*nz) + pv1*ny)/pv3;
    BC_A[1][3] = ( beta*c*c*nx - alpha*ubar*(u*nx*nx + v*nx*ny + u*nz*nz) + pv2*ny*nz)/pv3;
    BC_A[1][4] = (-rho*delta*c*c*ny + alpha*c*c*(w*nx - u*nz) + alpha*rho*ubar*ubar*ny)/(rho*c*c*pv3);

    BC_A[2][0] = nz;
    BC_A[2][1] = ( beta*c*c*ny - alpha*ubar*(v*nx*nx + w*ny*nz + v*ny*ny) + pv2*nx*nz)/pv3;
    BC_A[2][2] = (-beta*c*c*nx + alpha*ubar*(u*nx*nx + w*nx*nz + u*ny*ny) + pv2*ny*nz)/pv3;
    BC_A[2][3] = nz*ubar*(alpha*(u*ny - v*nx) + pv1*nz)/pv3;
    BC_A[2][4] = (-rho*delta*c*c*nz + alpha*c*c*(u*ny - v*nx) + alpha*rho*ubar*ubar*nz)/(rho*c*c*pv3);

    BC_A[3][0] = 0.0;
    BC_A[3][1] = pv4 * nx;
    BC_A[3][2] = pv4 * ny;
    BC_A[3][3] = pv4 * nz;
    BC_A[3][4] = (beta*c*c + alpha*ubar*(c_pre - ubar*zp))/(2.0*rho*c_pre*pv3);

    BC_A[4][0] = 0.0;
    BC_A[4][1] = -pv5 * nx;
    BC_A[4][2] = -pv5 * ny;
    BC_A[4][3] = -pv5 * nz;
    BC_A[4][4] = (beta*c*c - alpha*ubar*(c_pre + ubar*zp))/(2.0*rho*c_pre*pv3);
        
    // Determine how to set boundary conditions for edge
    switch (BCType) {
        // Euler Solid Wall
        case BC_TYPE_EULER_WALL:
            instance = BC_TYPE_EULER_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Euler Solid Wall boundary condition - 1");
                }
            }
            break;
        // Far Field
        case BC_TYPE_FAR_FIELD:
            instance = BC_TYPE_FAR_FIELD;
            // Supersonic inflow
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW;
            // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW;
            // Subsonic outflow
            else if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW;
            // Subsonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply boundary condition - 2");
            }
            break;
        // Inflow
        case BC_TYPE_INFLOW:
            instance = BC_TYPE_INFLOW;
            // Subsonic inflow
            if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUBSONIC_INFLOW;
            // Supersonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUPERSONIC_INFLOW;
            else {
                // Rescue Mode
                if (fabs(lamda1) < c)
                    instance = BC_TYPE_SUBSONIC_INFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_INFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Rescue Inflow boundary condition edge[%d] - 3", BEdgeID);
            }
            break;
        // Outflow
        case BC_TYPE_OUTFLOW:
            instance = BC_TYPE_OUTFLOW;
            // Subsonic outflow
            if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_SUBSONIC_OUTFLOW;
            // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            else {
                // Rescue Mode
                if (fabs(lamda1) < c)
                    instance = BC_TYPE_SUBSONIC_OUTFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_OUTFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Rescue Outflow boundary condition edge[%d] - 4", BEdgeID);
            }
            break;
        // Subsonic Inflow
        case BC_TYPE_SUBSONIC_INFLOW:
            instance = BC_TYPE_SUBSONIC_INFLOW;
            if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Subsonic Inflow boundary condition edge[%d] - 5", BEdgeID);
            }
            break;
        // Subsonic Outflow
        case BC_TYPE_SUBSONIC_OUTFLOW:
            instance = BC_TYPE_SUBSONIC_OUTFLOW;
            if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_SUBSONIC_OUTFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Subsonic Outflow boundary condition edge[%d] - 6", BEdgeID);
            }
            break;
        // Supersonic Inflow
        case BC_TYPE_SUPERSONIC_INFLOW:
            instance = BC_TYPE_SUPERSONIC_INFLOW;
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUPERSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Supersonic Inflow boundary condition edge[%d] - 7", BEdgeID);
            }
            break;
        // Supersonic Outflow
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Supersonic Outflow boundary condition edge[%d] - 8", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply boundary condition - 9");
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_TYPE_EULER_WALL:
            if (BCMethod == BC_METHOD_CHARCTERISTIC_WEAK) {
                // Compute the Normal Component of Velocity and Subtract From Velocity
                u_b_n   = u_i - nx * ubar_i;
                v_b_n   = v_i - ny * ubar_i;
                w_b_n   = w_i - nz * ubar_i;
                p_b_n   = p_i;
                rho_b_n = rho_i;
            } else {
                BC_L[0] = rho_i;
                BC_L[1] = u_i;
                BC_L[2] = v_i;
                BC_L[3] = w_i;
                BC_L[4] = p_i;   // Note Perturbation input

                // Compute the Characteristic variable in BC_B
                MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
                // Note for grid speed at negative of Grid speed
                BC_B[4] = 0.0;

                // Ignoring the "u_pre - c_pre" characteristics
                BC_A[4][0] = 0.0;
                BC_A[4][1] = nx;
                BC_A[4][2] = ny;
                BC_A[4][3] = nz;
                BC_A[4][4] = 0.0;

                // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
                GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
                // Ly = Pb
                Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
                // Pivot Ux = y
                Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);

                // Solution
                rho_b_n = BC_X[0];
                u_b_n   = BC_X[1];
                v_b_n   = BC_X[2];
                w_b_n   = BC_X[3];
                p_b_n   = BC_X[4];
                // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
                // This causes the flow to be allowed through the surface initially
                if (Iteration < ZPGIteration) {
                    if ((p_b_n + Gauge_Pressure) < 0.0) {
                        p_b_n   = p_i;
                        rho_b_n = rho_i;
                        u_b_n   = u_i - nx * ubar_i;
                        v_b_n   = v_i - ny * ubar_i;
                        w_b_n   = w_i - nz * ubar_i;
                    }
                }
            }
            break;
        case BC_TYPE_SUBSONIC_INFLOW:
            // Compute the Infinity Characteristic variable
            BC_L[0] = Inf_Rho;
            BC_L[1] = Inf_U;
            BC_L[2] = Inf_V;
            BC_L[3] = Inf_W;
            BC_L[4] = p_i; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Compute the Internal Characteristic Variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_X);
            
            // Update the B in forth Characteristic
            BC_B[3] = BC_X[3];
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_SUBSONIC_OUTFLOW:
            // Procedure validated in Mathematica
            // Compute the Internal Characteristic variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i;   // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Now Modify the last row of right eigenvector matrix to set pressure
            BC_A[4][0] = 0.0;
            BC_A[4][1] = 0.0;
            BC_A[4][2] = 0.0;
            BC_A[4][3] = 0.0;
            BC_A[4][4] = 1.0;
            
            // Set the last Characteristic to Outflow Pressure
            BC_B[4] = Outflow_Pressure;
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure; // Note Perturbation input
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;   // Note Perturbation input
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW:
            // Compute the Infinity Characteristic variable
            BC_L[0] = Inf_Rho;
            BC_L[1] = Inf_U;
            BC_L[2] = Inf_V;
            BC_L[3] = Inf_W;
            BC_L[4] = Inf_Pressure; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Compute the Internal Characteristic Variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i;          // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_X);
            
            // Update the B in forth Characteristic
            BC_B[3] = BC_X[3];
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW:
            // Compute the Internal Characteristic variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i;   // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Compute the Infinity Characteristic Variable
            BC_L[0] = Inf_Rho;
            BC_L[1] = Inf_U;
            BC_L[2] = Inf_V;
            BC_L[3] = Inf_W;
            BC_L[4] = Inf_Pressure; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_X);
            
            // Update the B in last Characteristic
            BC_B[4] = BC_X[4];
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure; // Note Perturbation input
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;   // Note Perturbation input
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD:
            p_b_n   = Inf_Pressure; // Note Perturbation input
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Precondition_Characteristic_BoundaryCondition: Unable to set BC Type");
            break;
    }

    if ((p_b_n + Gauge_Pressure) < 0.0) {
        warn("Compute_Precondition_Characteristic_BoundaryCondition: Negative Pressure Detected P_i: %e, P_b_n: %e", p_i, p_b_n);
        rvalue = 1;
    }

    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CON) {
        // Get the Total Energy
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        et_b_n = Material_Get_TotalEnergy(Q_B);
        
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = rho_b_n * u_b_n;
        Q_B[2] = rho_b_n * v_b_n;
        Q_B[3] = rho_b_n * w_b_n;
        Q_B[4] = rho_b_n * et_b_n;
        VariableType = ivVarType_Old;
    }
    
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_RUP) {
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
    }
    
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PUT) {
        // Get the Temperature
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        Q_B[4] = Material_Get_Temperature(Q_B);
        // Update the Q's for Ghost Node
        Q_B[0] = p_b_n; // Only Perturbation 
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        VariableType = ivVarType_Old;
    }
    
    // Primitive Variable Formulation Density Velocity Temperature
    if (VariableType == VARIABLE_RUT) {
        // Get the Temperature
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        Q_B[4] = Material_Get_Temperature(Q_B);
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        VariableType = ivVarType_Old;
    }
    
    return rvalue;
}

//------------------------------------------------------------------------------
//! Characteristic Based Boundary Condition: Numerically Solved
//! Note: Outward Pointing Boundary Normals
//! rvalue: 0 => Success, 1: Negative Pressure
//------------------------------------------------------------------------------
int Compute_Characteristic_BoundaryCondition(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    int i;
    double nx, ny, nz;
    double   rho,   u,   v,   w,   et,   p,   c,   T,   ubar,   ht,   q2,   mach;
    double rho_i, u_i, v_i, w_i, et_i, p_i, c_i, T_i, ubar_i, ht_i, q2_i, mach_i;
    double rho_b, u_b, v_b, w_b, et_b, p_b, c_b, T_b, ubar_b, ht_b, q2_b, mach_b;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    int instance = BC_TYPE_NONE;
    int rvalue = 0;
    
    // Initialization
    BC_Reset();
    
    // Default Initialization: Based on Variable Type Pressure will be Perturbation
    p_b_n   = Inf_Pressure;
    rho_b_n = Inf_Rho;
    u_b_n   = Inf_U;
    v_b_n   = Inf_V;
    w_b_n   = Inf_W;
    
    // Get area vector
    AreaVec.normalize();
    nx = AreaVec.vec[0];
    ny = AreaVec.vec[1];
    nz = AreaVec.vec[2];

    // Compute Equation of State
    // Physical node: Based on Variable Type Pressure will be Perturbation
    Material_Get_Face_Properties(Q_L, nx, ny, nz, rho_i, p_i, T_i, u_i, v_i, w_i, q2_i, c_i, mach_i, ubar_i, et_i, ht_i);
    
    // Ghost node: Based on Variable Type Pressure will be Perturbation
    Material_Get_Face_Properties(Q_R, nx, ny, nz, rho_b, p_b, T_b, u_b, v_b, w_b, q2_b, c_b, mach_b, ubar_b, et_b, ht_b);
    
    // Average State: Based on Variable Type Pressure will be Perturbation
    for (i = 0; i < NEQUATIONS; i++)
        Q_B[i] = 0.5*(Q_L[i] + Q_R[i]);
    Material_Get_Face_Properties(Q_B, nx, ny, nz, rho, p, T, u, v, w, q2, c, mach, ubar, et, ht);
    
    //==========================================================================
    // Compute Boundary Condition
    //==========================================================================
    // STEP 1:
    // Compute the Eigenvalues of the Boundary Condition
    lamda1 = ubar;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = ubar + c;
    lamda5 = ubar - c;
    
    // STEP 2:
    // Compute the Right Eigenvector
    BC_A[0][0] = nx;
    BC_A[0][1] = 0.0;
    BC_A[0][2] = nz;
    BC_A[0][3] = -ny;
    BC_A[0][4] = -nx/(c * c);

    BC_A[1][0] = ny;
    BC_A[1][1] = -nz;
    BC_A[1][2] = 0.0;
    BC_A[1][3] = nx;
    BC_A[1][4] = -ny/(c * c);

    BC_A[2][0] = nz;
    BC_A[2][1] = ny;
    BC_A[2][2] = -nx;
    BC_A[2][3] = 0.0;
    BC_A[2][4] = -nz/(c * c);

    BC_A[3][0] = 0.0;
    BC_A[3][1] = 0.5 * nx;
    BC_A[3][2] = 0.5 * ny;
    BC_A[3][3] = 0.5 * nz;
    BC_A[3][4] = 1.0/(2.0 * rho * c);
    
    BC_A[4][0] = 0.0;
    BC_A[4][1] = -0.5 * nx;
    BC_A[4][2] = -0.5 * ny;
    BC_A[4][3] = -0.5 * nz;
    BC_A[4][4] = 1.0/(2.0 * rho * c);
        
    // Determine how to set boundary conditions for edge
    switch (BCType) {
        // Euler Solid Wall
        case BC_TYPE_EULER_WALL:
            instance = BC_TYPE_EULER_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Euler Solid Wall boundary condition - 1");
                }
            }
            break;
        // Far Field
        case BC_TYPE_FAR_FIELD:
            instance = BC_TYPE_FAR_FIELD;
            // Supersonic inflow
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW;
            // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW;
            // Subsonic outflow
            else if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW;
            // Subsonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply boundary condition - 2");
            }
            break;
        // Inflow
        case BC_TYPE_INFLOW:
            instance = BC_TYPE_INFLOW;
            // Subsonic inflow
            if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUBSONIC_INFLOW;
            // Supersonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUPERSONIC_INFLOW;
            else {
                // Rescue Mode
                if (fabs(lamda1) < c)
                    instance = BC_TYPE_SUBSONIC_INFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_INFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Rescue Inflow boundary condition edge[%d] - 3", BEdgeID);
            }
            break;
        // Outflow
        case BC_TYPE_OUTFLOW:
            instance = BC_TYPE_OUTFLOW;
            // Subsonic outflow
            if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_SUBSONIC_OUTFLOW;
            // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            else {
                // Rescue Mode
                if (fabs(lamda1) < c)
                    instance = BC_TYPE_SUBSONIC_OUTFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_OUTFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Rescue Outflow boundary condition edge[%d] - 4", BEdgeID);
            }
            break;
        // Subsonic Inflow
        case BC_TYPE_SUBSONIC_INFLOW:
            instance = BC_TYPE_SUBSONIC_INFLOW;
            if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Subsonic Inflow boundary condition edge[%d] - 5", BEdgeID);
            }
            break;
        // Subsonic Outflow
        case BC_TYPE_SUBSONIC_OUTFLOW:
            instance = BC_TYPE_SUBSONIC_OUTFLOW;
            if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_SUBSONIC_OUTFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Subsonic Outflow boundary condition edge[%d] - 6", BEdgeID);
            }
            break;
        // Supersonic Inflow
        case BC_TYPE_SUPERSONIC_INFLOW:
            instance = BC_TYPE_SUPERSONIC_INFLOW;
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUPERSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Supersonic Inflow boundary condition edge[%d] - 7", BEdgeID);
            }
            break;
        // Supersonic Outflow
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Supersonic Outflow boundary condition edge[%d] - 8", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply boundary condition - 9");
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_TYPE_EULER_WALL:
            if (BCMethod == BC_METHOD_CHARCTERISTIC_WEAK) {
                // Compute the Normal Component of Velocity and Subtract From Velocity
                p_b_n   = p_i;
                rho_b_n = rho_i;
                u_b_n   = u_i - nx * ubar_i;
                v_b_n   = v_i - ny * ubar_i;
                w_b_n   = w_i - nz * ubar_i;
            } else {
                BC_L[0] = rho_i;
                BC_L[1] = u_i;
                BC_L[2] = v_i;
                BC_L[3] = w_i;
                BC_L[4] = p_i; // Note Perturbation input

                // Compute the Characteristic variable in BC_B
                MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
                // Note for grid speed at negative of Grid speed
                BC_B[4] = 0.0;

                // Ignoring the "u_pre - c_pre" characteristics
                BC_A[4][0] = 0.0;
                BC_A[4][1] = nx;
                BC_A[4][2] = ny;
                BC_A[4][3] = nz;
                BC_A[4][4] = 0.0;

                // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
                GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
                // Ly = Pb
                Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
                // Pivot Ux = y
                Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);

                // Solution
                rho_b_n = BC_X[0];
                u_b_n   = BC_X[1];
                v_b_n   = BC_X[2];
                w_b_n   = BC_X[3];
                p_b_n   = BC_X[4];
                // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
                // This causes the flow to be allowed through the surface initially
                if (Iteration < ZPGIteration) {
                    if ((p_b_n + Gauge_Pressure) < 0.0) {
                        p_b_n   = p_i;
                        rho_b_n = rho_i;
                        u_b_n   = u_i - nx * ubar_i;
                        v_b_n   = v_i - ny * ubar_i;
                        w_b_n   = w_i - nz * ubar_i;
                    }
                }
            }
            break;
        case BC_TYPE_SUBSONIC_INFLOW:
            // Compute the Infinity Characteristic variable
            BC_L[0] = Inf_Rho;
            BC_L[1] = Inf_U;
            BC_L[2] = Inf_V;
            BC_L[3] = Inf_W;
            BC_L[4] = p_i; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Compute the Internal Characteristic Variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_X);
            
            // Update the B in forth Characteristic
            BC_B[3] = BC_X[3];
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_SUBSONIC_OUTFLOW:
            // Procedure validated in Mathematica
            // Compute the Internal Characteristic variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i;   // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Now Modify the last row of right eigenvector matrix to set pressure
            BC_A[4][0] = 0.0;
            BC_A[4][1] = 0.0;
            BC_A[4][2] = 0.0;
            BC_A[4][3] = 0.0;
            BC_A[4][4] = 1.0;
            
            // Set the last Characteristic to Outflow Pressure
            BC_B[4] = Outflow_Pressure;
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure; // Note Perturbation input
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;   // Note Perturbation input
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW:
            // Compute the Infinity Characteristic variable
            BC_L[0] = Inf_Rho;
            BC_L[1] = Inf_U;
            BC_L[2] = Inf_V;
            BC_L[3] = Inf_W;
            BC_L[4] = Inf_Pressure; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Compute the Internal Characteristic Variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_X);
            
            // Update the B in forth Characteristic
            BC_B[3] = BC_X[3];
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW:
            // Compute the Internal Characteristic variable
            BC_L[0] = rho_i;
            BC_L[1] = u_i;
            BC_L[2] = v_i;
            BC_L[3] = w_i;
            BC_L[4] = p_i;   // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_B);
            
            // Compute the Infinity Characteristic Variable
            BC_L[0] = Inf_Rho;
            BC_L[1] = Inf_U;
            BC_L[2] = Inf_V;
            BC_L[3] = Inf_W;
            BC_L[4] = Inf_Pressure; // Note Perturbation input
            MC_Matrix_Mul_Vector(NEQUATIONS, NEQUATIONS, BC_A, BC_L, BC_X);
            
            // Update the B in last Characteristic
            BC_B[4] = BC_X[4];
            
            // Pivot Gaussain LU Decompose BC_A : Note BC_A is replace with LU
            GaussainLUDecompose(BC_A, BC_P, NEQUATIONS, 1);
            // Ly = Pb
            Solve_PivotForwardSubstitution(BC_A, BC_L, BC_B, BC_P, NEQUATIONS);
            // Pivot Ux = y
            Solve_PivotBackSubstitution(BC_A, BC_X, BC_L, BC_P, NEQUATIONS);
            
            // Solution
            rho_b_n = BC_X[0];
            u_b_n   = BC_X[1];
            v_b_n   = BC_X[2];
            w_b_n   = BC_X[3];
            p_b_n   = BC_X[4];
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure; // Note Perturbation input
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;   // Note Perturbation input
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD:
            p_b_n   = Inf_Pressure; // Note Perturbation input
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Precondition_Characteristic_BoundaryCondition: Unable to set BC Type");
            break;
    }
    
    if ((p_b_n + Gauge_Pressure) < 0.0) {
        warn("Compute_Precondition_Characteristic_BoundaryCondition: Negative Pressure Detected P_i: %e, P_b_n: %e", p_i, p_b_n);
        rvalue = 1;
    }

    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CON) {
        // Get the Total Energy
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation
        et_b_n = Material_Get_TotalEnergy(Q_B);
        
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = rho_b_n * u_b_n;
        Q_B[2] = rho_b_n * v_b_n;
        Q_B[3] = rho_b_n * w_b_n;
        Q_B[4] = rho_b_n * et_b_n;
        VariableType = ivVarType_Old;
    }
    
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_RUP) {
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n; // Only Perturbation 
    }
    
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PUT) {
        // Get the Temperature
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        Q_B[4] = Material_Get_Temperature(Q_B);
        // Update the Q's for Ghost Node
        Q_B[0] = p_b_n; // Only Perturbation 
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        VariableType = ivVarType_Old;
    }
    
    // Primitive Variable Formulation Density Velocity Temperature
    if (VariableType == VARIABLE_RUT) {
        // Get the Temperature
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        Q_B[4] = Material_Get_Temperature(Q_B);
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        VariableType = ivVarType_Old;
    }
    
    return rvalue;
}

//------------------------------------------------------------------------------
//! Characteristic Based Boundary Condition: Solved Directly
//! Note: Outward Pointing Boundary Normals
//! rvalue: 0 => Success, 1: Negative Pressure
//------------------------------------------------------------------------------
int Compute_Characteristic_BoundaryCondition_Direct(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    int i;
    double nx, ny, nz;
    double   rho,   u,   v,   w,   et,   p,   c,   T,   ubar,   ht,   q2,   mach;
    double rho_i, u_i, v_i, w_i, et_i, p_i, c_i, T_i, ubar_i, ht_i, q2_i, mach_i;
    double rho_b, u_b, v_b, w_b, et_b, p_b, c_b, T_b, ubar_b, ht_b, q2_b, mach_b;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    int instance = BC_TYPE_NONE;
    int rvalue = 0;
    
    // Default Initialization: Based on Variable Type Pressure will be Perturbation
    p_b_n   = Inf_Pressure;
    rho_b_n = Inf_Rho;
    u_b_n   = Inf_U;
    v_b_n   = Inf_V;
    w_b_n   = Inf_W;
    
    // Get area vector
    AreaVec.normalize();
    nx = AreaVec.vec[0];
    ny = AreaVec.vec[1];
    nz = AreaVec.vec[2];

    // Compute Equation of State
    // Physical node: Based on Variable Type Pressure will be Perturbation
    Material_Get_Face_Properties(Q_L, nx, ny, nz, rho_i, p_i, T_i, u_i, v_i, w_i, q2_i, c_i, mach_i, ubar_i, et_i, ht_i);
    
    // Ghost node: Based on Variable Type Pressure will be Perturbation
    Material_Get_Face_Properties(Q_R, nx, ny, nz, rho_b, p_b, T_b, u_b, v_b, w_b, q2_b, c_b, mach_b, ubar_b, et_b, ht_b);
    
    // Average State: Based on Variable Type Pressure will be Perturbation
    for (i = 0; i < NEQUATIONS; i++)
        Q_B[i] = 0.5*(Q_L[i] + Q_R[i]);
    Material_Get_Face_Properties(Q_B, nx, ny, nz, rho, p, T, u, v, w, q2, c, mach, ubar, et, ht);
    
    // Compute Eigenvalues
    lamda1 = ubar;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = lamda1 + c;
    lamda5 = lamda1 - c;

    // Determine how to set boundary conditions for edge
    switch (BCType) {
        // Euler Solid Wall
        case BC_TYPE_EULER_WALL:
            instance = BC_TYPE_EULER_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Compute_Characteristic_BoundaryCondition: Unable apply Euler Solid Wall boundary condition - 1");
                }
            }
            break;
        // Far Field
        case BC_TYPE_FAR_FIELD:
            instance = BC_TYPE_FAR_FIELD;
            // Supersonic inflow
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW;
            // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW;
            // Subsonic outflow
            else if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW;
            // Subsonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Compute_Characteristic_BoundaryCondition: Unable apply boundary condition - 2");
            }
            break;
        // Inflow
        case BC_TYPE_INFLOW:
            instance = BC_TYPE_INFLOW;
            // Subsonic inflow
            if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUBSONIC_INFLOW;
            // Supersonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUPERSONIC_INFLOW;
            else {
                // Rescue Mode
                if (fabs(lamda1) < c)
                    instance = BC_TYPE_SUBSONIC_INFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_INFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Rescue Inflow boundary condition edge[%d] - 3", BEdgeID);
            }
            break;
        // Outflow
        case BC_TYPE_OUTFLOW:
            instance = BC_TYPE_OUTFLOW;
            // Subsonic outflow
            if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_SUBSONIC_OUTFLOW;
            // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            else {
                // Rescue Mode
                if (fabs(lamda1) < c)
                    instance = BC_TYPE_SUBSONIC_OUTFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_OUTFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Rescue Outflow boundary condition edge[%d] - 4", BEdgeID);
            }
            break;
        // Subsonic Inflow
        case BC_TYPE_SUBSONIC_INFLOW:
            instance = BC_TYPE_SUBSONIC_INFLOW;
            if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Unable apply Subsonic Inflow boundary condition edge[%d] - 5", BEdgeID);
            }
            break;
        // Subsonic Outflow
        case BC_TYPE_SUBSONIC_OUTFLOW:
            instance = BC_TYPE_SUBSONIC_OUTFLOW;
            if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_SUBSONIC_OUTFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Unable apply Subsonic Outflow boundary condition edge[%d] - 6", BEdgeID);
            }
            break;
        // Supersonic Inflow
        case BC_TYPE_SUPERSONIC_INFLOW:
            instance = BC_TYPE_SUPERSONIC_INFLOW;
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUPERSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Unable apply Supersonic Inflow boundary condition edge[%d] - 7", BEdgeID);
            }
            break;
        // Supersonic Outflow
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_TYPE_SUPERSONIC_OUTFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Unable apply Supersonic Outflow boundary condition edge[%d] - 8", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Characteristic_BoundaryCondition: Unable apply boundary condition - 9");
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_TYPE_EULER_WALL:
            if (BCMethod == BC_METHOD_CHARCTERISTIC_WEAK) {
                // Compute the Normal Component of Velocity and Subtract From Velocity
                p_b_n   = p_i;
                rho_b_n = rho_i;
                u_b_n   = u_i - nx * ubar_i;
                v_b_n   = v_i - ny * ubar_i;
                w_b_n   = w_i - nz * ubar_i;
            } else {
                p_b_n   = p_i + rho * c * ubar_i;
                rho_b_n = rho_i + (p_b_n - p_i) / (c * c);
                u_b_n   = u_i - nx * (p_b_n - p_i) / (rho * c);
                v_b_n   = v_i - ny * (p_b_n - p_i) / (rho * c);
                w_b_n   = w_i - nz * (p_b_n - p_i) / (rho * c);
                // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
                // This causes the flow to be allowed through the surface initially
                if (Iteration < ZPGIteration) {
                    if (p_b_n < 0.0){
                        p_b_n   = p_i;
                        rho_b_n = rho_i;
                        u_b_n   = u_i - nx * ubar_i;
                        v_b_n   = v_i - ny * ubar_i;
                        w_b_n   = w_i - nz * ubar_i;
                    }
                }
            }
            break;
        case BC_TYPE_SUBSONIC_INFLOW:
            temp    = nx * (u_i - Inf_U) + ny * (v_i - Inf_V) + nz * (w_i - Inf_W);
            p_b_n   = p_i + 0.5 * rho * c * temp;
            rho_b_n = Inf_Rho + 0.5 * rho * temp / c;
            u_b_n   = Inf_U + 0.5 * nx*temp;
            v_b_n   = Inf_V + 0.5 * ny*temp;
            w_b_n   = Inf_W + 0.5 * nz*temp;
            break;
        case BC_TYPE_SUBSONIC_OUTFLOW:
            p_b_n   = Outflow_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c * c);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho * c);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho * c);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho * c);
            break;
        case BC_TYPE_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW:
            temp    = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
            p_b_n   = 0.5 * (Inf_Pressure + p_i - rho * c * temp);
            rho_b_n = Inf_Rho + (p_b_n - Inf_Pressure) / (c * c);
            u_b_n   = Inf_U + nx * (p_b_n - Inf_Pressure) / (rho * c);
            v_b_n   = Inf_V + ny * (p_b_n - Inf_Pressure) / (rho * c);
            w_b_n   = Inf_W + nz * (p_b_n - Inf_Pressure) / (rho * c);
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW:
            p_b_n   = Inf_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c * c);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho * c);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho * c);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho * c);
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD:
            p_b_n   = Inf_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Characteristic_BoundaryCondition: Unable to set BC Type");
            break;
    }
    
    if ((p_b_n + Gauge_Pressure) < 0.0) {
        warn("Compute_Characteristic_BoundaryCondition: Negative Pressure Detected P_i: %e, P_b_n: %e", p_i, p_b_n);
        rvalue = 1;
    }
    
    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CON) {
        // Get the Total Energy
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        et_b_n = Material_Get_TotalEnergy(Q_B);
        
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = rho_b_n * u_b_n;
        Q_B[2] = rho_b_n * v_b_n;
        Q_B[3] = rho_b_n * w_b_n;
        Q_B[4] = rho_b_n * et_b_n;
        VariableType = ivVarType_Old;
    }
    
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_RUP) {
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation
    }
    
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PUT) {
        // Get the Temperature
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        Q_B[4] = Material_Get_Temperature(Q_B);
        // Update the Q's for Ghost Node
        Q_B[0] = p_b_n; // Only Perturbation 
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        VariableType = ivVarType_Old;
    }
    
    // Primitive Variable Formulation Density Velocity Temperature
    if (VariableType == VARIABLE_RUT) {
        // Get the Temperature
        int ivVarType_Old = VariableType;
        VariableType = VARIABLE_RUP;
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;   // Only Perturbation 
        Q_B[4] = Material_Get_Temperature(Q_B);
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        VariableType = ivVarType_Old;
    }
    
    return rvalue;
}

//------------------------------------------------------------------------------
//! This Routine is created to assist Finite Difference Jacobian Computations
//! Note: Outward Pointing Boundary Normals
//------------------------------------------------------------------------------
void Apply_Characteristic_Boundary_Condition(int BEdgeID, int Iteration) {
    int rvalue;
    int node_L, node_R;
    Vector3D areavec;
    double Q_L[NEQUATIONS];
    double Q_R[NEQUATIONS];
    double Q_B[NEQUATIONS];
    
    node_L = bndEdge[BEdgeID].node[0];
    node_R = bndEdge[BEdgeID].node[1];

#ifdef DEBUG
assert(node_R > node_L);
#endif
    // Left Node - Physical Node
    Q_L[0] = Q1[node_L];
    Q_L[1] = Q2[node_L];
    Q_L[2] = Q3[node_L];
    Q_L[3] = Q4[node_L];
    Q_L[4] = Q5[node_L];

    // Right Node - Ghost Node
    Q_R[0] = Q1[node_R];
    Q_R[1] = Q2[node_R];
    Q_R[2] = Q3[node_R];
    Q_R[3] = Q4[node_R];
    Q_R[4] = Q5[node_R];

    // Get the Area Vector
    areavec = bndEdge[BEdgeID].areav;
    
    // Compute the Roe Flux for this edge
    rvalue = 0;
    switch (PrecondMethod) {
        case PRECOND_METHOD_NONE: // Roe
            rvalue = Compute_Characteristic_BoundaryCondition(Q_L, Q_R, areavec, BEdgeID, bndEdge[BEdgeID].type, Iteration, Q_B);
            break;
        case PRECOND_METHOD_LMFIX: // LMRoe
            rvalue = Compute_Characteristic_BoundaryCondition(Q_L, Q_R, areavec, BEdgeID, bndEdge[BEdgeID].type, Iteration, Q_B);
            break;
        case PRECOND_METHOD_THORNBER: // THORNBER
            rvalue = Compute_Characteristic_BoundaryCondition(Q_L, Q_R, areavec, BEdgeID, bndEdge[BEdgeID].type, Iteration, Q_B);
            break;
        case PRECOND_METHOD_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
            rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Turkel(Q_L, Q_R, areavec, BEdgeID, bndEdge[BEdgeID].type, Iteration, Q_B);
            break;
        case PRECOND_METHOD_ERIKSSON: // Roe Eriksson Pre-Conditioner
            rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Turkel(Q_L, Q_R, areavec, BEdgeID, bndEdge[BEdgeID].type, Iteration, Q_B);
            break;
        case PRECOND_METHOD_MERKEL: // Roe Merkel Pre-Conditioner
            rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Merkel(Q_L, Q_R, areavec, BEdgeID, bndEdge[BEdgeID].type, Iteration, Q_B);
            break;
        case PRECOND_METHOD_TURKEL: // Roe Turkel Pre-Conditioner
            rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Turkel(Q_L, Q_R, areavec, BEdgeID, bndEdge[BEdgeID].type, Iteration, Q_B);
            break;
        default:
            error("Apply_Characteristic_Boundary_Condition: Invalid Solver Precondition Scheme - %d -1", PrecondMethod);
            break;
    }
    
    if (rvalue != 0)
        error("Apply_Characteristic_Boundary_Condition: BndNode: %d,  GhostNode: %d -2", node_L, node_R);

    // Update the Boundary Value of Ghost Node
    Q1[node_R] = Q_B[0];
    Q2[node_R] = Q_B[1];
    Q3[node_R] = Q_B[2];
    Q4[node_R] = Q_B[3];
    Q5[node_R] = Q_B[4];
}

//------------------------------------------------------------------------------
//! Note: Outward Pointing Boundary Normals
//------------------------------------------------------------------------------
void Apply_Characteristic_Boundary_Condition(int Iteration) {
    int iBEdge, rvalue;
    int node_L, node_R;
    Vector3D areavec;
    double Q_L[NEQUATIONS];
    double Q_R[NEQUATIONS];
    double Q_B[NEQUATIONS];
    
    for (iBEdge = 0; iBEdge < nBEdge; iBEdge++) {
        node_L = bndEdge[iBEdge].node[0];
        node_R = bndEdge[iBEdge].node[1];
        
#ifdef DEBUG
    assert(node_R > node_L);
#endif
        // Left Node - Physical Node
        Q_L[0] = Q1[node_L];
        Q_L[1] = Q2[node_L];
        Q_L[2] = Q3[node_L];
        Q_L[3] = Q4[node_L];
        Q_L[4] = Q5[node_L];

        // Right Node - Ghost Node
        Q_R[0] = Q1[node_R];
        Q_R[1] = Q2[node_R];
        Q_R[2] = Q3[node_R];
        Q_R[3] = Q4[node_R];
        Q_R[4] = Q5[node_R];

        // Get the Area Vector
        areavec = bndEdge[iBEdge].areav;
        
        // Compute the Roe Flux for this edge
        rvalue = 0;
        switch (PrecondMethod) {
            case PRECOND_METHOD_NONE: // Roe
                rvalue = Compute_Characteristic_BoundaryCondition(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
                break;
            case PRECOND_METHOD_LMFIX: // LMRoe
                rvalue = Compute_Characteristic_BoundaryCondition(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
                break;
            case PRECOND_METHOD_THORNBER: // THORNBER
                rvalue = Compute_Characteristic_BoundaryCondition(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
                break;
            case PRECOND_METHOD_BTW: // Roe Briley Taylor Whitfield Pre-Conditioner
                rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Turkel(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
                break;
            case PRECOND_METHOD_ERIKSSON: // Roe Eriksson Pre-Conditioner
                rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Turkel(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
                break;
            case PRECOND_METHOD_MERKEL: // Roe Merkel Pre-Conditioner
                rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Merkel(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
                break;
            case PRECOND_METHOD_TURKEL: // Roe Turkel Pre-Conditioner
                rvalue = Compute_Precondition_Characteristic_BoundaryCondition_Turkel(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
                break;
            default:
                error("Apply_Characteristic_Boundary_Condition: Invalid Solver Precondition Scheme - %d -1", PrecondMethod);
                break;
        }
        
        if (rvalue != 0)
            error("Apply_Characteristic_Boundary_Condition: BndNode: %d,  GhostNode: %d -2", node_L, node_R);
        
        // Update the Boundary Value of Ghost Node
        Q1[node_R] = Q_B[0];
        Q2[node_R] = Q_B[1];
        Q3[node_R] = Q_B[2];
        Q4[node_R] = Q_B[3];
        Q5[node_R] = Q_B[4];
    }
}

//------------------------------------------------------------------------------
//! This Routine Apply Boundary Condition to Specified Edge
//------------------------------------------------------------------------------
void Apply_Boundary_Condition(int BEdgeID, int Iteration) {
    // Apply Characteristic Based Boundary Condition
    Apply_Characteristic_Boundary_Condition(BEdgeID, Iteration);
}

//------------------------------------------------------------------------------
//! This Routine Apply Boundary Condition to All Boundary Edges
//------------------------------------------------------------------------------
void Apply_Boundary_Condition(int Iteration) {
    // Apply Characteristic Based Boundary Condition
    Apply_Characteristic_Boundary_Condition(Iteration);
}

