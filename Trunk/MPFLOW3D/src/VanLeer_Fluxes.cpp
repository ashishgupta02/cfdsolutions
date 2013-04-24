/*******************************************************************************
 * File:        VanLeer_Fluxes.cpp
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
static int VanLeer_DB               = 0;
static double *VanLeer_fluxA        = NULL;
static double *VanLeer_flux_L       = NULL;
static double *VanLeer_flux_R       = NULL;
static double *VanLeer_Q_L          = NULL;
static double *VanLeer_Q_R          = NULL;
static double *VanLeer_dw           = NULL;
static double *VanLeer_dQ           = NULL;
static double **VanLeer_A           = NULL;
static double **VanLeer_Eigen       = NULL;
static double **VanLeer_M           = NULL;
static double **VanLeer_Minv        = NULL;
static double **VanLeer_P           = NULL;
static double **VanLeer_Pinv        = NULL;
static double **VanLeer_T           = NULL;
static double **VanLeer_Tinv        = NULL;

//------------------------------------------------------------------------------
//! Create VanLeer Scheme Data Structure
//------------------------------------------------------------------------------
void VanLeer_Init(void) {
    int i;
    
    // Check if VanLeer Data Structure is required
    if (VanLeer_DB == 0) {
        VanLeer_fluxA  = (double *)  malloc(NEQUATIONS*sizeof(double));
        VanLeer_flux_L = (double *)  malloc(NEQUATIONS*sizeof(double));
        VanLeer_flux_R = (double *)  malloc(NEQUATIONS*sizeof(double));
        VanLeer_Q_L    = (double *)  malloc(NEQUATIONS*sizeof(double));
        VanLeer_Q_R    = (double *)  malloc(NEQUATIONS*sizeof(double));
        VanLeer_dw     = (double *)  malloc(NEQUATIONS*sizeof(double));
        VanLeer_dQ     = (double *)  malloc(NEQUATIONS*sizeof(double));
        VanLeer_A      = (double **) malloc(NEQUATIONS*sizeof(double*));
        VanLeer_Eigen  = (double **) malloc(NEQUATIONS*sizeof(double*));
        VanLeer_M      = (double **) malloc(NEQUATIONS*sizeof(double*));
        VanLeer_Minv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        VanLeer_P      = (double **) malloc(NEQUATIONS*sizeof(double*));
        VanLeer_Pinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        VanLeer_T      = (double **) malloc(NEQUATIONS*sizeof(double*));
        VanLeer_Tinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        for (i = 0; i < NEQUATIONS; i++) {
            VanLeer_A[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            VanLeer_Eigen[i] = (double *) malloc(NEQUATIONS*sizeof(double));
            VanLeer_M[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            VanLeer_Minv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
            VanLeer_P[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            VanLeer_Pinv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
            VanLeer_T[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            VanLeer_Tinv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
        }
        VanLeer_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete VanLeer Scheme Data Structure
//------------------------------------------------------------------------------
void VanLeer_Finalize(void) {
    int i;
    
    // Free the Memory
    for (i = 0; i < NEQUATIONS; i++) {
        free(VanLeer_A[i]);
        free(VanLeer_Eigen[i]);
        free(VanLeer_M[i]);
        free(VanLeer_Minv[i]);
        free(VanLeer_P[i]);
        free(VanLeer_Pinv[i]);
        free(VanLeer_T[i]);
        free(VanLeer_Tinv[i]);
    }
    free(VanLeer_fluxA);
    free(VanLeer_flux_L);
    free(VanLeer_flux_R);
    free(VanLeer_Q_L);
    free(VanLeer_Q_R);
    free(VanLeer_dw);
    free(VanLeer_dQ);
    free(VanLeer_A);
    free(VanLeer_Eigen);
    free(VanLeer_M);
    free(VanLeer_Minv);
    free(VanLeer_P);
    free(VanLeer_Pinv);
    free(VanLeer_T);
    free(VanLeer_Tinv);
    VanLeer_DB = 0;
}

//------------------------------------------------------------------------------
//! Reset VanLeer Scheme Data Structure
//------------------------------------------------------------------------------
void VanLeer_Reset(void) {
    int i, j;
    
    if (VanLeer_DB == 0)
        VanLeer_Init();
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++) {
        VanLeer_fluxA[i]  = 0.0;
        VanLeer_flux_L[i] = 0.0;
        VanLeer_flux_R[i] = 0.0;
        VanLeer_Q_L[i]    = 0.0;
        VanLeer_Q_R[i]    = 0.0;
        VanLeer_dw[i]     = 0.0;
        VanLeer_dQ[i]     = 0.0;
        for (j = 0; j < NEQUATIONS; j++) {
            VanLeer_A[i][j]     = 0.0;
            VanLeer_Eigen[i][j] = 0.0;
            VanLeer_M[i][j]     = 0.0;
            VanLeer_Minv[i][j]  = 0.0;
            VanLeer_P[i][j]     = 0.0;
            VanLeer_Pinv[i][j]  = 0.0;
            VanLeer_T[i][j]     = 0.0;
            VanLeer_Tinv[i][j]  = 0.0;
        }
    }
}

//------------------------------------------------------------------------------
//! Compute VanLeer Flux
//------------------------------------------------------------------------------
void Compute_Flux_VanLeer(int node_L, int node_R, Vector3D areavec, double *Flux_VanLeer, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    
    double area, maxlambda;
    double nx, ny, nz;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_VanLeer[i] = 0.0;
    VanLeer_Reset();
    
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
                Compute_SecondOrderReconstructQ(node_L, node_R, VanLeer_Q_L, VanLeer_Q_R);
            } else {
                VanLeer_Q_L[0] = Q1[node_L];
                VanLeer_Q_L[1] = Q2[node_L];
                VanLeer_Q_L[2] = Q3[node_L];
                VanLeer_Q_L[3] = Q4[node_L];
                VanLeer_Q_L[4] = Q5[node_L];
                VanLeer_Q_R[0] = Q1[node_R];
                VanLeer_Q_R[1] = Q2[node_R];
                VanLeer_Q_R[2] = Q3[node_R];
                VanLeer_Q_R[3] = Q4[node_R];
                VanLeer_Q_R[4] = Q5[node_R];
            }
        } else { // Boundary and Ghost Node
            // Make Solution Second Order
            // Note: Boundary Residual cannot be made second order
            // because node_R is ghost node with no physical coordinates value
            // Hence boundary residual always remains first order
            // Ghost Edge always remains first order
            VanLeer_Q_L[0] = Q1[node_L];
            VanLeer_Q_L[1] = Q2[node_L];
            VanLeer_Q_L[2] = Q3[node_L];
            VanLeer_Q_L[3] = Q4[node_L];
            VanLeer_Q_L[4] = Q5[node_L];
            VanLeer_Q_R[0] = Q1[node_R];
            VanLeer_Q_R[1] = Q2[node_R];
            VanLeer_Q_R[2] = Q3[node_R];
            VanLeer_Q_R[3] = Q4[node_R];
            VanLeer_Q_R[4] = Q5[node_R];
        }
        
        // Compute Equation of State
        Material_Get_Face_Properties(VanLeer_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Material_Get_Face_Properties(VanLeer_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        p_L += Gauge_Pressure;
        p_R += Gauge_Pressure;
        
        double rho_avg, ubar_avg, p_avg, c_avg;
        rho_avg  = 0.5*(rho_L + rho_R);
        ubar_avg = 0.5*((u_L + u_R)*nx + (v_L + v_R)*ny + (w_L + w_R)*nz);
        p_avg    = 0.5*(p_L + p_R);
        // Get the Speed of Sound
        int ivVarType_Old = VariableType;
        double daQ[NEQUATIONS];
        daQ[0] = rho_avg;
        daQ[1] = 0.5*(u_L + u_R);
        daQ[2] = 0.5*(v_L + v_R);
        daQ[3] = 0.5*(w_L + w_R);
        daQ[4] = p_avg - Gauge_Pressure;
        VariableType = VARIABLE_RUP;
        c_avg = Material_Get_SpeedSound(daQ);
        VariableType = ivVarType_Old;
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(fabs(ubar_avg), MAX(fabs(ubar_avg + c_avg), fabs(ubar_avg - c_avg)));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
            maxlambda = MAX(fabs(ubar_avg), MAX(fabs(ubar_avg + c_avg), fabs(ubar_avg - c_avg)));
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
        
        double fluxp1, fluxp2, fluxp3, fluxp4, fluxp5;
        double fluxm1, fluxm2, fluxm3, fluxm4, fluxm5;
        double fmach;
        double ubl2a, ubr2a;
        
        // Construct Flux Plus
        fmach = ubar_L/c_L;
        ubl2a = -ubar_L + 2.0*c_L;
        
        if (ABS(fmach) < 1.0) {
            fluxp1 = 0.25*rho_L*c_L*(fmach + 1.0)*(fmach + 1.0);
            fluxp2 = fluxp1*(nx*ubl2a/Gamma + u_L);
            fluxp3 = fluxp1*(ny*ubl2a/Gamma + v_L);
            fluxp4 = fluxp1*(nz*ubl2a/Gamma + w_L);
            fluxp5 = fluxp1*((-(Gamma - 1.0)*ubar_L*ubar_L + 2.0*(Gamma - 1.0)*ubar_L*c_L + 2.0*c_L*c_L)/(Gamma*Gamma - 1.0) + 0.5*q2_L);
        } else if (fmach >= 1.0) {
            fluxp1 = rho_L*ubar_L;
            fluxp2 = rho_L*u_L*ubar_L + nx*p_L;
            fluxp3 = rho_L*v_L*ubar_L + ny*p_L;
            fluxp4 = rho_L*w_L*ubar_L + nz*p_L;
            fluxp5 = (rho_L*et_L + p_L)*ubar_L;
        } else {
            fluxp1 = 0.0;
            fluxp2 = 0.0;
            fluxp3 = 0.0;
            fluxp4 = 0.0;
            fluxp5 = 0.0;
        }
        
        // Construct Flux Minus
        fmach = ubar_R/c_R;
        ubr2a = -ubar_R - 2.0*c_R;
        
        if (ABS(fmach) < 1.0) {
            fluxm1 = -0.25*rho_R*c_R*(fmach - 1.0)*(fmach - 1.0);
            fluxm2 = fluxm1*(nx*ubr2a/Gamma + u_R);
            fluxm3 = fluxm1*(ny*ubr2a/Gamma + v_R);
            fluxm4 = fluxm1*(nz*ubr2a/Gamma + w_R);
            fluxm5 = fluxm1*((-(Gamma - 1.0)*ubar_R*ubar_R - 2.0*(Gamma - 1.0)*ubar_R*c_R + 2.0*c_R*c_R)/(Gamma*Gamma - 1.0) + 0.5*q2_R);
        } else if (fmach <= -1.0) {
            fluxm1 = rho_R*ubar_R;
            fluxm2 = rho_R*u_R*ubar_R + nx*p_R;
            fluxm3 = rho_R*v_R*ubar_R + ny*p_R;
            fluxm4 = rho_R*w_R*ubar_R + nz*p_R;
            fluxm5 = (rho_R*et_R + p_R)*ubar_R;
        } else {
            fluxm1 = 0.0;
            fluxm2 = 0.0;
            fluxm3 = 0.0;
            fluxm4 = 0.0;
            fluxm5 = 0.0;
        }
        
        // Compute the VanLeer Flux
        Flux_VanLeer[0] = (fluxp1 + fluxm1)*area;
        Flux_VanLeer[1] = (fluxp2 + fluxm2)*area;
        Flux_VanLeer[2] = (fluxp3 + fluxm3)*area;
        Flux_VanLeer[3] = (fluxp4 + fluxm4)*area;
        Flux_VanLeer[4] = (fluxp5 + fluxm5)*area;
    } else
        error("Compute_Flux_VanLeer: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Computes the VanLeer Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_VanLeer(int AddTime) {
    int i;
    int node_L, node_R;
    Vector3D areavec;
    double flux_vanleer[5];
    
    // Internal Edges
    for (i = 0; i < nEdge; i++) {
        // Get two nodes of edge
        node_L = intEdge[i].node[0];
        node_R = intEdge[i].node[1];
        
        // Get area vector
        areavec = intEdge[i].areav;
        
        // Compute the VanLeer Flux for this edge
        Compute_Flux_VanLeer(node_L, node_R, areavec, flux_vanleer, AddTime);
        
        // Compute for LHS
        Res1_Conv[node_L] += flux_vanleer[0];
        Res2_Conv[node_L] += flux_vanleer[1];
        Res3_Conv[node_L] += flux_vanleer[2];
        Res4_Conv[node_L] += flux_vanleer[3];
        Res5_Conv[node_L] += flux_vanleer[4];

        Res1_Conv[node_R] -= flux_vanleer[0];
        Res2_Conv[node_R] -= flux_vanleer[1];
        Res3_Conv[node_R] -= flux_vanleer[2];
        Res4_Conv[node_R] -= flux_vanleer[3];
        Res5_Conv[node_R] -= flux_vanleer[4];
    }

    // Boundary Edges
    for (i = 0; i < nBEdge; i++) {
        // Get two nodes of edge
        node_L = bndEdge[i].node[0];
        node_R = bndEdge[i].node[1];
        
        // Get area vector
        areavec = bndEdge[i].areav;
        
        // Compute the VanLeer Flux for this edge
        Compute_Flux_VanLeer(node_L, node_R, areavec, flux_vanleer, AddTime);
        
        // Compute for LHS
        Res1_Conv[node_L] += flux_vanleer[0];
        Res2_Conv[node_L] += flux_vanleer[1];
        Res3_Conv[node_L] += flux_vanleer[2];
        Res4_Conv[node_L] += flux_vanleer[3];
        Res5_Conv[node_L] += flux_vanleer[4];
    }
    
    // Precondition and Transform the Residual in Conservative or Primitive Form
    // Compute_Precondition_Residual_VanLeer(); // TODO
}

