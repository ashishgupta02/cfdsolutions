/*******************************************************************************
 * File:        HLLC_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    4
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
static double **HLLC_A           = NULL;
static double **HLLC_Eigen       = NULL;
static double **HLLC_M           = NULL;
static double **HLLC_Minv        = NULL;
static double **HLLC_P           = NULL;
static double **HLLC_Pinv        = NULL;
static double **HLLC_T           = NULL;
static double **HLLC_Tinv        = NULL;

//------------------------------------------------------------------------------
//! Create HLLC Scheme Data Structure
//------------------------------------------------------------------------------
void HLLC_Init(void) {
    int i;
    
    // Check if HLLC Data Structure is required
    if (HLLC_DB == 0) {
        HLLC_fluxA  = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_flux_L = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_flux_R = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_Q_L    = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_Q_R    = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_dw     = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_dQ     = (double *)  malloc(NEQUATIONS*sizeof(double));
        HLLC_A      = (double **) malloc(NEQUATIONS*sizeof(double*));
        HLLC_Eigen  = (double **) malloc(NEQUATIONS*sizeof(double*));
        HLLC_M      = (double **) malloc(NEQUATIONS*sizeof(double*));
        HLLC_Minv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        HLLC_P      = (double **) malloc(NEQUATIONS*sizeof(double*));
        HLLC_Pinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        HLLC_T      = (double **) malloc(NEQUATIONS*sizeof(double*));
        HLLC_Tinv   = (double **) malloc(NEQUATIONS*sizeof(double*));
        for (i = 0; i < NEQUATIONS; i++) {
            HLLC_A[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            HLLC_Eigen[i] = (double *) malloc(NEQUATIONS*sizeof(double));
            HLLC_M[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            HLLC_Minv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
            HLLC_P[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            HLLC_Pinv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
            HLLC_T[i]     = (double *) malloc(NEQUATIONS*sizeof(double));
            HLLC_Tinv[i]  = (double *) malloc(NEQUATIONS*sizeof(double));
        }
        HLLC_DB = 1;
    }
}

//------------------------------------------------------------------------------
//! Delete HLLC Scheme Data Structure
//------------------------------------------------------------------------------
void HLLC_Finalize(void) {
    int i;
    
    // Free the Memory
    for (i = 0; i < NEQUATIONS; i++) {
        free(HLLC_A[i]);
        free(HLLC_Eigen[i]);
        free(HLLC_M[i]);
        free(HLLC_Minv[i]);
        free(HLLC_P[i]);
        free(HLLC_Pinv[i]);
        free(HLLC_T[i]);
        free(HLLC_Tinv[i]);
    }
    free(HLLC_fluxA);
    free(HLLC_flux_L);
    free(HLLC_flux_R);
    free(HLLC_Q_L);
    free(HLLC_Q_R);
    free(HLLC_dw);
    free(HLLC_dQ);
    free(HLLC_A);
    free(HLLC_Eigen);
    free(HLLC_M);
    free(HLLC_Minv);
    free(HLLC_P);
    free(HLLC_Pinv);
    free(HLLC_T);
    free(HLLC_Tinv);
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
        for (j = 0; j < NEQUATIONS; j++) {
            HLLC_A[i][j]     = 0.0;
            HLLC_Eigen[i][j] = 0.0;
            HLLC_M[i][j]     = 0.0;
            HLLC_Minv[i][j]  = 0.0;
            HLLC_P[i][j]     = 0.0;
            HLLC_Pinv[i][j]  = 0.0;
            HLLC_T[i][j]     = 0.0;
            HLLC_Tinv[i][j]  = 0.0;
        }
    }
}

//------------------------------------------------------------------------------
//! Compute HLLC Flux
//------------------------------------------------------------------------------
void Compute_HLLCFlux(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, maxlambda;
    
    double area;
    double nx, ny, nz;
    
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
        Compute_EOS_Variables_Face(HLLC_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(HLLC_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        p_L += Gauge_Pressure;
        p_R += Gauge_Pressure;
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        phi   = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c     = (Gamma - 1.0)*ht - phi;
        c     = sqrt(c);
        ubar  = u*nx + v*ny + w*nz;
        
        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(fabs(ubar), MAX(fabs(ubar + c), fabs(ubar - c)));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
            maxlambda = MAX(fabs(ubar), MAX(fabs(ubar + c), fabs(ubar - c)));
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
        
        // Compute Min and Max Eigenvalues
        double eig5_L = ubar_L - c_L;
        double eig4_R = ubar_R + c_R;
        double eig4   = ubar + c;
        double eig5   = ubar - c;
        
        // Compute SL, SR and SM
        double SL = MIN(eig5_L, eig5);
        double SR = MAX(eig4_R, eig4);
        double SM = (rho_R*ubar_R*(SR - ubar_R) - rho_L*ubar_L*(SL - ubar_L) + p_L - p_R)/
                        (rho_R*(SR - ubar_R) - rho_L*(SL - ubar_L));
        
        double omega_L, rhostar_L, rhoustar_L, rhovstar_L, rhowstar_L, rhoetstar_L;
        double omega_R, rhostar_R, rhoustar_R, rhovstar_R, rhowstar_R, rhoetstar_R;
        double pstar, cnst1, cnst2;
        
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
            pstar   = p_L + rho_L*(ubar_L - SL)*(ubar_L - SM);
            omega_L = 1.0/(SL - SM);
            cnst1   = SL - ubar_L;
            cnst2   = pstar - p_L;
            
            rhostar_L   = omega_L*(cnst1*rho_L);
            rhoustar_L  = omega_L*(cnst1*rho_L*u_L + cnst2*nx);
            rhovstar_L  = omega_L*(cnst1*rho_L*v_L + cnst2*ny);
            rhowstar_L  = omega_L*(cnst1*rho_L*w_L + cnst2*nz);
            rhoetstar_L = omega_L*(cnst1*rho_L*et_L - p_L*ubar_L + pstar*SM);
            
            Flux_HLLC[0] = rhostar_L*SM;
            Flux_HLLC[1] = rhoustar_L*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_L*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_L*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_L + pstar)*SM;        
        } else if ((SM <= 0.0) && (SR > 0.0)) {
            pstar   = p_R + rho_R*(ubar_R - SR)*(ubar_R - SM);
            omega_R = 1.0/(SR - SM);
            cnst1   = SR - ubar_R;
            cnst2   = pstar - p_R;
            
            rhostar_R   = omega_R*(cnst1*rho_R);
            rhoustar_R  = omega_R*(cnst1*rho_R*u_R + cnst2*nx);
            rhovstar_R  = omega_R*(cnst1*rho_R*v_R + cnst2*ny);
            rhowstar_R  = omega_R*(cnst1*rho_R*w_R + cnst2*nz);
            rhoetstar_R = omega_R*(cnst1*rho_R*et_R - p_R*ubar_R + pstar*SM);
            
            Flux_HLLC[0] = rhostar_R*SM;
            Flux_HLLC[1] = rhoustar_R*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_R*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_R*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_R + pstar)*SM;
        } else
            error("Compute_HLLCFlux: Exception Anomaly");
        
        // Compute the HLLC Flux
        Flux_HLLC[0] *= area;
        Flux_HLLC[1] *= area;
        Flux_HLLC[2] *= area;
        Flux_HLLC[3] *= area;
        Flux_HLLC[4] *= area;
    } else
        error("Compute_HLLCFlux: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute HLLC Flux with LMFix
//------------------------------------------------------------------------------
void Compute_HLLCFlux_LMFix(int node_L, int node_R, Vector3D areavec, double *Flux_HLLC, int AddTime) {
    int i;
    double rho_L, u_L, v_L, w_L, et_L, p_L, c_L, ht_L, ubar_L, q2_L, T_L, mach_L;
    double rho_R, u_R, v_R, w_R, et_R, p_R, c_R, ht_R, ubar_R, q2_R, T_R, mach_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, maxlambda;
    
    double area;
    double nx, ny, nz;
    
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
        Compute_EOS_Variables_Face(HLLC_Q_L, nx, ny, nz, rho_L, p_L, T_L, u_L, v_L, w_L, q2_L, c_L, mach_L, ubar_L, et_L, ht_L);
        Compute_EOS_Variables_Face(HLLC_Q_R, nx, ny, nz, rho_R, p_R, T_R, u_R, v_R, w_R, q2_R, c_R, mach_R, ubar_R, et_R, ht_R);
        p_L += Gauge_Pressure;
        p_R += Gauge_Pressure;
        
        // ROE AVERAGE VARIABLES
        rho   = sqrt(rho_R * rho_L);
        sigma = rho/(rho_L + rho);
        u     = u_L  + sigma*(u_R  - u_L);
        v     = v_L  + sigma*(v_R  - v_L);
        w     = w_L  + sigma*(w_R  - w_L);
        ht    = ht_L + sigma*(ht_R - ht_L);
        phi   = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c     = (Gamma - 1.0)*ht - phi;
        c     = sqrt(c);
        ubar  = u*nx + v*ny + w*nz;

        // Add to Time only if Required
        if (AddTime == TRUE) {
            // Time Computation : Add Max Eigenvalue
            maxlambda = MAX(fabs(ubar), MAX(fabs(ubar + c), fabs(ubar - c)));
            DeltaT[node_L] += area*maxlambda;
            if (node_R < nNode)
                DeltaT[node_R] += area*maxlambda;
        }
        
        // Computations for Variable Relaxation Residual Smoothing
        if (ResidualSmoothMethod == RESIDUAL_SMOOTH_METHOD_IMPLICIT) {
            maxlambda = MAX(fabs(ubar), MAX(fabs(ubar + c), fabs(ubar - c)));
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
        
        // Compute Min and Max Eigenvalues
        double eig5_L = ubar_L - c_L;
        double eig4_R = ubar_R + c_R;
        double eig4   = ubar + c;
        double eig5   = ubar - c;
        
        // Compute SL, SR and SM
        double SL = MIN(eig5_L, eig5);
        double SR = MAX(eig4_R, eig4);
        double SM = (rho_R*ubar_R*(SR - ubar_R) - rho_L*ubar_L*(SL - ubar_L) + p_L - p_R)/
                        (rho_R*(SR - ubar_R) - rho_L*(SL - ubar_L));
        
        double omega_L, rhostar_L, rhoustar_L, rhovstar_L, rhowstar_L, rhoetstar_L;
        double omega_R, rhostar_R, rhoustar_R, rhovstar_R, rhowstar_R, rhoetstar_R;
        double pstar, cnst1, cnst2;
        
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
            pstar   = p_L + rho_L*(ubar_L - SL)*(ubar_L - SM);
            omega_L = 1.0/(SL - SM);
            cnst1   = SL - ubar_L;
            cnst2   = pstar - p_L;
            
            rhostar_L   = omega_L*(cnst1*rho_L);
            rhoustar_L  = omega_L*(cnst1*rho_L*u_L + cnst2*nx);
            rhovstar_L  = omega_L*(cnst1*rho_L*v_L + cnst2*ny);
            rhowstar_L  = omega_L*(cnst1*rho_L*w_L + cnst2*nz);
            rhoetstar_L = omega_L*(cnst1*rho_L*et_L  - p_L*ubar_L + pstar*SM);
            
            Flux_HLLC[0] = rhostar_L*SM;
            Flux_HLLC[1] = rhoustar_L*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_L*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_L*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_L + pstar)*SM;        
        } else if ((SM <= 0.0) && (SR > 0.0)) {
            pstar   = p_R + rho_R*(ubar_R - SR)*(ubar_R - SM);
            omega_R = 1.0/(SR - SM);
            cnst1   = SR - ubar_R;
            cnst2   = pstar - p_R;
            
            rhostar_R   = omega_R*(cnst1*rho_R);
            rhoustar_R  = omega_R*(cnst1*rho_R*u_R + cnst2*nx);
            rhovstar_R  = omega_R*(cnst1*rho_R*v_R + cnst2*ny);
            rhowstar_R  = omega_R*(cnst1*rho_R*w_R + cnst2*nz);
            rhoetstar_R = omega_R*(cnst1*rho_R*et_R  - p_R*ubar_R + pstar*SM);
            
            Flux_HLLC[0] = rhostar_R*SM;
            Flux_HLLC[1] = rhoustar_R*SM + pstar*nx;
            Flux_HLLC[2] = rhovstar_R*SM + pstar*ny;
            Flux_HLLC[3] = rhowstar_R*SM + pstar*nz;
            Flux_HLLC[4] = (rhoetstar_R + pstar)*SM;
        } else
            error("Compute_HLLCFlux: Exception Anomaly");
        
        // Compute the HLLC Flux
        Flux_HLLC[0] *= area;
        Flux_HLLC[1] *= area;
        Flux_HLLC[2] *= area;
        Flux_HLLC[3] *= area;
        Flux_HLLC[4] *= area;
    } else
        error("Compute_HLLCFlux: Invalid Node - %d", node_L);
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
        Compute_HLLCFlux(node_L, node_R, areavec, flux_hllc, AddTime);
        
        // Compute for LHS
        Res1[node_L] += flux_hllc[0];
        Res2[node_L] += flux_hllc[1];
        Res3[node_L] += flux_hllc[2];
        Res4[node_L] += flux_hllc[3];
        Res5[node_L] += flux_hllc[4];

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
        Compute_HLLCFlux(node_L, node_R, areavec, flux_hllc, AddTime);
        
        // Compute for LHS
        Res1[node_L] += flux_hllc[0];
        Res2[node_L] += flux_hllc[1];
        Res3[node_L] += flux_hllc[2];
        Res4[node_L] += flux_hllc[3];
        Res5[node_L] += flux_hllc[4];
    }
    
    // Precondition and Transform the Residual in Conservative or Primitive Form
    // Compute_Precondition_Residual_HLLC(); // TODO
}
