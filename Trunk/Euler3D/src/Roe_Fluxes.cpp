/*******************************************************************************
 * File:        Roe_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    3
 * Reference:   A Low-Mach Number Fix for Roe's Approximate Riemann Solver
 *              Felix Rieper, Journal of Computational Physics, 230 (2011)
 ******************************************************************************/

#ifdef DEBUG
#include <assert.h>
#endif

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Solver.h"
#include "Vector3D.h"

// Static Variable for Speed Up
static int Roe_DB               = 0;
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
//! Compute Euler Flux
//------------------------------------------------------------------------------
void Compute_EulerFlux(double *Q_Node, Vector3D areavec, double *Flux_Euler) {
    double nx, ny, nz;
    double rho, u, v, w, rhoet, p, ht, ubar;
    
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    rho   = Q_Node[0];
    u     = Q_Node[1] / rho;
    v     = Q_Node[2] / rho;
    w     = Q_Node[3] / rho;
    rhoet = Q_Node[4];
    p     = (Gamma - 1.0)*(rhoet - 0.5*rho*(u*u + v*v + w*w));
    ht    = (rhoet + p)/rho;
    ubar  = u*nx + v*ny + w*nz;
    
    // Compute Flux
    Flux_Euler[0] = rho*ubar;
    Flux_Euler[1] =   u*Flux_Euler[0] + p*nx;
    Flux_Euler[2] =   v*Flux_Euler[0] + p*ny;
    Flux_Euler[3] =   w*Flux_Euler[0] + p*nz;
    Flux_Euler[4] =  ht*Flux_Euler[0];
}

//------------------------------------------------------------------------------
//! Compute Roe Flux
//------------------------------------------------------------------------------
void Compute_RoeFlux(int node_L, int node_R, Vector3D areavec, double *Flux_Roe) {
    int i, j, k;
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, c_L, ht_L, ubar_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, c_R, ht_R, ubar_R;
    double rho, u, v, w, ht, c, phi, ubar, ubar1, ubar2;
    double Mach, Vmag, fix, sigma;
    Vector3D V, Vn1, Vn2;
    
    double area;
    double nx, ny, nz;
    
    // Initialization
    for (i = 0; i < NEQUATIONS; i++)
        Flux_Roe[i] = 0.0;
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
            if (Order == 2) {
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
        
        // Left Node
        rho_L   = Roe_Q_L[0];
        u_L     = Roe_Q_L[1] / rho_L;
        v_L     = Roe_Q_L[2] / rho_L;
        w_L     = Roe_Q_L[3] / rho_L;
        rhoet_L = Roe_Q_L[4];
        p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L));
        ht_L    = (rhoet_L + p_L)/rho_L;
        c_L     = sqrt((Gamma * p_L) / rho_L);
        ubar_L  = u_L*nx + v_L*ny + w_L*nz;

        // Compute flux_L
        Roe_flux_L[0] = rho_L*ubar_L;
        Roe_flux_L[1] =   u_L*Roe_flux_L[0] + p_L*nx;
        Roe_flux_L[2] =   v_L*Roe_flux_L[0] + p_L*ny;
        Roe_flux_L[3] =   w_L*Roe_flux_L[0] + p_L*nz;
        Roe_flux_L[4] =  ht_L*Roe_flux_L[0];
        
        // Right Node
        rho_R   = Roe_Q_R[0];
        u_R     = Roe_Q_R[1] / rho_R;
        v_R     = Roe_Q_R[2] / rho_R;
        w_R     = Roe_Q_R[3] / rho_R;
        rhoet_R = Roe_Q_R[4];
        p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R));
        ht_R    = (rhoet_R + p_R)/rho_R;
        c_R     = sqrt((Gamma * p_R) / rho_R);
        ubar_R  = u_R*nx + v_R*ny + w_R*nz;
        
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
        c     = (Gamma - 1.0)*ht - phi;
        c     = sqrt(c);
        ubar  = u*nx + v*ny + w*nz;

        // Do if Low Mach Number Fix is Requested
        if (LMRoeFix == 1) {
            // Compute the Normal to Normal of Area vector using Velocity
            V.vec[0] = u;
            V.vec[1] = v;
            V.vec[2] = w;
            Vmag = V.magnitude();
            V.normalize();
            
            Vn1 = V%areavec;
            Vn1.normalize();
            ubar1 = Vmag*(Vn1*V);
            
            Vn2 = Vn1%areavec;
            Vn2.normalize();
            ubar2 = Vmag*(Vn2*V);
            
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

        Roe_M[4][0] = phi/(Gamma - 1.0);
        Roe_M[4][1] = rho * u;
        Roe_M[4][2] = rho * v;
        Roe_M[4][3] = rho * w;
        Roe_M[4][4] = 1.0/(Gamma - 1.0);

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

        Roe_Minv[4][0] = phi;
        Roe_Minv[4][1] = -u * (Gamma - 1.0);
        Roe_Minv[4][2] = -v * (Gamma - 1.0);
        Roe_Minv[4][3] = -w * (Gamma - 1.0);
        Roe_Minv[4][4] = (Gamma - 1.0);

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
        Roe_Pinv[3][1] = 0.5 * nx * fix;
        Roe_Pinv[3][2] = 0.5 * ny * fix;
        Roe_Pinv[3][3] = 0.5 * nz * fix;
        Roe_Pinv[3][4] = 1.0/(2.0 * rho * c);

        Roe_Pinv[4][0] = 0.0;
        Roe_Pinv[4][1] = -0.5 * nx * fix;
        Roe_Pinv[4][2] = -0.5 * ny * fix;
        Roe_Pinv[4][3] = -0.5 * nz * fix;
        Roe_Pinv[4][4] = 1.0/(2.0 * rho * c);

        // START: Computing 
        // if Roe:      A*dQ = M*P*|Lambda|*Pinv*Minv*dQ
        // if LMRoeFix: A*dw = M*P*|Lambda|*Pinv*dw
        // Note: w = non-conservative {rho, v, u, w, p}
        // Where: dQ = Minv*dw
        
        // Calculate T = M*P
        MC_Matrix_Mul_Matrix(5, 5, Roe_M, Roe_P, Roe_T);
        
        // ROE: Calculate Tinv = Pinv*Minv
        if (LMRoeFix != 1)
            MC_Matrix_Mul_Matrix(5, 5, Roe_Pinv, Roe_Minv, Roe_Tinv);
        
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
        
        // Apply Low Mach Fix if Requested
        if (LMRoeFix == 1) {
            // Temporary = Tinv = EigenMatrix*Pinv
            MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Pinv, Roe_Tinv);

            // Get Matrix A = M*P*|Lambda|*Pinv
            MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_Tinv, Roe_A);

            // Compute dw
            Roe_dw[0] = rho_R - rho_L;
            Roe_dw[1] = u_R   - u_L;
            Roe_dw[2] = v_R   - v_L;
            Roe_dw[3] = w_R   - w_L;
            Roe_dw[4] = p_R   - p_L;

            // Compute |A|*dw
            MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dw, Roe_fluxA);
            // END: Computing A*dw = M*P*|Lambda|*Pinv*dw
        } else { // ROE
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
            Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
            Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
            Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
            Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
            Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];

            // Compute |A|*dQ
            MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
            // END: Computing A*dQ = M*P*|Lambda|*Pinv*Minv*dQ
        }
        
        // Compute the Roe Flux
        Flux_Roe[0] = 0.5 * (Roe_flux_L[0] + Roe_flux_R[0] - Roe_fluxA[0]) * area;
        Flux_Roe[1] = 0.5 * (Roe_flux_L[1] + Roe_flux_R[1] - Roe_fluxA[1]) * area;
        Flux_Roe[2] = 0.5 * (Roe_flux_L[2] + Roe_flux_R[2] - Roe_fluxA[2]) * area;
        Flux_Roe[3] = 0.5 * (Roe_flux_L[3] + Roe_flux_R[3] - Roe_fluxA[3]) * area;
        Flux_Roe[4] = 0.5 * (Roe_flux_L[4] + Roe_flux_R[4] - Roe_fluxA[4]) * area;
    } else
        error("Compute_RoeFlux: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Averaged Variables
//------------------------------------------------------------------------------
void Compute_RoeVariables(double *Q_L, double *Q_R, double *Q_Roe) {
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, ht_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, ht_R;
    double rho, u, v, w, ht;
    double sigma;

    // Left Node
    rho_L   = Q_L[0];
    u_L     = Q_L[1] / rho_L;
    v_L     = Q_L[2] / rho_L;
    w_L     = Q_L[3] / rho_L;
    rhoet_L = Q_L[4];
    p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L));
    ht_L    = (rhoet_L + p_L)/rho_L;

    // Right Node
    rho_R   = Q_R[0];
    u_R     = Q_R[1] / rho_R;
    v_R     = Q_R[2] / rho_R;
    w_R     = Q_R[3] / rho_R;
    rhoet_R = Q_R[4];
    p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R));
    ht_R    = (rhoet_R + p_R)/rho_R;
    
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
}

//------------------------------------------------------------------------------
//! Computes the Roe Flux Residual for all Edges Internal and Boundary
//------------------------------------------------------------------------------
void Compute_Residual_Roe(void) {
    int i;
    int node_L, node_R;
    Vector3D areavec;
    double flux_roe[5];
    
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
        Compute_RoeFlux(node_L, node_R, areavec, flux_roe);
        
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
        Compute_RoeFlux(node_L, node_R, areavec, flux_roe);
        
        // Compute for LHS
        Res1[node_L] += flux_roe[0];
        Res2[node_L] += flux_roe[1];
        Res3[node_L] += flux_roe[2];
        Res4[node_L] += flux_roe[3];
        Res5[node_L] += flux_roe[4];
    }
}

