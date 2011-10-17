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
#include "Solver.h"

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
//! Compute Roe Flux Jacobian A
//! Note this function should not be called from Compute_RoeFlux
//------------------------------------------------------------------------------
void Compute_RoeAJacobian(double *Q_L, double *Q_R, Vector3D areavec, double **AJacobian_Roe) {
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, c_L, ht_L, ubar_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, c_R, ht_R, ubar_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, area, nx, ny, nz;
    
    // Initialization
    Roe_Reset();
    
    // Get area vector
    area = areavec.magnitude();
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];

    // Left Node
    rho_L   = Q_L[0];
    u_L     = Q_L[1] / rho_L;
    v_L     = Q_L[2] / rho_L;
    w_L     = Q_L[3] / rho_L;
    rhoet_L = Q_L[4];
    p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L)) + Gauge_Pressure;
    ht_L    = (rhoet_L + p_L)/rho_L;
    c_L     = sqrt((Gamma * p_L) / rho_L);
    ubar_L  = u_L*nx + v_L*ny + w_L*nz;

    // Right Node
    rho_R   = Q_R[0];
    u_R     = Q_R[1] / rho_R;
    v_R     = Q_R[2] / rho_R;
    w_R     = Q_R[3] / rho_R;
    rhoet_R = Q_R[4];
    p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R)) + Gauge_Pressure;
    ht_R    = (rhoet_R + p_R)/rho_R;
    c_R     = sqrt((Gamma * p_R) / rho_R);
    ubar_R  = u_R*nx + v_R*ny + w_R*nz;

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
//! Compute Roe Flux
//------------------------------------------------------------------------------
void Compute_RoeFlux(int node_L, int node_R, Vector3D areavec, double *Flux_Roe) {
    int i, j, k, maxcount;
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, c_L, ht_L, ubar_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, c_R, ht_R, ubar_R;
    double rho, u, v, w, ht, c, phi, ubar, ubar1, ubar2;
    double Mach, nmax, fix, fixUn, sigma;
    Vector3D Vn1, Vn2;
    
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
        p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L)) + Gauge_Pressure;
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
        p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R)) + Gauge_Pressure;
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
        
        // Apply Low Mach Fix if Requested
        if (LMRoeFix == 1) {
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
//! Compute Roe Flux with Weiss and Smith Pre-Conditioner 
//------------------------------------------------------------------------------
void Compute_RoeFlux_WeissSmith_Old(int node_L, int node_R, Vector3D areavec, double *Flux_Roe) {
    int i;
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, c_L, ht_L, ubar_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, c_R, ht_R, ubar_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma;
    
    double area;
    double nx, ny, nz;
    double eps, alpha, beta, pubar, pc, Ur, Mstar, Cstar, deltaUbar;
    double deltaU, deltaP, mach;
    double dtmp;
    
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
        p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L)) + Gauge_Pressure;
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
        p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R)) + Gauge_Pressure;
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
        for (i = 0; i < NEQUATIONS; i++)
            Roe_fluxA[i] = dtmp*(Roe_Q_R[i] - Roe_Q_L[i]);
        
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
        for (i = 0; i < NEQUATIONS; i++)
            Flux_Roe[i] = 0.5 * (Roe_flux_L[i] + Roe_flux_R[i] - Roe_fluxA[i]) * area;
    } else
        error("Compute_RoeFlux: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Weiss and Smith Pre-Conditioner 
//------------------------------------------------------------------------------
void Compute_RoeFlux_WeissSmith(int node_L, int node_R, Vector3D areavec, double *Flux_Roe) {
    int i;
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, c_L, ht_L, ubar_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, c_R, ht_R, ubar_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, area, nx, ny, nz;
    double lambda1, lambda2, lambda3, lambda4, lambda5;
    
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
        p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L)) + Gauge_Pressure;
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
        p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R)) + Gauge_Pressure;
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
        phi  = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c    = (Gamma - 1.0)*ht - phi;
        c    = sqrt(c);
        ubar = u*nx + v*ny + w*nz;
        
        // Compute dQ
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
            
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
        for (i = 0; i < NEQUATIONS; i++)
            Flux_Roe[i] = 0.5 * (Roe_flux_L[i] + Roe_flux_R[i] - Roe_fluxA[i]) * area;
    } else
        error("Compute_RoeFlux: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Cecile Voizat Pre-Conditioner
//! This doesnot work because of preconditioner is not multiplied
//------------------------------------------------------------------------------
void Compute_RoeFlux_CecileVoizat_Old(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss) {
    int i;
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, c_L, ht_L, ubar_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, c_R, ht_R, ubar_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, area, nx, ny, nz;
    
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
        p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L)) + Gauge_Pressure;
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
        p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R)) + Gauge_Pressure;
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
        phi  = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c    = (Gamma - 1.0)*ht - phi;
        c    = sqrt(c);
        ubar = u*nx + v*ny + w*nz;
        
        // Compute dQ
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
            
        double alpha, beta, zeta, q2;
        
        // Set the precondition parameters
        q2    = u*u + v*v + w*w;
        sigma = MIN(1.0, MAX(Inf_Mach, sqrt(u*u + v*v + w*w)/c));
        alpha = ubar*ubar*(sigma*sigma*sigma*sigma - 2.0*sigma*sigma + 1.0) + 4.0*c*c*sigma*sigma;
        alpha = sqrt(alpha);
        beta  = (sigma*sigma - 1.0)*ubar + alpha;
        zeta  = (1.0 - sigma*sigma)*ubar + alpha;
        
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(0.5*(ubar + sigma*sigma*ubar + alpha));
        Roe_Eigen[4][4] = fabs(0.5*(ubar + sigma*sigma*ubar - alpha));
        
        // Compute T
        Roe_T[0][0] = -rho*nx/Gamma;
        Roe_T[0][1] = -rho*ny/Gamma;
        Roe_T[0][2] = -rho*nz/Gamma;
        Roe_T[0][3] = 2.0*rho/zeta;
        Roe_T[0][4] = 2.0*rho/beta;
        
        Roe_T[1][0] = -u*rho*nx/Gamma;
        Roe_T[1][1] = -rho*(u*ny + Gamma*nz)/Gamma;
        Roe_T[1][2] = rho*(ny - u*nz/Gamma);
        Roe_T[1][3] = rho*(nx + 2.0*u/zeta);
        Roe_T[1][4] = rho*(-nx + 2.0*u/beta);
        
        Roe_T[2][0] = rho*(-v*nx/Gamma + nz);
        Roe_T[2][1] = -v*rho*ny/Gamma;
        Roe_T[2][2] = -rho*(Gamma*nx + v*nz)/Gamma;
        Roe_T[2][3] = rho*(ny + 2.0*v/zeta);
        Roe_T[2][4] = rho*(-ny + 2.0*v/beta);
        
        Roe_T[3][0] = -rho*(w*nx + Gamma*ny)/Gamma;
        Roe_T[3][1] = rho*(nx - w*ny/Gamma);
        Roe_T[3][2] = -w*rho*nz/Gamma;
        Roe_T[3][3] = rho*(nz+ 2.0*w/zeta);
        Roe_T[3][4] = rho*(-nz + 2.0*w/beta);
        
        Roe_T[4][0] = -w*rho*ny + v*rho*nz - rho*nx*q2/(2.0*Gamma);
        Roe_T[4][1] =  w*rho*nx - u*rho*nz - rho*ny*q2/(2.0*Gamma);
        Roe_T[4][2] = -v*rho*nx + u*rho*ny - rho*nz*q2/(2.0*Gamma);
        Roe_T[4][3] = rho*(ubar + (2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*zeta));
        Roe_T[4][4] = rho*(-ubar + (2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*beta));
        
        // Compute Tinv
        Roe_Tinv[0][0] = (2.0*c*c*(w*ny - v*nz) + Gamma*nx*(-2.0*c*c + (Gamma - 1.0)*q2))/(2.0*c*c*rho);
        Roe_Tinv[0][1] = -(u*(Gamma - 1.0)*Gamma*nx)/(c*c*rho);
        Roe_Tinv[0][2] = ((-v*(Gamma - 1.0)*Gamma*nx)/(c*c) + nz)/rho;
        Roe_Tinv[0][3] = -((w*(Gamma - 1.0)*Gamma*nx)/(c*c) + ny)/rho;
        Roe_Tinv[0][4] = (Gamma - 1.0)*Gamma*nx/(c*c*rho);
        
        Roe_Tinv[1][0] = (2.0*c*c*(u*nz - w*nx) + Gamma*ny*(-2.0*c*c + (Gamma - 1.0)*q2))/(2.0*c*c*rho);
        Roe_Tinv[1][1] = -((u*(Gamma - 1.0)*Gamma*ny)/(c*c) + nz)/rho;
        Roe_Tinv[1][2] = -(v*(Gamma - 1.0)*Gamma*ny)/(c*c*rho);
        Roe_Tinv[1][3] = ((-w*(Gamma - 1.0)*Gamma*ny)/(c*c) + nx)/rho;
        Roe_Tinv[1][4] = (Gamma - 1.0)*Gamma*ny/(c*c*rho);
        
        Roe_Tinv[2][0] = (2.0*c*c*(v*nx - u*ny) + Gamma*nz*(-2.0*c*c + (Gamma - 1.0)*q2))/(2.0*c*c*rho);
        Roe_Tinv[2][1] = ((-u*(Gamma - 1.0)*Gamma*nz)/(c*c) + ny)/rho;
        Roe_Tinv[2][2] = -((v*(Gamma - 1.0)*Gamma*nz)/(c*c) + nx)/rho;
        Roe_Tinv[2][3] = -(w*(Gamma - 1.0)*Gamma*nz)/(c*c*rho);
        Roe_Tinv[2][4] = (Gamma - 1.0)*Gamma*nz/(c*c*rho);
        
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
        
        // Tmp => P = |Lambda|*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);

        // Get Matrix |A| = T*|Lambda|*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);

        // Compute |A|*dQ
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // Compute the Flux_Roe_Conv and Flux_Roe_Diss
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = 0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux_CecileVoizat: Invalid Node - %d", node_L);
}

//------------------------------------------------------------------------------
//! Compute Roe Flux with Cecile Voizat Pre-Conditioner
//------------------------------------------------------------------------------
void Compute_RoeFlux_CecileVoizat(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss) {
    int i;
    double rho_L, u_L, v_L, w_L, rhoet_L, p_L, c_L, ht_L, ubar_L;
    double rho_R, u_R, v_R, w_R, rhoet_R, p_R, c_R, ht_R, ubar_R;
    double rho, u, v, w, ht, c, phi, ubar;
    double sigma, area, nx, ny, nz;
    
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
        p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L)) + Gauge_Pressure;
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
        p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R)) + Gauge_Pressure;
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
        phi  = 0.5*(Gamma - 1.0)*(u*u + v*v + w*w);
        c    = (Gamma - 1.0)*ht - phi;
        c    = sqrt(c);
        ubar = u*nx + v*ny + w*nz;
        
        // Compute dQ
        Roe_dQ[0] = Roe_Q_R[0] - Roe_Q_L[0];
        Roe_dQ[1] = Roe_Q_R[1] - Roe_Q_L[1];
        Roe_dQ[2] = Roe_Q_R[2] - Roe_Q_L[2];
        Roe_dQ[3] = Roe_Q_R[3] - Roe_Q_L[3];
        Roe_dQ[4] = Roe_Q_R[4] - Roe_Q_L[4];
        
        int nid;
        double alpha, beta, zeta, q2;
        double lc, lp, dp_L, dp_R, dp, lmach_L, lmach_R, lmach;
        
        // Set the precondition parameters
        q2   = u*u + v*v + w*w;
        //sigma = MIN(1.0, MAX(Inf_Mach, sqrt(q2)/c));
        dp_L = dp_R = 0.0;
        lmach_L = lmach_R = 0.0;
        for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L+1]; i++) {
            nid     = crs_JA_Node2Node[i];
            lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]) + Gauge_Pressure;
            lc      = sqrt((Gamma * lp) / Q1[nid]);
            lmach_L = MAX(lmach_L, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
            dp_L    = MAX(dp_L, fabs(p_L - lp));
        }
        dp    = dp_L;
        lmach = lmach_L;
        // Check if Right Node is not a Ghost Boundary Node
        if (node_R < nNode) {
            for (i = crs_IA_Node2Node[node_R]; i < crs_IA_Node2Node[node_R+1]; i++) {
                nid  = crs_JA_Node2Node[i];
                lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]) + Gauge_Pressure;
                lc      = sqrt((Gamma * lp) / Q1[nid]);
                lmach_R = MAX(lmach_R, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
                dp_R    = MAX(dp_R, fabs(p_R - lp));
            }
            dp    = 0.5*(dp_L + dp_R);
            lmach = 0.5*(lmach_L + lmach_R);
        }
        sigma = MIN(1.0, MAX(MAX(MAX(1.0e-5, sqrt(q2)/c), lmach), dp/(rho*c*c)));
        if (sqrt(q2)/c > 0.3)
            sigma = 1.0;
        
        alpha = ubar*ubar*(sigma*sigma*sigma*sigma - 2.0*sigma*sigma + 1.0) + 4.0*c*c*sigma*sigma;
        alpha = sqrt(alpha);
        beta  = (sigma*sigma - 1.0)*ubar + alpha;
        zeta  = (1.0 - sigma*sigma)*ubar + alpha;
        
        // Compute the Eigenvalues of the Precondition Dissipation Term |Lambda|
        Roe_Eigen[0][0] = fabs(ubar);
        Roe_Eigen[1][1] = Roe_Eigen[0][0];
        Roe_Eigen[2][2] = Roe_Eigen[0][0];
        Roe_Eigen[3][3] = fabs(0.5*(ubar + sigma*sigma*ubar + alpha));
        Roe_Eigen[4][4] = fabs(0.5*(ubar + sigma*sigma*ubar - alpha));
        
        // Compute T = M3.P
        Roe_T[0][0] = -rho*nx/Gamma;
        Roe_T[0][1] = -rho*ny/Gamma;
        Roe_T[0][2] = -rho*nz/Gamma;
        Roe_T[0][3] = 2.0*rho*sigma*sigma/zeta;
        Roe_T[0][4] = 2.0*rho*sigma*sigma/beta;
        
        Roe_T[1][0] = -u*rho*nx/Gamma;
        Roe_T[1][1] = -rho*(u*ny + Gamma*nz)/Gamma;
        Roe_T[1][2] = rho*(ny - u*nz/Gamma);
        Roe_T[1][3] = rho*(nx + 2.0*u*sigma*sigma/zeta);
        Roe_T[1][4] = rho*(-nx + 2.0*u*sigma*sigma/beta);
        
        Roe_T[2][0] = rho*(-v*nx/Gamma + nz);
        Roe_T[2][1] = -v*rho*ny/Gamma;
        Roe_T[2][2] = -rho*(Gamma*nx + v*nz)/Gamma;
        Roe_T[2][3] = rho*(ny + 2.0*v*sigma*sigma/zeta);
        Roe_T[2][4] = rho*(-ny + 2.0*v*sigma*sigma/beta);
        
        Roe_T[3][0] = -rho*(w*nx + Gamma*ny)/Gamma;
        Roe_T[3][1] = rho*(nx - w*ny/Gamma);
        Roe_T[3][2] = -w*rho*nz/Gamma;
        Roe_T[3][3] = rho*(nz+ 2.0*w*sigma*sigma/zeta);
        Roe_T[3][4] = rho*(-nz + 2.0*w*sigma*sigma/beta);
        
        Roe_T[4][0] = -w*rho*ny + v*rho*nz - rho*nx*q2/(2.0*Gamma);
        Roe_T[4][1] =  w*rho*nx - u*rho*nz - rho*ny*q2/(2.0*Gamma);
        Roe_T[4][2] = -v*rho*nx + u*rho*ny - rho*nz*q2/(2.0*Gamma);
        Roe_T[4][3] = rho*(ubar + sigma*sigma*(2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*zeta));
        Roe_T[4][4] = rho*(-ubar + sigma*sigma*(2.0*c*c + (Gamma - 1.0)*q2)/((Gamma - 1.0)*beta));
        
        // Compute Tinv = Pinv.M3inv
        Roe_Tinv[0][0] = (2.0*c*c*(w*ny - v*nz) + Gamma*nx*(-2.0*c*c + (Gamma - 1.0)*q2))/(2.0*c*c*rho);
        Roe_Tinv[0][1] = -(u*(Gamma - 1.0)*Gamma*nx)/(c*c*rho);
        Roe_Tinv[0][2] = ((-v*(Gamma - 1.0)*Gamma*nx)/(c*c) + nz)/rho;
        Roe_Tinv[0][3] = -((w*(Gamma - 1.0)*Gamma*nx)/(c*c) + ny)/rho;
        Roe_Tinv[0][4] = (Gamma - 1.0)*Gamma*nx/(c*c*rho);
        
        Roe_Tinv[1][0] = (2.0*c*c*(u*nz - w*nx) + Gamma*ny*(-2.0*c*c + (Gamma - 1.0)*q2))/(2.0*c*c*rho);
        Roe_Tinv[1][1] = -((u*(Gamma - 1.0)*Gamma*ny)/(c*c) + nz)/rho;
        Roe_Tinv[1][2] = -(v*(Gamma - 1.0)*Gamma*ny)/(c*c*rho);
        Roe_Tinv[1][3] = ((-w*(Gamma - 1.0)*Gamma*ny)/(c*c) + nx)/rho;
        Roe_Tinv[1][4] = (Gamma - 1.0)*Gamma*ny/(c*c*rho);
        
        Roe_Tinv[2][0] = (2.0*c*c*(v*nx - u*ny) + Gamma*nz*(-2.0*c*c + (Gamma - 1.0)*q2))/(2.0*c*c*rho);
        Roe_Tinv[2][1] = ((-u*(Gamma - 1.0)*Gamma*nz)/(c*c) + ny)/rho;
        Roe_Tinv[2][2] = -((v*(Gamma - 1.0)*Gamma*nz)/(c*c) + nx)/rho;
        Roe_Tinv[2][3] = -(w*(Gamma - 1.0)*Gamma*nz)/(c*c*rho);
        Roe_Tinv[2][4] = (Gamma - 1.0)*Gamma*nz/(c*c*rho);
        
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
        
        // Tmp => P = |Lambda|*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_Eigen, Roe_Tinv, Roe_P);

        // Get Matrix |A| = T*|Lambda|*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Roe_T, Roe_P, Roe_A);

        // Compute |A|*dQ
        MC_Matrix_Mul_Vector(5, 5, Roe_A, Roe_dQ, Roe_fluxA);
        
        // Compute the Flux_Roe_Conv and Flux_Roe_Diss
        for (i = 0; i < NEQUATIONS; i++) {
            Flux_Roe_Conv[i] = 0.5*(Roe_flux_L[i] + Roe_flux_R[i])*area;
            Flux_Roe_Diss[i] = 0.5*Roe_fluxA[i]*area;
        }
    } else
        error("Compute_RoeFlux_CecileVoizat: Invalid Node - %d", node_L);
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
    p_L     = (Gamma - 1.0)*(rhoet_L - 0.5*rho_L*(u_L*u_L + v_L*v_L + w_L*w_L)) + Gauge_Pressure;
    ht_L    = (rhoet_L + p_L)/rho_L;

    // Right Node
    rho_R   = Q_R[0];
    u_R     = Q_R[1] / rho_R;
    v_R     = Q_R[2] / rho_R;
    w_R     = Q_R[3] / rho_R;
    rhoet_R = Q_R[4];
    p_R     = (Gamma - 1.0)*(rhoet_R - 0.5*rho_R*(u_R*u_R + v_R*v_R + w_R*w_R)) + Gauge_Pressure;
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
        switch (SolverScheme) {
            case SOLVER_SCHEME_ROE: // Roe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe);
                break;
            case SOLVER_SCHEME_LMROE: // LMRoe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe);
                break;
            case SOLVER_SCHEME_ROE_WS: // Roe Weiss Smith Precondition
                Compute_RoeFlux_WeissSmith(node_L, node_R, areavec, flux_roe);
                break;
            case SOLVER_SCHEME_ROE_CV: // Roe Cecile Voizat Precondition
                Compute_RoeFlux_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss);
                break;
            default:
                error("Compute_Residual_Roe: Invalid Solver Scheme - %d -1", SolverScheme);
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
        if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
            // Compute for LHS
            Res1_Diss[node_L] -= flux_roe_diss[0];
            Res2_Diss[node_L] -= flux_roe_diss[1];
            Res3_Diss[node_L] -= flux_roe_diss[2];
            Res4_Diss[node_L] -= flux_roe_diss[3];
            Res5_Diss[node_L] -= flux_roe_diss[4];

            Res1_Diss[node_R] += flux_roe_diss[0];
            Res2_Diss[node_R] += flux_roe_diss[1];
            Res3_Diss[node_R] += flux_roe_diss[2];
            Res4_Diss[node_R] += flux_roe_diss[3];
            Res5_Diss[node_R] += flux_roe_diss[4];
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
        switch (SolverScheme) {
            case SOLVER_SCHEME_ROE: // Roe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe);
                break;
            case SOLVER_SCHEME_LMROE: // LMRoe
                Compute_RoeFlux(node_L, node_R, areavec, flux_roe);
                break;
            case SOLVER_SCHEME_ROE_WS: // Roe Weiss Smith Precondition
                Compute_RoeFlux_WeissSmith(node_L, node_R, areavec, flux_roe);
                break;
            case SOLVER_SCHEME_ROE_CV: // Roe Cecile Voizat Precondition
                Compute_RoeFlux_CecileVoizat(node_L, node_R, areavec, flux_roe, flux_roe_diss);
                break;
            default:
                error("Compute_Residual_Roe: Invalid Solver Scheme - %d -2", SolverScheme);
                break;
        }
        
        // Compute for LHS
        Res1[node_L] += flux_roe[0];
        Res2[node_L] += flux_roe[1];
        Res3[node_L] += flux_roe[2];
        Res4[node_L] += flux_roe[3];
        Res5[node_L] += flux_roe[4];
        
        // Dissipation Term
        if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
            // Compute for LHS
            Res1_Diss[node_L] -= flux_roe_diss[0];
            Res2_Diss[node_L] -= flux_roe_diss[1];
            Res3_Diss[node_L] -= flux_roe_diss[2];
            Res4_Diss[node_L] -= flux_roe_diss[3];
            Res5_Diss[node_L] -= flux_roe_diss[4];
        }
    }
    
    // Multiply by Precondition Matrix and Construct the Flux
    if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
        int j, nid;
        double rho, u, v, w, c, ht, rhoet, p, q2, sigma, reta;
        double lp, lc, dp, lmach;
        for (i = 0; i < nNode; i++) {
            rho   = Q1[i];
            u     = Q2[i]/rho;
            v     = Q3[i]/rho;
            w     = Q4[i]/rho;
            rhoet = Q5[i];
            p     = (Gamma - 1.0)*(rhoet - 0.5*rho*(u*u + v*v + w*w)) + Gauge_Pressure;
            ht    = (rhoet + p)/rho;
            c     = sqrt((Gamma * p) / rho);
            q2    = u*u + v*v + w*w;
            //sigma = MAX(MAX(MIN(1.0, Ref_Mach), sqrt(u*u + v*v + w*w)/c), 1.0e-5);
            dp    = 0.0;
            lmach = 0.0;
            for (j = crs_IA_Node2Node[i]; j < crs_IA_Node2Node[i+1]; j++) {
                nid     = crs_JA_Node2Node[j];
                lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]) + Gauge_Pressure;
                lc      = sqrt((Gamma * lp) / Q1[nid]);
                lmach = MAX(lmach, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
                dp    = MAX(dp, fabs(p - lp));
            }
            
             sigma = MIN(1.0, MAX(MAX(MAX(1.0e-5, sqrt(q2)/c), lmach), dp/(rho*c*c)));
             if (sqrt(q2)/c > 0.3)
                 sigma  = 1.0;
            //sigma = MIN(1.0, MAX(1.0e-6, sqrt(u*u + v*v + w*w)/c));
            //sigma = sqrt(u*u + v*v + w*w)/c;
            
            // Compute the Precondition in Conservative form
            reta = (Gamma - 1.0)*(sigma*sigma - 1.0)/(c*c);

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

            Roe_P[4][0] = 0.5*q2*((sigma*sigma - 1.0) + 0.5*q2*reta);
            Roe_P[4][1] = -u*((sigma*sigma - 1.0) + 0.5*q2*reta);
            Roe_P[4][2] = -v*((sigma*sigma - 1.0) + 0.5*q2*reta);
            Roe_P[4][3] = -w*((sigma*sigma - 1.0) + 0.5*q2*reta);
            Roe_P[4][4] = sigma*sigma + 0.5*q2*reta;

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
            MC_Matrix_Mul_Vector(5, 5, Roe_P, flux_roe_diss, flux_roe);

            // Compute the Preconditioned Roe Flux Residual
            Res1[i] = flux_roe[0] + Res1_Diss[i];
            Res2[i] = flux_roe[1] + Res2_Diss[i];
            Res3[i] = flux_roe[2] + Res3_Diss[i];
            Res4[i] = flux_roe[3] + Res4_Diss[i];
            Res5[i] = flux_roe[4] + Res5_Diss[i];
        }
    }
}

