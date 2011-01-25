/*******************************************************************************
 * File:        Roe_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <assert.h>

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Compute_Residual(void) {
    int i, j, k;
    int node_L, node_R;
    double rho_L, u_L, v_L, w_L, et_L, e_L, p_L, c_L, h_L, ht_L;
    double rho_R, u_R, v_R, w_R, et_R, e_R, p_R, c_R, h_R, ht_R;
    double rho, u, v, w, ht, c, phi;
    double temp;
    
    Vector3D areavec;
    double   area;
    double nx, ny, nz;
    double f_vec[5];
    double g_vec[5];
    double h_vec[5];
    double fluxA[5];
    double flux_L[5];
    double flux_R[5];
    double *dQ;
    double **A;
    double **Eigen;
    double **M;
    double **Minv;
    double **P;
    double **Pinv;
    double **T;
    double **Tinv;

    M     = (double **) malloc(5*sizeof(double*));
    Minv  = (double **) malloc(5*sizeof(double*));
    P     = (double **) malloc(5*sizeof(double*));
    Pinv  = (double **) malloc(5*sizeof(double*));
    T     = (double **) malloc(5*sizeof(double*));
    Tinv  = (double **) malloc(5*sizeof(double*));
    A     = (double **) malloc(5*sizeof(double*));
    Eigen = (double **) malloc(5*sizeof(double*));
    dQ    = (double *) malloc(5*sizeof(double));
    for (i = 0; i < 5; i++) {
        M[i]     = (double *) malloc(5*sizeof(double));
        Minv[i]  = (double *) malloc(5*sizeof(double));
        P[i]     = (double *) malloc(5*sizeof(double));
        Pinv[i]  = (double *) malloc(5*sizeof(double));
        T[i]     = (double *) malloc(5*sizeof(double));
        Tinv[i]  = (double *) malloc(5*sizeof(double));
        A[i]     = (double *) malloc(5*sizeof(double));
        Eigen[i] = (double *) malloc(5*sizeof(double));
    }

    // Initialialize
    for (i = 0; i < 5; i++) {
        fluxA[i]  = 0.0;
        flux_L[i] = 0.0;
        flux_R[i] = 0.0;
        f_vec[i]  = 0.0;
        g_vec[i]  = 0.0;
        h_vec[i]  = 0.0;
        dQ[i]     = 0.0;
        for (j = 0; j < 5; j++) {
            A[i][j]     = 0.0;
            M[i][j]     = 0.0;
            Minv[i][j]  = 0.0;
            P[i][j]     = 0.0;
            Pinv[i][j]  = 0.0;
            T[i][j]     = 0.0;
            Tinv[i][j]  = 0.0;
            Eigen[i][j] = 0.0;
        }
    }
    
    // Internal Edges
    for (i = 0; i < nEdge; i++) {
        // Get two nodes of edge
        node_L = intEdge[i].node[0];
        node_R = intEdge[i].node[1];

        assert(node_R > node_L);

        // Get area vector
        areavec = intEdge[i].areav;
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Left Node
        rho_L = Q1[node_L];
        u_L   = Q2[node_L] / rho_L;
        v_L   = Q3[node_L] / rho_L;
        w_L   = Q4[node_L] / rho_L;
        et_L  = Q5[node_L] / rho_L;
        e_L   = et_L - 0.5 * (u_L * u_L + v_L * v_L + w_L * w_L);
        p_L   = (Gamma - 1.0) * rho_L * e_L;
        c_L   = sqrt((Gamma * p_L) / rho_L);
        h_L   = (c_L * c_L) / (Gamma - 1.0);
        ht_L  = h_L + 0.5 * (u_L * u_L + v_L * v_L + w_L * w_L);

        // F Vector
        f_vec[0] = rho_L * u_L;
        f_vec[1] = rho_L * u_L * u_L + p_L;
        f_vec[2] = rho_L * u_L * v_L;
        f_vec[3] = rho_L * u_L * w_L;
        f_vec[4] = rho_L * ht_L * u_L;

        // G Vector
        g_vec[0] = rho_L * v_L;
        g_vec[1] = rho_L * v_L * u_L;
        g_vec[2] = rho_L * v_L * v_L + p_L;
        g_vec[3] = rho_L * v_L * w_L;
        g_vec[4] = rho_L * ht_L * v_L;

        // H Vector
        h_vec[0] = rho_L * w_L;
        h_vec[1] = rho_L * w_L * u_L;
        h_vec[2] = rho_L * w_L * v_L;
        h_vec[3] = rho_L * w_L * w_L + p_L;
        h_vec[4] = rho_L * ht_L * w_L;

        // Accumulate to flux_L
        flux_L[0] = f_vec[0] * nx + g_vec[0] * ny + h_vec[0] * nz;
        flux_L[1] = f_vec[1] * nx + g_vec[1] * ny + h_vec[1] * nz;
        flux_L[2] = f_vec[2] * nx + g_vec[2] * ny + h_vec[2] * nz;
        flux_L[3] = f_vec[3] * nx + g_vec[3] * ny + h_vec[3] * nz;
        flux_L[4] = f_vec[4] * nx + g_vec[4] * ny + h_vec[4] * nz;

        // Right Node
        rho_R = Q1[node_R];
        u_R   = Q2[node_R] / rho_R;
        v_R   = Q3[node_R] / rho_R;
        w_R   = Q4[node_R] / rho_R;
        et_R  = Q5[node_R] / rho_R;
        e_R   = et_R - 0.5 * (u_R * u_R + v_R * v_R + w_R * w_R);
        p_R   = (Gamma - 1.0) * rho_R * e_R;
        c_R   = sqrt((Gamma * p_R) / rho_R);
        h_R   = (c_R * c_R) / (Gamma - 1.0);
        ht_R  = h_R + 0.5 * (u_R * u_R + v_R * v_R + w_R * w_R);

        // F Vector
        f_vec[0] = rho_R * u_R;
        f_vec[1] = rho_R * u_R * u_R + p_R;
        f_vec[2] = rho_R * u_R * v_R;
        f_vec[3] = rho_R * u_R * w_R;
        f_vec[4] = rho_R * ht_R * u_R;

        // G Vector
        g_vec[0] = rho_R * v_R;
        g_vec[1] = rho_R * v_R * u_R;
        g_vec[2] = rho_R * v_R * v_R + p_R;
        g_vec[3] = rho_R * v_R * w_R;
        g_vec[4] = rho_R * ht_R * v_R;

        // H Vector
        h_vec[0] = rho_R * w_R;
        h_vec[1] = rho_R * w_R * u_R;
        h_vec[2] = rho_R * w_R * v_R;
        h_vec[3] = rho_R * w_R * w_R + p_R;
        h_vec[4] = rho_R * ht_R * w_R;

        // Accumulate to flux_R
        flux_R[0] = f_vec[0] * nx + g_vec[0] * ny + h_vec[0] * nz;
        flux_R[1] = f_vec[1] * nx + g_vec[1] * ny + h_vec[1] * nz;
        flux_R[2] = f_vec[2] * nx + g_vec[2] * ny + h_vec[2] * nz;
        flux_R[3] = f_vec[3] * nx + g_vec[3] * ny + h_vec[3] * nz;
        flux_R[4] = f_vec[4] * nx + g_vec[4] * ny + h_vec[4] * nz;
        
        // ROE AVERAGE VARIABLES
        temp = sqrt(rho_R) + sqrt(rho_L);
        rho  = sqrt(rho_R * rho_L);
        u    = (u_L * sqrt(rho_L)  + u_R * sqrt(rho_R)) / temp;
        v    = (v_L * sqrt(rho_L)  + v_R * sqrt(rho_R)) / temp;
        w    = (w_L * sqrt(rho_L)  + w_R * sqrt(rho_R)) / temp;
        ht   = (ht_L * sqrt(rho_L) + ht_R * sqrt(rho_R)) / temp;
        c    = (Gamma - 1.0) * (ht - 0.5 * (u * u + v * v + w * w));
        c    = sqrt(c);
        phi  = 0.5 * (Gamma - 1.0)*(u * u + v * v + w * w);

        // M
        M[0][0] = 1.0;
        M[0][1] = 0.0;
        M[0][2] = 0.0;
        M[0][3] = 0.0;
        M[0][4] = 0.0;

        M[1][0] = u;
        M[1][1] = rho;
        M[1][2] = 0.0;
        M[1][3] = 0.0;
        M[1][4] = 0.0;

        M[2][0] = v;
        M[2][1] = 0.0;
        M[2][2] = rho;
        M[2][3] = 0.0;
        M[2][4] = 0.0;

        M[3][0] = w;
        M[3][1] = 0.0;
        M[3][2] = 0.0;
        M[3][3] = rho;
        M[3][4] = 0.0;

        M[4][0] = phi/(Gamma - 1.0);
        M[4][1] = rho * u;
        M[4][2] = rho * v;
        M[4][3] = rho * w;
        M[4][4] = 1.0/(Gamma - 1.0);

        // Minv
        Minv[0][0] = 1.0;
        Minv[0][1] = 0.0;
        Minv[0][2] = 0.0;
        Minv[0][3] = 0.0;
        Minv[0][4] = 0.0;

        Minv[1][0] = -u/rho;
        Minv[1][1] = 1.0/rho;
        Minv[1][2] = 0.0;
        Minv[1][3] = 0.0;
        Minv[1][4] = 0.0;

        Minv[2][0] = -v/rho;
        Minv[2][1] = 0.0;
        Minv[2][2] = 1.0/rho;
        Minv[2][3] = 0.0;
        Minv[2][4] = 0.0;

        Minv[3][0] = -w/rho;
        Minv[3][1] = 0.0;
        Minv[3][2] = 0.0;
        Minv[3][3] = 1.0/rho;
        Minv[3][4] = 0.0;

        Minv[4][0] = phi;
        Minv[4][1] = -u * (Gamma - 1.0);
        Minv[4][2] = -v * (Gamma - 1.0);
        Minv[4][3] = -w * (Gamma - 1.0);
        Minv[4][4] = (Gamma - 1.0);

        // P
        P[0][0] = nx;
        P[0][1] = ny;
        P[0][2] = nz;
        P[0][3] = rho/c;
        P[0][4] = rho/c;

        P[1][0] = 0.0;
        P[1][1] = -nz;
        P[1][2] = ny;
        P[1][3] = nx;
        P[1][4] = -nx;

        P[2][0] = nz;
        P[2][1] = 0.0;
        P[2][2] = -nx;
        P[2][3] = ny;
        P[2][4] = -ny;

        P[3][0] = -ny;
        P[3][1] = nx;
        P[3][2] = 0.0;
        P[3][3] = nz;
        P[3][4] = -nz;

        P[4][0] = 0.0;
        P[4][1] = 0.0;
        P[4][2] = 0.0;
        P[4][3] = rho * c;
        P[4][4] = rho * c;

        // Pinv
        Pinv[0][0] = nx;
        Pinv[0][1] = 0.0;
        Pinv[0][2] = nz;
        Pinv[0][3] = -ny;
        Pinv[0][4] = -nx/(c * c);

        Pinv[1][0] = ny;
        Pinv[1][1] = -nz;
        Pinv[1][2] = 0.0;
        Pinv[1][3] = nx;
        Pinv[1][4] = -ny/(c * c);

        Pinv[2][0] = nz;
        Pinv[2][1] = ny;
        Pinv[2][2] = -nx;
        Pinv[2][3] = 0.0;
        Pinv[2][4] = -nz/(c * c);

        Pinv[3][0] = 0.0;
        Pinv[3][1] = 0.5 * nx;
        Pinv[3][2] = 0.5 * ny;
        Pinv[3][3] = 0.5 * nz;
        Pinv[3][4] = 1.0/(2.0 * rho * c);

        Pinv[4][0] = 0.0;
        Pinv[4][1] = -0.5 * nx;
        Pinv[4][2] = -0.5 * ny;
        Pinv[4][3] = -0.5 * nz;
        Pinv[4][4] = 1.0/(2.0 * rho * c);

        // Calculate T = M*P
        MC_Matrix_Mul_Matrix(5, 5, M, P, T);
        // Calculate Tinv = Pinv*Minv
        MC_Matrix_Mul_Matrix(5, 5, Pinv, Minv, Tinv);
        
        // EigenMatrix
        for (j = 0; j < 5; j++)
            for (k = 0; k < 5; k++)
                Eigen[j][k] = 0.0;
        
        Eigen[0][0] = fabs(u * nx + v * ny + w * nz);
        Eigen[1][1] = fabs(u * nx + v * ny + w * nz);
        Eigen[2][2] = fabs(u * nx + v * ny + w * nz);
        Eigen[3][3] = fabs((u * nx + v * ny + w * nz) + c);
        Eigen[4][4] = fabs((u * nx + v * ny + w * nz) - c);

        // EigenMatrix*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Eigen, Tinv, A);

        // Tinv = EigenMatrix*Tinv
        for (j = 0; j < 5; j++) {
            for (k = 0; k < 5; k++) {
                Tinv[j][k] = A[j][k];
                A[j][k] = 0.0;
            }
        }
        
        // Get Matrix A
        MC_Matrix_Mul_Matrix(5, 5, T, Tinv, A);

        // Set up matrix vector multiplication
        dQ[0] = dQ[1] = dQ[2] = dQ[3] = dQ[4] = 0.0;
        dQ[0] = Q1[node_R] - Q1[node_L];
        dQ[1] = Q2[node_R] - Q2[node_L];
        dQ[2] = Q3[node_R] - Q3[node_L];
        dQ[3] = Q4[node_R] - Q4[node_L];
        dQ[4] = Q5[node_R] - Q5[node_L];

        fluxA[0] = fluxA[1] = fluxA[2] = fluxA[3] = fluxA[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, A, dQ, fluxA);
        
        // Compute for LHS
        Res1[node_L] += 0.5 * (flux_L[0] + flux_R[0] - fluxA[0]) * area;
        Res2[node_L] += 0.5 * (flux_L[1] + flux_R[1] - fluxA[1]) * area;
        Res3[node_L] += 0.5 * (flux_L[2] + flux_R[2] - fluxA[2]) * area;
        Res4[node_L] += 0.5 * (flux_L[3] + flux_R[3] - fluxA[3]) * area;
        Res5[node_L] += 0.5 * (flux_L[4] + flux_R[4] - fluxA[4]) * area;

        Res1[node_R] -= 0.5 * (flux_L[0] + flux_R[0] - fluxA[0]) * area;
        Res2[node_R] -= 0.5 * (flux_L[1] + flux_R[1] - fluxA[1]) * area;
        Res3[node_R] -= 0.5 * (flux_L[2] + flux_R[2] - fluxA[2]) * area;
        Res4[node_R] -= 0.5 * (flux_L[3] + flux_R[3] - fluxA[3]) * area;
        Res5[node_R] -= 0.5 * (flux_L[4] + flux_R[4] - fluxA[4]) * area;
    }
    
    // Initialialize
    for (i = 0; i < 5; i++) {
        fluxA[i]  = 0.0;
        flux_L[i] = 0.0;
        flux_R[i] = 0.0;
        f_vec[i]  = 0.0;
        g_vec[i]  = 0.0;
        h_vec[i]  = 0.0;
        dQ[i]     = 0.0;
        for (j = 0; j < 5; j++) {
            A[i][j]     = 0.0;
            M[i][j]     = 0.0;
            Minv[i][j]  = 0.0;
            P[i][j]     = 0.0;
            Pinv[i][j]  = 0.0;
            T[i][j]     = 0.0;
            Tinv[i][j]  = 0.0;
            Eigen[i][j] = 0.0;
        }
    }

    // Boundary Edges
    for (i = 0; i < nBEdge; i++) {
        // Get two nodes of edge
        node_L = bndEdge[i].node[0];
        node_R = bndEdge[i].node[1];

        assert(node_R > node_L);

        // Get area vector
        areavec = bndEdge[i].areav;
        area = areavec.magnitude();
        areavec.normalize();
        nx = areavec.vec[0];
        ny = areavec.vec[1];
        nz = areavec.vec[2];
        
        // Left Node
        rho_L = Q1[node_L];
        u_L   = Q2[node_L] / rho_L;
        v_L   = Q3[node_L] / rho_L;
        w_L   = Q4[node_L] / rho_L;
        et_L  = Q5[node_L] / rho_L;
        e_L   = et_L - 0.5 * (u_L * u_L + v_L * v_L + w_L * w_L);
        p_L   = (Gamma - 1.0) * rho_L * e_L;
        c_L   = sqrt((Gamma * p_L) / rho_L);
        h_L   = (c_L * c_L) / (Gamma - 1.0);
        ht_L  = h_L + 0.5 * (u_L * u_L + v_L * v_L + w_L * w_L);

        // F Vector
        f_vec[0] = rho_L * u_L;
        f_vec[1] = rho_L * u_L * u_L + p_L;
        f_vec[2] = rho_L * u_L * v_L;
        f_vec[3] = rho_L * u_L * w_L;
        f_vec[4] = rho_L * ht_L * u_L;

        // G Vector
        g_vec[0] = rho_L * v_L;
        g_vec[1] = rho_L * v_L * u_L;
        g_vec[2] = rho_L * v_L * v_L + p_L;
        g_vec[3] = rho_L * v_L * w_L;
        g_vec[4] = rho_L * ht_L * v_L;

        // H Vector
        h_vec[0] = rho_L * w_L;
        h_vec[1] = rho_L * w_L * u_L;
        h_vec[2] = rho_L * w_L * v_L;
        h_vec[3] = rho_L * w_L * w_L + p_L;
        h_vec[4] = rho_L * ht_L * w_L;

        // flux_L
        flux_L[0] = f_vec[0] * nx + g_vec[0] * ny + h_vec[0] * nz;
        flux_L[1] = f_vec[1] * nx + g_vec[1] * ny + h_vec[1] * nz;
        flux_L[2] = f_vec[2] * nx + g_vec[2] * ny + h_vec[2] * nz;
        flux_L[3] = f_vec[3] * nx + g_vec[3] * ny + h_vec[3] * nz;
        flux_L[4] = f_vec[4] * nx + g_vec[4] * ny + h_vec[4] * nz;

        // Right Node
        rho_R = Q1[node_R];
        u_R   = Q2[node_R] / rho_R;
        v_R   = Q3[node_R] / rho_R;
        w_R   = Q4[node_R] / rho_R;
        et_R  = Q5[node_R] / rho_R;
        e_R   = et_R - 0.5 * (u_R * u_R + v_R * v_R + w_R * w_R);
        p_R   = (Gamma - 1.0) * rho_R * e_R;
        c_R   = sqrt((Gamma * p_R) / rho_R);
        h_R   = (c_R * c_R) / (Gamma - 1.0);
        ht_R  = h_R + 0.5 * (u_R * u_R + v_R * v_R + w_R * w_R);

        // F Vector
        f_vec[0] = rho_R * u_R;
        f_vec[1] = rho_R * u_R * u_R + p_R;
        f_vec[2] = rho_R * u_R * v_R;
        f_vec[3] = rho_R * u_R * w_R;
        f_vec[4] = rho_R * ht_R * u_R;

        // G Vector
        g_vec[0] = rho_R * v_R;
        g_vec[1] = rho_R * v_R * u_R;
        g_vec[2] = rho_R * v_R * v_R + p_R;
        g_vec[3] = rho_R * v_R * w_R;
        g_vec[4] = rho_R * ht_R * v_R;

        // H Vector
        h_vec[0] = rho_R * w_R;
        h_vec[1] = rho_R * w_R * u_R;
        h_vec[2] = rho_R * w_R * v_R;
        h_vec[3] = rho_R * w_R * w_R + p_R;
        h_vec[4] = rho_R * ht_R * w_R;

        // flux_R
        flux_R[0] = f_vec[0] * nx + g_vec[0] * ny + h_vec[0] * nz;
        flux_R[1] = f_vec[1] * nx + g_vec[1] * ny + h_vec[1] * nz;
        flux_R[2] = f_vec[2] * nx + g_vec[2] * ny + h_vec[2] * nz;
        flux_R[3] = f_vec[3] * nx + g_vec[3] * ny + h_vec[3] * nz;
        flux_R[4] = f_vec[4] * nx + g_vec[4] * ny + h_vec[4] * nz;

        // ROE AVERAGE VARIABLES
        temp = sqrt(rho_R) + sqrt(rho_L);
        rho  = sqrt(rho_R * rho_L);
        u    = (u_L * sqrt(rho_L) + u_R * sqrt(rho_R)) / temp;
        v    = (v_L * sqrt(rho_L) + v_R * sqrt(rho_R)) / temp;
        w    = (w_L * sqrt(rho_L) + w_R * sqrt(rho_R)) / temp;
        ht   = (ht_L * sqrt(rho_L) + ht_R * sqrt(rho_R)) / temp;
        c    = (Gamma - 1.0) * (ht - 0.5 * (u * u + v * v + w * w));
        c    = sqrt(c);
        phi  = 0.5 * (Gamma - 1.0)*(u * u + v * v + w * w);

        // M
        M[0][0] = 1.0;
        M[0][1] = 0.0;
        M[0][2] = 0.0;
        M[0][3] = 0.0;
        M[0][4] = 0.0;

        M[1][0] = u;
        M[1][1] = rho;
        M[1][2] = 0.0;
        M[1][3] = 0.0;
        M[1][4] = 0.0;

        M[2][0] = v;
        M[2][1] = 0.0;
        M[2][2] = rho;
        M[2][3] = 0.0;
        M[2][4] = 0.0;

        M[3][0] = w;
        M[3][1] = 0.0;
        M[3][2] = 0.0;
        M[3][3] = rho;
        M[3][4] = 0.0;

        M[4][0] = phi/(Gamma - 1.0);
        M[4][1] = rho * u;
        M[4][2] = rho * v;
        M[4][3] = rho * w;
        M[4][4] = 1.0/(Gamma - 1.0);

        // Minv
        Minv[0][0] = 1.0;
        Minv[0][1] = 0.0;
        Minv[0][2] = 0.0;
        Minv[0][3] = 0.0;
        Minv[0][4] = 0.0;

        Minv[1][0] = -u/rho;
        Minv[1][1] = 1.0/rho;
        Minv[1][2] = 0.0;
        Minv[1][3] = 0.0;
        Minv[1][4] = 0.0;

        Minv[2][0] = -v/rho;
        Minv[2][1] = 0.0;
        Minv[2][2] = 1.0/rho;
        Minv[2][3] = 0.0;
        Minv[2][4] = 0.0;

        Minv[3][0] = -w/rho;
        Minv[3][1] = 0.0;
        Minv[3][2] = 0.0;
        Minv[3][3] = 1.0/rho;
        Minv[3][4] = 0.0;

        Minv[4][0] = phi;
        Minv[4][1] = -u * (Gamma - 1.0);
        Minv[4][2] = -v * (Gamma - 1.0);
        Minv[4][3] = -w * (Gamma - 1.0);
        Minv[4][4] = (Gamma - 1.0);

        // P
        P[0][0] = nx;
        P[0][1] = ny;
        P[0][2] = nz;
        P[0][3] = rho/c;
        P[0][4] = rho/c;

        P[1][0] = 0.0;
        P[1][1] = -nz;
        P[1][2] = ny;
        P[1][3] = nx;
        P[1][4] = -nx;

        P[2][0] = nz;
        P[2][1] = 0.0;
        P[2][2] = -nx;
        P[2][3] = ny;
        P[2][4] = -ny;

        P[3][0] = -ny;
        P[3][1] = nx;
        P[3][2] = 0.0;
        P[3][3] = nz;
        P[3][4] = -nz;

        P[4][0] = 0.0;
        P[4][1] = 0.0;
        P[4][2] = 0.0;
        P[4][3] = rho * c;
        P[4][4] = rho * c;

        // Pinv
        Pinv[0][0] = nx;
        Pinv[0][1] = 0.0;
        Pinv[0][2] = nz;
        Pinv[0][3] = -ny;
        Pinv[0][4] = -nx/(c * c);

        Pinv[1][0] = ny;
        Pinv[1][1] = -nz;
        Pinv[1][2] = 0.0;
        Pinv[1][3] = nx;
        Pinv[1][4] = -ny/(c * c);

        Pinv[2][0] = nz;
        Pinv[2][1] = ny;
        Pinv[2][2] = -nx;
        Pinv[2][3] = 0.0;
        Pinv[2][4] = -nz/(c * c);

        Pinv[3][0] = 0.0;
        Pinv[3][1] = 0.5 * nx;
        Pinv[3][2] = 0.5 * ny;
        Pinv[3][3] = 0.5 * nz;
        Pinv[3][4] = 1.0/(2.0 * rho * c);

        Pinv[4][0] = 0.0;
        Pinv[4][1] = -0.5 * nx;
        Pinv[4][2] = -0.5 * ny;
        Pinv[4][3] = -0.5 * nz;
        Pinv[4][4] = 1.0/(2.0 * rho * c);

        // Calculate T = M*P
        MC_Matrix_Mul_Matrix(5, 5, M, P, T);
        // Calculate Tinv = Pinv*Minv
        MC_Matrix_Mul_Matrix(5, 5, Pinv, Minv, Tinv);

        // EigenMatrix
        for (j = 0; j < 5; j++)
            for (k = 0; k < 5; k++)
                Eigen[j][k] = 0.0;

        Eigen[0][0] = fabs(u * nx + v * ny + w * nz);
        Eigen[1][1] = fabs(u * nx + v * ny + w * nz);
        Eigen[2][2] = fabs(u * nx + v * ny + w * nz);
        Eigen[3][3] = fabs((u * nx + v * ny + w * nz) + c);
        Eigen[4][4] = fabs((u * nx + v * ny + w * nz) - c);
        
        // EigenMatrix*Tinv
        MC_Matrix_Mul_Matrix(5, 5, Eigen, Tinv, A);

        // Tinv = EigenMatrix*Tinv
        for (j = 0; j < 5; j++) {
            for (k = 0; k < 5; k++) {
                Tinv[j][k] = A[j][k];
                A[j][k] = 0.0;
            }
        }

        // Get Matrix A
        MC_Matrix_Mul_Matrix(5, 5, T, Tinv, A);
        
        // Set up matrix vector multiplication
        dQ[0] = dQ[1] = dQ[2] = dQ[3] = dQ[4] = 0.0;
        dQ[0] = Q1[node_R] - Q1[node_L];
        dQ[1] = Q2[node_R] - Q2[node_L];
        dQ[2] = Q3[node_R] - Q3[node_L];
        dQ[3] = Q4[node_R] - Q4[node_L];
        dQ[4] = Q5[node_R] - Q5[node_L];

        fluxA[0] = fluxA[1] = fluxA[2] = fluxA[3] = fluxA[4] = 0.0;
        MC_Matrix_Mul_Vector(5, 5, A, dQ, fluxA);
        
        // Compute LHS for Boundary Nodes
        Res1[node_L] += 0.5 * (flux_L[0] + flux_R[0] - fluxA[0]) * area;
        Res2[node_L] += 0.5 * (flux_L[1] + flux_R[1] - fluxA[1]) * area;
        Res3[node_L] += 0.5 * (flux_L[2] + flux_R[2] - fluxA[2]) * area;
        Res4[node_L] += 0.5 * (flux_L[3] + flux_R[3] - fluxA[3]) * area;
        Res5[node_L] += 0.5 * (flux_L[4] + flux_R[4] - fluxA[4]) * area;
    }

    // Free the Memory
    for (i = 0; i < 5; i++) {
        free(M[i]);
        free(Minv[i]);
        free(P[i]);
        free(Pinv[i]);
        free(T[i]);
        free(Tinv[i]);
        free(A[i]);
        free(Eigen[i]);
    }
    free(dQ);
    free(M);
    free(Minv);
    free(P);
    free(Pinv);
    free(T);
    free(Tinv);
    free(A);
    free(Eigen);
}

