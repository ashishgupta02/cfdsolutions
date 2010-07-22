/*******************************************************************************
 * File:        Roe_Fluxes.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

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
    int node_L, node_R;
    double rho_L, u_L, v_L, w_L, et_L, e_L, p_L, c_L, h_L, ht_L;
    double rho_R, u_R, v_R, w_R, et_R, e_R, p_R, c_R, h_R, ht_R;
    double rho, u, v, w, ht, c, phi;
    double temp;
    
    double dQ[5];
    double flux_L[5];
    double flux_R[5];
    
    Vector3D areavec;
    double   area;

    double f_vec[5];
    double g_vec[5];
    double h_vec[5];
    
    int i, j;
    double fluxA[5];
    double **MP;
    double **MPinv;
    double **A;

    MP    = (double **) malloc(5*sizeof(double*));
    MPinv = (double **) malloc(5*sizeof(double*));
    A     = (double **) malloc(5*sizeof(double*));
    for (i = 0; i < 5; i++) {
        MP[i]    = (double *) malloc(5*sizeof(double));
        MPinv[i] = (double *) malloc(5*sizeof(double));
        A[i]     = (double *) malloc(5*sizeof(double));
    }

    for (i = 0; i < 5; i++) {
        fluxA[i] = 0.0;
        flux_L[i] = 0.0;
        flux_R[i] = 0.0;
        f_vec[i] = 0.0;
        g_vec[i] = 0.0;
        h_vec[i] = 0.0;

        for (j = 0; j < 5; j++) {
            A[i][j] = 0.0;
            MP[i][j] = 0.0;
            MPinv[i][j] = 0.0;
        }
    }

    // Internal Edges
    for (i = 0; i < nEdge; i++) {
        // Get two nodes of edge
        node_L = int_edge_info[i].node[0];
        node_R = int_edge_info[i].node[1];

        // Get area vector
        areavec = int_edge_info[i].areav;
        area = areavec.magnitude();
        areavec.normalize();

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
        flux_L[0] = f_vec[0] * areavec.vec[0] + g_vec[0] * areavec.vec[1] + h_vec[0] * areavec.vec[2];
        flux_L[1] = f_vec[1] * areavec.vec[0] + g_vec[1] * areavec.vec[1] + h_vec[1] * areavec.vec[2];
        flux_L[2] = f_vec[2] * areavec.vec[0] + g_vec[2] * areavec.vec[1] + h_vec[2] * areavec.vec[2];
        flux_L[3] = f_vec[3] * areavec.vec[0] + g_vec[3] * areavec.vec[1] + h_vec[3] * areavec.vec[2];
        flux_L[4] = f_vec[4] * areavec.vec[0] + g_vec[4] * areavec.vec[1] + h_vec[4] * areavec.vec[2];

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
        flux_R[0] = f_vec[0] * areavec.vec[0] + g_vec[0] * areavec.vec[1] + h_vec[0] * areavec.vec[2];
        flux_R[1] = f_vec[1] * areavec.vec[0] + g_vec[1] * areavec.vec[1] + h_vec[1] * areavec.vec[2];
        flux_R[2] = f_vec[2] * areavec.vec[0] + g_vec[2] * areavec.vec[1] + h_vec[2] * areavec.vec[2];
        flux_R[3] = f_vec[3] * areavec.vec[0] + g_vec[3] * areavec.vec[1] + h_vec[3] * areavec.vec[2];
        flux_R[4] = f_vec[4] * areavec.vec[0] + g_vec[4] * areavec.vec[1] + h_vec[4] * areavec.vec[2];
        
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

        // MP
        MP[0][0] = areavec.vec[0];
        MP[0][1] = areavec.vec[1];
        MP[0][2] = areavec.vec[2];
        MP[0][3] = rho / c;
        MP[0][4] = rho / c;

        MP[1][0] =  areavec.vec[0] * u;
        MP[1][1] = -areavec.vec[2] * rho + areavec.vec[1] * u;
        MP[1][2] =  areavec.vec[1] * rho + areavec.vec[2] * u;
        MP[1][3] =  areavec.vec[0] * rho + (rho * u / c);
        MP[1][4] = -areavec.vec[0] * rho + (rho * u / c);

        MP[2][0] =  areavec.vec[2] * rho + areavec.vec[0] * v;
        MP[2][1] =  areavec.vec[1] * v;
        MP[2][2] = -areavec.vec[0] * rho + areavec.vec[2] * v;
        MP[2][3] =  areavec.vec[1] * rho + (rho * v / c);
        MP[2][4] = -areavec.vec[1] * rho + (rho * v / c);

        MP[3][0] = -areavec.vec[1] * rho + areavec.vec[0] * w;
        MP[3][1] =  areavec.vec[0] * rho + areavec.vec[1] * w;
        MP[3][2] =  areavec.vec[2] * w;
        MP[3][3] =  areavec.vec[2] * rho + (rho * w / c);
        MP[3][4] = -areavec.vec[2] * rho + (rho * w / c);


        MP[4][0] = (areavec.vec[0] * phi) / (Gamma - 1.0) 
                  + areavec.vec[2] * rho * v - areavec.vec[1] * rho * w;
        MP[4][1] = (areavec.vec[1] * phi) / (Gamma - 1.0)
                  - areavec.vec[2] * rho * u + areavec.vec[0] * rho * w;
        MP[4][2] = (areavec.vec[2] * phi) / (Gamma - 1.0) 
                  + areavec.vec[1] * rho * u - areavec.vec[0] * rho * v;
        MP[4][3] = (c * rho) / (Gamma - 1.0) + (phi * rho) / (c * (Gamma - 1.0)) 
                  + areavec.vec[0] * rho * u + areavec.vec[1] * rho * v + areavec.vec[2] * rho * w;
        MP[4][4] = (c * rho) / (Gamma - 1.0) + (phi * rho) / (c * (Gamma - 1.0)) 
                  - areavec.vec[0] * rho * u - areavec.vec[1] * rho * v - areavec.vec[2] * rho * w;

        // MP inverse
        MPinv[0][0] = areavec.vec[0] - areavec.vec[0] * phi / (c * c)
                    - (areavec.vec[2] * v) / rho + (areavec.vec[1] * w) / rho;
        MPinv[0][1] = (Gamma - 1.0) * areavec.vec[0] * u / (c * c);
        MPinv[0][2] = (areavec.vec[2] / rho) + (Gamma - 1.0) * areavec.vec[0] * v / (c * c);
        MPinv[0][3] = (-areavec.vec[1] / rho) + (Gamma - 1.0) * areavec.vec[0] * w / (c * c);
        MPinv[0][4] = -(Gamma - 1.0) * areavec.vec[0] / (c * c);

        MPinv[1][0] = areavec.vec[1] - areavec.vec[1] * phi / (c * c)
                     + areavec.vec[2] * u / rho - areavec.vec[0] * w / rho;
        MPinv[1][1] = -areavec.vec[2] / rho + (Gamma - 1.0) * areavec.vec[1] * u / (c * c);
        MPinv[1][2] = (Gamma - 1.0) * areavec.vec[1] * v / (c * c);
        MPinv[1][3] = areavec.vec[0] / rho + (Gamma - 1.0) * areavec.vec[1] * w / (c * c);
        MPinv[1][4] = -(Gamma - 1.0) * areavec.vec[1] / (c * c);

        MPinv[2][0] = areavec.vec[2] - areavec.vec[2] * phi / (c * c)
                     - areavec.vec[1] * u / rho + areavec.vec[0] * v / rho;
        MPinv[2][1] = areavec.vec[1] / rho + (Gamma - 1.0) * areavec.vec[2] * u / (c * c);
        MPinv[2][2] = -areavec.vec[0] / rho + (Gamma - 1.0) * areavec.vec[2] * v / (c * c);
        MPinv[2][3] = (Gamma - 1.0) * areavec.vec[2] * w / (c * c);
        MPinv[2][4] = -(Gamma - 1.0) * areavec.vec[2] / (c * c);

        MPinv[3][0] = phi / (2.0 * c * rho) - areavec.vec[0] * u / (2.0 * rho)
                     - areavec.vec[1] * v / (2.0 * rho) - areavec.vec[2] * w / (2.0 * rho);
        MPinv[3][1] = areavec.vec[0] / (2.0 * rho) - (Gamma - 1.0) * u / (2.0 * c * rho);
        MPinv[3][2] = areavec.vec[1] / (2.0 * rho) - (Gamma - 1.0) * v / (2.0 * c * rho);
        MPinv[3][3] = areavec.vec[2] / (2.0 * rho) - (Gamma - 1.0) * w / (2.0 * c * rho);
        MPinv[3][4] = (Gamma - 1.0) / (2.0 * c * rho);

        MPinv[4][0] = phi / (2.0 * c * rho) + (areavec.vec[0] * u) / (2.0 * rho)
                     + (areavec.vec[1] * v) / (2.0 * rho) + (areavec.vec[2] * w) / (2.0 * rho);
        MPinv[4][1] = -areavec.vec[0] / (2.0 * rho) - (Gamma - 1.0) * u / (2.0 * c * rho);
        MPinv[4][2] = -areavec.vec[1] / (2.0 * rho) - (Gamma - 1.0) * v / (2.0 * c * rho);
        MPinv[4][3] = -areavec.vec[2] / (2.0 * rho) - (Gamma - 1.0) * w / (2.0 * c * rho);
        MPinv[4][4] = (Gamma - 1.0) / (2.0 * c * rho);

        // Multiply Diagonal of T=MP by absolute values of Eigenmatrix
        MP[0][0] = MP[0][0] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2]);
        MP[1][1] = MP[1][1] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2]);
        MP[2][2] = MP[2][2] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2]);
        MP[3][3] = MP[3][3] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2] + c);
        MP[4][4] = MP[4][4] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2] - c);
        
        MC_Matrix_Mul_Matrix(5, 5, MP, MPinv, A);

        // Set up matrix vector multiplication
        dQ[0] = Q1[node_R] - Q1[node_L];
        dQ[1] = Q2[node_R] - Q2[node_L];
        dQ[2] = Q3[node_R] - Q3[node_L];
        dQ[3] = Q4[node_R] - Q4[node_L];
        dQ[4] = Q5[node_R] - Q5[node_L];

        MC_Matrix_Mul_Vector(5, 5, A, dQ, fluxA);

        // Compute for RHS
        Res1[node_L] -= 0.5 * (flux_L[0] + flux_R[0] - fluxA[0]) * area;
        Res2[node_L] -= 0.5 * (flux_L[1] + flux_R[1] - fluxA[1]) * area;
        Res3[node_L] -= 0.5 * (flux_L[2] + flux_R[2] - fluxA[2]) * area;
        Res4[node_L] -= 0.5 * (flux_L[3] + flux_R[3] - fluxA[3]) * area;
        Res5[node_L] -= 0.5 * (flux_L[4] + flux_R[4] - fluxA[4]) * area;

        Res1[node_R] += 0.5 * (flux_L[0] + flux_R[0] - fluxA[0]) * area;
        Res2[node_R] += 0.5 * (flux_L[1] + flux_R[1] - fluxA[1]) * area;
        Res3[node_R] += 0.5 * (flux_L[2] + flux_R[2] - fluxA[2]) * area;
        Res4[node_R] += 0.5 * (flux_L[3] + flux_R[3] - fluxA[3]) * area;
        Res5[node_R] += 0.5 * (flux_L[4] + flux_R[4] - fluxA[4]) * area;
    }

    for (i = 0; i < 5; i++) {
        flux_L[i] = 0.0;
        flux_R[i] = 0.0;
        fluxA[i]  = 0.0;
    }
    
    // Boundary Edges
    for (i = 0; i < nBEdge; i++) {
        // Get two nodes of edge
        node_L = bndry_edge_info[i].node[0];
        node_R = bndry_edge_info[i].node[1];

        // Get area vector
        areavec = bndry_edge_info[i].areav;
        area = areavec.magnitude();
        areavec.normalize();
        
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
        flux_L[0] = f_vec[0] * areavec.vec[0] + g_vec[0] * areavec.vec[1] + h_vec[0] * areavec.vec[2];
        flux_L[1] = f_vec[1] * areavec.vec[0] + g_vec[1] * areavec.vec[1] + h_vec[1] * areavec.vec[2];
        flux_L[2] = f_vec[2] * areavec.vec[0] + g_vec[2] * areavec.vec[1] + h_vec[2] * areavec.vec[2];
        flux_L[3] = f_vec[3] * areavec.vec[0] + g_vec[3] * areavec.vec[1] + h_vec[3] * areavec.vec[2];
        flux_L[4] = f_vec[4] * areavec.vec[0] + g_vec[4] * areavec.vec[1] + h_vec[4] * areavec.vec[2];

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
        ht_R  = h_R + 0.5 * (u_R * u_R + v_R * v_R + w_R*w_R);

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
        flux_R[0] = f_vec[0] * areavec.vec[0] + g_vec[0] * areavec.vec[1] + h_vec[0] * areavec.vec[2];
        flux_R[1] = f_vec[1] * areavec.vec[0] + g_vec[1] * areavec.vec[1] + h_vec[1] * areavec.vec[2];
        flux_R[2] = f_vec[2] * areavec.vec[0] + g_vec[2] * areavec.vec[1] + h_vec[2] * areavec.vec[2];
        flux_R[3] = f_vec[3] * areavec.vec[0] + g_vec[3] * areavec.vec[1] + h_vec[3] * areavec.vec[2];
        flux_R[4] = f_vec[4] * areavec.vec[0] + g_vec[4] * areavec.vec[1] + h_vec[4] * areavec.vec[2];

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

        // Compute MP
        MP[0][0] = areavec.vec[0];
        MP[0][1] = areavec.vec[1];
        MP[0][2] = areavec.vec[2];
        MP[0][3] = rho / c;
        MP[0][4] = rho / c;

        MP[1][0] =  areavec.vec[0] * u;
        MP[1][1] = -areavec.vec[2] * rho + areavec.vec[1] * u;
        MP[1][2] =  areavec.vec[1] * rho + areavec.vec[2] * u;
        MP[1][3] =  areavec.vec[0] * rho + (rho * u / c);
        MP[1][4] = -areavec.vec[0] * rho + (rho * u / c);

        MP[2][0] =  areavec.vec[2] * rho + areavec.vec[0] * v;
        MP[2][1] =  areavec.vec[1] * v;
        MP[2][2] = -areavec.vec[0] * rho + areavec.vec[2] * v;
        MP[2][3] =  areavec.vec[1] * rho + (rho * v / c);
        MP[2][4] = -areavec.vec[1] * rho + (rho * v / c);

        MP[3][0] = -areavec.vec[1] * rho + areavec.vec[0] * w;
        MP[3][1] =  areavec.vec[0] * rho + areavec.vec[1] * w;
        MP[3][2] =  areavec.vec[2] * w;
        MP[3][3] =  areavec.vec[2] * rho + (rho * w / c);
        MP[3][4] = -areavec.vec[2] * rho + (rho * w / c);

        MP[4][0] = (areavec.vec[0] * phi) / (Gamma - 1.0)
                  + areavec.vec[2] * rho * v - areavec.vec[1] * rho * w;
        MP[4][1] = (areavec.vec[1] * phi) / (Gamma - 1.0)
                  - areavec.vec[2] * rho * u + areavec.vec[0] * rho * w;
        MP[4][2] = (areavec.vec[2] * phi) / (Gamma - 1.0)
                  + areavec.vec[1] * rho * u - areavec.vec[0] * rho * v;
        MP[4][3] = (c * rho) / (Gamma - 1.0) + (phi * rho) / (c * (Gamma - 1.0))
                  + areavec.vec[0] * rho * u + areavec.vec[1] * rho * v + areavec.vec[2] * rho * w;
        MP[4][4] = (c * rho) / (Gamma - 1.0) + (phi * rho) / (c * (Gamma - 1.0))
                  - areavec.vec[0] * rho * u - areavec.vec[1] * rho * v - areavec.vec[2] * rho * w;

        // Compute MP inverse
        MPinv[0][0] = areavec.vec[0] - areavec.vec[0] * phi / (c * c)
                     - (areavec.vec[2] * v) / rho + (areavec.vec[1] * w) / rho;
        MPinv[0][1] = (Gamma - 1.0) * areavec.vec[0] * u / (c * c);
        MPinv[0][2] = (areavec.vec[2] / rho) + (Gamma - 1.0) * areavec.vec[0] * v / (c * c);
        MPinv[0][3] = (-areavec.vec[1] / rho) + (Gamma - 1.0) * areavec.vec[0] * w / (c * c);
        MPinv[0][4] = -(Gamma - 1.0) * areavec.vec[0] / (c * c);

        MPinv[1][0] = areavec.vec[1] - areavec.vec[1] * phi / (c * c)
                     + areavec.vec[2] * u / rho - areavec.vec[0] * w / rho;
        MPinv[1][1] = -areavec.vec[2] / rho + (Gamma - 1.0) * areavec.vec[1] * u / (c * c);
        MPinv[1][2] = (Gamma - 1.0) * areavec.vec[1] * v / (c * c);
        MPinv[1][3] = areavec.vec[0] / rho + (Gamma - 1.0) * areavec.vec[1] * w / (c * c);
        MPinv[1][4] = -(Gamma - 1.0) * areavec.vec[1] / (c * c);

        MPinv[2][0] = areavec.vec[2] - areavec.vec[2] * phi / (c * c)
                     - areavec.vec[1] * u / rho + areavec.vec[0] * v / rho;
        MPinv[2][1] = areavec.vec[1] / rho + (Gamma - 1.0) * areavec.vec[2] * u / (c * c);
        MPinv[2][2] = -areavec.vec[0] / rho + (Gamma - 1.0) * areavec.vec[2] * v / (c * c);
        MPinv[2][3] = (Gamma - 1.0) * areavec.vec[2] * w / (c * c);
        MPinv[2][4] = -(Gamma - 1.0) * areavec.vec[2] / (c * c);

        MPinv[3][0] = phi / (2.0 * c * rho) - areavec.vec[0] * u / (2.0 * rho)
                     - areavec.vec[1] * v / (2.0 * rho) - areavec.vec[2] * w / (2.0 * rho);
        MPinv[3][1] = areavec.vec[0] / (2.0 * rho) - (Gamma - 1.0) * u / (2.0 * c * rho);
        MPinv[3][2] = areavec.vec[1] / (2.0 * rho) - (Gamma - 1.0) * v / (2.0 * c * rho);
        MPinv[3][3] = areavec.vec[2] / (2.0 * rho) - (Gamma - 1.0) * w / (2.0 * c * rho);
        MPinv[3][4] = (Gamma - 1.0) / (2.0 * c * rho);

        MPinv[4][0] = phi / (2.0 * c * rho) + (areavec.vec[0] * u) / (2.0 * rho)
                     + (areavec.vec[1] * v) / (2.0 * rho) + (areavec.vec[2] * w) / (2.0 * rho);
        MPinv[4][1] = -areavec.vec[0] / (2.0 * rho) - (Gamma - 1.0) * u / (2.0 * c * rho);
        MPinv[4][2] = -areavec.vec[1] / (2.0 * rho) - (Gamma - 1.0) * v / (2.0 * c * rho);
        MPinv[4][3] = -areavec.vec[2] / (2.0 * rho) - (Gamma - 1.0) * w / (2.0 * c * rho);
        MPinv[4][4] = (Gamma - 1.0) / (2.0 * c * rho);

        //Multiply Diagonal of T=MP by absolute values of Eigenmatrix
        MP[0][0] = MP[0][0] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2]);
        MP[1][1] = MP[1][1] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2]);
        MP[2][2] = MP[2][2] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2]);
        MP[3][3] = MP[3][3] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2] + c);
        MP[4][4] = MP[4][4] * fabs(u * areavec.vec[0] + v * areavec.vec[1] + w * areavec.vec[2] - c);
        
        MC_Matrix_Mul_Matrix(5, 5, MP, MPinv, A);

        // Set up matrix vector multiplication
        dQ[0] = Q1[node_R] - Q1[node_L];
        dQ[1] = Q2[node_R] - Q2[node_L];
        dQ[2] = Q3[node_R] - Q3[node_L];
        dQ[3] = Q4[node_R] - Q4[node_L];
        dQ[4] = Q5[node_R] - Q5[node_L];
        
        MC_Matrix_Mul_Vector(5, 5, A, dQ, fluxA);

        // Compute RHS for Boundary Nodes
        Res1[node_L] -= 0.5 * (flux_L[0] + flux_R[0] - fluxA[0]) * area;
        Res2[node_L] -= 0.5 * (flux_L[1] + flux_R[1] - fluxA[1]) * area;
        Res3[node_L] -= 0.5 * (flux_L[2] + flux_R[2] - fluxA[2]) * area;
        Res4[node_L] -= 0.5 * (flux_L[3] + flux_R[3] - fluxA[3]) * area;
        Res5[node_L] -= 0.5 * (flux_L[4] + flux_R[4] - fluxA[4]) * area;
    }
}

