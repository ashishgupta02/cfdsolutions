/*******************************************************************************
 * File:        BC.cpp
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
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Initialize_Boundary_Condition() {
    // Get the Boundary Type
    for (int i = 0; i < nBEdge; i++)
        bndEdge[i].type = bndType[bndEdge[i].tag - 1];
}

//------------------------------------------------------------------------------
//! This Routine is created to assist Finite Difference Jacobian Computations
//! Note: Outward Pointing Boundary Normals
//------------------------------------------------------------------------------
void Compute_Characteristic_BoundaryCondition(int BEdgeID, int Iteration) {
    int node_P, node_G;
    Vector3D normal;
    double nx, ny, nz;
    double rho_b, u_b, v_b, w_b, et_b, e_b, p_b, c_b;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double rho_i, u_i, v_i, w_i, et_i, e_i, p_i, c_i;
    double rho_0, u_0, v_0, w_0, et_0, e_0, p_0, c_0, ubar_0;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    double eps, alpha, beta, Ur, mach;
    double dtmp;
    int instance = BC_NONE;
    
    // Default Initialization
    p_b_n   = Inf_Pressure + Gauge_Pressure;
    rho_b_n = Inf_Rho;
    u_b_n   = Inf_U;
    v_b_n   = Inf_V;
    w_b_n   = Inf_W;

    node_P = bndEdge[BEdgeID].node[0];
    node_G = bndEdge[BEdgeID].node[1];

#ifdef DEBUG
    assert(node_G > node_P);
#endif

    // Get area vector
    normal = bndEdge[BEdgeID].areav;
    normal.normalize();
    nx = normal.vec[0];
    ny = normal.vec[1];
    nz = normal.vec[2];

    // Physical node
    rho_i = Q1[node_P];
    u_i   = Q2[node_P] / rho_i;
    v_i   = Q3[node_P] / rho_i;
    w_i   = Q4[node_P] / rho_i;
    et_i  = Q5[node_P] / rho_i;
    e_i   = et_i - 0.5 * (u_i * u_i + v_i * v_i + w_i * w_i);
    p_i   = rho_i * e_i * (Gamma - 1.0) + Gauge_Pressure;
    c_i   = sqrt((Gamma * p_i) / rho_i);

    // Ghost node
    rho_b = Q1[node_G];
    u_b   = Q2[node_G] / rho_b;
    v_b   = Q3[node_G] / rho_b;
    w_b   = Q4[node_G] / rho_b;
    et_b  = Q5[node_G] / rho_b;
    e_b   = et_b - 0.5 * (u_b * u_b + v_b * v_b + w_b * w_b);
    p_b   = rho_b * e_b * (Gamma - 1.0) + Gauge_Pressure;
    c_b   = sqrt((Gamma * p_b) / rho_b);

    // Compute the average
    rho_0  = 0.5 * (Q1[node_P] + Q1[node_G]);
    u_0    = 0.5 * (Q2[node_P] + Q2[node_G])/rho_0;
    v_0    = 0.5 * (Q3[node_P] + Q3[node_G])/rho_0;
    w_0    = 0.5 * (Q4[node_P] + Q4[node_G])/rho_0;
    et_0   = 0.5 * (Q5[node_P] + Q5[node_G])/rho_0;
    e_0    = et_0 - 0.5 * (u_0 * u_0 + v_0 * v_0 + w_0 * w_0);
    p_0    = rho_0 * e_0 * (Gamma - 1.0) + Gauge_Pressure;
    c_0    = sqrt((Gamma * p_0) / rho_0);
    ubar_0 = u_0 * nx + v_0 * ny + w_0 * nz;

    // Get eigenvalues
    // Compute Weiss Smith Precondition Variables
    if (SolverScheme == SOLVER_SCHEME_ROE_WS) {
        // Compute Precondition Variables
        dtmp = sqrt(u_0 * u_0 + v_0 * v_0 + w_0 * w_0);
        mach = dtmp/c_0;
        eps  = MIN(1.0, MAX(1e-10, mach*mach));

        // Compute beta for Ideal Gas
        beta   = 1.0/(c_0*c_0);

        // Compute Ur for Ideal Gas
        if (dtmp < eps*c_0)
            Ur = eps*c_0;
        else if ((eps*c_0 < dtmp) && (dtmp < c_0))
            Ur = dtmp;
        else
            Ur = c_0;

        // Compute alpha
        alpha  = 0.5*(1.0 - beta*Ur*Ur);

        // Compute U'
        ubar_0 = ubar_0*(1.0 - alpha);

        // Compute C'
        c_0 = sqrt(alpha*alpha*ubar_0*ubar_0 + Ur*Ur);
    }
    // Finally Compute Pre/Non-Condition Eigenvalues
    lamda1 = ubar_0;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = lamda1 + c_0;
    lamda5 = lamda1 - c_0;

    // Determine how to set boundary conditions for edge
    switch (bndEdge[BEdgeID].type) {
        // Solid Wall
        case BC_SOLID_WALL:
            instance = BC_SOLID_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c_0);
                    info("Nodes: %d %d", node_P, node_G);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Apply_Boundary_Condition: Unable apply Solid Wall boundary condition edge[%d] - 1", BEdgeID);
                }
            }
            break;
        // Free Stream
        case BC_FREE_STREAM:
            instance = BC_FREE_STREAM;
            // Supersonic inflow
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_SUPERSONIC_INFLOW;
                // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_SUPERSONIC_OUTFLOW;
                // Subsonic outflow
            else if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_SUBSONIC_OUTFLOW;
                // Subsonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Apply_Boundary_Condition: Unable apply boundary condition edge[%d] - 2", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c_0);
            info("Nodes: %d %d", node_P, node_G);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Apply_Boundary_Condition: Unable apply boundary condition edge[%d] - 3", BEdgeID);
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_SOLID_WALL:
            temp    = rho_0 * c_0* (nx * u_i + ny * v_i + nz * w_i);
            p_b_n = p_i + temp;
            // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
            // This causes the flow to be allowed through the surface initially
            if (Iteration < ZPGIteration) {
                if (p_b_n < 0.0)
                    p_b_n = p_i;
            }   
            rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            break;
        case BC_SUBSONIC_OUTFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            break;
        case BC_SUBSONIC_INFLOW:
            temp    = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
            p_b_n   = 0.5 * ((Inf_Pressure + Gauge_Pressure) + p_i - rho_0 * c_0 * temp);
            rho_b_n = Inf_Rho + (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (c_0 * c_0);
            u_b_n   = Inf_U + nx * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            v_b_n   = Inf_V + ny * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            w_b_n   = Inf_W + nz * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            break;
        case BC_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_FREE_STREAM:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Apply_Boundary_Condition: Unable to set BC Type edge[%d]", BEdgeID);
            break;
    }

    if (p_b_n < 0.0)
        error("Negative Pressure: Node: %d Ghost: %d, P_i: %e, P_b_n: %e", node_P, node_G, p_i, p_b_n);

    et_b_n = (p_b_n - Gauge_Pressure)/ ((Gamma - 1.0) * rho_b_n) + 0.5 * (u_b_n * u_b_n + v_b_n * v_b_n + w_b_n * w_b_n);

    // Update the Q's for Ghost Node
    Q1[node_G] = rho_b_n;
    Q2[node_G] = rho_b_n * u_b_n;
    Q3[node_G] = rho_b_n * v_b_n;
    Q4[node_G] = rho_b_n * w_b_n;
    Q5[node_G] = rho_b_n * et_b_n;
}

//------------------------------------------------------------------------------
//! This Routine is created to assist Finite Difference Jacobian Computations
//! Note: Outward Pointing Boundary Normals
//! rvalue: 0 => Success, 1: Negative Pressure
//------------------------------------------------------------------------------
int Compute_Characteristic_BoundaryCondition_New(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    double nx, ny, nz;
    double rho_b, u_b, v_b, w_b, et_b, e_b, p_b, c_b;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double rho_i, u_i, v_i, w_i, et_i, e_i, p_i, c_i;
    double rho_0, u_0, v_0, w_0, et_0, e_0, p_0, c_0, ubar_0;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    double eps, alpha, beta, Ur, mach, sigma;
    double ubar_pre, c_pre;
    double c1, c2, c3;
    double dtmp;
    int instance = BC_NONE;
    int rvalue = 0;
    int node_L, node_R;
    
    // Default Initialization
    p_b_n   = Inf_Pressure + Gauge_Pressure;
    rho_b_n = Inf_Rho;
    u_b_n   = Inf_U;
    v_b_n   = Inf_V;
    w_b_n   = Inf_W;
    
    // Get area vector
    AreaVec.normalize();
    nx = AreaVec.vec[0];
    ny = AreaVec.vec[1];
    nz = AreaVec.vec[2];

    // Physical node
    rho_i = Q_L[0];
    u_i   = Q_L[1] / rho_i;
    v_i   = Q_L[2] / rho_i;
    w_i   = Q_L[3] / rho_i;
    et_i  = Q_L[4] / rho_i;
    e_i   = et_i - 0.5 * (u_i * u_i + v_i * v_i + w_i * w_i);
    p_i   = rho_i * e_i * (Gamma - 1.0) + Gauge_Pressure;
    c_i   = sqrt((Gamma * p_i) / rho_i);

    // Ghost node
    rho_b = Q_R[0];
    u_b   = Q_R[1] / rho_b;
    v_b   = Q_R[2] / rho_b;
    w_b   = Q_R[3] / rho_b;
    et_b  = Q_R[4] / rho_b;
    e_b   = et_b - 0.5 * (u_b * u_b + v_b * v_b + w_b * w_b);
    p_b   = rho_b * e_b * (Gamma - 1.0) + Gauge_Pressure;
    c_b   = sqrt((Gamma * p_b) / rho_b);

    // Compute the average
    rho_0 = 0.5 * (Q_L[0] + Q_R[0]);
    u_0   = 0.5 * (Q_L[1] + Q_R[1])/rho_0;
    v_0   = 0.5 * (Q_L[2] + Q_R[2])/rho_0;
    w_0   = 0.5 * (Q_L[3] + Q_R[3])/rho_0;
    et_0  = 0.5 * (Q_L[4] + Q_R[4])/rho_0;
    e_0   = et_0 - 0.5 * (u_0 * u_0 + v_0 * v_0 + w_0 * w_0);
    p_0   = rho_0 * e_0 * (Gamma - 1.0) + Gauge_Pressure;
    c_0   = sqrt((Gamma * p_0) / rho_0);
    ubar_0 = u_0 * nx + v_0 * ny + w_0 * nz;
    
    // Get eigenvalues
    // Compute Weiss Smith Precondition Variables
    if (SolverScheme == SOLVER_SCHEME_ROE_WS) {
        // Compute Precondition Variables
        dtmp = sqrt(u_0 * u_0 + v_0 * v_0 + w_0 * w_0);
        mach = dtmp/c_0;
        eps  = MIN(1.0, MAX(1e-10, mach*mach));

        // Compute beta for Ideal Gas
        beta   = 1.0/(c_0*c_0);

        // Compute Ur for Ideal Gas
        if (dtmp < eps*c_0)
            Ur = eps*c_0;
        else if ((eps*c_0 < dtmp) && (dtmp < c_0))
            Ur = dtmp;
        else
            Ur = c_0;

        // Compute alpha
        alpha  = 0.5*(1.0 - beta*Ur*Ur);

        // Compute ubar_pre
        ubar_pre = ubar_0*(1.0 - alpha);

        // Compute c_pre
        c_pre = sqrt(alpha*alpha*ubar_pre*ubar_pre + Ur*Ur);
    }
    
    // Compute Cecile Voizat Precondition Variables
    if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
        // Get two nodes of edge
        node_L = bndEdge[BEdgeID].node[0];
        node_R = bndEdge[BEdgeID].node[1];
        
        double q2, lc, lp, dp_L, dp, lmach_L, lmach;
        int i, nid;
        
        q2 = u_0*u_0 + v_0*v_0 + w_0*w_0;
        dp_L    = 0.0;
        lmach_L = 0.0;
        for (i = crs_IA_Node2Node[node_L]; i < crs_IA_Node2Node[node_L+1]; i++) {
            nid     = crs_JA_Node2Node[i];
            lp      = (Gamma - 1.0)*(Q5[nid] - 0.5*(Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/Q1[nid]) + Gauge_Pressure;
            lc      = sqrt((Gamma * lp) / Q1[nid]);
            lmach_L = MAX(lmach_L, sqrt((Q2[nid]*Q2[nid] + Q3[nid]*Q3[nid] + Q4[nid]*Q4[nid])/(Q1[nid]*Q1[nid]))/lc);
            dp_L    = MAX(dp_L, fabs(p_i - lp));
        }
        dp    = dp_L;
        lmach = lmach_L;
        
        // Compute Precondition Variables
        dtmp  = sqrt(q2);
        mach  = dtmp/c_0;
        sigma = MIN(1.0, MAX(MAX(MAX(1.0e-5, mach), lmach), dp/(rho_0*c_0*c_0)));
        if (mach > 0.3)
            sigma = 1.0;
        //sigma = MIN(1.0, MAX(Inf_Mach, mach));
        
        // alpha, beta, c1, c2, c3
        alpha = ubar_0*ubar_0*(sigma*sigma*sigma*sigma - 2.0*sigma*sigma + 1.0) + 4.0*c_0*c_0*sigma*sigma;
        alpha = sqrt(alpha);
        beta  = (sigma*sigma - 1.0)*ubar_0;
        c1    = 1.0/(rho_0*alpha);
        c2    = (alpha - beta)/(2.0*alpha);
        c3    = (alpha + beta)/(2.0*alpha);
        
        // Compute c_pre
        c_pre = 0.5*sqrt(4.0*c_0*c_0*sigma*sigma + (1.0 - sigma*sigma)*(1.0 - sigma*sigma)*ubar_0*ubar_0);

        // Compute ubar_pre
        ubar_pre = 0.5*(1.0 + sigma*sigma)*ubar_0;
    }
    
    // No Precondition
    if (SolverScheme == SOLVER_SCHEME_ROE) {
        ubar_pre = ubar_0;
        c_pre    = c_0;
    }
    
    // Finally Compute Pre/Non-Condition Eigenvalues
    lamda1 = ubar_pre;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = lamda1 + c_pre;
    lamda5 = lamda1 - c_pre;

    // Determine how to set boundary conditions for edge
    switch (BCType) {
        // Solid Wall
        case BC_SOLID_WALL:
            instance = BC_SOLID_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c_0);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Compute_Characteristic_BoundaryCondition: Unable apply Solid Wall boundary condition - 1");
                }
            }
            break;
        // Free Stream
        case BC_FREE_STREAM:
            instance = BC_FREE_STREAM;
            // Supersonic inflow
            if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                instance = BC_SUPERSONIC_INFLOW;
                // Supersonic outflow
            else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                instance = BC_SUPERSONIC_OUTFLOW;
                // Subsonic outflow
            else if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_SUBSONIC_OUTFLOW;
                // Subsonic inflow
            else if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Compute_Characteristic_BoundaryCondition: Unable apply boundary condition - 2");
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c_0);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Characteristic_BoundaryCondition: Unable apply boundary condition - 3");
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_SOLID_WALL:
            // No Preconditioning 
            if (SolverScheme == SOLVER_SCHEME_ROE) {
                temp    = rho_0 * c_0* (nx * u_i + ny * v_i + nz * w_i);
                p_b_n = p_i + temp;
                // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
                // This causes the flow to be allowed through the surface initially
                if (Iteration < ZPGIteration) {
                    if (p_b_n < 0.0)
                        p_b_n = p_i;
                }   
                rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
                u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
                v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
                w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            }
            
            // Weiss Smith Precondition
            if (SolverScheme == SOLVER_SCHEME_ROE_WS) {
                // TODO
            }
            
            // Cecile Voizat Precondition
            if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
                temp = (c2/c1)*(nx * u_i + ny * v_i + nz * w_i);
                p_b_n   = p_i + temp;
                // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
                // This causes the flow to be allowed through the surface initially
                if (Iteration < ZPGIteration) {
                    if (p_b_n < 0.0)
                        p_b_n = p_i;
                }
                rho_b_n = rho_i * pow((p_b_n/p_i), (1.0/Gamma));
                u_b_n   = u_i - nx * (c1/c2)*(p_b_n - p_i);
                v_b_n   = v_i - ny * (c1/c2)*(p_b_n - p_i);
                w_b_n   = w_i - nz * (c1/c2)*(p_b_n - p_i);
            }
            break;
        case BC_SUBSONIC_OUTFLOW:
            // No Preconditioning 
            if (SolverScheme == SOLVER_SCHEME_ROE) {
                p_b_n   = Inf_Pressure + Gauge_Pressure;
                rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
                u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
                v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
                w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            }
            
            // Weiss Smith Precondition
            if (SolverScheme == SOLVER_SCHEME_ROE_WS) {
                // TODO
            }
            
            // Cecile Voizat Precondition
            if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
                temp  = nx * (u_i - Inf_U) + ny * (v_i - Inf_V) + nz * (w_i - Inf_W);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c2*(Inf_Pressure + Gauge_Pressure) + c3*p_i) + c2*c3*temp);
                rho_b_n = rho_i * pow((p_b_n/p_i), (1.0/Gamma));
                u_b_n   = u_i - nx * (c1/c2)*(p_b_n - p_i);
                v_b_n   = v_i - ny * (c1/c2)*(p_b_n - p_i);
                w_b_n   = w_i - nz * (c1/c2)*(p_b_n - p_i);
            }
            break;
        case BC_SUBSONIC_INFLOW:
            // No Preconditioning 
            if (SolverScheme == SOLVER_SCHEME_ROE) {
                temp    = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
                p_b_n   = 0.5 * ((Inf_Pressure + Gauge_Pressure) + p_i - rho_0 * c_0 * temp);
                rho_b_n = Inf_Rho + (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (c_0 * c_0);
                u_b_n   = Inf_U + nx * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
                v_b_n   = Inf_V + ny * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
                w_b_n   = Inf_W + nz * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            }
            
            // Weiss Smith Precondition
            if (SolverScheme == SOLVER_SCHEME_ROE_WS) {
                // TODO
            }
            
            // Cecile Voizat Precondition
            if (SolverScheme == SOLVER_SCHEME_ROE_CV) {
                temp  = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c3*(Inf_Pressure + Gauge_Pressure) + c2*p_i) - c2*c3*temp);
                rho_b_n = Inf_Rho * pow((p_b_n/(Inf_Pressure + Gauge_Pressure)), (1.0/Gamma));
                u_b_n   = Inf_U + nx * (c1/c2)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                v_b_n   = Inf_V + ny * (c1/c2)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                w_b_n   = Inf_W + nz * (c1/c2)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
            }
            break;
        case BC_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        case BC_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_FREE_STREAM:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Characteristic_BoundaryCondition: Unable to set BC Type");
            break;
    }

    if (p_b_n < 0.0) {
        warn("Compute_Characteristic_BoundaryCondition: Negative Pressure Detected P_i: %e, P_b_n: %e", p_i, p_b_n);
        rvalue = 1;
    }

    et_b_n = (p_b_n - Gauge_Pressure) / ((Gamma - 1.0) * rho_b_n) + 0.5 * (u_b_n * u_b_n + v_b_n * v_b_n + w_b_n * w_b_n);

    // Update the Q's for Ghost Node
    Q_B[0] = rho_b_n;
    Q_B[1] = rho_b_n * u_b_n;
    Q_B[2] = rho_b_n * v_b_n;
    Q_B[3] = rho_b_n * w_b_n;
    Q_B[4] = rho_b_n * et_b_n;
    
    return rvalue;
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
        
        rvalue = Compute_Characteristic_BoundaryCondition_New(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
        if (rvalue != 0)
            error("Apply_Characteristic_Boundary_Condition: BndNode: %d,  GhostNode: %d", node_L, node_R);
        
        // Update the Boundary Value of Ghost Node
        Q1[node_R] = Q_B[0];
        Q2[node_R] = Q_B[1];
        Q3[node_R] = Q_B[2];
        Q4[node_R] = Q_B[3];
        Q5[node_R] = Q_B[4];
    }
}

//------------------------------------------------------------------------------
//! Note: Outward Pointing Boundary Normals
//------------------------------------------------------------------------------
void Compute_Pressure_BoundaryCondition(int BEdgeID, int Iteration) {
    
}

//------------------------------------------------------------------------------
//! Note: Outward Pointing Boundary Normals
//------------------------------------------------------------------------------
void Apply_Pressure_Boundary_Condition(int Iteration) {
    
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
void Apply_Boundary_Condition(int Iteration) {
    switch (SolverBCScheme) {
        case SOLVER_BC_SCHEME_CHARCTERISTIC: // Characteristic Boundary Condition
            Apply_Characteristic_Boundary_Condition(Iteration);
            break;
        case SOLVER_BC_SCHEME_PRESSURE: // Pressure Based Boundary Condition
            Apply_Pressure_Boundary_Condition(Iteration);
            break;
        default:
            error("Apply_Boundary_Condition: Invalid Boundary Condition Scheme - %d", SolverBCScheme);
            break;
    }
}

