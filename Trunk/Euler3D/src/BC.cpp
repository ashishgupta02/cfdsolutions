/*******************************************************************************
 * File:        BC.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <assert.h>

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
//!
//------------------------------------------------------------------------------
void Apply_Boundary_Condition() {
    int pnode;
    int ghostnode;
    Vector3D normal;
    double nx, ny, nz;
    double rho_b, u_b, v_b, w_b, et_b, e_b, p_b, c_b;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double rho_i, u_i, v_i, w_i, et_i, e_i, p_i, c_i;
    double rho_0, u_0, v_0, w_0, et_0, e_0, p_0, c_0;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    int instance = BC_NONE;

    // Default Initialization
    p_b_n   = Inf_Pressure;
    rho_b_n = Inf_Rho;
    u_b_n   = Inf_U;
    v_b_n   = Inf_V;
    w_b_n   = Inf_W;

    for (int i = 0; i < nBEdge; i++) {
        pnode     = bndEdge[i].node[0];
        ghostnode = bndEdge[i].node[1];

        assert(ghostnode > pnode);
        
        // Get area vector
        normal = bndEdge[i].areav;
        normal.normalize();
        nx = normal.vec[0];
        ny = normal.vec[1];
        nz = normal.vec[2];
        
        // Physical node
        rho_i = Q1[pnode];
        u_i   = Q2[pnode] / rho_i;
        v_i   = Q3[pnode] / rho_i;
        w_i   = Q4[pnode] / rho_i;
        et_i  = Q5[pnode] / rho_i;
        e_i   = et_i - 0.5 * (u_i * u_i + v_i * v_i + w_i * w_i);
        p_i   = rho_i * e_i * (Gamma - 1.0);
        c_i   = sqrt((Gamma * p_i) / rho_i);
        
        // Ghost node
        rho_b = Q1[ghostnode];
        u_b   = Q2[ghostnode] / rho_b;
        v_b   = Q3[ghostnode] / rho_b;
        w_b   = Q4[ghostnode] / rho_b;
        et_b  = Q5[ghostnode] / rho_b;
        e_b   = et_b - 0.5 * (u_b * u_b + v_b * v_b + w_b * w_b);
        p_b   = rho_b * e_b * (Gamma - 1.0);
        c_b   = sqrt((Gamma * p_b) / rho_b);

        rho_0 = 0.5 * (rho_i + rho_b);
        u_0   = 0.5 * (u_i + u_b);
        v_0   = 0.5 * (v_i + v_b);
        w_0   = 0.5 * (w_i + w_b);
        et_0  = 0.5 * (et_i + et_b);
        e_0   = 0.5 * (e_i + e_b);
        p_0   = rho_0 * e_0 *(Gamma - 1.0);
        c_0   = sqrt(Gamma * p_0/rho_0);
        
        
        // Get eigenvalues
        lamda1 = u_i * nx + v_i * ny + w_i * nz;
        lamda2 = lamda1;
        lamda3 = lamda1;
        lamda4 = lamda1 + c_i;
        lamda5 = lamda1 - c_i;
        
        // Determine how to set boundary conditions for edge
        switch (bndEdge[i].type) {
            // Solid Wall
            case BC_SOLID_WALL:
                instance = BC_SOLID_WALL;
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lamda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Nodes: %d %d", pnode, ghostnode);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    error("Apply_Boundary_Condition: Unable apply Solid Wall boundary condition edge[%d]", i);
                }
                break;
            // Free Stream
            case BC_FREE_STREAM:
                instance = BC_FREE_STREAM;
                if (fabs(lamda1)/ c_0 >= 1.0) {
                    // Supersonic inflow
                    if ((lamda1 <= 0.0) && (lamda4 <= 0.0) && (lamda5 <= 0.0))
                        instance = BC_SUPERSONIC_INFLOW;
                    // Supersonic outflow
                    else if ((lamda1 >= 0.0) && (lamda4 >= 0.0) && (lamda5 >= 0.0))
                        instance = BC_SUPERSONIC_OUTFLOW;
                    else {
                        info(" Supersonic Lamda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                                lamda3, lamda4, lamda5);
                        info("Nodes: %d %d", pnode, ghostnode);
                        info("Normal: %lf %lf %lf ", nx, ny, nz);
                        info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                        error("Apply_Boundary_Condition: Unable apply boundary condition edge[%d]", i);
                    }
                } else {
                    // Subsonic outflow
                    if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                        instance = BC_SUBSONIC_OUTFLOW;
                    // Subsonic inflow
                    else if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                        instance = BC_SUBSONIC_INFLOW;
                    else {
                        info(" Subsonic Lamda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                                lamda3, lamda4, lamda5);
                        info("Nodes: %d %d", pnode, ghostnode);
                        info("Normal: %lf %lf %lf ", nx, ny, nz);
                        info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                        error("Apply_Boundary_Condition: Unable apply boundary condition edge[%d]", i);
                    }
                }
                break;
            default:
                info("Lamda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Nodes: %d %d", pnode, ghostnode);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Apply_Boundary_Condition: Unable apply boundary condition edge[%d]", i);
                break;
        }

        // Apply Boundary Conditions
        switch (instance) {
            case BC_SOLID_WALL:
                temp    = nx * u_i + ny * v_i + nz * w_i;
                p_b_n   = p_i + rho_0 * c_0 * temp;
                rho_b_n = rho_i + (p_b - p_i) / (c_0 * c_0);
                u_b_n   = u_i - nx * (p_b - p_i) / (rho_0 * c_0);
                v_b_n   = v_i - ny * (p_b - p_i) / (rho_0 * c_0);
                w_b_n   = w_i - nz * (p_b - p_i) / (rho_0 * c_0);
                break;
            case BC_SUBSONIC_OUTFLOW:
                p_b_n   = Inf_Pressure;
                rho_b_n = rho_i + (p_b - p_i) / (c_0 * c_0);
                u_b_n   = u_i - nx * (p_b - p_i) / (rho_0 * c_0);
                v_b_n   = v_i - ny * (p_b - p_i) / (rho_0 * c_0);
                w_b_n   = w_i - nz * (p_b - p_i) / (rho_0 * c_0);
                break;
            case BC_SUBSONIC_INFLOW:
                temp    =   nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
                p_b_n   = 0.5 * (Inf_Pressure + p_i - rho_0 * c_0 * temp);
                rho_b_n = Inf_Rho + (p_b - Inf_Pressure) / (c_0 * c_0);
                u_b_n   = Inf_U + nx * (p_b - Inf_Pressure) / (rho_0 * c_0);
                v_b_n   = Inf_V + ny * (p_b - Inf_Pressure) / (rho_0 * c_0);
                w_b_n   = Inf_W + nz * (p_b - Inf_Pressure) / (rho_0 * c_0);
                break;
            case BC_SUPERSONIC_INFLOW:
                p_b_n   = Inf_Pressure;
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
                p_b_n   = Inf_Pressure;
                rho_b_n = Inf_Rho;
                u_b_n   = Inf_U;
                v_b_n   = Inf_V;
                w_b_n   = Inf_W;
                break;
            default:
                error("Apply_Boundary_Condition: Unable to set BC Type edge[%d]", i);
                break;
        }
        
        et_b_n = p_b_n / ((Gamma - 1.0) * rho_b_n) + 0.5 * (u_b_n * u_b_n + v_b_n * v_b_n + w_b_n * w_b_n);

        // Update the Q's for Ghostnode
        Q1[ghostnode] = rho_b_n;
        Q2[ghostnode] = rho_b_n * u_b_n;
        Q3[ghostnode] = rho_b_n * v_b_n;
        Q4[ghostnode] = rho_b_n * w_b_n;
        Q5[ghostnode] = rho_b_n * et_b_n;
    }
}

