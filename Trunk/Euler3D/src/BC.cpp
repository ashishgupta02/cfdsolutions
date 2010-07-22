/*******************************************************************************
 * File:        BC.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

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
        bndry_edge_info[i].type = bndType[bndry_edge_info[i].tag - 1];
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Apply_Boundary_Condition() {
    int pnode;
    int ghostnode;
    Vector3D normal;
    double p_b, p_i;
    double rho_b, u_b, v_b, w_b, et_b, e_b, c_b;
    double rho_i, u_i, v_i, w_i, et_i, e_i, c_i;
    double c_0, rho_0;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    double mach;
    int instance = -1;
    
    for (int i = 0; i < nBEdge; i++) {
        pnode     = bndry_edge_info[i].node[0];
        ghostnode = bndry_edge_info[i].node[1];

        // Get area vector
        normal = bndry_edge_info[i].areav;
        normal.normalize();
        
        // Physical node
        rho_i = Q1[pnode];
        u_i   = Q2[pnode] / rho_i;
        v_i   = Q3[pnode] / rho_i;
        w_i   = Q4[pnode] / rho_i;
        et_i  = Q5[pnode] / rho_i;
        e_i   = et_i - 0.5 * (u_i * u_i + v_i * v_i + w_i * w_i);
        p_i   = rho_i * e_i * (Gamma - 1.0);
        c_i   = sqrt(Gamma * p_i / rho_i);
        
        // Ghost node
        rho_b = Q1[ghostnode];
        u_b   = Q2[ghostnode] / rho_b;
        v_b   = Q3[ghostnode] / rho_b;
        w_b   = Q4[ghostnode] / rho_b;
        et_b  = Q5[ghostnode] / rho_b;
        e_b   = et_b - 0.5 * (u_b * u_b + v_b * v_b + w_b * w_b);
        p_b   = rho_b * e_b * (Gamma - 1.0);
        c_b   = sqrt(Gamma * p_b / rho_b);
        
        c_0   = 0.5 * (c_i   + c_b);
        rho_0 = 0.5 * (rho_i + rho_b);
        
        // Get eigenvalues
        lamda1 = (u_i + u_b) * normal.vec[0] + (v_i + v_b) * normal.vec[1] + (w_i + w_b) * normal.vec[2];
        lamda1 = 0.5*lamda1;
        lamda2 = lamda1;
        lamda3 = lamda1;
        lamda4 = lamda1 + c_0;
        lamda5 = lamda1 - c_0;
        
        // Determine how to set boundary conditions for edge
        switch (bndry_edge_info[i].type) {
            // Solid Wall
            case BC_SOLID_WALL:
                instance = BC_SOLID_WALL;
                if ((lamda5 >= 0.0) || (lamda4 <= 0.0)) {
                    info("Lamda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Normal: %lf %lf %lf ", normal.vec[0], normal.vec[1], normal.vec[2]);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    error("BC: Unable apply Solid Wall boundary condition edge[%d]", i);
                }
                break;
            // Free Stream
            case BC_FREE_STREAM:
                mach = fabs(lamda1/c_0);
                if (mach < 1.0) {
                    // Subsonic outflow
                    if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                        instance = BC_SUBSONIC_OUTFLOW;
                    // Subsonic inflow
                    else if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                        instance = BC_SUBSONIC_INFLOW;
                    else {
                        info("Subsonic Lamda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                                lamda3, lamda4, lamda5);
                        info("Mach %lg Nodes: %d %d", mach, pnode, ghostnode);
                        info("Normal: %lf %lf %lf ", normal.vec[0], normal.vec[1], normal.vec[2]);
                        info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                        error("BC: Unable apply boundary condition edge[%d]", i);
                    }
                } else {
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
                        info("Normal: %lf %lf %lf ", normal.vec[0], normal.vec[1], normal.vec[2]);
                        info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                        error("BC: Unable apply boundary condition edge[%d]", i);
                    }
                }
                break;
        }

        // Apply Boundary Conditions
        switch (instance) {
            case BC_SOLID_WALL:
                temp  = normal.vec[0] * u_i + normal.vec[1] * v_i + normal.vec[2] * w_i;
                p_b   = p_i + rho_0 * c_0 * temp;
                rho_b = rho_i + (p_b - p_i) / (c_0 * c_0);
                u_b   = u_i - normal.vec[0] * (p_b - p_i) / (rho_0 * c_0);
                v_b   = v_i - normal.vec[1] * (p_b - p_i) / (rho_0 * c_0);
                w_b   = w_i - normal.vec[2] * (p_b - p_i) / (rho_0 * c_0);
                break;
            case BC_SUBSONIC_OUTFLOW:
                p_b   = Inf_Pressure;
                rho_b = rho_i + (p_b - p_i) / (c_0 * c_0);
                u_b   = u_i - normal.vec[0] * (p_b - p_i) / (rho_0 * c_0);
                v_b   = v_i - normal.vec[1] * (p_b - p_i) / (rho_0 * c_0);
                w_b   = w_i - normal.vec[2] * (p_b - p_i) / (rho_0 * c_0);
                break;
            case BC_SUBSONIC_INFLOW:
                temp  =   normal.vec[0] * (Inf_U - u_i)
                        + normal.vec[1] * (Inf_V - v_i)
                        + normal.vec[2] * (Inf_W - w_i);
                p_b   = 0.5 * (Inf_Pressure + p_i - rho_0 * c_0 * temp);
                rho_b = Inf_Rho + (p_b - Inf_Pressure) / (c_0 * c_0);
                u_b   = Inf_U + normal.vec[0] * (p_b - Inf_Pressure) / (rho_0 * c_0);
                v_b   = Inf_V + normal.vec[1] * (p_b - Inf_Pressure) / (rho_0 * c_0);
                w_b   = Inf_W + normal.vec[2] * (p_b - Inf_Pressure) / (rho_0 * c_0);
                break;
            case BC_SUPERSONIC_INFLOW:
                p_b   = Inf_Pressure;
                rho_b = Inf_Rho;
                u_b   = Inf_U;
                v_b   = Inf_V;
                w_b   = Inf_W;
                break;
            case BC_SUPERSONIC_OUTFLOW:
                p_b   = p_i;
                rho_b = rho_i;
                u_b   = u_i;
                v_b   = v_i;
                w_b   = w_i;
                break;
        }
        
        et_b = p_b / ((Gamma - 1.0) * rho_b) + 0.5 * (u_b * u_b + v_b * v_b + w_b * w_b);

        // Read in these values to Q for the Ghostnode under consideration
        Q1[ghostnode] = rho_b;
        Q2[ghostnode] = rho_b * u_b;
        Q3[ghostnode] = rho_b * v_b;
        Q4[ghostnode] = rho_b * w_b;
        Q5[ghostnode] = rho_b * et_b;
    }
}

