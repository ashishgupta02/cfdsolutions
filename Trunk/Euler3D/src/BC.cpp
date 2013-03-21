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
#include "Material.h"
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
    double Q[NEQUATIONS];
    double rho_i, u_i, v_i, w_i, et_i, p_i, c_i, T_i, ubar_i, ht_i, q2_i, mach_i;
    double rho_b, u_b, v_b, w_b, et_b, p_b, c_b, T_b, ubar_b, ht_b, q2_b, mach_b;
    double rho_0, u_0, v_0, w_0, et_0, p_0, c_0, T_0, ubar_0, ht_0, q2_0, mach_0;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    int instance = BC_TYPE_NONE;
    
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

    // Compute Equation of State
    // Physical node
    Q[0] = Q1[node_P];
    Q[1] = Q2[node_P];
    Q[2] = Q3[node_P];
    Q[3] = Q4[node_P];
    Q[4] = Q5[node_P];
    Compute_EOS_Variables_Face(Q, nx, ny, nz, rho_i, p_i, T_i, u_i, v_i, w_i, q2_i, c_i, mach_i, ubar_i, et_i, ht_i);
    p_i += Gauge_Pressure;
    // Ghost node
    Q[0] = Q1[node_G];
    Q[1] = Q2[node_G];
    Q[2] = Q3[node_G];
    Q[3] = Q4[node_G];
    Q[4] = Q5[node_G];
    Compute_EOS_Variables_Face(Q, nx, ny, nz, rho_b, p_b, T_b, u_b, v_b, w_b, q2_b, c_b, mach_b, ubar_b, et_b, ht_b);
    p_b += Gauge_Pressure;
    // Average State
    Q[0] = 0.5*(Q1[node_P] + Q1[node_G]);
    Q[1] = 0.5*(Q2[node_P] + Q2[node_G]);
    Q[2] = 0.5*(Q3[node_P] + Q3[node_G]);
    Q[3] = 0.5*(Q4[node_P] + Q4[node_G]);
    Q[4] = 0.5*(Q5[node_P] + Q5[node_G]);
    Compute_EOS_Variables_Face(Q, nx, ny, nz, rho_0, p_0, T_0, u_0, v_0, w_0, q2_0, c_0, mach_0, ubar_0, et_0, ht_0);
    p_0 += Gauge_Pressure;
    
    // Finally Compute Eigenvalues
    lamda1 = ubar_0;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = lamda1 + c_0;
    lamda5 = lamda1 - c_0;

    // Determine how to set boundary conditions for edge
    switch (bndEdge[BEdgeID].type) {
        // Euler Solid Wall
        case BC_TYPE_EULER_WALL:
            instance = BC_TYPE_EULER_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c_0);
                    info("Nodes: %d %d", node_P, node_G);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Compute_Characteristic_BoundaryCondition: Unable apply Euler Solid Wall boundary condition edge[%d] - 1", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Compute_Characteristic_BoundaryCondition: Unable apply boundary condition edge[%d] - 2", BEdgeID);
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_INFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_INFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_OUTFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_OUTFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Rescue Outflow boundary condition edge[%d] - 4", BEdgeID);
            }
            break;
        // Subsonic Inflow
        case BC_TYPE_SUBSONIC_INFLOW:
            instance = BC_TYPE_SUBSONIC_INFLOW;
            // Subsonic inflow
            if ((lamda1 <= 0.0) && (lamda4 >= 0.0) && (lamda5 <= 0.0))
                instance = BC_TYPE_SUBSONIC_INFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Unable apply Subsonic Inflow boundary condition edge[%d] - 5", BEdgeID);
            }
            break;
        // Subsonic Outflow
        case BC_TYPE_SUBSONIC_OUTFLOW:
            instance = BC_TYPE_SUBSONIC_OUTFLOW;
            // Subsonic outflow
            if ((lamda1 >= 0.0) && (lamda5 <= 0.0) && (lamda4 >= 0.0))
                instance = BC_TYPE_SUBSONIC_OUTFLOW;
            else {
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Nodes: %d %d", node_P, node_G);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition: Unable apply Supersonic Outflow boundary condition edge[%d] - 8", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c_0);
            info("Nodes: %d %d", node_P, node_G);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Characteristic_BoundaryCondition: Unable apply boundary condition edge[%d] - 9", BEdgeID);
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_TYPE_EULER_WALL:
            p_b_n = p_i + rho_0 * c_0 * ubar_i;
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
        case BC_TYPE_SUBSONIC_INFLOW:
            break;
        case BC_TYPE_SUBSONIC_OUTFLOW:
            p_b_n   = Outflow_Pressure + Gauge_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            break;
        case BC_TYPE_SUPERSONIC_INFLOW:
            break;
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW:
            temp    = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
            p_b_n   = 0.5 * (Inf_Pressure + Gauge_Pressure + p_i - rho_0 * c_0 * temp);
            rho_b_n = Inf_Rho + (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (c_0 * c_0);
            u_b_n   = Inf_U + nx * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            v_b_n   = Inf_V + ny * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            w_b_n   = Inf_W + nz * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
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
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Characteristic_BoundaryCondition: Unable to set BC Type edge[%d]", BEdgeID);
            break;
    }

    if (p_b_n < 0.0)
        error("Compute_Characteristic_BoundaryCondition: Negative Pressure: Node: %d Ghost: %d, P_i: %e, P_b_n: %e", node_P, node_G, p_i, p_b_n);
    else
        p_b_n -= Gauge_Pressure; // Only Perturbation

    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CONSERVATIVE) {
        et_b_n = Get_TotalEnergy(rho_b_n, p_b_n, u_b_n, v_b_n, w_b_n);
        // Update the Q's for Ghost Node
        Q1[node_G] = rho_b_n;
        Q2[node_G] = rho_b_n * u_b_n;
        Q3[node_G] = rho_b_n * v_b_n;
        Q4[node_G] = rho_b_n * w_b_n;
        Q5[node_G] = rho_b_n * et_b_n;
    }
    
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PRIMITIVE_PUT) {
        // Update the Q's for Ghost Node
        Q1[node_G] = p_b_n;
        Q2[node_G] = u_b_n;
        Q3[node_G] = v_b_n;
        Q4[node_G] = w_b_n;
        Q5[node_G] = Get_Temperature(rho_b_n, p_b_n);
    }
    
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_PRIMITIVE_RUP) {
        // Update the Q's for Ghost Node
        Q1[node_G] = rho_b_n;
        Q2[node_G] = u_b_n;
        Q3[node_G] = v_b_n;
        Q4[node_G] = w_b_n;
        Q5[node_G] = p_b_n;
    }
}

//------------------------------------------------------------------------------
//! This Routine is created to assist Finite Difference Jacobian Computations
//! Note: Outward Pointing Boundary Normals
//! rvalue: 0 => Success, 1: Negative Pressure
//------------------------------------------------------------------------------
int Compute_Characteristic_BoundaryCondition_New(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    int i;
    double nx, ny, nz;
    double rho_i, u_i, v_i, w_i, et_i, p_i, c_i, T_i, ubar_i, ht_i, q2_i, mach_i;
    double rho_b, u_b, v_b, w_b, et_b, p_b, c_b, T_b, ubar_b, ht_b, q2_b, mach_b;
    double rho_0, u_0, v_0, w_0, et_0, p_0, c_0, T_0, ubar_0, ht_0, q2_0, mach_0;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    int instance = BC_TYPE_NONE;
    int rvalue = 0;
    
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

    // Compute Equation of State
    // Physical node
    Compute_EOS_Variables_Face(Q_L, nx, ny, nz, rho_i, p_i, T_i, u_i, v_i, w_i, q2_i, c_i, mach_i, ubar_i, et_i, ht_i);
    p_i += Gauge_Pressure;
    // Ghost node
    Compute_EOS_Variables_Face(Q_R, nx, ny, nz, rho_b, p_b, T_b, u_b, v_b, w_b, q2_b, c_b, mach_b, ubar_b, et_b, ht_b);
    p_b += Gauge_Pressure;
    // Average State
    for (i = 0; i < NEQUATIONS; i++)
        Q_B[i] = 0.5*(Q_L[i] + Q_R[i]);
    Compute_EOS_Variables_Face(Q_B, nx, ny, nz, rho_0, p_0, T_0, u_0, v_0, w_0, q2_0, c_0, mach_0, ubar_0, et_0, ht_0);
    p_0 += Gauge_Pressure;
    
    // Compute Eigenvalues
    lamda1 = ubar_0;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = lamda1 + c_0;
    lamda5 = lamda1 - c_0;

    // Determine how to set boundary conditions for edge
    switch (BCType) {
        // Euler Solid Wall
        case BC_TYPE_EULER_WALL:
            instance = BC_TYPE_EULER_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c_0);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Compute_Characteristic_BoundaryCondition_New: Unable apply Euler Solid Wall boundary condition - 1");
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Compute_Characteristic_BoundaryCondition_New: Unable apply boundary condition - 2");
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_INFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_INFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition_New: Rescue Inflow boundary condition edge[%d] - 3", BEdgeID);
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_OUTFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_OUTFLOW;
                
//                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
//                        lamda3, lamda4, lamda5);
//                info("Mach: %lf", fabs(lamda1) / c_0);
//                info("Normal: %lf %lf %lf ", nx, ny, nz);
//                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
//                warn("Compute_Characteristic_BoundaryCondition_New: Rescue Outflow boundary condition edge[%d] - 4", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition_New: Unable apply Subsonic Inflow boundary condition edge[%d] - 5", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition_New: Unable apply Subsonic Outflow boundary condition edge[%d] - 6", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition_New: Unable apply Supersonic Inflow boundary condition edge[%d] - 7", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_BoundaryCondition_New: Unable apply Supersonic Outflow boundary condition edge[%d] - 8", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c_0);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Characteristic_BoundaryCondition_New: Unable apply boundary condition - 9");
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_TYPE_EULER_WALL:
            p_b_n = p_i + rho_0 * c_0 * ubar_i;
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
        case BC_TYPE_SUBSONIC_INFLOW:
            break;
        case BC_TYPE_SUBSONIC_OUTFLOW:
            p_b_n   = Outflow_Pressure + Gauge_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            break;
        case BC_TYPE_SUPERSONIC_INFLOW:
            break;
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            p_b_n   = p_i;
            rho_b_n = rho_i;
            u_b_n   = u_i;
            v_b_n   = v_i;
            w_b_n   = w_i;
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW:
            temp    = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
            p_b_n   = 0.5 * (Inf_Pressure + Gauge_Pressure + p_i - rho_0 * c_0 * temp);
            rho_b_n = Inf_Rho + (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (c_0 * c_0);
            u_b_n   = Inf_U + nx * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            v_b_n   = Inf_V + ny * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            w_b_n   = Inf_W + nz * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
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
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Characteristic_BoundaryCondition_New: Unable to set BC Type");
            break;
    }

    if (p_b_n < 0.0) {
        warn("Compute_Characteristic_BoundaryCondition_New: Negative Pressure Detected P_i: %e, P_b_n: %e", p_i, p_b_n);
        rvalue = 1;
    } else
        p_b_n -= Gauge_Pressure; // Only Perturbation

    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CONSERVATIVE) {
        et_b_n = Get_TotalEnergy(rho_b_n, p_b_n, u_b_n, v_b_n, w_b_n);
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = rho_b_n * u_b_n;
        Q_B[2] = rho_b_n * v_b_n;
        Q_B[3] = rho_b_n * w_b_n;
        Q_B[4] = rho_b_n * et_b_n;
    }
    
     // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PRIMITIVE_PUT) {
        // Update the Q's for Ghost Node
        Q_B[0] = p_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = Get_Temperature(rho_b_n, p_b_n);
    }
    
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_PRIMITIVE_RUP) {
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;
    }
    
    return rvalue;
}

//------------------------------------------------------------------------------
//! This Routine is created to assist Finite Difference Jacobian Computations
//! Note: Outward Pointing Boundary Normals
//! rvalue: 0 => Success, 1: Negative Pressure
//------------------------------------------------------------------------------
int Compute_Precondition_Characteristic_BoundaryCondition(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    int i;
    double nx, ny, nz;
    double rho_i, u_i, v_i, w_i, et_i, p_i, c_i, T_i, ubar_i, ht_i, q2_i, mach_i;
    double rho_b, u_b, v_b, w_b, et_b, p_b, c_b, T_b, ubar_b, ht_b, q2_b, mach_b;
    double rho_0, u_0, v_0, w_0, et_0, p_0, c_0, T_0, ubar_0, ht_0, q2_0, mach_0;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    double eps, alpha, beta, omega, Ur, mach, sigma;
    double ubar_pre, c_pre;
    double c1, c2, c3;
    double dtmp;
    int instance = BC_TYPE_NONE;
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

    // Compute Equation of State
    // Physical node
    Compute_EOS_Variables_Face(Q_L, nx, ny, nz, rho_i, p_i, T_i, u_i, v_i, w_i, q2_i, c_i, mach_i, ubar_i, et_i, ht_i);
    p_i += Gauge_Pressure;
    // Ghost node
    Compute_EOS_Variables_Face(Q_R, nx, ny, nz, rho_b, p_b, T_b, u_b, v_b, w_b, q2_b, c_b, mach_b, ubar_b, et_b, ht_b);
    p_b += Gauge_Pressure;
    // Average State
    for (i = 0; i < NEQUATIONS; i++)
        Q_B[i] = 0.5*(Q_L[i] + Q_R[i]);
    Compute_EOS_Variables_Face(Q_B, nx, ny, nz, rho_0, p_0, T_0, u_0, v_0, w_0, q2_0, c_0, mach_0, ubar_0, et_0, ht_0);
    p_0 += Gauge_Pressure;
    
    // Get eigenvalues
    // Compute Weiss Smith Precondition Variables
    if (PrecondMethod == PRECOND_METHOD_ROE_WS) {
        // Compute Precondition Variables
        mach = mach_0;
        eps  = MIN(1.0, MAX(1e-10, mach*mach));

        // Compute beta for Ideal Gas
        beta   = 1.0/(c_0*c_0);
        
        dtmp = sqrt(q2_0);
        
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
    if (PrecondMethod == PRECOND_METHOD_ROE_CV || PrecondMethod == PRECOND_METHOD_ROE_CV_ORIGINAL) {
        // Get two nodes of edge
        node_L = bndEdge[BEdgeID].node[0];
        node_R = bndEdge[BEdgeID].node[1];
        
        int i, nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double dp_max, mach_max;
        
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
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    mach_max = MAX(mach_max, lmach);
                    dp_max   = MAX(dp_max, fabs(p_i - lp));
                }
            } else {
                mach_max = MAX(mach_i, mach_b);
                dp_max   = fabs(p_i - p_b);
            }
        }
        
        // Compute Precondition Variables
        mach  = mach_0;
        sigma = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            mach = MAX(mach, mach_max);
            if (mach < Ref_Mach)
                sigma = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
            else
                sigma = 0.5;
            sigma = sigma*(sqrt(mach*Ref_Mach)+ mach);
//            sigma = sqrt(mach*mach + MAX(0.005*Ref_Mach*Ref_Mach, 0.05*Ref_Mach*mach));
//            if (sigma < Ref_Mach)
//                sigma = sqrt(Ref_Mach*sigma);
            sigma = MIN(1.0, sigma);
//            //if (Inf_Mach < 0.5) {
//                sigma = MIN(sqrt(1.0e-8/fabs(ubar_0)), Inf_Mach);
//                // Check if Precondition Smoother is Requested
//                //if (PrecondSmooth != 0)
//                    sigma = MIN(1.0, MAX(MAX(MAX(sigma, mach), mach_max), 2.0*dp_max/(rho_0*c_0*c_0)));
//                //else
//                //    sigma = MIN(1.0, MAX(sigma, mach));
//            //sigma = MAX(mach, dp_max/(rho_0*c_0*c_0));
//            //sigma = MAX(1.0e-6, MIN(mach_max, dp_max/(Inf_Mach*Inf_Mach*rho_0*c_0*c_0)));
//                if (sigma > Inf_Mach)
//                    sigma = Inf_Mach;
//            //}
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            sigma = MIN(1.0, Ref_Mach);
        
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, sigma);
        MaxPrecondSigma = MAX(MaxPrecondSigma, sigma);
        
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
    
    // Compute Briley Taylor Whitfield Precondition Variables
    if (PrecondMethod == PRECOND_METHOD_ROE_BTW || PrecondMethod == PRECOND_METHOD_ROE_BTW_ORIGINAL) {
        // Get two nodes of edge
        node_L = bndEdge[BEdgeID].node[0];
        node_R = bndEdge[BEdgeID].node[1];
        
        int i, nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double beta_m, beta_p;
        double dp_max, mach_max;
        
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
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    mach_max = MAX(mach_max, lmach);
                    dp_max   = MAX(dp_max, fabs(p_i - lp));
                }
            } else {
                mach_max = MAX(mach_i, mach_b);
                dp_max   = fabs(p_i - p_b);
            }
        }
        
        // Compute Precondition Variables
        mach = mach_0;
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            beta = MAX(mach, mach_max);
            if (beta < Ref_Mach) {
                beta = MAX(5.0e-5, MAX(sqrt(1.0e-11/sqrt(q2_0)), sqrt(Ref_Mach*beta)));
                beta = MIN(beta, Ref_Mach);
            }
            beta = MIN(1.0, beta);
//            beta = MAX(beta, sqrt(1.0e-12/sqrt(q2_0)));
//            beta = MAX(beta, 2.0*dp_max/(rho_0*c_0*c_0));
//            beta = MAX(MAX(beta, 1.0e-5), 0.01*Ref_Mach);
//            beta = MIN(1.0, MIN(beta, Ref_Mach));
//            if (beta > 0.7) {
//                beta = 1.0;
//            } else {
//                beta = 2.0*beta*beta;
//                beta = beta/(1.0 - beta);
//                beta = sqrt(beta/(1 + (Gamma - 1.0)*beta));
//            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, beta);
        MaxPrecondSigma = MAX(MaxPrecondSigma, beta);
        beta   = beta*beta;
        beta_m = 0.5*(1.0 - beta);
        beta_p = 0.5*(1.0 + beta);
        sigma  = sqrt(ubar_0*ubar_0*beta_m*beta_m + beta*c_0*c_0);
        
        // alpha, c1, c2, c3
        alpha = 0.5*ubar_0*beta_m/sigma;
        omega = 2.0*beta_m*rho_0*ubar_0/(c_0*c_0*beta);
        c1    = 0.5/(rho_0*sigma);
        c2    = (0.5 + alpha);
        c3    = (0.5 - alpha);
        
        // Compute c_pre
        c_pre = sigma;

        // Compute ubar_pre
        ubar_pre = ubar_0*beta_p;
    }
    
    // Compute Eriksson Precondition Variables
    if (PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON || PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON_ORIGINAL) {
        // Get two nodes of edge
        node_L = bndEdge[BEdgeID].node[0];
        node_R = bndEdge[BEdgeID].node[1];
        
        int i, nid;
        double lQ[NEQUATIONS];
        double lrho, lu, lv, lw, lp, lT, lc, lht, let, lq2, lmach;
        double beta_m, beta_p;
        double dp_max, mach_max;
        
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
                    Compute_EOS_Variables_ControlVolume(lQ, lrho, lp, lT, lu, lv, lw, lq2, lc, lmach, let, lht);

                    // Compute Smoothing Parameters
                    mach_max = MAX(mach_max, lmach);
                    dp_max   = MAX(dp_max, fabs(p_i - lp));
                }
            } else {
                mach_max = MAX(mach_i, mach_b);
                dp_max   = fabs(p_i - p_b);
            }
        }
        
        // Compute Precondition Variables
        mach = mach_0;
        beta = 1.0; // Corresponds to no pre-conditioning
        if (PrecondType == PRECOND_TYPE_LOCAL) {
            mach = MAX(mach, mach_max);
//            if (beta < Ref_Mach){
//                beta = sqrt(Ref_Mach*mach);
//                beta = sqrt(0.5*(MAX(Ref_Mach, mach)*mach + mach*mach));
            //beta = sqrt(mach*mach + MAX(0.005*Ref_Mach*Ref_Mach, 0.05*Ref_Mach*mach));
            if (mach < Ref_Mach)
                beta = 0.5*(mach*mach)/(Ref_Mach*Ref_Mach) - mach/(Ref_Mach) + 1.0;
            else
                beta = 0.5;
            beta = beta*(sqrt(mach*Ref_Mach)+ mach);
//            }
            beta = MIN(1.0, beta);
//            beta = MAX(beta, sqrt(1.0e-12/sqrt(q2_0)));
//            beta = MAX(beta, 2.0*dp_max/(rho_0*c_0*c_0));
//            beta = MAX(MAX(beta, 1.0e-5), 0.01*Ref_Mach);
//            beta = MIN(1.0, MIN(beta, Ref_Mach));
//            if (beta > 0.7) {
//                beta = 1.0;
//            } else {
//                beta = 2.0*beta*beta;
//                beta = beta/(1.0 - beta);
//                beta = sqrt(beta/(1 + (Gamma - 1.0)*beta));
//            }
        }
        if (PrecondType == PRECOND_TYPE_GLOBAL)
            beta = MIN(1.0, Ref_Mach);
        
        // Min and Max Precondition Variable
        MinPrecondSigma = MIN(MinPrecondSigma, beta);
        MaxPrecondSigma = MAX(MaxPrecondSigma, beta);
        beta   = beta*beta;
        beta_m = 0.5*(1.0 - beta);
        beta_p = 0.5*(1.0 + beta);
        alpha  = 2.0*sqrt(ubar_0*ubar_0*beta_m*beta_m + beta*c_0*c_0);
        
        // c1, c2, c3
        c1 = 1.0/(rho_0*alpha);
        c2 = 0.5 + (beta_m*ubar_0)/alpha;
        c3 = 0.5 + (beta_m*ubar_0)/alpha;
        
        // Compute c_pre
        c_pre = 0.5*alpha;

        // Compute ubar_pre
        ubar_pre = ubar_0*beta_p;
    }
    
    // No Precondition
    if (PrecondMethod == PRECOND_METHOD_NONE || PrecondMethod == PRECOND_METHOD_ROE_LMFIX) {
        ubar_pre = ubar_0;
        c_pre    = c_0;
    }
    
    // Finally Compute Pre/Non-Condition Eigenvalues
    lamda1 = ubar_0;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = ubar_pre + c_pre;
    lamda5 = ubar_pre - c_pre;

    // Determine how to set boundary conditions for edge
    switch (BCType) {
        // Euler Solid Wall
        case BC_TYPE_EULER_WALL:
            instance = BC_TYPE_EULER_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c_0);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_INFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_INFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_OUTFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_OUTFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply Supersonic Outflow boundary condition edge[%d] - 8", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c_0);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Precondition_Characteristic_BoundaryCondition: Unable apply boundary condition - 9");
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_TYPE_EULER_WALL:
            // No Preconditioning 
            if (PrecondMethod == PRECOND_METHOD_NONE || PrecondMethod == PRECOND_METHOD_ROE_LMFIX) {
                p_b_n = p_i + rho_0 * c_0 * ubar_i;
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
            if (PrecondMethod == PRECOND_METHOD_ROE_WS) {
                // TODO
            }
            
            // Cecile Voizat Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_CV || PrecondMethod == PRECOND_METHOD_ROE_CV_ORIGINAL) {
                p_b_n   = p_i + (c2/c1)*ubar_i;
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
            
            // Briley Taylor Whitfield Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_BTW || PrecondMethod == PRECOND_METHOD_ROE_BTW_ORIGINAL) {
                p_b_n   = p_i + (c2/c1)*ubar_i;
                // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
                // This causes the flow to be allowed through the surface initially
                if (Iteration < ZPGIteration) {
                    if (p_b_n < 0.0)
                        p_b_n = p_i;
                }
                rho_b_n = rho_i  + (1.0/(c_0*c_0*beta) - omega*c1/c2)*(p_b_n - p_i);
                u_b_n   = u_i - nx * (c1/c2)*(p_b_n - p_i);
                v_b_n   = v_i - ny * (c1/c2)*(p_b_n - p_i);
                w_b_n   = w_i - nz * (c1/c2)*(p_b_n - p_i);
            }
            
            // Eriksson Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON || PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON_ORIGINAL) {
                p_b_n   = p_i + (c2/c1)*ubar_i;
                // Apply Zero Pressure Gradient (ZPG) BC for Harsh initialization
                // This causes the flow to be allowed through the surface initially
                if (Iteration < ZPGIteration) {
                    if (p_b_n < 0.0)
                        p_b_n = p_i;
                }
                rho_b_n = rho_i  + (p_b_n - p_i)/(c_0 * c_0);
                u_b_n   = u_i - nx * (c1/c2)*(p_b_n - p_i);
                v_b_n   = v_i - ny * (c1/c2)*(p_b_n - p_i);
                w_b_n   = w_i - nz * (c1/c2)*(p_b_n - p_i);
            }
            break;
        case BC_TYPE_SUBSONIC_INFLOW:
            break;
        case BC_TYPE_SUBSONIC_OUTFLOW:
            break;
        case BC_TYPE_SUPERSONIC_INFLOW:
            break;
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW:
            // No Preconditioning 
            if (PrecondMethod == PRECOND_METHOD_NONE || PrecondMethod == PRECOND_METHOD_ROE_LMFIX) {
                p_b_n   = Inf_Pressure + Gauge_Pressure;
                rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
                u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
                v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
                w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            }
            
            // Weiss Smith Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_WS) {
                // TODO
            }
            
            // Cecile Voizat Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_CV || PrecondMethod == PRECOND_METHOD_ROE_CV_ORIGINAL) {
                temp  = nx * (u_i - Inf_U) + ny * (v_i - Inf_V) + nz * (w_i - Inf_W);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c2*(Inf_Pressure + Gauge_Pressure) + c3*p_i) + c2*c3*temp);
                rho_b_n = rho_i * pow((p_b_n/p_i), (1.0/Gamma));
                u_b_n   = u_i - nx * (c1/c2)*(p_b_n - p_i);
                v_b_n   = v_i - ny * (c1/c2)*(p_b_n - p_i);
                w_b_n   = w_i - nz * (c1/c2)*(p_b_n - p_i);
            }
            
            // Briley Taylor Whitfield Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_BTW || PrecondMethod == PRECOND_METHOD_ROE_BTW_ORIGINAL) {
                temp  = nx * (u_i - Inf_U) + ny * (v_i - Inf_V) + nz * (w_i - Inf_W);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c2*(Inf_Pressure + Gauge_Pressure) + c3*p_i) + c2*c3*temp);
                rho_b_n = rho_i + (1.0/(c_0*c_0*beta) - omega*c1/c2)*(p_b_n - p_i);
                u_b_n   = u_i - nx * (c1/c2)*(p_b_n - p_i);
                v_b_n   = v_i - ny * (c1/c2)*(p_b_n - p_i);
                w_b_n   = w_i - nz * (c1/c2)*(p_b_n - p_i);
            }
            
            // Eriksson Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON || PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON_ORIGINAL) {
                temp  = nx * (u_i - Inf_U) + ny * (v_i - Inf_V) + nz * (w_i - Inf_W);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c2*(Inf_Pressure + Gauge_Pressure) + c3*p_i) + c2*c3*temp);
                rho_b_n = rho_i + (p_b_n - p_i)/(c_0 * c_0);
                u_b_n   = u_i - nx * (c1/c2)*(p_b_n - p_i);
                v_b_n   = v_i - ny * (c1/c2)*(p_b_n - p_i);
                w_b_n   = w_i - nz * (c1/c2)*(p_b_n - p_i);
            }
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW:
            // No Preconditioning 
            if (PrecondMethod == PRECOND_METHOD_NONE || PrecondMethod == PRECOND_METHOD_ROE_LMFIX) {
                temp    = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
                p_b_n   = 0.5 * (Inf_Pressure + Gauge_Pressure + p_i - rho_0 * c_0 * temp);
                rho_b_n = Inf_Rho + (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (c_0 * c_0);
                u_b_n   = Inf_U + nx * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
                v_b_n   = Inf_V + ny * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
                w_b_n   = Inf_W + nz * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            }
            
            // Weiss Smith Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_WS) {
                // TODO
            }
            
            // Cecile Voizat Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_CV || PrecondMethod == PRECOND_METHOD_ROE_CV_ORIGINAL) {
                temp  = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c2*(Inf_Pressure + Gauge_Pressure) + c3*p_i) - c2*c3*temp);
                rho_b_n = Inf_Rho * pow((p_b_n/(Inf_Pressure + Gauge_Pressure)), (1.0/Gamma));
                u_b_n   = Inf_U + nx*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                v_b_n   = Inf_V + ny*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                w_b_n   = Inf_W + nz*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
            }
            
            // Briley Taylor Whitfield Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_BTW || PrecondMethod == PRECOND_METHOD_ROE_BTW_ORIGINAL) {
                temp  = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c2*(Inf_Pressure + Gauge_Pressure) + c3*p_i) - c2*c3*temp);
                rho_b_n = Inf_Rho + (omega*c1/c3 + 1.0/(c_0*c_0*beta))*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                u_b_n   = Inf_U + nx*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                v_b_n   = Inf_V + ny*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                w_b_n   = Inf_W + nz*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
            }
            
            // Eriksson Precondition
            if (PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON || PrecondMethod == PRECOND_METHOD_ROE_ERIKSSON_ORIGINAL) {
                temp  = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
                p_b_n = (1.0/(c1*(c2+c3)))*(c1*(c2*(Inf_Pressure + Gauge_Pressure) + c3*p_i) - c2*c3*temp);
                rho_b_n = Inf_Rho + (p_b_n - (Inf_Pressure + Gauge_Pressure))/(c_0 * c_0);
                u_b_n   = Inf_U + nx*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                v_b_n   = Inf_V + ny*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
                w_b_n   = Inf_W + nz*(c1/c3)*(p_b_n - (Inf_Pressure + Gauge_Pressure));
            }
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
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
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Precondition_Characteristic_BoundaryCondition: Unable to set BC Type");
            break;
    }

    if (p_b_n < 0.0) {
        warn("Compute_Precondition_Characteristic_BoundaryCondition: Negative Pressure Detected P_i: %e, P_b_n: %e", p_i, p_b_n);
        rvalue = 1;
    } else
        p_b_n -= Gauge_Pressure; // Only Perturbation

    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CONSERVATIVE) {
        et_b_n = Get_TotalEnergy(rho_b_n, p_b_n, u_b_n, v_b_n, w_b_n);
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = rho_b_n * u_b_n;
        Q_B[2] = rho_b_n * v_b_n;
        Q_B[3] = rho_b_n * w_b_n;
        Q_B[4] = rho_b_n * et_b_n;
    }
    
     // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PRIMITIVE_PUT) {
        // Update the Q's for Ghost Node
        Q_B[0] = p_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = Get_Temperature(rho_b_n, p_b_n);
    }
    
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_PRIMITIVE_RUP) {
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;
    }
    
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
        
        if (PrecondMethod == PRECOND_METHOD_NONE)
            rvalue = Compute_Characteristic_BoundaryCondition_New(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
        else
            rvalue = Compute_Precondition_Characteristic_BoundaryCondition(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
        
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
//! Euler Solid Wall Boundary Condition is Pressure Based Hence Weak Formulation
//! Note: Outward Pointing Boundary Normals
//------------------------------------------------------------------------------
int Compute_Characteristic_Weak_BoundaryCondition(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], 
        Vector3D AreaVec, int BEdgeID, int BCType, int Iteration, double Q_B[NEQUATIONS]) {
    int i;
    double nx, ny, nz;
    double rho_i, u_i, v_i, w_i, et_i, p_i, c_i, T_i, ubar_i, ht_i, q2_i, mach_i;
    double rho_b, u_b, v_b, w_b, et_b, p_b, c_b, T_b, ubar_b, ht_b, q2_b, mach_b;
    double rho_0, u_0, v_0, w_0, et_0, p_0, c_0, T_0, ubar_0, ht_0, q2_0, mach_0;
    double rho_b_n, u_b_n, v_b_n, w_b_n, et_b_n, p_b_n;
    double lamda1, lamda2, lamda3, lamda4, lamda5;
    double temp;
    int instance = BC_TYPE_NONE;
    int rvalue = 0;
    
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

    // Compute Equation of State
    // Physical node
    Compute_EOS_Variables_Face(Q_L, nx, ny, nz, rho_i, p_i, T_i, u_i, v_i, w_i, q2_i, c_i, mach_i, ubar_i, et_i, ht_i);
    p_i += Gauge_Pressure;
    // Ghost node
    Compute_EOS_Variables_Face(Q_R, nx, ny, nz, rho_b, p_b, T_b, u_b, v_b, w_b, q2_b, c_b, mach_b, ubar_b, et_b, ht_b);
    p_b += Gauge_Pressure;
    // Average State
    for (i = 0; i < NEQUATIONS; i++)
        Q_B[i] = 0.5*(Q_L[i] + Q_R[i]);
    Compute_EOS_Variables_Face(Q_B, nx, ny, nz, rho_0, p_0, T_0, u_0, v_0, w_0, q2_0, c_0, mach_0, ubar_0, et_0, ht_0);
    p_0 += Gauge_Pressure;
    
    // Finally Compute Eigenvalues
    lamda1 = ubar_0;
    lamda2 = lamda1;
    lamda3 = lamda1;
    lamda4 = lamda1 + c_0;
    lamda5 = lamda1 - c_0;

    // Determine how to set boundary conditions for edge
    switch (BCType) {
        // Euler Solid Wall
        case BC_TYPE_EULER_WALL:
            instance = BC_TYPE_EULER_WALL;
            if (Iteration >= ZPGIteration) {
                if ((lamda5 > 0.0) || (lamda4 < 0.0)) {
                    info("Wall Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                            lamda3, lamda4, lamda5);
                    info("Mach: %lf", fabs(lamda1) / c_0);
                    info("Normal: %lf %lf %lf ", nx, ny, nz);
                    info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                    warn("Compute_Characteristic_Weak_BoundaryCondition: Unable apply Euler Solid Wall boundary condition - 1");
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                error("Compute_Characteristic_Weak_BoundaryCondition: Unable apply boundary condition - 2");
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_INFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_INFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_Weak_BoundaryCondition: Rescue Inflow boundary condition edge[%d] - 3", BEdgeID);
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
                if (fabs(lamda1) < c_0)
                    instance = BC_TYPE_SUBSONIC_OUTFLOW;
                else
                    instance = BC_TYPE_SUPERSONIC_OUTFLOW;
                
                info("Lambda: %lf %lf %lf %10.5e %lf", lamda1, lamda2,
                        lamda3, lamda4, lamda5);
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_Weak_BoundaryCondition: Rescue Outflow boundary condition edge[%d] - 4", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_Weak_BoundaryCondition: Unable apply Subsonic Inflow boundary condition edge[%d] - 5", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_Weak_BoundaryCondition: Unable apply Subsonic Outflow boundary condition edge[%d] - 6", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_Weak_BoundaryCondition: Unable apply Supersonic Inflow boundary condition edge[%d] - 7", BEdgeID);
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
                info("Mach: %lf", fabs(lamda1) / c_0);
                info("Normal: %lf %lf %lf ", nx, ny, nz);
                info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
                warn("Compute_Characteristic_Weak_BoundaryCondition: Unable apply Supersonic Outflow boundary condition edge[%d] - 8", BEdgeID);
            }
            break;
        default:
            info("Lambda: %lf %lf %lf %lf %lf", lamda1, lamda2,
                    lamda3, lamda4, lamda5);
            info("Mach: %lf", fabs(lamda1) / c_0);
            info("Normal: %lf %lf %lf ", nx, ny, nz);
            info("Velocity: %lf %lf %lf ", u_i, v_i, w_i);
            error("Compute_Characteristic_Weak_BoundaryCondition: Unable apply boundary condition - 9");
            break;
    }

    // Apply Boundary Conditions
    switch (instance) {
        case BC_TYPE_EULER_WALL:
            // Compute the Normal Component of Velocity and Subtract From Velocity
            u_b_n   = u_i - nx * ubar_i;
            v_b_n   = v_i - ny * ubar_i;
            w_b_n   = w_i - nz * ubar_i;
            p_b_n   = p_i;
            rho_b_n = rho_i;
            break;
        case BC_TYPE_SUBSONIC_INFLOW:
            break;
        case BC_TYPE_SUBSONIC_OUTFLOW:
            break;
        case BC_TYPE_SUPERSONIC_INFLOW:
            break;
        case BC_TYPE_SUPERSONIC_OUTFLOW:
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_OUTFLOW:
            // According to Characteristics
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = rho_i + (p_b_n - p_i) / (c_0 * c_0);
            u_b_n   = u_i - nx * (p_b_n - p_i) / (rho_0 * c_0);
            v_b_n   = v_i - ny * (p_b_n - p_i) / (rho_0 * c_0);
            w_b_n   = w_i - nz * (p_b_n - p_i) / (rho_0 * c_0);
            break;
        case BC_TYPE_FAR_FIELD_SUBSONIC_INFLOW:
            // According to Characteristics
            temp    = nx * (Inf_U - u_i) + ny * (Inf_V - v_i) + nz * (Inf_W - w_i);
            p_b_n   = 0.5 * (Inf_Pressure + Gauge_Pressure + p_i - rho_0 * c_0 * temp);
            rho_b_n = Inf_Rho + (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (c_0 * c_0);
            u_b_n   = Inf_U + nx * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            v_b_n   = Inf_V + ny * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            w_b_n   = Inf_W + nz * (p_b_n - (Inf_Pressure + Gauge_Pressure)) / (rho_0 * c_0);
            break;
        case BC_TYPE_FAR_FIELD_SUPERSONIC_INFLOW:
            p_b_n   = Inf_Pressure + Gauge_Pressure;
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
            p_b_n   = Inf_Pressure + Gauge_Pressure;
            rho_b_n = Inf_Rho;
            u_b_n   = Inf_U;
            v_b_n   = Inf_V;
            w_b_n   = Inf_W;
            break;
        default:
            error("Compute_Characteristic_Weak_BoundaryCondition: Unable to set BC Type");
            break;
    }

    if (p_b_n < 0.0) {
        warn("Compute_Characteristic_Weak_BoundaryCondition: Negative Pressure Detected P_i: %e, P_b_n: %e", p_i, p_b_n);
        rvalue = 1;
    } else
        p_b_n -= Gauge_Pressure; // Only Perturbation

    // Conservative Variable Formulation
    if (VariableType == VARIABLE_CONSERVATIVE) {
        et_b_n = Get_TotalEnergy(rho_b_n, p_b_n, u_b_n, v_b_n, w_b_n);
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = rho_b_n * u_b_n;
        Q_B[2] = rho_b_n * v_b_n;
        Q_B[3] = rho_b_n * w_b_n;
        Q_B[4] = rho_b_n * et_b_n;
    }
    
    // Primitive Variable Formulation Pressure Velocity Temperature
    if (VariableType == VARIABLE_PRIMITIVE_PUT) {
        // Update the Q's for Ghost Node
        Q_B[0] = p_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = Get_Temperature(rho_b_n, p_b_n);
    }
    
    // Primitive Variable Formulation Density Velocity Pressure
    if (VariableType == VARIABLE_PRIMITIVE_RUP) {
        // Update the Q's for Ghost Node
        Q_B[0] = rho_b_n;
        Q_B[1] = u_b_n;
        Q_B[2] = v_b_n;
        Q_B[3] = w_b_n;
        Q_B[4] = p_b_n;
    }
    
    return rvalue;
}

//------------------------------------------------------------------------------
//! Note: Outward Pointing Boundary Normals
//------------------------------------------------------------------------------
void Apply_Characteristic_Weak_Boundary_Condition(int Iteration) {
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
        
        rvalue = Compute_Characteristic_Weak_BoundaryCondition(Q_L, Q_R, areavec, iBEdge, bndEdge[iBEdge].type, Iteration, Q_B);
        if (rvalue != 0)
            error("Apply_Characteristic_Weak_Boundary_Condition: BndNode: %d,  GhostNode: %d", node_L, node_R);
        
        // Update the Boundary Value of Ghost Node
        Q1[node_R] = Q_B[0];
        Q2[node_R] = Q_B[1];
        Q3[node_R] = Q_B[2];
        Q4[node_R] = Q_B[3];
        Q5[node_R] = Q_B[4];
    }
}


//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
void Apply_Boundary_Condition(int Iteration) {
    switch (BCMethod) {
        case BC_METHOD_CHARCTERISTIC: // Characteristic Boundary Condition
            Apply_Characteristic_Boundary_Condition(Iteration);
            break;
        case BC_METHOD_CHARCTERISTIC_WEAK: // Characteristic Weak Boundary Condition
            Apply_Characteristic_Weak_Boundary_Condition(Iteration);
            break;
        default:
            error("Apply_Boundary_Condition: Invalid Boundary Condition Method - %d", BCMethod);
            break;
    }
}

