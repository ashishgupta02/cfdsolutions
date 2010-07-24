/*******************************************************************************
 * File:        Time_Step.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "Commons.h"
#include "Solver.h"
#include <assert.h>

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Compute_DeltaT(int Iteration) {
    int nodeL, nodeR;
    double denom = 0.0;
    
    Vector3D areavec;
    double   area;

    double rho_L, u_L, v_L, w_L, et_L, e_L, p_L, c_L;
    double rho_R, u_R, v_R, w_R, et_R, e_R, p_R, c_R;

    double lamda4, lamda5;
    double max_lamda;
    double rho_avg, p_avg, e_avg, u_avg, v_avg, w_avg, c_avg;

    // Internal Edges
    for (int i = 0; i < nEdge; i++) {
        // Get two nodes of edge
        nodeL = int_edge_info[i].node[0];
        nodeR = int_edge_info[i].node[1];

        assert(nodeR > nodeL);
        
        // Get area vector
        areavec = int_edge_info[i].areav;
        area  = areavec.magnitude();
        areavec.normalize();
        
        // Left Node
        rho_L = Q1[nodeL];
        u_L   = Q2[nodeL] / rho_L;
        v_L   = Q3[nodeL] / rho_L;
        w_L   = Q4[nodeL] / rho_L;
        et_L  = Q5[nodeL] / rho_L;
        e_L   = et_L - 0.5 * (u_L * u_L + v_L * v_L + w_L * w_L);
        p_L   = (Gamma - 1.0) * rho_L * e_L;
        c_L   = sqrt((Gamma * p_L) / rho_L);

        // Right Node
        rho_R = Q1[nodeR];
        u_R   = Q2[nodeR] / rho_R;
        v_R   = Q3[nodeR] / rho_R;
        w_R   = Q4[nodeR] / rho_R;
        et_R  = Q5[nodeR] / rho_R;
        e_R   = et_R - 0.5 * (u_R * u_R + v_R * v_R + w_R * w_R);
        p_R   = (Gamma - 1.0) * rho_R * e_R;
        c_R   = sqrt((Gamma * p_R) / rho_R);

        // Get averaged variables
        rho_avg = 0.5*(rho_L + rho_R);
        e_avg   = 0.5*(e_L   + e_R);
        u_avg   = 0.5*(u_L   + u_R);
        v_avg   = 0.5*(v_L   + v_R);
        w_avg   = 0.5*(w_L   + w_R);
        p_avg   = 0.5*(p_L   + p_R);
        c_avg   = 0.5*(c_L   + c_R);
        
        // Find maximum eigenvalue for nodeL
        lamda4 = fabs((u_avg * areavec.vec[0] + v_avg * areavec.vec[1] + w_avg * areavec.vec[2]) + c_avg);
        lamda5 = fabs((u_avg * areavec.vec[0] + v_avg * areavec.vec[1] + w_avg * areavec.vec[2]) - c_avg);
        max_lamda = MAX(lamda4, lamda5);

        denom = (max_lamda * area);

        // Add to Each Node
        DeltaT[nodeL] += denom;
        DeltaT[nodeR] += denom;
    }

    // Boundary Edges
    for (int i = 0; i < nBEdge; i++) {
        // Get two nodes of edge
        nodeL = bndry_edge_info[i].node[0];
        nodeR = bndry_edge_info[i].node[1];

        assert(nodeR > nodeL);
        
        // Get area vector
        areavec = bndry_edge_info[i].areav;
        area    = areavec.magnitude();
        areavec.normalize();

        // Note we only need to accumulate for the left node
        // Left Node
        rho_L = Q1[nodeL];
        u_L   = Q2[nodeL] / rho_L;
        v_L   = Q3[nodeL] / rho_L;
        w_L   = Q4[nodeL] / rho_L;
        et_L  = Q5[nodeL] / rho_L;
        e_L   = et_L - 0.5 * (u_L * u_L + v_L * v_L + w_L * w_L);
        p_L   = (Gamma - 1.0) * rho_L * e_L;
        c_L   = sqrt((Gamma * p_L) / rho_L);

        // Right Node
        rho_R = Q1[nodeR];
        u_R   = Q2[nodeR] / rho_R;
        v_R   = Q3[nodeR] / rho_R;
        w_R   = Q4[nodeR] / rho_R;
        et_R  = Q5[nodeR] / rho_R;
        e_R   = et_R - 0.5 * (u_R * u_R + v_R * v_R + w_R * w_R);
        p_R   = (Gamma - 1.0) * rho_R * e_R;
        c_R   = sqrt((Gamma * p_R) / rho_R);

        // Get averaged variables
        rho_avg = 0.5*(rho_L + rho_R);
        e_avg   = 0.5*(e_L   + e_R);
        u_avg   = 0.5*(u_L   + u_R);
        v_avg   = 0.5*(v_L   + v_R);
        w_avg   = 0.5*(w_L   + w_R);
        p_avg   = 0.5*(p_L   + p_R);
        c_avg   = 0.5*(c_L   + c_R);

        // Find maximum eigenvalue for nodeL
        lamda4 = fabs((u_avg * areavec.vec[0] + v_avg * areavec.vec[1] + w_avg * areavec.vec[2]) + c_avg);
        lamda5 = fabs((u_avg * areavec.vec[0] + v_avg * areavec.vec[1] + w_avg * areavec.vec[2]) - c_avg);
        max_lamda = MAX(lamda4, lamda5);

        denom = (max_lamda * area);
        
        // Add to Left Node Only
        DeltaT[nodeL] += denom;
    }
    
    // Finally Compute the Local Time Stepping
    if (CFL_Ramp > 1) {
        if (Iteration < CFL_Ramp)
            CFL = CFL_MIN + (CFL_MAX - CFL_MIN)*(((double)Iteration)/((double)(CFL_Ramp-1)));
        else
            CFL = CFL_MAX;
    } else
        CFL = CFL_MAX;

//    double dt_max = DBL_MIN;
//    double dt_min = DBL_MAX;
    for (int i = 0; i < nNode; i++) {
        DeltaT[i] = cVolume[i] * CFL / DeltaT[i];
//        dt_min = MIN(dt_min, DeltaT[i]);
//        dt_max = MAX(dt_max, DeltaT[i]);
    }
//    printf("CFL %10.5e DeltaT Min: %10.5e Max: %10.5e \n", CFL, dt_min, dt_max);
}

