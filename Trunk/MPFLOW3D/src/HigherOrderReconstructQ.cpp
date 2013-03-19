/*******************************************************************************
 * File:        HigherOrderReconstructQ.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! Compute Second Order Recontruction of Q
//! Limiting Option and Pressure Corrections are applied
//------------------------------------------------------------------------------
void Compute_SecondOrderReconstructQ(int node_L, int node_R, double *Q_L, double *Q_R) {
    double Phi_L[5], Phi_R[5];
    double Rx_L, Ry_L, Rz_L;
    double Rx_R, Ry_R, Rz_R;
    double frac;

    // Initialize
    frac = 0.0;
    if (LimiterOrder == LIMITER_ORDER_SECOND)
        frac = 0.5;
    
    // Check if Limiting is applied
    if (LimiterMethod != LIMITER_METHOD_NONE) {
        // Left
        Phi_L[0] = Limiter_Phi1[node_L];
        Phi_L[1] = Limiter_Phi2[node_L];
        Phi_L[2] = Limiter_Phi3[node_L];
        Phi_L[3] = Limiter_Phi4[node_L];
        Phi_L[4] = Limiter_Phi5[node_L];
        // Right
        Phi_R[0] = Limiter_Phi1[node_R];
        Phi_R[1] = Limiter_Phi2[node_R];
        Phi_R[2] = Limiter_Phi3[node_R];
        Phi_R[3] = Limiter_Phi4[node_R];
        Phi_R[4] = Limiter_Phi5[node_R];
    } else {
        Phi_L[0] = Phi_L[1] = Phi_L[2] = Phi_L[3] = Phi_L[4] = 1.0;
        Phi_R[0] = Phi_R[1] = Phi_R[2] = Phi_R[3] = Phi_R[4] = 1.0;
    }

    // Compute vector R
    Rx_L = 0.5*(coordXYZ[PHY_DIM * node_R + 0] - coordXYZ[PHY_DIM * node_L + 0]);
    Ry_L = 0.5*(coordXYZ[PHY_DIM * node_R + 1] - coordXYZ[PHY_DIM * node_L + 1]);
    Rz_L = 0.5*(coordXYZ[PHY_DIM * node_R + 2] - coordXYZ[PHY_DIM * node_L + 2]);
    // Compute Left Q's
    Q_L[0]  = Q1[node_L] + Phi_L[0]*(frac*(Q1[node_R] - Q1[node_L]) + (1.0 - frac)*(Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L));
    Q_L[1]  = Q2[node_L] + Phi_L[1]*(frac*(Q2[node_R] - Q2[node_L]) + (1.0 - frac)*(Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L));
    Q_L[2]  = Q3[node_L] + Phi_L[2]*(frac*(Q3[node_R] - Q3[node_L]) + (1.0 - frac)*(Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L));
    Q_L[3]  = Q4[node_L] + Phi_L[3]*(frac*(Q4[node_R] - Q4[node_L]) + (1.0 - frac)*(Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L));
    Q_L[4]  = Q5[node_L] + Phi_L[4]*(frac*(Q5[node_R] - Q5[node_L]) + (1.0 - frac)*(Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L));
    
    // Compute vector R
    Rx_R = - Rx_L;
    Ry_R = - Ry_L;
    Rz_R = - Rz_L;
    // Finally Compute Right Q's
    Q_R[0]  = Q1[node_R] + Phi_R[0]*(frac*(Q1[node_L] - Q1[node_R]) + (1.0 - frac)*(Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R));
    Q_R[1]  = Q2[node_R] + Phi_R[1]*(frac*(Q2[node_L] - Q2[node_R]) + (1.0 - frac)*(Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R));
    Q_R[2]  = Q3[node_R] + Phi_R[2]*(frac*(Q3[node_L] - Q3[node_R]) + (1.0 - frac)*(Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R));
    Q_R[3]  = Q4[node_R] + Phi_R[3]*(frac*(Q4[node_L] - Q4[node_R]) + (1.0 - frac)*(Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R));
    Q_R[4]  = Q5[node_R] + Phi_R[4]*(frac*(Q5[node_L] - Q5[node_R]) + (1.0 - frac)*(Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R));
}

