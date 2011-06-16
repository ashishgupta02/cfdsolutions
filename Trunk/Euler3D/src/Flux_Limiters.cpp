/*******************************************************************************
 * File:        Flux_Limiters.cpp
 * Author:      Ashish Gupta
 * Revision:    2
 ******************************************************************************/

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! Flux Limiter: Barth and Jespersen
//! Limiter = 1
//------------------------------------------------------------------------------
void Compute_Flux_Limiter_Barth_Jespersen(int node_L, int node_R, double *Phi_L, double *Phi_R) {
    int i, n, nodeid;
    double Rx_L, Ry_L, Rz_L;
    double Rx_R, Ry_R, Rz_R;
    double Q_L[5], Q_R[5], Qtmp[5], Qmin[5], Qmax[5];

    // Initialize Flux Limiters
    Phi_L[0] = Phi_L[1] = Phi_L[2] = Phi_L[3] = Phi_L[4] = 1.0;
    Phi_R[0] = Phi_R[1] = Phi_R[2] = Phi_R[3] = Phi_R[4] = 1.0;

    // Compute Limiter for Left Node
    Qmin[0] = Qmax[0] = Q1[node_L];
    Qmin[1] = Qmax[1] = Q2[node_L];
    Qmin[2] = Qmax[2] = Q3[node_L];
    Qmin[3] = Qmax[3] = Q4[node_L];
    Qmin[4] = Qmax[4] = Q5[node_L];
    // Compute Qmin and Qmax around Left Node
    for (n = crs_IA_Node2Node[node_L]; n < crs_IA_Node2Node[node_L + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Get Minimum
        Qmin[0] = MIN(Qmin[0], Q1[nodeid]);
        Qmin[1] = MIN(Qmin[1], Q2[nodeid]);
        Qmin[2] = MIN(Qmin[2], Q3[nodeid]);
        Qmin[3] = MIN(Qmin[3], Q4[nodeid]);
        Qmin[4] = MIN(Qmin[4], Q5[nodeid]);
        // Get Maximum
        Qmax[0] = MAX(Qmax[0], Q1[nodeid]);
        Qmax[1] = MAX(Qmax[1], Q2[nodeid]);
        Qmax[2] = MAX(Qmax[2], Q3[nodeid]);
        Qmax[3] = MAX(Qmax[3], Q4[nodeid]);
        Qmax[4] = MAX(Qmax[4], Q5[nodeid]);
    }
    // Compute the Minimum Phi_L for all edges around Left Node
    Qtmp[0] = Q1[node_L];
    Qtmp[1] = Q2[node_L];
    Qtmp[2] = Q3[node_L];
    Qtmp[3] = Q4[node_L];
    Qtmp[4] = Q5[node_L];
    for (n = crs_IA_Node2Node[node_L]; n < crs_IA_Node2Node[node_L + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Compute Vector R for this edge node_L-->nodeid
        Rx_L = 0.5*(coordXYZ[PHY_DIM * nodeid + 0] - coordXYZ[PHY_DIM * node_L + 0]);
        Ry_L = 0.5*(coordXYZ[PHY_DIM * nodeid + 1] - coordXYZ[PHY_DIM * node_L + 1]);
        Rz_L = 0.5*(coordXYZ[PHY_DIM * nodeid + 2] - coordXYZ[PHY_DIM * node_L + 2]);
        // Compute second order Left Q's for this edge node_L-->nodeid without limiter
        Q_L[0]  = Q1[node_L] + Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L;
        Q_L[1]  = Q2[node_L] + Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L;
        Q_L[2]  = Q3[node_L] + Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L;
        Q_L[3]  = Q4[node_L] + Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L;
        Q_L[4]  = Q5[node_L] + Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L;
        // Compute the limiters for Left Q's
        for (i = 0; i < 5; i++) {
            if ((Q_L[i] - Qtmp[i]) > 0.0)
                Phi_L[i] = MIN(Phi_L[i], MIN(1.0, (Qmax[i] - Qtmp[i]) / (Q_L[i] - Qtmp[i])));
            else if ((Q_L[i] - Qtmp[i]) < 0.0)
                Phi_L[i] = MIN(Phi_L[i], MIN(1.0, (Qmin[i] - Qtmp[i]) / (Q_L[i] - Qtmp[i])));
            else
                Phi_L[i] = MIN(Phi_L[i], 1.0);
        }
    }

    // Compute Limiter for Right Node
    Qmin[0] = Qmax[0] = Q1[node_R];
    Qmin[1] = Qmax[1] = Q2[node_R];
    Qmin[2] = Qmax[2] = Q3[node_R];
    Qmin[3] = Qmax[3] = Q4[node_R];
    Qmin[4] = Qmax[4] = Q5[node_R];
    // Compute Qmin and Qmax around Right Node
    for (n = crs_IA_Node2Node[node_R]; n < crs_IA_Node2Node[node_R + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Get the minimum
        Qmin[0] = MIN(Qmin[0], Q1[nodeid]);
        Qmin[1] = MIN(Qmin[1], Q2[nodeid]);
        Qmin[2] = MIN(Qmin[2], Q3[nodeid]);
        Qmin[3] = MIN(Qmin[3], Q4[nodeid]);
        Qmin[4] = MIN(Qmin[4], Q5[nodeid]);
        // Get the maximum
        Qmax[0] = MAX(Qmax[0], Q1[nodeid]);
        Qmax[1] = MAX(Qmax[1], Q2[nodeid]);
        Qmax[2] = MAX(Qmax[2], Q3[nodeid]);
        Qmax[3] = MAX(Qmax[3], Q4[nodeid]);
        Qmax[4] = MAX(Qmax[4], Q5[nodeid]);
    }
    // Compute the Minimum Phi_R for all edges around Right Node
    Qtmp[0] = Q1[node_R];
    Qtmp[1] = Q2[node_R];
    Qtmp[2] = Q3[node_R];
    Qtmp[3] = Q4[node_R];
    Qtmp[4] = Q5[node_R];
    for (n = crs_IA_Node2Node[node_R]; n < crs_IA_Node2Node[node_R + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Compute Vector R for this edge node_R-->nodeid
        Rx_R = 0.5*(coordXYZ[PHY_DIM * nodeid + 0] - coordXYZ[PHY_DIM * node_R + 0]);
        Ry_R = 0.5*(coordXYZ[PHY_DIM * nodeid + 1] - coordXYZ[PHY_DIM * node_R + 1]);
        Rz_R = 0.5*(coordXYZ[PHY_DIM * nodeid + 2] - coordXYZ[PHY_DIM * node_R + 2]);
        // Compute second order Right Q's for this edge node_R-->nodeid without limiter
        Q_R[0]  = Q1[node_R] + Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R;
        Q_R[1]  = Q2[node_R] + Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R;
        Q_R[2]  = Q3[node_R] + Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R;
        Q_R[3]  = Q4[node_R] + Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R;
        Q_R[4]  = Q5[node_R] + Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R;
        // Compute the limiters for Left Q's
        for (i = 0; i < 5; i++) {
            if ((Q_R[i] - Qtmp[i]) > 0.0)
                Phi_R[i] = MIN(Phi_R[i], MIN(1.0, (Qmax[i] - Qtmp[i]) / (Q_R[i] - Qtmp[i])));
            else if ((Q_R[i] - Qtmp[i]) < 0.0)
                Phi_R[i] = MIN(Phi_R[i], MIN(1.0, (Qmin[i] - Qtmp[i]) / (Q_R[i] - Qtmp[i])));
            else
                Phi_R[i] = MIN(Phi_R[i], 1.0);
        }
    }
}

//------------------------------------------------------------------------------
//! Flux Limiter: Venkatakrishnan
//! Limiter = 2
//------------------------------------------------------------------------------
void Compute_Flux_Limiter_Venkatakrishnan(int node_L, int node_R, double *Phi_L, double *Phi_R) {
    int i, n, nodeid;
    double temp, K_threshold;
    double eps2, deltaP, deltaM;
    double Rx_L, Ry_L, Rz_L;
    double Rx_R, Ry_R, Rz_R;
    double Q_L[5], Q_R[5], Qtmp[5], Qmin[5], Qmax[5];

    // Initialize Flux Limiters
    Phi_L[0] = Phi_L[1] = Phi_L[2] = Phi_L[3] = Phi_L[4] = 1.0;
    Phi_R[0] = Phi_R[1] = Phi_R[2] = Phi_R[3] = Phi_R[4] = 1.0;
    K_threshold = 1.0;

    // Compute Limiter for Left Node
    Qmin[0] = Qmax[0] = Q1[node_L];
    Qmin[1] = Qmax[1] = Q2[node_L];
    Qmin[2] = Qmax[2] = Q3[node_L];
    Qmin[3] = Qmax[3] = Q4[node_L];
    Qmin[4] = Qmax[4] = Q5[node_L];

    eps2 = pow(K_threshold, 3.0);
    eps2 = eps2*cVolume[node_L];
    
    // Compute Qmin and Qmax around Left Node
    for (n = crs_IA_Node2Node[node_L]; n < crs_IA_Node2Node[node_L + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Get Minimum
        Qmin[0] = MIN(Qmin[0], Q1[nodeid]);
        Qmin[1] = MIN(Qmin[1], Q2[nodeid]);
        Qmin[2] = MIN(Qmin[2], Q3[nodeid]);
        Qmin[3] = MIN(Qmin[3], Q4[nodeid]);
        Qmin[4] = MIN(Qmin[4], Q5[nodeid]);
        // Get Maximum
        Qmax[0] = MAX(Qmax[0], Q1[nodeid]);
        Qmax[1] = MAX(Qmax[1], Q2[nodeid]);
        Qmax[2] = MAX(Qmax[2], Q3[nodeid]);
        Qmax[3] = MAX(Qmax[3], Q4[nodeid]);
        Qmax[4] = MAX(Qmax[4], Q5[nodeid]);
    }
    // Compute the Minimum Phi_L for all edges around Left Node
    Qtmp[0] = Q1[node_L];
    Qtmp[1] = Q2[node_L];
    Qtmp[2] = Q3[node_L];
    Qtmp[3] = Q4[node_L];
    Qtmp[4] = Q5[node_L];
    for (n = crs_IA_Node2Node[node_L]; n < crs_IA_Node2Node[node_L + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Compute Vector R for this edge node_L-->nodeid
        Rx_L = 0.5*(coordXYZ[PHY_DIM * nodeid + 0] - coordXYZ[PHY_DIM * node_L + 0]);
        Ry_L = 0.5*(coordXYZ[PHY_DIM * nodeid + 1] - coordXYZ[PHY_DIM * node_L + 1]);
        Rz_L = 0.5*(coordXYZ[PHY_DIM * nodeid + 2] - coordXYZ[PHY_DIM * node_L + 2]);
        // Compute second order Left Q's for this edge node_L-->nodeid without limiter
        Q_L[0]  = Q1[node_L] + Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L;
        Q_L[1]  = Q2[node_L] + Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L;
        Q_L[2]  = Q3[node_L] + Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L;
        Q_L[3]  = Q4[node_L] + Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L;
        Q_L[4]  = Q5[node_L] + Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L;
        // Compute the limiters for Left Q's
        for (i = 0; i < 5; i++) {
            if ((Q_L[i] - Qtmp[i]) > 0.0) {
                deltaP = Qmax[i] - Qtmp[i];
                deltaM = Q_L[i] - Qtmp[i];
                temp   = (deltaP*deltaP + eps2) + 2.0*(deltaM*deltaP);
                temp   = temp/(deltaP*deltaP + 2.0*deltaM*deltaM + deltaM*deltaP + eps2);
                Phi_L[i] = MIN(Phi_L[i], temp);
            } else if ((Q_L[i] - Qtmp[i]) < 0.0) {
                deltaP = Qmin[i] - Qtmp[i];
                deltaM = Q_L[i] - Qtmp[i];
                temp   = (deltaP*deltaP + eps2) + 2.0*(deltaM*deltaP);
                temp   = temp/(deltaP*deltaP + 2.0*deltaM*deltaM + deltaM*deltaP + eps2);
                Phi_L[i] = MIN(Phi_L[i], temp);
            } else
                Phi_L[i] = MIN(Phi_L[i], 1.0);
        }
    }

    // Compute Limiter for Right Node
    Qmin[0] = Qmax[0] = Q1[node_R];
    Qmin[1] = Qmax[1] = Q2[node_R];
    Qmin[2] = Qmax[2] = Q3[node_R];
    Qmin[3] = Qmax[3] = Q4[node_R];
    Qmin[4] = Qmax[4] = Q5[node_R];

    eps2 = pow(K_threshold, 3.0);
    eps2 = eps2*cVolume[node_R];
    
    // Compute Qmin and Qmax around Right Node
    for (n = crs_IA_Node2Node[node_R]; n < crs_IA_Node2Node[node_R + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Get the minimum
        Qmin[0] = MIN(Qmin[0], Q1[nodeid]);
        Qmin[1] = MIN(Qmin[1], Q2[nodeid]);
        Qmin[2] = MIN(Qmin[2], Q3[nodeid]);
        Qmin[3] = MIN(Qmin[3], Q4[nodeid]);
        Qmin[4] = MIN(Qmin[4], Q5[nodeid]);
        // Get the maximum
        Qmax[0] = MAX(Qmax[0], Q1[nodeid]);
        Qmax[1] = MAX(Qmax[1], Q2[nodeid]);
        Qmax[2] = MAX(Qmax[2], Q3[nodeid]);
        Qmax[3] = MAX(Qmax[3], Q4[nodeid]);
        Qmax[4] = MAX(Qmax[4], Q5[nodeid]);
    }
    // Compute the Minimum Phi_R for all edges around Right Node
    Qtmp[0] = Q1[node_R];
    Qtmp[1] = Q2[node_R];
    Qtmp[2] = Q3[node_R];
    Qtmp[3] = Q4[node_R];
    Qtmp[4] = Q5[node_R];
    for (n = crs_IA_Node2Node[node_R]; n < crs_IA_Node2Node[node_R + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Compute Vector R for this edge node_R-->nodeid
        Rx_R = 0.5*(coordXYZ[PHY_DIM * nodeid + 0] - coordXYZ[PHY_DIM * node_R + 0]);
        Ry_R = 0.5*(coordXYZ[PHY_DIM * nodeid + 1] - coordXYZ[PHY_DIM * node_R + 1]);
        Rz_R = 0.5*(coordXYZ[PHY_DIM * nodeid + 2] - coordXYZ[PHY_DIM * node_R + 2]);
        // Compute second order Right Q's for this edge node_R-->nodeid without limiter
        Q_R[0]  = Q1[node_R] + Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R;
        Q_R[1]  = Q2[node_R] + Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R;
        Q_R[2]  = Q3[node_R] + Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R;
        Q_R[3]  = Q4[node_R] + Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R;
        Q_R[4]  = Q5[node_R] + Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R;
        // Compute the limiters for Left Q's
        for (i = 0; i < 5; i++) {
            if ((Q_R[i] - Qtmp[i]) > 0.0) {
                deltaP = Qmax[i] - Qtmp[i];
                deltaM = Q_R[i] - Qtmp[i];
                temp   = (deltaP*deltaP + eps2) + 2.0*(deltaM*deltaP);
                temp   = temp/(deltaP*deltaP + 2.0*deltaM*deltaM + deltaM*deltaP + eps2);
                Phi_R[i] = MIN(Phi_R[i], temp);
            } else if ((Q_R[i] - Qtmp[i]) < 0.0) {
                deltaP = Qmin[i] - Qtmp[i];
                deltaM = Q_R[i] - Qtmp[i];
                temp   = (deltaP*deltaP + eps2) + 2.0*(deltaM*deltaP);
                temp   = temp/(deltaP*deltaP + 2.0*deltaM*deltaM + deltaM*deltaP + eps2);
                Phi_R[i] = MIN(Phi_R[i], temp);
            } else
                Phi_R[i] = MIN(Phi_R[i], 1.0);
        }
    }
}

