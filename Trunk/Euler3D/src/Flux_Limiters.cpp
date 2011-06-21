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
//! Higher Order Limiter: Barth and Jespersen
//! Limiter = 1
//------------------------------------------------------------------------------
void Compute_Limiter_Barth_Jespersen(int node_L, int node_R, double *Phi_L, double *Phi_R) {
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
        if (LimiterOrder == 2) {
            Q_L[0]  = Q1[node_L] + 0.5*(Q1[nodeid] - Q1[node_L]) + 0.5*(Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L);
            Q_L[1]  = Q2[node_L] + 0.5*(Q2[nodeid] - Q2[node_L]) + 0.5*(Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L);
            Q_L[2]  = Q3[node_L] + 0.5*(Q3[nodeid] - Q3[node_L]) + 0.5*(Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L);
            Q_L[3]  = Q4[node_L] + 0.5*(Q4[nodeid] - Q4[node_L]) + 0.5*(Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L);
            Q_L[4]  = Q5[node_L] + 0.5*(Q5[nodeid] - Q5[node_L]) + 0.5*(Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L);
        } else {
            Q_L[0]  = Q1[node_L] + Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L;
            Q_L[1]  = Q2[node_L] + Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L;
            Q_L[2]  = Q3[node_L] + Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L;
            Q_L[3]  = Q4[node_L] + Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L;
            Q_L[4]  = Q5[node_L] + Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L;
        }
        
        // Compute the limiters for Left Q's
        for (i = 0; i < 5; i++) {
            if ((Q_L[i] - Qtmp[i]) > 0.0)
                Phi_L[i] = MIN(Phi_L[i], MIN(1.0, (Qmax[i] - Qtmp[i]) / (Q_L[i] - Qtmp[i])));
            else if ((Q_L[i] - Qtmp[i]) < 0.0)
                Phi_L[i] = MIN(Phi_L[i], MIN(1.0, (Qmin[i] - Qtmp[i]) / (Q_L[i] - Qtmp[i])));
            else
                Phi_L[i] = MIN(Phi_L[i], 1.0);
            Phi_L[i] = MAX(0.0, Phi_L[i]);
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
        if (LimiterOrder == 2) {
            Q_R[0]  = Q1[node_R] + 0.5*(Q1[nodeid] - Q1[node_R]) + 0.5*(Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R);
            Q_R[1]  = Q2[node_R] + 0.5*(Q2[nodeid] - Q2[node_R]) + 0.5*(Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R);
            Q_R[2]  = Q3[node_R] + 0.5*(Q3[nodeid] - Q3[node_R]) + 0.5*(Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R);
            Q_R[3]  = Q4[node_R] + 0.5*(Q4[nodeid] - Q4[node_R]) + 0.5*(Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R);
            Q_R[4]  = Q5[node_R] + 0.5*(Q5[nodeid] - Q5[node_R]) + 0.5*(Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R);
        } else {
            Q_R[0]  = Q1[node_R] + Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R;
            Q_R[1]  = Q2[node_R] + Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R;
            Q_R[2]  = Q3[node_R] + Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R;
            Q_R[3]  = Q4[node_R] + Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R;
            Q_R[4]  = Q5[node_R] + Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R;
        }
        
        // Compute the limiters for Left Q's
        for (i = 0; i < 5; i++) {
            if ((Q_R[i] - Qtmp[i]) > 0.0)
                Phi_R[i] = MIN(Phi_R[i], MIN(1.0, (Qmax[i] - Qtmp[i]) / (Q_R[i] - Qtmp[i])));
            else if ((Q_R[i] - Qtmp[i]) < 0.0)
                Phi_R[i] = MIN(Phi_R[i], MIN(1.0, (Qmin[i] - Qtmp[i]) / (Q_R[i] - Qtmp[i])));
            else
                Phi_R[i] = MIN(Phi_R[i], 1.0);
            Phi_R[i] = MAX(0.0, Phi_R[i]);
        }
    }

    // Check if Smooth Limiter is Reqested - All Phi's for All variable will be same
    if (LimiterSmooth == 1) {
        // Get the Minimum Phi
        for (i = 1; i < 5; i++) {
            Phi_L[0] = MIN(Phi_L[0], Phi_L[i]);
            Phi_R[0] = MIN(Phi_R[0], Phi_R[i]);
        }
        // Apply the Minimum Phi to all Variables
        for (i = 1; i < 5; i++) {
            Phi_L[i] = Phi_L[0];
            Phi_R[i] = Phi_R[0];
        }
    }

    // Apply Pressure Correction To Limiters
    Compute_Limiter_PressureCorrection(node_L, node_R, Phi_L, Phi_R);
}

//------------------------------------------------------------------------------
//! Higher Order Limiter: Venkatakrishnan
//! Limiter = 2
//------------------------------------------------------------------------------
void Compute_Limiter_Venkatakrishnan(int node_L, int node_R, double *Phi_L, double *Phi_R) {
    int i, n, nodeid;
    double temp;
    double eps2, deltaP, deltaM;
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

    // Compute length scale as a diameter of a sphere that has the same volume
    // as the control volume
    eps2 = pow(Venkat_KThreshold, 3.0);
    eps2 = 6.0*eps2*cVolume[node_L]/M_PI;
    
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
        if (LimiterOrder == 2) {
            Q_L[0]  = Q1[node_L] + 0.5*(Q1[nodeid] - Q1[node_L]) + 0.5*(Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L);
            Q_L[1]  = Q2[node_L] + 0.5*(Q2[nodeid] - Q2[node_L]) + 0.5*(Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L);
            Q_L[2]  = Q3[node_L] + 0.5*(Q3[nodeid] - Q3[node_L]) + 0.5*(Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L);
            Q_L[3]  = Q4[node_L] + 0.5*(Q4[nodeid] - Q4[node_L]) + 0.5*(Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L);
            Q_L[4]  = Q5[node_L] + 0.5*(Q5[nodeid] - Q5[node_L]) + 0.5*(Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L);
        } else {
            Q_L[0]  = Q1[node_L] + Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L;
            Q_L[1]  = Q2[node_L] + Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L;
            Q_L[2]  = Q3[node_L] + Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L;
            Q_L[3]  = Q4[node_L] + Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L;
            Q_L[4]  = Q5[node_L] + Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L;
        }
        
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

            Phi_L[i] = MAX(0.0, Phi_L[i]);
        }
    }

    // Compute Limiter for Right Node
    Qmin[0] = Qmax[0] = Q1[node_R];
    Qmin[1] = Qmax[1] = Q2[node_R];
    Qmin[2] = Qmax[2] = Q3[node_R];
    Qmin[3] = Qmax[3] = Q4[node_R];
    Qmin[4] = Qmax[4] = Q5[node_R];

    // Compute length scale as a diameter of a sphere that has the same volume
    // as the control volume
    eps2 = pow(Venkat_KThreshold, 3.0);
    eps2 = 6.0*eps2*cVolume[node_R]/M_PI;

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
        if (LimiterOrder == 2) {
            Q_R[0]  = Q1[node_R] + 0.5*(Q1[nodeid] - Q1[node_R]) + 0.5*(Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R);
            Q_R[1]  = Q2[node_R] + 0.5*(Q2[nodeid] - Q2[node_R]) + 0.5*(Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R);
            Q_R[2]  = Q3[node_R] + 0.5*(Q3[nodeid] - Q3[node_R]) + 0.5*(Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R);
            Q_R[3]  = Q4[node_R] + 0.5*(Q4[nodeid] - Q4[node_R]) + 0.5*(Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R);
            Q_R[4]  = Q5[node_R] + 0.5*(Q5[nodeid] - Q5[node_R]) + 0.5*(Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R);
        } else {
            Q_R[0]  = Q1[node_R] + Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R;
            Q_R[1]  = Q2[node_R] + Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R;
            Q_R[2]  = Q3[node_R] + Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R;
            Q_R[3]  = Q4[node_R] + Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R;
            Q_R[4]  = Q5[node_R] + Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R;
        }
        
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

            Phi_R[i] = MAX(0.0, Phi_R[i]);
        }
    }

    // Check if Smooth Limiter is Reqested - All Phi's for All variable will be same
    if (LimiterSmooth == 1) {
        // Get the Minimum Phi
        for (i = 1; i < 5; i++) {
            Phi_L[0] = MIN(Phi_L[0], Phi_L[i]);
            Phi_R[0] = MIN(Phi_R[0], Phi_R[i]);
        }
        // Apply the Minimum Phi to all Variables
        for (i = 1; i < 5; i++) {
            Phi_L[i] = Phi_L[0];
            Phi_R[i] = Phi_R[0];
        }
    }

    // Apply Pressure Correction To Limiters
    Compute_Limiter_PressureCorrection(node_L, node_R, Phi_L, Phi_R);
}

//------------------------------------------------------------------------------
//! Higher Order Limiter: Pressure Correction
//! If Pressure becomes negative set Limiters = 0.0, Hence First Order Only
//------------------------------------------------------------------------------
void Compute_Limiter_PressureCorrection(int node_L, int node_R, double *Phi_L, double *Phi_R) {
    int i;
    double Rx_L, Ry_L, Rz_L;
    double Rx_R, Ry_R, Rz_R;
    double Q_L[5], Q_R[5];
    double P_L, P_R;
    
    // Compute Vector R for this edge node_L-->node_R
    Rx_L = 0.5*(coordXYZ[PHY_DIM * node_R + 0] - coordXYZ[PHY_DIM * node_L + 0]);
    Ry_L = 0.5*(coordXYZ[PHY_DIM * node_R + 1] - coordXYZ[PHY_DIM * node_L + 1]);
    Rz_L = 0.5*(coordXYZ[PHY_DIM * node_R + 2] - coordXYZ[PHY_DIM * node_L + 2]);

    // Compute second order Left Q's for this edge node_L-->node_R without limiter
    if (LimiterOrder == 2) {
        Q_L[0]  = Q1[node_L] + Phi_L[0]*(0.5*(Q1[node_R] - Q1[node_L]) + 0.5*(Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L));
        Q_L[1]  = Q2[node_L] + Phi_L[1]*(0.5*(Q2[node_R] - Q2[node_L]) + 0.5*(Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L));
        Q_L[2]  = Q3[node_L] + Phi_L[2]*(0.5*(Q3[node_R] - Q3[node_L]) + 0.5*(Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L));
        Q_L[3]  = Q4[node_L] + Phi_L[3]*(0.5*(Q4[node_R] - Q4[node_L]) + 0.5*(Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L));
        Q_L[4]  = Q5[node_L] + Phi_L[4]*(0.5*(Q5[node_R] - Q5[node_L]) + 0.5*(Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L));
    } else {
        Q_L[0]  = Q1[node_L] + Phi_L[0]*(Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L);
        Q_L[1]  = Q2[node_L] + Phi_L[1]*(Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L);
        Q_L[2]  = Q3[node_L] + Phi_L[2]*(Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L);
        Q_L[3]  = Q4[node_L] + Phi_L[3]*(Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L);
        Q_L[4]  = Q5[node_L] + Phi_L[4]*(Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L);
    }

    // Compute Vector R for this edge node_R-->nodeid
    Rx_R = - Rx_L;
    Ry_R = - Ry_L;
    Rz_R = - Rz_L;

    // Compute second order Right Q's for this edge node_R-->node_L without limiter
    if (LimiterOrder == 2) {
        Q_R[0]  = Q1[node_R] + Phi_R[0]*(0.5*(Q1[node_L] - Q1[node_R]) + 0.5*(Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R));
        Q_R[1]  = Q2[node_R] + Phi_R[1]*(0.5*(Q2[node_L] - Q2[node_R]) + 0.5*(Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R));
        Q_R[2]  = Q3[node_R] + Phi_R[2]*(0.5*(Q3[node_L] - Q3[node_R]) + 0.5*(Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R));
        Q_R[3]  = Q4[node_R] + Phi_R[3]*(0.5*(Q4[node_L] - Q4[node_R]) + 0.5*(Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R));
        Q_R[4]  = Q5[node_R] + Phi_R[4]*(0.5*(Q5[node_L] - Q5[node_R]) + 0.5*(Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R));
    } else {
        Q_R[0]  = Q1[node_R] + Phi_R[0]*(Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R);
        Q_R[1]  = Q2[node_R] + Phi_R[1]*(Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R);
        Q_R[2]  = Q3[node_R] + Phi_R[2]*(Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R);
        Q_R[3]  = Q4[node_R] + Phi_R[3]*(Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R);
        Q_R[4]  = Q5[node_R] + Phi_R[4]*(Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R);
    }

    // Apply Pressure Correction on Limiter
    // Compute Roe Variables
    if (SolverScheme == SOLVER_ROE) {
        double Q_Roe[5], P_Roe;
        Compute_RoeVariables(Q_L, Q_R, Q_Roe);
        P_Roe = (Gamma - 1.0)*(Q_Roe[4] - (0.5/Q_Roe[0])*(Q_Roe[1]*Q_Roe[1] + Q_Roe[2]*Q_Roe[2] + Q_Roe[3]*Q_Roe[3]));
        if (P_Roe < 0.0) {
            for (i = 0; i < 5; i++) {
                Phi_L[i] = 0.0;
                Phi_R[i] = 0.0;
            }
        }
    }
    // Left
    P_L = (Gamma - 1.0)*(Q_L[4] - (0.5/Q_L[0])*(Q_L[1]*Q_L[1] + Q_L[2]*Q_L[2] + Q_L[3]*Q_L[3]));
    if (P_L < 0.0) {
        for (i = 0; i < 5; i++)
            Phi_L[i] = 0.0;
    }
    // Right
    P_R = (Gamma - 1.0)*(Q_R[4] - (0.5/Q_R[0])*(Q_R[1]*Q_R[1] + Q_R[2]*Q_R[2] + Q_R[3]*Q_R[3]));
    if (P_R < 0.0) {
        for (i = 0; i < 5; i++)
            Phi_R[i] = 0.0;
    }
}

