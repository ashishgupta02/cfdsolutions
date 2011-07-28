/*******************************************************************************
 * File:        Limiters.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Vector3D.h"
#include "MC.h"
#include "Commons.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! Higher Order Limiter of Node
//------------------------------------------------------------------------------
void Compute_Limiter(void) {
    int i, n, nodeid, iNode;
    double Phi_C[5];

    // Initialize Phi's to 1.0
    for (iNode = 0; iNode < nNode; iNode++) {
        Limiter_Phi1[iNode] = 1.0;
        Limiter_Phi2[iNode] = 1.0;
        Limiter_Phi3[iNode] = 1.0;
        Limiter_Phi4[iNode] = 1.0;
        Limiter_Phi5[iNode] = 1.0;
    }
    
    // Compute the Limiter
    for (iNode = 0; iNode < nNode; iNode++) {
        // Initialize Flux Limiters
        Phi_C[0] = Phi_C[1] = Phi_C[2] = Phi_C[3] = Phi_C[4] = 1.0;

        // Compute Flux Limiter
        switch (Limiter) {
            case 1: // Barth and Jespersen
                Compute_Limiter_Barth_Jespersen(iNode, Phi_C);
                break;
            case 2: // Venkatakrishnan
                Compute_Limiter_Venkatakrishnan(iNode, Phi_C);
                break;
        }

        // Check if Smooth Limiter is Reqested - All Phi_C's for All variable will be same
        if (LimiterSmooth == 1) {
            // Get the Minimum Phi_C
            for (i = 1; i < 5; i++)
                Phi_C[0] = MIN(Phi_C[0], Phi_C[i]);

            // Apply the Minimum Phi_C to all Variables
            for (i = 1; i < 5; i++)
                Phi_C[i] = Phi_C[0];
        }

        // Apply Pressure Correction To Limiters
        Compute_Limiter_PressureCorrection(iNode, Phi_C);
        
        // Set the Limiter Values
        Limiter_Phi1[iNode] = Phi_C[0];
        Limiter_Phi2[iNode] = Phi_C[1];
        Limiter_Phi3[iNode] = Phi_C[2];
        Limiter_Phi4[iNode] = Phi_C[3];
        Limiter_Phi5[iNode] = Phi_C[4];
    }

    // Apply Roe Pressure Correction to Limiters
    if ((SolverScheme == SOLVER_SCHEME_ROE) || (SolverScheme == SOLVER_SCHEME_LMROE)) {
        double Phi_N[5];
        
        for (iNode = 0; iNode < nNode; iNode++) {
            // Get the Limiter Values
            Phi_C[0] = Limiter_Phi1[iNode];
            Phi_C[1] = Limiter_Phi2[iNode];
            Phi_C[2] = Limiter_Phi3[iNode];
            Phi_C[3] = Limiter_Phi4[iNode];
            Phi_C[4] = Limiter_Phi5[iNode];

            // Loop over the neighbour nodes
            for (n = crs_IA_Node2Node[iNode]; n < crs_IA_Node2Node[iNode + 1]; n++) {
                nodeid = crs_JA_Node2Node[n];
                // Get the Limiter Value for Neighbour Node
                Phi_N[0] = Limiter_Phi1[nodeid];
                Phi_N[1] = Limiter_Phi2[nodeid];
                Phi_N[2] = Limiter_Phi3[nodeid];
                Phi_N[3] = Limiter_Phi4[nodeid];
                Phi_N[4] = Limiter_Phi5[nodeid];
                
                // Compute Roe Pressure Correction
                Compute_Limiter_RoePressureCorrection(iNode, nodeid, Phi_C, Phi_N);

                // Update the Limiter Values
                Limiter_Phi1[iNode] = Phi_C[0];
                Limiter_Phi2[iNode] = Phi_C[1];
                Limiter_Phi3[iNode] = Phi_C[2];
                Limiter_Phi4[iNode] = Phi_C[3];
                Limiter_Phi5[iNode] = Phi_C[4];

                Limiter_Phi1[nodeid] = Phi_N[0];
                Limiter_Phi2[nodeid] = Phi_N[1];
                Limiter_Phi3[nodeid] = Phi_N[2];
                Limiter_Phi4[nodeid] = Phi_N[3];
                Limiter_Phi5[nodeid] = Phi_N[4];
            }
        }
    }
}

//------------------------------------------------------------------------------
//! Higher Order Limiter : Barth and Jespersen
//! Limiter = 1
//------------------------------------------------------------------------------
void Compute_Limiter_Barth_Jespersen(int node_C, double *Phi_C) {
    int i, n, nodeid;
    double Rx_C, Ry_C, Rz_C;
    double Q_C[5], Qtmp[5], Qmin[5], Qmax[5], frac;

    frac = 0.0;
    if (LimiterOrder == 2)
        frac = 0.5;
    
    // Initialize Flux Limiters
    Phi_C[0] = Phi_C[1] = Phi_C[2] = Phi_C[3] = Phi_C[4] = 1.0;

    // Compute Limiter for Current Node
    Qmin[0] = Qmax[0] = Q1[node_C];
    Qmin[1] = Qmax[1] = Q2[node_C];
    Qmin[2] = Qmax[2] = Q3[node_C];
    Qmin[3] = Qmax[3] = Q4[node_C];
    Qmin[4] = Qmax[4] = Q5[node_C];
    // Compute Qmin and Qmax around Candidate Node
    for (n = crs_IA_Node2Node[node_C]; n < crs_IA_Node2Node[node_C + 1]; n++) {
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
    // Compute the Minimum Phi_C for all edges around Candidate Node
    Qtmp[0] = Q1[node_C];
    Qtmp[1] = Q2[node_C];
    Qtmp[2] = Q3[node_C];
    Qtmp[3] = Q4[node_C];
    Qtmp[4] = Q5[node_C];
    for (n = crs_IA_Node2Node[node_C]; n < crs_IA_Node2Node[node_C + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Compute Vector R for this edge node_C-->nodeid
        Rx_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 0] - coordXYZ[PHY_DIM * node_C + 0]);
        Ry_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 1] - coordXYZ[PHY_DIM * node_C + 1]);
        Rz_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 2] - coordXYZ[PHY_DIM * node_C + 2]);
        // Compute second order Q's for this edge node_C-->nodeid without limiter
        Q_C[0]  = Q1[node_C] + frac*(Q1[nodeid] - Q1[node_C]) + (1.0 - frac)*(Q1x[node_C]*Rx_C + Q1y[node_C]*Ry_C + Q1z[node_C]*Rz_C);
        Q_C[1]  = Q2[node_C] + frac*(Q2[nodeid] - Q2[node_C]) + (1.0 - frac)*(Q2x[node_C]*Rx_C + Q2y[node_C]*Ry_C + Q2z[node_C]*Rz_C);
        Q_C[2]  = Q3[node_C] + frac*(Q3[nodeid] - Q3[node_C]) + (1.0 - frac)*(Q3x[node_C]*Rx_C + Q3y[node_C]*Ry_C + Q3z[node_C]*Rz_C);
        Q_C[3]  = Q4[node_C] + frac*(Q4[nodeid] - Q4[node_C]) + (1.0 - frac)*(Q4x[node_C]*Rx_C + Q4y[node_C]*Ry_C + Q4z[node_C]*Rz_C);
        Q_C[4]  = Q5[node_C] + frac*(Q5[nodeid] - Q5[node_C]) + (1.0 - frac)*(Q5x[node_C]*Rx_C + Q5y[node_C]*Ry_C + Q5z[node_C]*Rz_C);

        // Compute the limiters for Candidate Node Q's
        for (i = 0; i < 5; i++) {
            if ((Q_C[i] - Qtmp[i]) > 0.0)
                Phi_C[i] = MIN(Phi_C[i], MIN(1.0, (Qmax[i] - Qtmp[i]) / (Q_C[i] - Qtmp[i])));
            else if ((Q_C[i] - Qtmp[i]) < 0.0)
                Phi_C[i] = MIN(Phi_C[i], MIN(1.0, (Qmin[i] - Qtmp[i]) / (Q_C[i] - Qtmp[i])));
            else
                Phi_C[i] = MIN(Phi_C[i], 1.0);
            Phi_C[i] = MAX(0.0, Phi_C[i]);
        }
    }
}

//------------------------------------------------------------------------------
//! Higher Order Limiter : Venkatakrishnan
//! Limiter = 2
//------------------------------------------------------------------------------
void Compute_Limiter_Venkatakrishnan(int node_C, double *Phi_C) {
    int i, n, nodeid;
    double temp;
    double eps2, deltaP, deltaM, frac;
    double Rx_C, Ry_C, Rz_C;
    double Q_C[5], Qtmp[5], Qmin[5], Qmax[5];

    frac = 0.0;
    if (LimiterOrder == 2)
        frac = 0.5;
    
    // Initialize Flux Limiters
    Phi_C[0] = Phi_C[1] = Phi_C[2] = Phi_C[3] = Phi_C[4] = 1.0;

    // Compute Limiter for Candidate Node
    Qmin[0] = Qmax[0] = Q1[node_C];
    Qmin[1] = Qmax[1] = Q2[node_C];
    Qmin[2] = Qmax[2] = Q3[node_C];
    Qmin[3] = Qmax[3] = Q4[node_C];
    Qmin[4] = Qmax[4] = Q5[node_C];

    // Compute length scale as a diameter of a sphere that has the same volume
    // as the control volume
    eps2 = pow(Venkat_KThreshold, 3.0);
    eps2 = 6.0*eps2*cVolume[node_C]/M_PI;

    // Compute Qmin and Qmax around Candidate Node
    for (n = crs_IA_Node2Node[node_C]; n < crs_IA_Node2Node[node_C + 1]; n++) {
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
    // Compute the Minimum Phi_C for all edges around Candidate Node
    Qtmp[0] = Q1[node_C];
    Qtmp[1] = Q2[node_C];
    Qtmp[2] = Q3[node_C];
    Qtmp[3] = Q4[node_C];
    Qtmp[4] = Q5[node_C];
    for (n = crs_IA_Node2Node[node_C]; n < crs_IA_Node2Node[node_C + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Compute Vector R for this edge node_C-->nodeid
        Rx_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 0] - coordXYZ[PHY_DIM * node_C + 0]);
        Ry_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 1] - coordXYZ[PHY_DIM * node_C + 1]);
        Rz_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 2] - coordXYZ[PHY_DIM * node_C + 2]);

        // Compute second order Q's for this edge node_C-->nodeid without limiter
        Q_C[0]  = Q1[node_C] + frac*(Q1[nodeid] - Q1[node_C]) + (1.0 - frac)*(Q1x[node_C]*Rx_C + Q1y[node_C]*Ry_C + Q1z[node_C]*Rz_C);
        Q_C[1]  = Q2[node_C] + frac*(Q2[nodeid] - Q2[node_C]) + (1.0 - frac)*(Q2x[node_C]*Rx_C + Q2y[node_C]*Ry_C + Q2z[node_C]*Rz_C);
        Q_C[2]  = Q3[node_C] + frac*(Q3[nodeid] - Q3[node_C]) + (1.0 - frac)*(Q3x[node_C]*Rx_C + Q3y[node_C]*Ry_C + Q3z[node_C]*Rz_C);
        Q_C[3]  = Q4[node_C] + frac*(Q4[nodeid] - Q4[node_C]) + (1.0 - frac)*(Q4x[node_C]*Rx_C + Q4y[node_C]*Ry_C + Q4z[node_C]*Rz_C);
        Q_C[4]  = Q5[node_C] + frac*(Q5[nodeid] - Q5[node_C]) + (1.0 - frac)*(Q5x[node_C]*Rx_C + Q5y[node_C]*Ry_C + Q5z[node_C]*Rz_C);

        // Compute the limiters for Candidate Node Q's
        for (i = 0; i < 5; i++) {
            if ((Q_C[i] - Qtmp[i]) > 0.0) {
                deltaP = Qmax[i] - Qtmp[i];
                deltaM = Q_C[i] - Qtmp[i];
                temp   = (deltaP*deltaP + eps2) + 2.0*(deltaM*deltaP);
                temp   = temp/(deltaP*deltaP + 2.0*deltaM*deltaM + deltaM*deltaP + eps2);
                Phi_C[i] = MIN(Phi_C[i], temp);
            } else if ((Q_C[i] - Qtmp[i]) < 0.0) {
                deltaP = Qmin[i] - Qtmp[i];
                deltaM = Q_C[i] - Qtmp[i];
                temp   = (deltaP*deltaP + eps2) + 2.0*(deltaM*deltaP);
                temp   = temp/(deltaP*deltaP + 2.0*deltaM*deltaM + deltaM*deltaP + eps2);
                Phi_C[i] = MIN(Phi_C[i], temp);
            } else
                Phi_C[i] = MIN(Phi_C[i], 1.0);

            Phi_C[i] = MAX(0.0, Phi_C[i]);
        }
    }
}

//------------------------------------------------------------------------------
//! Higher Order Limiter : Pressure Correction
//! If Pressure becomes negative set Limiters = 0.0, Hence First Order Only
//------------------------------------------------------------------------------
void Compute_Limiter_PressureCorrection(int node_C, double *Phi_C) {
    int i, n, nodeid;
    double Rx_C, Ry_C, Rz_C;
    double Q_C[5];
    double P_C, frac;

    frac = 0.0;
    if (LimiterOrder == 2)
        frac = 0.5;
    
    for (n = crs_IA_Node2Node[node_C]; n < crs_IA_Node2Node[node_C + 1]; n++) {
        nodeid = crs_JA_Node2Node[n];
        // Compute Vector R for this edge node_C-->nodeid
        Rx_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 0] - coordXYZ[PHY_DIM * node_C + 0]);
        Ry_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 1] - coordXYZ[PHY_DIM * node_C + 1]);
        Rz_C = 0.5*(coordXYZ[PHY_DIM * nodeid + 2] - coordXYZ[PHY_DIM * node_C + 2]);

        // Compute second order Left Q's for this edge node_C-->nodeid without limiter
        Q_C[0]  = Q1[node_C] + Phi_C[0]*(frac*(Q1[nodeid] - Q1[node_C]) + (1.0 - frac)*(Q1x[node_C]*Rx_C + Q1y[node_C]*Ry_C + Q1z[node_C]*Rz_C));
        Q_C[1]  = Q2[node_C] + Phi_C[1]*(frac*(Q2[nodeid] - Q2[node_C]) + (1.0 - frac)*(Q2x[node_C]*Rx_C + Q2y[node_C]*Ry_C + Q2z[node_C]*Rz_C));
        Q_C[2]  = Q3[node_C] + Phi_C[2]*(frac*(Q3[nodeid] - Q3[node_C]) + (1.0 - frac)*(Q3x[node_C]*Rx_C + Q3y[node_C]*Ry_C + Q3z[node_C]*Rz_C));
        Q_C[3]  = Q4[node_C] + Phi_C[3]*(frac*(Q4[nodeid] - Q4[node_C]) + (1.0 - frac)*(Q4x[node_C]*Rx_C + Q4y[node_C]*Ry_C + Q4z[node_C]*Rz_C));
        Q_C[4]  = Q5[node_C] + Phi_C[4]*(frac*(Q5[nodeid] - Q5[node_C]) + (1.0 - frac)*(Q5x[node_C]*Rx_C + Q5y[node_C]*Ry_C + Q5z[node_C]*Rz_C));

        // Compute the Pressure
        P_C = (Gamma - 1.0)*(Q_C[4] - (0.5/Q_C[0])*(Q_C[1]*Q_C[1] + Q_C[2]*Q_C[2] + Q_C[3]*Q_C[3]));
        if (P_C < 0.0) {
            for (i = 0; i < 5; i++)
                Phi_C[i] = 0.0;
        }
    }
}


//------------------------------------------------------------------------------
//! Higher Order Limiter : Roe Pressure Correction
//! If Pressure becomes negative set Limiters = 0.0, Hence First Order Only
//------------------------------------------------------------------------------
void Compute_Limiter_RoePressureCorrection(int node_L, int node_R, double *Phi_L, double *Phi_R) {
    int i;
    double Rx_L, Ry_L, Rz_L;
    double Rx_R, Ry_R, Rz_R;
    double Q_L[5], Q_R[5], Q_Roe[5];
    double P_Roe, frac;

    frac = 0.0;
    if (LimiterOrder == 2)
        frac = 0.5;

    // Compute Vector R for this edge node_L-->node_R
    Rx_L = 0.5*(coordXYZ[PHY_DIM * node_R + 0] - coordXYZ[PHY_DIM * node_L + 0]);
    Ry_L = 0.5*(coordXYZ[PHY_DIM * node_R + 1] - coordXYZ[PHY_DIM * node_L + 1]);
    Rz_L = 0.5*(coordXYZ[PHY_DIM * node_R + 2] - coordXYZ[PHY_DIM * node_L + 2]);

    // Compute second order Left Q's for this edge node_L-->node_R without limiter
    Q_L[0]  = Q1[node_L] + Phi_L[0]*(frac*(Q1[node_R] - Q1[node_L]) + (1.0 - frac)*(Q1x[node_L]*Rx_L + Q1y[node_L]*Ry_L + Q1z[node_L]*Rz_L));
    Q_L[1]  = Q2[node_L] + Phi_L[1]*(frac*(Q2[node_R] - Q2[node_L]) + (1.0 - frac)*(Q2x[node_L]*Rx_L + Q2y[node_L]*Ry_L + Q2z[node_L]*Rz_L));
    Q_L[2]  = Q3[node_L] + Phi_L[2]*(frac*(Q3[node_R] - Q3[node_L]) + (1.0 - frac)*(Q3x[node_L]*Rx_L + Q3y[node_L]*Ry_L + Q3z[node_L]*Rz_L));
    Q_L[3]  = Q4[node_L] + Phi_L[3]*(frac*(Q4[node_R] - Q4[node_L]) + (1.0 - frac)*(Q4x[node_L]*Rx_L + Q4y[node_L]*Ry_L + Q4z[node_L]*Rz_L));
    Q_L[4]  = Q5[node_L] + Phi_L[4]*(frac*(Q5[node_R] - Q5[node_L]) + (1.0 - frac)*(Q5x[node_L]*Rx_L + Q5y[node_L]*Ry_L + Q5z[node_L]*Rz_L));

    // Compute Vector R for this edge node_R-->nodeid
    Rx_R = - Rx_L;
    Ry_R = - Ry_L;
    Rz_R = - Rz_L;

    // Compute second order Right Q's for this edge node_R-->node_L without limiter
    Q_R[0]  = Q1[node_R] + Phi_R[0]*(frac*(Q1[node_L] - Q1[node_R]) + (1.0 - frac)*(Q1x[node_R]*Rx_R + Q1y[node_R]*Ry_R + Q1z[node_R]*Rz_R));
    Q_R[1]  = Q2[node_R] + Phi_R[1]*(frac*(Q2[node_L] - Q2[node_R]) + (1.0 - frac)*(Q2x[node_R]*Rx_R + Q2y[node_R]*Ry_R + Q2z[node_R]*Rz_R));
    Q_R[2]  = Q3[node_R] + Phi_R[2]*(frac*(Q3[node_L] - Q3[node_R]) + (1.0 - frac)*(Q3x[node_R]*Rx_R + Q3y[node_R]*Ry_R + Q3z[node_R]*Rz_R));
    Q_R[3]  = Q4[node_R] + Phi_R[3]*(frac*(Q4[node_L] - Q4[node_R]) + (1.0 - frac)*(Q4x[node_R]*Rx_R + Q4y[node_R]*Ry_R + Q4z[node_R]*Rz_R));
    Q_R[4]  = Q5[node_R] + Phi_R[4]*(frac*(Q5[node_L] - Q5[node_R]) + (1.0 - frac)*(Q5x[node_R]*Rx_R + Q5y[node_R]*Ry_R + Q5z[node_R]*Rz_R));

    // Apply Pressure Correction on Limiter
    // Compute Roe Variables
    Compute_RoeVariables(Q_L, Q_R, Q_Roe);
    P_Roe = (Gamma - 1.0)*(Q_Roe[4] - (0.5/Q_Roe[0])*(Q_Roe[1]*Q_Roe[1] + Q_Roe[2]*Q_Roe[2] + Q_Roe[3]*Q_Roe[3]));
    if (P_Roe < 0.0) {
        for (i = 0; i < 5; i++) {
            Phi_L[i] = 0.0;
            Phi_R[i] = 0.0;
        }
    }
}

