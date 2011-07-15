/*******************************************************************************
 * File:        Time_Step.cpp
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
void Compute_DeltaT(int Iteration) {
    int i, nodeL, nodeR;
    double denom = 0.0;
    
    Vector3D areavec;
    double   area;

    double Ubar_avg, p_avg, c_avg,Q_avg[5];
    double lamda4, lamda5;
    double max_lamda;
    

    // Internal Edges
    for (i = 0; i < nEdge; i++) {
        // Get two nodes of edge
        nodeL = intEdge[i].node[0];
        nodeR = intEdge[i].node[1];

#ifdef DEBUG
        assert(nodeR > nodeL);
#endif
        
        // Get area vector
        areavec = intEdge[i].areav;
        area  = areavec.magnitude();
        areavec.normalize();

        // Get the average
        Q_avg[0] = 0.5*(Q1[nodeL] + Q1[nodeR]);
        Q_avg[1] = 0.5*(Q2[nodeL] + Q2[nodeR]);
        Q_avg[2] = 0.5*(Q3[nodeL] + Q3[nodeR]);
        Q_avg[3] = 0.5*(Q4[nodeL] + Q4[nodeR]);
        Q_avg[4] = 0.5*(Q5[nodeL] + Q5[nodeR]);
        Ubar_avg = (Q_avg[1]*areavec.vec[0] + Q_avg[2]*areavec.vec[1] + Q_avg[3]*areavec.vec[2])/Q_avg[0];
        p_avg    = (Gamma - 1.0)*(Q_avg[4] - 0.5*((Q_avg[1]*Q_avg[1] + Q_avg[2]*Q_avg[2] + Q_avg[3]*Q_avg[3])/Q_avg[0]));
        c_avg    = sqrt(Gamma*p_avg/Q_avg[0]);

        // Find maximum eigenvalue
        lamda4 = Ubar_avg + c_avg;
        lamda5 = Ubar_avg - c_avg;
        max_lamda = MAX(fabs(lamda4), fabs(lamda5));

        denom = (max_lamda * area);

        // Add to Each Node
        DeltaT[nodeL] += denom;
        DeltaT[nodeR] += denom;
    }

    // Boundary Edges
    for (i = 0; i < nBEdge; i++) {
        // Get two nodes of edge
        nodeL = bndEdge[i].node[0];
        nodeR = bndEdge[i].node[1];

#ifdef DEBUG
        assert(nodeR > nodeL);
#endif
        
        // Get area vector
        areavec = bndEdge[i].areav;
        area    = areavec.magnitude();
        areavec.normalize();

        // Note we only need to accumulate for the left node
        // Get the average
        Q_avg[0] = 0.5*(Q1[nodeL] + Q1[nodeR]);
        Q_avg[1] = 0.5*(Q2[nodeL] + Q2[nodeR]);
        Q_avg[2] = 0.5*(Q3[nodeL] + Q3[nodeR]);
        Q_avg[3] = 0.5*(Q4[nodeL] + Q4[nodeR]);
        Q_avg[4] = 0.5*(Q5[nodeL] + Q5[nodeR]);
        Ubar_avg = (Q_avg[1]*areavec.vec[0] + Q_avg[2]*areavec.vec[1] + Q_avg[3]*areavec.vec[2])/Q_avg[0];
        p_avg    = (Gamma - 1.0)*(Q_avg[4] - 0.5*((Q_avg[1]*Q_avg[1] + Q_avg[2]*Q_avg[2] + Q_avg[3]*Q_avg[3])/Q_avg[0]));
        c_avg    = sqrt(Gamma*p_avg/Q_avg[0]);

        // Find maximum eigenvalue
        lamda4 = Ubar_avg + c_avg;
        lamda5 = Ubar_avg - c_avg;
        max_lamda = MAX(fabs(lamda4), fabs(lamda5));

        denom = (max_lamda * area);
        
        // Add to Left Node Only
        DeltaT[nodeL] += denom;
    }
    
    // Finally Compute the Local Time Stepping
    if ((CFL_Ramp > 1) && (CFL_MAX > CFL_MIN)) {
        if (Iteration < CFL_Ramp)
            CFL = CFL_MIN + (CFL_MAX - CFL_MIN)*(((double)Iteration)/((double)(CFL_Ramp-1)));
        else
            CFL = CFL_MAX;
    } else
        CFL = CFL_MAX;

    for (i = 0; i < nNode; i++)
        DeltaT[i] = cVolume[i] * CFL / DeltaT[i];

    // Check if Global Time Stepping is Required
    if (TimeAccuracy == 1) {
        denom = DeltaT[0];
        for (i = 0; i < nNode; i++)
            denom = MIN(denom, DeltaT[i]);
        for (i = 0; i < nNode; i++)
            DeltaT[i] = denom;
    }
}

