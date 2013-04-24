/*******************************************************************************
 * File:        Test.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifdef DEBUG
#include <assert.h>
#endif

// Custom header files
#include "Trim_Utils.h"
#include "RestartIO.h"
#include "MeshIO.h"
#include "Commons.h"
#include "Solver.h"
#include "MC.h"
#include "Gradient.h"
#include "CompressibleUtils.h"
#include "Material.h"
#include "DebugSolver.h"
#include "Residual_Smoothing.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Unsteady_Initialization(void) {
    
    int iNode;
    double x0, y0, z0;
    double r, r2, beta, fac1, fac2, fac3;
    double rho0, u0, v0, w0, p0, T0, s0, c0;
    double T, u, v, w, p, rho, et;
    double daQ[NEQUATIONS];
    int dvVariableType_Old;
    
    // This domain is 0 < x < xmax; -5 < y < 5; translated in z
    // An isentropic perturbation is being added
    // This is based on a modification of Yee's test case
    // and comes to us from AEDC (Bobby Nichols) via Roy at UAB
    // This is good for the compressible and variable Mach regimes
    x0   = 5.0; // x location of vortex at time t = 0
    y0   = 0.0; // y location of vortex at time t = 0
    z0   = 0.0;  // z location of vortex at time t = 0
    beta = 5.0; // Vortex strength
    rho0 = Inf_Rho;
    p0   = Inf_Pressure + Gauge_Pressure;
    T0   = Inf_Temperature;
    c0   = Inf_SpeedSound;
    u0   = Inf_U;
    v0   = Inf_V;
    w0   = Inf_W;
    s0   = p0/pow(rho0, Gamma);
    
    // Initialize the Physical Nodes
    for (iNode = 0; iNode < nNode; iNode++) {
        r2 = (coordXYZ[3*iNode] - x0)*(coordXYZ[3*iNode] - x0) + (coordXYZ[3*iNode + 1] - y0)*(coordXYZ[3*iNode + 1] - y0);
        r = sqrt(r2);
        
        fac1 = beta / (2.0 * M_PI);
        fac2 = -0.5 * (Gamma - 1.0) / Gamma * fac1*fac1;
        fac3 = exp((1.0 - r2) / 2.0);
        
        T = T0 + fac2*fac3*fac3;
        u = u0 - (fac1 * fac3)*(coordXYZ[3*iNode + 1] - y0);
        v = v0 + (fac1 * fac3)*(coordXYZ[3*iNode] - x0);
        w = w0;
        p = pow(pow(T/(Gamma), Gamma)/s0, 1.0/(Gamma - 1.0));
        
        // Compute Density and Total Energy
        dvVariableType_Old = VariableType;
        VariableType = VARIABLE_PUT;
        daQ[0] = p;
        daQ[1] = u;
        daQ[2] = v;
        daQ[3] = w;
        daQ[4] = T;
        rho = Material_Get_Density(daQ);
        et  = Material_Get_TotalEnergy(daQ);
        VariableType = dvVariableType_Old;

        Q1[iNode] = rho;
        Q2[iNode] = rho*u;
        Q3[iNode] = rho*v;
        Q4[iNode] = rho*w;
        Q5[iNode] = rho*et;
    }
    
    // Initialize the Ghost Nodes
    for (int iBEdge = 0; iBEdge < nBEdge; iBEdge++) {
        Q1[bndEdge[iBEdge].node[1]] = Q1[bndEdge[iBEdge].node[0]];
        Q2[bndEdge[iBEdge].node[1]] = Q2[bndEdge[iBEdge].node[0]];
        Q3[bndEdge[iBEdge].node[1]] = Q3[bndEdge[iBEdge].node[0]];
        Q4[bndEdge[iBEdge].node[1]] = Q4[bndEdge[iBEdge].node[0]];
        Q5[bndEdge[iBEdge].node[1]] = Q5[bndEdge[iBEdge].node[0]];
    }
}

