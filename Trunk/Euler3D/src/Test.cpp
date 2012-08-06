/*******************************************************************************
 * File:        Test.cpp
 * Author:      Ashish Gupta
 * Revision:    4
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

//------------------------------------------------------------------------------
//! Compute Transformed Residual in Conservative to Primitive Form
//------------------------------------------------------------------------------
void Transform_Residual_ConservativeToPrimitive(void) {
    int i;
    double Q[5];
    double flux_roe[5];
    double flux_roe_diss[5];
    double rho, u, v, w, p, T, c, ht, et, q2, mach;
    double dtmp;
    double **Roe_P;
    
    // Transform the Residual in Conservative Form to Primitive Form
    for (i = 0; i < nNode; i++) {
        // Get the Variables
        Q[0] = Q1[i];
        Q[1] = Q2[i];
        Q[2] = Q3[i];
        Q[3] = Q4[i];
        Q[4] = Q5[i];

        // Compute Equation of State
        Compute_EOS_Variables_ControlVolume(Q, rho, p, T, u, v, w, q2, c, mach, et, ht);

        // Compute the Transformation Matrix Based on Variable Type and Non Dimensionalization Used
        switch (NonDimensional_Type) {
            case NONDIMENSIONAL_GENERIC:
                switch (Variable_Type) {
                    case VARIABLE_PRIMITIVE_PUT:
                        dtmp = 0.5*(Gamma - 1.0)*q2;
                        Roe_P[0][0] = dtmp;
                        Roe_P[0][1] = (1.0 - Gamma)*u;
                        Roe_P[0][2] = (1.0 - Gamma)*v;
                        Roe_P[0][3] = (1.0 - Gamma)*w;
                        Roe_P[0][4] = (Gamma - 1.0);

                        Roe_P[1][0] = -u/rho;
                        Roe_P[1][1] = 1.0/rho;
                        Roe_P[1][2] = 0.0;
                        Roe_P[1][3] = 0.0;
                        Roe_P[1][4] = 0.0;

                        Roe_P[2][0] = -v/rho;
                        Roe_P[2][1] = 0.0;
                        Roe_P[2][2] = 1.0/rho;
                        Roe_P[2][3] = 0.0;
                        Roe_P[2][4] = 0.0;

                        Roe_P[3][0] = -w/rho;
                        Roe_P[3][1] = 0.0;
                        Roe_P[3][2] = 0.0;
                        Roe_P[3][3] = 1.0/rho;
                        Roe_P[3][4] = 0.0;

                        Roe_P[4][0] = (dtmp*Gamma - T)/rho;
                        Roe_P[4][1] = (1.0 - Gamma)*Gamma*u/rho;
                        Roe_P[4][2] = (1.0 - Gamma)*Gamma*v/rho;
                        Roe_P[4][3] = (1.0 - Gamma)*Gamma*w/rho;
                        Roe_P[4][4] = (1.0 - Gamma)*Gamma/rho;
                        break;
                    case VARIABLE_PRIMITIVE_RUP:
                        break;
                    default:
                        error("Transform_Residual_ConservativeToPrimitive: Undefined Variable Type - %d - Error-1", Variable_Type);
                        break;
                }
                break;
            case NONDIMENSIONAL_BTW:
                    break;
            case NONDIMENSIONAL_LMROE:
                break;
            default:
                error("Transform_Residual_ConservativeToPrimitive: Undefined Non-Dimensional Type - %d", NonDimensional_Type);
                break;
        }
        // Get the Conservative Residual
        flux_roe[0] = Res1[i];
        flux_roe[1] = Res2[i];
        flux_roe[2] = Res3[i];
        flux_roe[3] = Res4[i];
        flux_roe[4] = Res5[i];
        MC_Matrix_Mul_Vector(5, 5, Roe_P, flux_roe, flux_roe_diss);

        // Update the Residual in Primitive Form
        Res1[i] = flux_roe_diss[0];
        Res2[i] = flux_roe_diss[1];
        Res3[i] = flux_roe_diss[2];
        Res4[i] = flux_roe_diss[3];
        Res5[i] = flux_roe_diss[4]; 
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Unsteady_Initialization(void) {
    
    int iNode;
    double x0, y0, z0;
    double r, r2, beta, fac1, fac2, fac3;
    double rho0, u0, v0, w0, p0, T0, s0, c0;
    double T, u, v, w, p, rho;
    
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
    p0   = Inf_Pressure;
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
        rho = Get_Rho(p, T);

        Q1[iNode] = rho;
        Q2[iNode] = rho*u;
        Q3[iNode] = rho*v;
        Q4[iNode] = rho*w;
        Q5[iNode] = rho*Get_TotalEnergy(rho, p, u, v, w);
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

