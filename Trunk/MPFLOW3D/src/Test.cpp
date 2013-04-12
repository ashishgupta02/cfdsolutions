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

//------------------------------------------------------------------------------
//! Computes Implicit Linearized Residual System of Equation (SOE): Roe Approximate
//! Note: Jacobian are computed using first order Q's
//------------------------------------------------------------------------------
void Compute_Implicit_Linearized_Residual_SOE_Roe_Approximate(int AddTime, int Iteration) {
    int i, j, k, iNode, iEdge, ibEdge;
    int node_L, node_R;
    int idgn, idgnL, idgnR, ofdgnL, ofdgnR;
    double area, tmp;
    double **dFdL;
    double **dFdR;
    double **ARoe;
    Vector3D areavec;
    
    // Create the Helper Matrix
    dFdL = (double **) malloc(NEQUATIONS*sizeof(double*));
    dFdR = (double **) malloc(NEQUATIONS*sizeof(double*));
    ARoe = (double **) malloc(NEQUATIONS*sizeof(double*));
    for (i = 0; i < NEQUATIONS; i++) {
        dFdL[i] = (double *) malloc(NEQUATIONS*sizeof(double));
        dFdR[i] = (double *) malloc(NEQUATIONS*sizeof(double));
        ARoe[i] = (double *) malloc(NEQUATIONS*sizeof(double));
    }
    
    // Initialize the CRS Matrix
    for (i = 0; i < SolverBlockMatrix.DIM; i++) {
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                SolverBlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Internal Edges
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
        // Get two nodes of edge
        node_L = intEdge[iEdge].node[0];
        node_R = intEdge[iEdge].node[1];
        
        // Get area vector
        areavec = intEdge[iEdge].areav;
        area    = areavec.magnitude();
        
        // Get the diagonal Locations
        idgnL = SolverBlockMatrix.IAU[node_L];
        idgnR = SolverBlockMatrix.IAU[node_R];
        
        // Get the Off-Diagonal Locations
        // node_L: ofdgnL-> node_R;
        // node_R: ofdgnR-> node_L;
        ofdgnL = -1;
        for (i = SolverBlockMatrix.IA[node_L]; i < SolverBlockMatrix.IA[node_L+1]; i++) {
            if (SolverBlockMatrix.JA[i] == node_R) {
                ofdgnL = i;
                break;
            }
        }
        ofdgnR = -1;
        for (i = SolverBlockMatrix.IA[node_R]; i < SolverBlockMatrix.IA[node_R+1]; i++) {
            if (SolverBlockMatrix.JA[i] == node_L) {
                ofdgnR = i;
                break;
            }
        }
        
        // Initialize the Helper Matrix
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                dFdL[i][j] = 0.0;
                dFdR[i][j] = 0.0;
                ARoe[i][j] = 0.0;
            }
        }
        
        // Compute dFdL
        Compute_Flux_Jacobian_Euler_Convective(node_L, areavec, dFdL);
        
        // Compute dFdR
        Compute_Flux_Jacobian_Euler_Convective(node_R, areavec, dFdR);
        
        // Compute Dissipation Matrix Roe
        Compute_Dissipation_Matrix_Roe(node_L, node_R, areavec, ARoe);
        
        // Finally Compute the dFlux/dQ_L and dFlux/dQ_R
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                dFdL[i][j] = 0.5*(dFdL[i][j] + ARoe[i][j])*area;
                dFdR[i][j] = 0.5*(dFdR[i][j] - ARoe[i][j])*area;
            }
        }
        
        // Update the Diagonal and Off-Diagonal Terms
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                // Diagonal
                SolverBlockMatrix.A[idgnL][i][j] += dFdL[i][j];
                SolverBlockMatrix.A[idgnR][i][j] -= dFdR[i][j];
                // Off-Diagonal
                SolverBlockMatrix.A[ofdgnL][i][j] += dFdR[i][j];
                SolverBlockMatrix.A[ofdgnR][i][j] -= dFdL[i][j];
            }
        }
    }
    
    // Boundary Edges
    for (ibEdge = 0; ibEdge < nBEdge; ibEdge++) {
        // Get two nodes of edge
        node_L = bndEdge[ibEdge].node[0];
        node_R = bndEdge[ibEdge].node[1];
        
        // Get area vector
        areavec = bndEdge[ibEdge].areav;
        area    = areavec.magnitude();
        
        // Get the diagonal Locations - Only Physical
        // No diagonal and off-diagonal locations exists for Ghost Nodes
        idgnL = SolverBlockMatrix.IAU[node_L];
        
        // Initialize the Helper Matrix - Only Physical
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                dFdL[i][j] = 0.0;
                ARoe[i][j] = 0.0;
            }
        }
        
        // Compute dFdL
        Compute_Flux_Jacobian_Euler_Convective(node_L, areavec, dFdL);
        
        // Compute Dissipation Matrix Roe
        Compute_Dissipation_Matrix_Roe(node_L, node_R, areavec, ARoe);
        
        
        // Finally Compute the dFlux/dQ_L
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                dFdL[i][j] = 0.5*(dFdL[i][j] + ARoe[i][j])*area;
        }
        
        // Update the Diagonal Term of Physical Node Only
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++)
                SolverBlockMatrix.A[idgnL][i][j] += dFdL[i][j];
        }
    }
    
    // Copy the Residuals to Block Matrix which is B
    // And Copy I/DeltaT
    for (iNode = 0; iNode < nNode; iNode++) {
        // Get the LHS
        SolverBlockMatrix.B[iNode][0] = -Res1[iNode];
        SolverBlockMatrix.B[iNode][1] = -Res2[iNode];
        SolverBlockMatrix.B[iNode][2] = -Res3[iNode];
        SolverBlockMatrix.B[iNode][3] = -Res4[iNode];
        SolverBlockMatrix.B[iNode][4] = -Res5[iNode];
        
        if (AddTime == 1) {
            // Initialize the Helper Matrix
            for (i = 0; i < NEQUATIONS; i++)
                for (j = 0; j < NEQUATIONS; j++)
                    ARoe[i][j] = 0.0;

            // Compute the Diagonal Term: Cinv
            Compute_Roe_Transformed_Precondition_Matrix(iNode, 1, ARoe);

            // Scale the Diagonal Term
            tmp = cVolume[iNode]/DeltaT[iNode];
            for (i = 0; i < NEQUATIONS; i++)
                for (j = 0; j < NEQUATIONS; j++)
                    ARoe[i][j] *= tmp;
            
            // Get the diagonal location
            idgn = SolverBlockMatrix.IAU[iNode];
            
            // Add the Diagonal Term
            for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
                for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                    if (k == j)
                        SolverBlockMatrix.A[idgn][j][k] += ARoe[j][k];
            }
        }
    }
    
    // Delete the Helper Matrix
    for (i = 0; i < NEQUATIONS; i++) {
        free(ARoe[i]);
        free(dFdL[i]);
        free(dFdR[i]);
    }
    free(ARoe);
    free(dFdL);
    free(dFdR);
}

//------------------------------------------------------------------------------
//! Solver in Steady State Mode with Explicit Mode and RK4 Time Integration
//------------------------------------------------------------------------------
int Solver_Steady_Explicit_RungeKutta4(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    int SaveRMS     = 0;
    double *WQ;
    double *WDeltaT, *Dummy;
    double phi[4];
    double eta;
    double scale_RMS[5], scale_RMS_Res;
    
    // Check if Euler or Runge-Kutta Scheme
    WQ      = NULL;
    WDeltaT = NULL;
    Dummy   = NULL;
    
    // Create the helper arrays for Runge-Kutta Method
    // WQ
    WQ = new double[NEQUATIONS*nNode];
    
    // WDeltaT : Needed because DeltaT computations is inside Residual Computation
    WDeltaT = new double[nNode];
    
    // Helper Variables
    eta  = 0.0;
    
    // Runge-Kutta 4 Coefficients
    phi[0] = 0.25;
    phi[1] = 0.333333333333333333;
    phi[2] = 0.5;
    phi[3] = 1.0;
    
    // Save Solver Parameters
    SaveOrder   = SolverOrder;
    SaveLimiter = LimiterMethod;
    
    for (int iter = RestartIteration; iter < SolverNIteration; iter++) {
        SolverIteration = iter;
        
        // Reset Residuals and DeltaT
        for (int i = 0; i < nNode; i++) {
            Res1[i]      = 0.0;
            Res2[i]      = 0.0;
            Res3[i]      = 0.0;
            Res4[i]      = 0.0;
            Res5[i]      = 0.0;
            Res1_Diss[i] = 0.0;
            Res2_Diss[i] = 0.0;
            Res3_Diss[i] = 0.0;
            Res4_Diss[i] = 0.0;
            Res5_Diss[i] = 0.0;
            DeltaT[i]    = 0.0;
        }
        
        // Min and Max EigenValue Value
        MinEigenLamda1   = DBL_MAX;
        MaxEigenLamda1   = DBL_MIN;
        MinEigenLamda4   = DBL_MAX;
        MaxEigenLamda4   = DBL_MIN;
        MinEigenLamda5   = DBL_MAX;
        MaxEigenLamda5   = DBL_MIN;
        
        // Min and Max Time
        MinDeltaT = DBL_MAX;
        MaxDeltaT = DBL_MIN;
        
        // Min and Max Precondition Variable
        if (PrecondMethod != PRECOND_METHOD_NONE) {
            for (int i = 0; i < nNode; i++)
                PrecondSigma[i] = 0.0;
            MinPrecondSigma = DBL_MAX;
            MaxPrecondSigma = DBL_MIN;
        }
        
        // Compute Far Field Conditions with Mach Ramping
        Material_Set_InfinityCondition(iter);
        
        // Check if First Order Iterations are Required
        SolverOrder = SaveOrder;
        if (FirstOrderNIteration > iter)
            SolverOrder = SOLVER_ORDER_FIRST;

        // Set the Start and End of Limiter Iterations
        LimiterMethod = LIMITER_METHOD_NONE;
        if (LimiterStartSolverIteration < LimiterEndSolverIteration) {
            if ((LimiterStartSolverIteration <= iter+1) && (LimiterEndSolverIteration > iter))
                LimiterMethod = SaveLimiter;
        }

        // Compute Least Square Gradient -- Unweighted
        if (SolverOrder == SOLVER_ORDER_SECOND) {
            Compute_Least_Square_Gradient(0); // Need to be fixed: replace zero with enum type
            if (LimiterMethod != LIMITER_METHOD_NONE)
                Compute_Limiter();
        }
        
        // Apply boundary conditions
        Apply_Boundary_Condition(iter);

        // Reset the Residual Smoothing Variables: Required before Residual Computation
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
            Residual_Smoothing_Reset();
        
        // Compute Residuals
        AddTime = TRUE;
        Compute_Residual(AddTime);
        
        // Compute Local Time Stepping
        Compute_DeltaT(iter);

        // Smooth the Residual: DeltaT Computation is required
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
            Residual_Smoothing();
        
        // Compute RMS
        RMS[0] = RMS[1] = RMS[2] = RMS[3] = RMS[4] = 0.0;
        for (int i = 0; i < nNode; i++) {
            RMS[0] += (Res1[i]+ Res1_Diss[i])*(Res1[i]+ Res1_Diss[i]);
            RMS[1] += (Res2[i]+ Res2_Diss[i])*(Res2[i]+ Res2_Diss[i]);
            RMS[2] += (Res3[i]+ Res3_Diss[i])*(Res3[i]+ Res3_Diss[i]);
            RMS[3] += (Res4[i]+ Res4_Diss[i])*(Res4[i]+ Res4_Diss[i]);
            RMS[4] += (Res5[i]+ Res5_Diss[i])*(Res5[i]+ Res5_Diss[i]);
        }
        RMS_Res = RMS[0] + RMS[1] + RMS[2] + RMS[3] + RMS[4];
        RMS_Res = sqrt(RMS_Res/(5.0 * (double)nNode));
        RMS[0] = sqrt(RMS[0]/(double)nNode);
        RMS[1] = sqrt(RMS[1]/(double)nNode);
        RMS[2] = sqrt(RMS[2]/(double)nNode);
        RMS[3] = sqrt(RMS[3]/(double)nNode);
        RMS[4] = sqrt(RMS[4]/(double)nNode);
        
        // Normalize the RMS
        if (SaveRMS == 0) {
            scale_RMS_Res = 0.0;
            for (int i = 0; i < NEQUATIONS; i++) {
                scale_RMS[i] = RMS[i];
                scale_RMS_Res += RMS[i];
            }
            scale_RMS_Res /= NEQUATIONS;
            SaveRMS = 1;
        }
        for (int i = 0; i < NEQUATIONS; i++)
            RMS[i] /= scale_RMS[i];
        RMS_Res /= scale_RMS_Res;
        
        // Write RMS
        RMS_Writer(iter+1, RMS);
        
        // Check for Residual NAN
        if (isnan(RMS_Res)) {
            info("Solve_Explicit: NAN Encountered ! - Abort");
            iter = SolverNIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }
        
        // Runge-Kutta 4 Time Integration Method
        for (int i = 0; i < nNode; i++) {
            // W0 = Wn
            WQ[NEQUATIONS*i + 0] = Q1[i];
            WQ[NEQUATIONS*i + 1] = Q2[i];
            WQ[NEQUATIONS*i + 2] = Q3[i];
            WQ[NEQUATIONS*i + 3] = Q4[i];
            WQ[NEQUATIONS*i + 4] = Q5[i];
            WDeltaT[i] = DeltaT[i];

            // W1 = W0 - phi*dT*RES0/vol
            eta   = WDeltaT[i]/cVolume[i];
            Q1[i] = WQ[NEQUATIONS*i + 0] - phi[0]*eta*Res1[i];
            Q2[i] = WQ[NEQUATIONS*i + 1] - phi[0]*eta*Res2[i];
            Q3[i] = WQ[NEQUATIONS*i + 2] - phi[0]*eta*Res3[i];
            Q4[i] = WQ[NEQUATIONS*i + 3] - phi[0]*eta*Res4[i];
            Q5[i] = WQ[NEQUATIONS*i + 4] - phi[0]*eta*Res5[i];
            
            // Reset the Residual
            Res1[i]      = 0.0;
            Res2[i]      = 0.0;
            Res3[i]      = 0.0;
            Res4[i]      = 0.0;
            Res5[i]      = 0.0;
            Res1_Diss[i] = 0.0;
            Res2_Diss[i] = 0.0;
            Res3_Diss[i] = 0.0;
            Res4_Diss[i] = 0.0;
            Res5_Diss[i] = 0.0;
            DeltaT[i]    = 0.0;
        }

        // Compute RES1
        // Compute Least Square Gradient -- Unweighted
        if (SolverOrder == SOLVER_ORDER_SECOND) {
            Compute_Least_Square_Gradient(0);
            if (LimiterMethod != LIMITER_METHOD_NONE)
                Compute_Limiter();
        }

        // Apply boundary conditions
        Apply_Boundary_Condition(iter);

        // Reset the Residual Smoothing Variables: Required before Residual Computation
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
            Residual_Smoothing_Reset();

        // Compute Residuals
        AddTime = FALSE;
        Compute_Residual(AddTime);

        // Smooth the Residual: DeltaT Computation is required
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE) {
            // Set the original Time
            Dummy  = DeltaT;
            DeltaT = WDeltaT;
            Residual_Smoothing();
            // Swith back
            DeltaT = Dummy;
            Dummy  = NULL;
        }

        for (int i = 0; i < nNode; i++) {
            // W2 = W0 - phi2*dT*RES1/vol
            eta   = WDeltaT[i]/cVolume[i];
            Q1[i] = WQ[NEQUATIONS*i + 0] - phi[1]*eta*Res1[i];
            Q2[i] = WQ[NEQUATIONS*i + 1] - phi[1]*eta*Res2[i];
            Q3[i] = WQ[NEQUATIONS*i + 2] - phi[1]*eta*Res3[i];
            Q4[i] = WQ[NEQUATIONS*i + 3] - phi[1]*eta*Res4[i];
            Q5[i] = WQ[NEQUATIONS*i + 4] - phi[1]*eta*Res5[i];

            // Reset the Residual
            Res1[i]      = 0.0;
            Res2[i]      = 0.0;
            Res3[i]      = 0.0;
            Res4[i]      = 0.0;
            Res5[i]      = 0.0;
            Res1_Diss[i] = 0.0;
            Res2_Diss[i] = 0.0;
            Res3_Diss[i] = 0.0;
            Res4_Diss[i] = 0.0;
            Res5_Diss[i] = 0.0;
            DeltaT[i]    = 0.0;
        }

        // Compute RES2
        // Compute Least Square Gradient -- Unweighted
        if (SolverOrder == SOLVER_ORDER_SECOND) {
            Compute_Least_Square_Gradient(0);
            if (LimiterMethod != LIMITER_METHOD_NONE)
                Compute_Limiter();
        }

        // Apply boundary conditions
        Apply_Boundary_Condition(iter);

        // Reset the Residual Smoothing Variables: Required before Residual Computation
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
            Residual_Smoothing_Reset();

        // Compute Residuals
        AddTime = FALSE;
        Compute_Residual(AddTime);

        // Smooth the Residual: DeltaT Computation is required
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE) {
            // Set the original Time
            Dummy  = DeltaT;
            DeltaT = WDeltaT;
            Residual_Smoothing();
            // Swith back
            DeltaT = Dummy;
            Dummy  = NULL;
        }

        for (int i = 0; i < nNode; i++) {
            // W3 = W0 - phi3*dT*RES2/vol
            eta   = WDeltaT[i]/cVolume[i];
            Q1[i] = WQ[NEQUATIONS*i + 0] - phi[2]*eta*Res1[i];
            Q2[i] = WQ[NEQUATIONS*i + 1] - phi[2]*eta*Res2[i];
            Q3[i] = WQ[NEQUATIONS*i + 2] - phi[2]*eta*Res3[i];
            Q4[i] = WQ[NEQUATIONS*i + 3] - phi[2]*eta*Res4[i];
            Q5[i] = WQ[NEQUATIONS*i + 4] - phi[2]*eta*Res5[i];

            // Reset the Residual
            Res1[i]      = 0.0;
            Res2[i]      = 0.0;
            Res3[i]      = 0.0;
            Res4[i]      = 0.0;
            Res5[i]      = 0.0;
            Res1_Diss[i] = 0.0;
            Res2_Diss[i] = 0.0;
            Res3_Diss[i] = 0.0;
            Res4_Diss[i] = 0.0;
            Res5_Diss[i] = 0.0;
            DeltaT[i]    = 0.0;
        }

        // Compute RES3
        // Compute Least Square Gradient -- Unweighted
        if (SolverOrder == SOLVER_ORDER_SECOND) {
            Compute_Least_Square_Gradient(0);
            if (LimiterMethod != LIMITER_METHOD_NONE)
                Compute_Limiter();
        }

        // Apply boundary conditions
        Apply_Boundary_Condition(iter);

        // Reset the Residual Smoothing Variables: Required before Residual Computation
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
            Residual_Smoothing_Reset();

        // Compute Residuals
        AddTime = FALSE;
        Compute_Residual(AddTime);

        // Smooth the Residual: DeltaT Computation is required
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE) {
            // Set the original Time
            Dummy  = DeltaT;
            DeltaT = WDeltaT;
            Residual_Smoothing();
            // Swith back
            DeltaT = Dummy;
            Dummy  = NULL;
        }

        for (int i = 0; i < nNode; i++) {
            // W4 = W0 - phi4*dT*RES3/vol
            eta   = WDeltaT[i]/cVolume[i];
            Q1[i] = WQ[NEQUATIONS*i + 0] - phi[3]*eta*Res1[i];
            Q2[i] = WQ[NEQUATIONS*i + 1] - phi[3]*eta*Res2[i];
            Q3[i] = WQ[NEQUATIONS*i + 2] - phi[3]*eta*Res3[i];
            Q4[i] = WQ[NEQUATIONS*i + 3] - phi[3]*eta*Res4[i];
            Q5[i] = WQ[NEQUATIONS*i + 4] - phi[3]*eta*Res5[i];
        }

        // Precondition Variable Normalize
        if (PrecondMethod != PRECOND_METHOD_NONE) {
            for (int i = 0; i < nNode; i++)
                PrecondSigma[i] /= (crs_IA_Node2Node[i+1] - crs_IA_Node2Node[i]);
        }
        
        // Check Cyclic Restart is Requested
        RestartIteration = iter+1;
        Check_Restart(iter);
        
        if ((RMS_Res < (DBL_EPSILON*10.0))|| ((iter+1) == SolverNIteration)) {
            iter = SolverNIteration + 1;
        }
    }
    
    // Check if Solution Restart is Requested
    if (RestartOutput && CheckNAN != 1)
        Restart_Writer(RestartOutputFilename, 1);

    // Debug the NAN
    if (CheckNAN)
        DebugNAN();

    // Free Memory
    // Wo
    if (WQ != NULL)
        delete[] WQ;
    
    // WDeltaT
    if (WDeltaT != NULL)
        delete[] WDeltaT;

    WQ = NULL;
    WDeltaT = NULL;
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

