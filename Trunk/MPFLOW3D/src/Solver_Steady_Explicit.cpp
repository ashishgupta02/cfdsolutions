/*******************************************************************************
 * File:        Solver_Steady_Explicit.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "RestartIO.h"
#include "MeshIO.h"
#include "Commons.h"
#include "Solver.h"
#include "MC.h"
#include "Gradient.h"
#include "CompressibleUtils.h"
#include "Residual_Smoothing.h"
#include "Material.h"
#include "DebugSolver.h"

//------------------------------------------------------------------------------
//! Solver in Steady State Mode with Explicit Mode and RungeKutta Time Integration
//------------------------------------------------------------------------------
int Solver_Steady_Explicit_RungeKutta(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    int SaveRMS     = 0;
    int RKStage     = 0;
    double eta      = 0.0;
    double *WQ      = NULL;
    double *WDeltaT = NULL;
    double *Dummy   = NULL;
    double phi[5];
    double scale_RMS[NEQUATIONS], scale_RMS_Res;
    
    // Initialize
    for (int i = 0; i < NEQUATIONS; i++)
        scale_RMS[i] = 1.0;
    scale_RMS_Res = 1.0;
    
    // Create the helper arrays for Runge-Kutta Method
    // WQ
    WQ = new double[NEQUATIONS*nNode];
    
    // WDeltaT : Needed because DeltaT computations is inside Residual Computation
    WDeltaT = new double[nNode];
    
    // Select the Explicit Time Integration Method Type
    switch (TimeIntegrationMethod) {
        case TIME_INTEGRATION_METHOD_EXPLICIT_EULER:
            RKStage = 1;
            phi[0]  = 1.0;
            phi[1]  = 0.0;
            phi[2]  = 0.0;
            phi[3]  = 0.0;
            phi[4]  = 0.0;
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK3:
            RKStage = 3;
            phi[0]  = 2.0/3.0;
            phi[1]  = 2.0/3.0;
            phi[2]  = 1.0;
            phi[3]  = 0.0;
            phi[4]  = 0.0;
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK4:
            RKStage = 4;
            phi[0]  = 1.0/4.0;
            phi[1]  = 1.0/3.0;
            phi[2]  = 1.0/2.0;
            phi[3]  = 1.0;
            phi[4]  = 0.0;
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK5:
            RKStage = 5;
            phi[0]  = 1.0/4.0;
            phi[1]  = 1.0/6.0;
            phi[2]  = 3.0/8.0;
            phi[3]  = 1.0/2.0;
            phi[4]  = 1.0;
            break;
        default:
            error("Solver_Steady_Explicit_RungeKutta: Invalid Time Integration Method - %d", TimeIntegrationMethod);
            break;
    }
    
    // Save Solver Parameters
    SaveOrder   = SolverOrder;
    SaveLimiter = LimiterMethod;
    
    for (int iter = RestartIteration; iter < SolverNIteration; iter++) {
        // Reset Residuals and DeltaT
        for (int i = 0; i < nNode; i++) {
            Res1_Conv[i] = 0.0;
            Res2_Conv[i] = 0.0;
            Res3_Conv[i] = 0.0;
            Res4_Conv[i] = 0.0;
            Res5_Conv[i] = 0.0;
            Res1_Diss[i] = 0.0;
            Res2_Diss[i] = 0.0;
            Res3_Diss[i] = 0.0;
            Res4_Diss[i] = 0.0;
            Res5_Diss[i] = 0.0;
            DeltaT[i]    = 0.0;
        }
        
        // Min and Max EigenValue Value
        MinEigenLamda1 = DBL_MAX;
        MaxEigenLamda1 = DBL_MIN;
        MinEigenLamda4 = DBL_MAX;
        MaxEigenLamda4 = DBL_MIN;
        MinEigenLamda5 = DBL_MAX;
        MaxEigenLamda5 = DBL_MIN;
        
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
            RMS[0] += (Res1_Conv[i] + Res1_Diss[i])*(Res1_Conv[i] + Res1_Diss[i]);
            RMS[1] += (Res2_Conv[i] + Res2_Diss[i])*(Res2_Conv[i] + Res2_Diss[i]);
            RMS[2] += (Res3_Conv[i] + Res3_Diss[i])*(Res3_Conv[i] + Res3_Diss[i]);
            RMS[3] += (Res4_Conv[i] + Res4_Diss[i])*(Res4_Conv[i] + Res4_Diss[i]);
            RMS[4] += (Res5_Conv[i] + Res5_Diss[i])*(Res5_Conv[i] + Res5_Diss[i]);
        }
        RMS_Res = RMS[0] + RMS[1] + RMS[2] + RMS[3] + RMS[4];
        RMS_Res = sqrt(RMS_Res/(5.0 * (double)nNode));
        for (int i = 0; i < NEQUATIONS; i++)
            RMS[i] = sqrt(RMS[i]/(double)nNode);
        
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
            info("Solver_Steady_Explicit_RungeKutta: NAN Encountered ! - Abort");
            iter = SolverNIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }
        
        // Start Runge-Kutta Stages Computation
        for (int iRKS = 0; iRKS < RKStage; iRKS++) {
            // Stage n
            for (int i = 0; i < nNode; i++) {
                
                // W(n) = W(n-1)
                WQ[NEQUATIONS*i + 0] = Q1[i];
                WQ[NEQUATIONS*i + 1] = Q2[i];
                WQ[NEQUATIONS*i + 2] = Q3[i];
                WQ[NEQUATIONS*i + 3] = Q4[i];
                WQ[NEQUATIONS*i + 4] = Q5[i];
                if (iRKS == 0)
                    WDeltaT[i] = DeltaT[i];

                // W(n+1) = W(n) - phi*dT*RES(n)/vol
                eta   = WDeltaT[i]/cVolume[i];
                Q1[i] = WQ[NEQUATIONS*i + 0] - phi[iRKS]*eta*(Res1_Conv[i] + Res1_Diss[i]);
                Q2[i] = WQ[NEQUATIONS*i + 1] - phi[iRKS]*eta*(Res2_Conv[i] + Res2_Diss[i]);
                Q3[i] = WQ[NEQUATIONS*i + 2] - phi[iRKS]*eta*(Res3_Conv[i] + Res3_Diss[i]);
                Q4[i] = WQ[NEQUATIONS*i + 3] - phi[iRKS]*eta*(Res4_Conv[i] + Res4_Diss[i]);
                Q5[i] = WQ[NEQUATIONS*i + 4] - phi[iRKS]*eta*(Res5_Conv[i] + Res5_Diss[i]);
                
                if (iRKS < (RKStage - 1)) {
                    // Reset the Residual
                    Res1_Conv[i] = 0.0;
                    Res2_Conv[i] = 0.0;
                    Res3_Conv[i] = 0.0;
                    Res4_Conv[i] = 0.0;
                    Res5_Conv[i] = 0.0;
                    Res1_Diss[i] = 0.0;
                    Res2_Diss[i] = 0.0;
                    Res3_Diss[i] = 0.0;
                    Res4_Diss[i] = 0.0;
                    Res5_Diss[i] = 0.0;
                    DeltaT[i]    = 0.0;
                }
            }
        
            if (iRKS < (RKStage - 1)) {
                // Compute RESn
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
            }
        }

        // Precondition Variable Normalize
        if (PrecondMethod != PRECOND_METHOD_NONE) {
            for (int i = 0; i < nNode; i++)
                PrecondSigma[i] /= RKStage*(crs_IA_Node2Node[i+1] - crs_IA_Node2Node[i]);
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

    WQ      = NULL;
    WDeltaT = NULL;
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Steady State Mode with Explicit Mode
//  and RungeKutta Martinelli Jameson Time Integration
//------------------------------------------------------------------------------
int Solver_Steady_Explicit_RungeKutta_Martinelli_Jameson(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    int SaveRMS     = 0;
    int RKStage     = 0;
    double eta      = 0.0;
    double *WQ      = NULL;
    double *Rd_Old  = NULL;
    double *WDeltaT = NULL;
    double *Dummy   = NULL;
    double phi[5], psi[5];
    double scale_RMS[NEQUATIONS], scale_RMS_Res;
    
    // Initialize
    for (int i = 0; i < NEQUATIONS; i++)
        scale_RMS[i] = 1.0;
    scale_RMS_Res = 1.0;
    
    // Create the helper arrays for Runge-Kutta Method
    // WQ and Residuals
    WQ     = new double[NEQUATIONS*nNode];
    Rd_Old = new double[NEQUATIONS*nNode];
    
    // WDeltaT : Needed because DeltaT computations is inside Residual Computation
    WDeltaT = new double[nNode];
    
    // Runge-Kutta 5 Coefficients
    RKStage = 5;
    phi[0]  = 1.0/4.0;
    phi[1]  = 1.0/6.0;
    phi[2]  = 3.0/8.0;
    phi[3]  = 1.0/2.0;
    phi[4]  = 1.0;
    
    psi[0]  = 1.0;
    psi[1]  = 0.0;
    psi[2]  = 14.0/25.0;
    psi[3]  = 0.0;
    psi[4]  = 11.0/25.0;
    
    // Save Solver Parameters
    SaveOrder   = SolverOrder;
    SaveLimiter = LimiterMethod;
    
    for (int iter = RestartIteration; iter < SolverNIteration; iter++) {
        // Reset Residuals and DeltaT
        for (int i = 0; i < nNode; i++) {
            Res1_Conv[i] = 0.0;
            Res2_Conv[i] = 0.0;
            Res3_Conv[i] = 0.0;
            Res4_Conv[i] = 0.0;
            Res5_Conv[i] = 0.0;
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
            RMS[0] += (Res1_Conv[i] + Res1_Diss[i])*(Res1_Conv[i] + Res1_Diss[i]);
            RMS[1] += (Res2_Conv[i] + Res2_Diss[i])*(Res2_Conv[i] + Res2_Diss[i]);
            RMS[2] += (Res3_Conv[i] + Res3_Diss[i])*(Res3_Conv[i] + Res3_Diss[i]);
            RMS[3] += (Res4_Conv[i] + Res4_Diss[i])*(Res4_Conv[i] + Res4_Diss[i]);
            RMS[4] += (Res5_Conv[i] + Res5_Diss[i])*(Res5_Conv[i] + Res5_Diss[i]);
        }
        RMS_Res = RMS[0] + RMS[1] + RMS[2] + RMS[3] + RMS[4];
        RMS_Res = sqrt(RMS_Res/(5.0 * (double)nNode));
        for (int i = 0; i < NEQUATIONS; i++)
            RMS[i] = sqrt(RMS[i]/(double)nNode);
        
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
            info("Solver_Steady_Explicit_RungeKutta_Martinelli_Jameson: NAN Encountered ! - Abort");
            iter = SolverNIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }
        
        // Save the Old Residuals
        for (int i = 0; i < nNode; i++) {
            Rd_Old[NEQUATIONS*i + 0] = Res1_Diss[i];
            Rd_Old[NEQUATIONS*i + 1] = Res2_Diss[i];
            Rd_Old[NEQUATIONS*i + 2] = Res3_Diss[i];
            Rd_Old[NEQUATIONS*i + 3] = Res4_Diss[i];
            Rd_Old[NEQUATIONS*i + 4] = Res5_Diss[i];
        }
        // Start Runge-Kutta Stages Computation
        for (int iRKS = 0; iRKS < RKStage; iRKS++) {
            // Stage n
            for (int i = 0; i < nNode; i++) {
                
                // W(n) = W(n-1)
                WQ[NEQUATIONS*i + 0] = Q1[i];
                WQ[NEQUATIONS*i + 1] = Q2[i];
                WQ[NEQUATIONS*i + 2] = Q3[i];
                WQ[NEQUATIONS*i + 3] = Q4[i];
                WQ[NEQUATIONS*i + 4] = Q5[i];
                if (iRKS == 0) {
                    WDeltaT[i] = DeltaT[i];
                    
                    // W(n+1) = W(n) - dT*(phi*Rc(n) + psi*Rd(n))/vol
                    eta   = WDeltaT[i]/cVolume[i];
                    Q1[i] = WQ[NEQUATIONS*i + 0] - phi[iRKS]*eta*(Res1_Conv[i] + psi[iRKS]*Rd_Old[NEQUATIONS*i + 0]);
                    Q2[i] = WQ[NEQUATIONS*i + 1] - phi[iRKS]*eta*(Res2_Conv[i] + psi[iRKS]*Rd_Old[NEQUATIONS*i + 1]);
                    Q3[i] = WQ[NEQUATIONS*i + 2] - phi[iRKS]*eta*(Res3_Conv[i] + psi[iRKS]*Rd_Old[NEQUATIONS*i + 2]);
                    Q4[i] = WQ[NEQUATIONS*i + 3] - phi[iRKS]*eta*(Res4_Conv[i] + psi[iRKS]*Rd_Old[NEQUATIONS*i + 3]);
                    Q5[i] = WQ[NEQUATIONS*i + 4] - phi[iRKS]*eta*(Res5_Conv[i] + psi[iRKS]*Rd_Old[NEQUATIONS*i + 4]);
                } else {
                    // W(n+1) = W(n) - dT*(phi*Rc(n) + psi*Rd(n) + (1-psi)*Rd(n-1))/vol
                    eta   = WDeltaT[i]/cVolume[i];
                    Q1[i] = WQ[NEQUATIONS*i + 0] - phi[iRKS]*eta*(Res1_Conv[i] + psi[iRKS]*Res1_Diss[i] + (1.0 - psi[iRKS])*Rd_Old[NEQUATIONS*i + 0]);
                    Q2[i] = WQ[NEQUATIONS*i + 1] - phi[iRKS]*eta*(Res2_Conv[i] + psi[iRKS]*Res2_Diss[i] + (1.0 - psi[iRKS])*Rd_Old[NEQUATIONS*i + 1]);
                    Q3[i] = WQ[NEQUATIONS*i + 2] - phi[iRKS]*eta*(Res3_Conv[i] + psi[iRKS]*Res3_Diss[i] + (1.0 - psi[iRKS])*Rd_Old[NEQUATIONS*i + 2]);
                    Q4[i] = WQ[NEQUATIONS*i + 3] - phi[iRKS]*eta*(Res4_Conv[i] + psi[iRKS]*Res4_Diss[i] + (1.0 - psi[iRKS])*Rd_Old[NEQUATIONS*i + 3]);
                    Q5[i] = WQ[NEQUATIONS*i + 4] - phi[iRKS]*eta*(Res5_Conv[i] + psi[iRKS]*Res5_Diss[i] + (1.0 - psi[iRKS])*Rd_Old[NEQUATIONS*i + 4]);
                    
                    Rd_Old[NEQUATIONS*i + 0] = Res1_Diss[i];
                    Rd_Old[NEQUATIONS*i + 1] = Res2_Diss[i];
                    Rd_Old[NEQUATIONS*i + 2] = Res3_Diss[i];
                    Rd_Old[NEQUATIONS*i + 3] = Res4_Diss[i];
                    Rd_Old[NEQUATIONS*i + 4] = Res5_Diss[i];
                    
                }
                if (iRKS < (RKStage - 1)) {
                    // Reset the Residual
                    Res1_Conv[i] = 0.0;
                    Res2_Conv[i] = 0.0;
                    Res3_Conv[i] = 0.0;
                    Res4_Conv[i] = 0.0;
                    Res5_Conv[i] = 0.0;
                    Res1_Diss[i] = 0.0;
                    Res2_Diss[i] = 0.0;
                    Res3_Diss[i] = 0.0;
                    Res4_Diss[i] = 0.0;
                    Res5_Diss[i] = 0.0;
                    DeltaT[i]    = 0.0;
                }
            }
        
            if (iRKS < (RKStage - 1)) {
                // Compute RESn
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
            }
        }

        // Precondition Variable Normalize
        if (PrecondMethod != PRECOND_METHOD_NONE) {
            for (int i = 0; i < nNode; i++)
                PrecondSigma[i] /= RKStage*(crs_IA_Node2Node[i+1] - crs_IA_Node2Node[i]);
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
    
    // Rd_Old
    if (Rd_Old != NULL)
        delete[] Rd_Old;
        
    // WDeltaT
    if (WDeltaT != NULL)
        delete[] WDeltaT;

    WQ      = NULL;
    Rd_Old  = NULL;
    WDeltaT = NULL;
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Explicit Steady State Mode
//------------------------------------------------------------------------------
int Solver_Steady_Explicit(void) {
    int rvalue = EXIT_FAILURE;
    
    // Select the Explicit Time Integration Method Type
    switch (TimeIntegrationMethod) {
        case TIME_INTEGRATION_METHOD_EXPLICIT_EULER:
            rvalue = Solver_Steady_Explicit_RungeKutta();
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK3:
            rvalue = Solver_Steady_Explicit_RungeKutta();
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK4:
            rvalue = Solver_Steady_Explicit_RungeKutta();
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK5:
            rvalue = Solver_Steady_Explicit_RungeKutta();
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK5_MJ:
            rvalue = Solver_Steady_Explicit_RungeKutta_Martinelli_Jameson();
            break;
        default:
            error("Solver_Steady_Explicit: Invalid Time Integration Method - %d", TimeIntegrationMethod);
            break;
    }
    
    return rvalue;
}

