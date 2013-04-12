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
//! Solver in Steady State Mode with Explicit Mode and Euler Time Integration
//------------------------------------------------------------------------------
int Solver_Steady_Explicit_Euler(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    int SaveRMS     = 0;
    double eta      = 0.0;
    double scale_RMS[5], scale_RMS_Res;
    
    // Initialize
    for (int i = 0; i < NEQUATIONS; i++)
        scale_RMS[i] = 1.0;
    scale_RMS_Res = 1.0;
    
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
            RMS[0] += Res1[i]*Res1[i];
            RMS[1] += Res2[i]*Res2[i];
            RMS[2] += Res3[i]*Res3[i];
            RMS[3] += Res4[i]*Res4[i];
            RMS[4] += Res5[i]*Res5[i];
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

        // Explicit Euler Time Integration Method
        // Update Conservative or Primitive Variables
        for (int i = 0; i < nNode; i++) {
            eta    = DeltaT[i]/cVolume[i];
            Q1[i] -= eta * Res1[i];
            Q2[i] -= eta * Res2[i];
            Q3[i] -= eta * Res3[i];
            Q4[i] -= eta * Res4[i];
            Q5[i] -= eta * Res5[i]; 
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
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Steady State Mode with Explicit Mode and RK4 Time Integration
//------------------------------------------------------------------------------
int Solver_Steady_Explicit_RK4(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    int SaveRMS     = 0;
    double *W01, *W02, *W03, *W04, *W05;
    double *WDeltaT, *Dummy;
    double phi1, phi2, phi3, phi4;
    double eta;
    double scale_RMS[5], scale_RMS_Res;
    
    // Initialize
    for (int i = 0; i < NEQUATIONS; i++)
        scale_RMS[i] = 1.0;
    scale_RMS_Res = 1.0;
    
    // Check if Euler or Runge-Kutta Scheme
    W01 = W02 = W03 = W04 = W05 = NULL;
    WDeltaT = NULL;
    Dummy   = NULL;
    
    // Create the helper arrays for Runge-Kutta Method
    // Wo
    W01 = new double[nNode];
    W02 = new double[nNode];
    W03 = new double[nNode];
    W04 = new double[nNode];
    W05 = new double[nNode];

    // WDeltaT : Needed because DeltaT computations is inside Residual Computation
    WDeltaT = new double[nNode];
    
    // Helper Variables
    eta  = 0.0;
    
    // Runge-Kutta 4 Coefficients
    phi1 = 0.25;
    phi2 = 0.333333333333333333;
    phi3 = 0.5;
    phi4 = 1.0;
    
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
            RMS[0] += Res1[i]*Res1[i];
            RMS[1] += Res2[i]*Res2[i];
            RMS[2] += Res3[i]*Res3[i];
            RMS[3] += Res4[i]*Res4[i];
            RMS[4] += Res5[i]*Res5[i];
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
            W01[i]     = Q1[i];
            W02[i]     = Q2[i];
            W03[i]     = Q3[i];
            W04[i]     = Q4[i];
            W05[i]     = Q5[i];
            WDeltaT[i] = DeltaT[i];

            // W1 = W0 - phi*dT*RES0/vol
            eta   = WDeltaT[i]/cVolume[i];
            Q1[i] = W01[i] - phi1*eta*Res1[i];
            Q2[i] = W02[i] - phi1*eta*Res2[i];
            Q3[i] = W03[i] - phi1*eta*Res3[i];
            Q4[i] = W04[i] - phi1*eta*Res4[i];
            Q5[i] = W05[i] - phi1*eta*Res5[i];

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
            Q1[i] = W01[i] - phi2*eta*Res1[i];
            Q2[i] = W02[i] - phi2*eta*Res2[i];
            Q3[i] = W03[i] - phi2*eta*Res3[i];
            Q4[i] = W04[i] - phi2*eta*Res4[i];
            Q5[i] = W05[i] - phi2*eta*Res5[i];

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
            Q1[i] = W01[i] - phi3*eta*Res1[i];
            Q2[i] = W02[i] - phi3*eta*Res2[i];
            Q3[i] = W03[i] - phi3*eta*Res3[i];
            Q4[i] = W04[i] - phi3*eta*Res4[i];
            Q5[i] = W05[i] - phi3*eta*Res5[i];

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
            Q1[i] = W01[i] - phi4*eta*Res1[i];
            Q2[i] = W02[i] - phi4*eta*Res2[i];
            Q3[i] = W03[i] - phi4*eta*Res3[i];
            Q4[i] = W04[i] - phi4*eta*Res4[i];
            Q5[i] = W05[i] - phi4*eta*Res5[i];
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
    if (W01 != NULL)
        delete[] W01;
    if (W02 != NULL)
        delete[] W02;
    if (W03 != NULL)
        delete[] W03;
    if (W04 != NULL)
        delete[] W04;
    if (W05 != NULL)
        delete[] W05;

    // WDeltaT
    if (WDeltaT != NULL)
        delete[] WDeltaT;

    W01 = W02 = W03 = W04 = W05 = NULL;
    WDeltaT = NULL;
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Steady State Mode with Explicit Mode and RK5 Time Integration
//------------------------------------------------------------------------------
int Solver_Steady_Explicit_RK5(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    int SaveRMS     = 0;
    double *W01, *W02, *W03, *W04, *W05;
    double *WDeltaT, *Dummy;
    double phi1, phi2, phi3, phi4, phi5;
    double psi1, psi2, psi3, psi4, psi5;
    double eta;
    double scale_RMS[5], scale_RMS_Res;
    
    // Initialize
    for (int i = 0; i < NEQUATIONS; i++)
        scale_RMS[i] = 1.0;
    scale_RMS_Res = 1.0;
    
    // Runge-Kutta Scheme
    W01 = W02 = W03 = W04 = W05 = NULL;
    WDeltaT = NULL;
    Dummy   = NULL;
    
    // Create the helper arrays for Runge-Kutta Method
    // Wo
    W01 = new double[nNode];
    W02 = new double[nNode];
    W03 = new double[nNode];
    W04 = new double[nNode];
    W05 = new double[nNode];

    // WDeltaT : Needed because DeltaT computations is inside Residual Computation
    WDeltaT = new double[nNode];
    
    // Helper Variables
    eta  = 0.0;
    
    // Runge-Kutta 5 Coefficients
    phi1 = 1.0/4.0;
    phi2 = 1.0/6.0;
    phi3 = 3.0/8.0;
    phi4 = 1.0/2.0;
    phi5 = 1.0;
    
    psi1 = 1.0;
    psi2 = 0.0;
    psi3 = 14.0/25.0;
    psi4 = 0.0;
    psi5 = 11.0/25.0;
    
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
            RMS[0] += Res1[i]*Res1[i];
            RMS[1] += Res2[i]*Res2[i];
            RMS[2] += Res3[i]*Res3[i];
            RMS[3] += Res4[i]*Res4[i];
            RMS[4] += Res5[i]*Res5[i];
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
            W01[i]     = Q1[i];
            W02[i]     = Q2[i];
            W03[i]     = Q3[i];
            W04[i]     = Q4[i];
            W05[i]     = Q5[i];
            WDeltaT[i] = DeltaT[i];

            // W1 = W0 - phi*dT*RES0/vol
            eta   = WDeltaT[i]/cVolume[i];
            Q1[i] = W01[i] - phi1*eta*Res1[i];
            Q2[i] = W02[i] - phi1*eta*Res2[i];
            Q3[i] = W03[i] - phi1*eta*Res3[i];
            Q4[i] = W04[i] - phi1*eta*Res4[i];
            Q5[i] = W05[i] - phi1*eta*Res5[i];

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
            Q1[i] = W01[i] - phi2*eta*Res1[i];
            Q2[i] = W02[i] - phi2*eta*Res2[i];
            Q3[i] = W03[i] - phi2*eta*Res3[i];
            Q4[i] = W04[i] - phi2*eta*Res4[i];
            Q5[i] = W05[i] - phi2*eta*Res5[i];

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
            Q1[i] = W01[i] - phi3*eta*Res1[i];
            Q2[i] = W02[i] - phi3*eta*Res2[i];
            Q3[i] = W03[i] - phi3*eta*Res3[i];
            Q4[i] = W04[i] - phi3*eta*Res4[i];
            Q5[i] = W05[i] - phi3*eta*Res5[i];

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
            Q1[i] = W01[i] - phi4*eta*Res1[i];
            Q2[i] = W02[i] - phi4*eta*Res2[i];
            Q3[i] = W03[i] - phi4*eta*Res3[i];
            Q4[i] = W04[i] - phi4*eta*Res4[i];
            Q5[i] = W05[i] - phi4*eta*Res5[i];
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
    if (W01 != NULL)
        delete[] W01;
    if (W02 != NULL)
        delete[] W02;
    if (W03 != NULL)
        delete[] W03;
    if (W04 != NULL)
        delete[] W04;
    if (W05 != NULL)
        delete[] W05;

    // WDeltaT
    if (WDeltaT != NULL)
        delete[] WDeltaT;

    W01 = W02 = W03 = W04 = W05 = NULL;
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
            rvalue = Solver_Steady_Explicit_Euler();
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK4:
            rvalue = Solver_Steady_Explicit_RK4();
            break;
        case TIME_INTEGRATION_METHOD_EXPLICIT_RK5:
            rvalue = Solver_Steady_Explicit_RK5();
            break;
        default:
            error("Solver_Steady_Explicit: Invalid Time Integration Method - %d", TimeIntegrationMethod);
            break;
    }
    
    return rvalue;
}

