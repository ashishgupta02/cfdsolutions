/*******************************************************************************
 * File:        Solver_Steady_Implicit.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "MC.h"
#include "Commons.h"
#include "DebugSolver.h"
#include "RestartIO.h"
#include "MeshIO.h"
#include "Solver.h"
#include "Gradient.h"
#include "CompressibleUtils.h"
#include "Residual_Smoothing.h"
#include "Material.h"
#include "Jacobian.h"


//------------------------------------------------------------------------------
//! Solver in Steady State Mode with Implicit Mode and Euler-Newton Time Integration
//------------------------------------------------------------------------------
int Solver_Steady_Implicit_Euler(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Implicit Steady State Mode
//------------------------------------------------------------------------------
int Solver_Steady_Implicit(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    double lrms     = 0.0;
    double RMS_Res;
    
    // Save Solver Parameters
    SaveOrder   = SolverOrder;
    SaveLimiter = LimiterMethod;
    
    // Initialize the Jacobian Data Structure
    Jacobian_Init();
    
    // Limit the Solution in Valid Range
    Material_Limit_Solution();
    
    // Smooth the Solution if Required
    if (SolutionSmooth == TRUE)
        CogSolver.Smooth_Solution();
    
    // Compute the Thermodynamic Properties of All Node (Physical + Ghosts)
    CogSolver.Compute_Bulk_Extended_Properties(0, nNode+nBNode);
    
    // Pseduo-Time Loop
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
        
        // Compute the Thermodynamic Properties of Changed Physical Nodes
        // Update is needed to aid the calculation of Jacobians of First Order
        CogSolver.Compute_Bulk_Extended_Properties(0, nNode);
        
        // Compute Least Square Gradient -- Unweighted
        if (SolverOrder == SOLVER_ORDER_SECOND) {
            Compute_Least_Square_Gradient(0);
            if (LimiterMethod != LIMITER_METHOD_NONE)
                Compute_Limiter();
        }
        
        // Apply boundary conditions
        Apply_Boundary_Condition(iter);
        
        // Update the Thermodynamic Properties of Boundary Ghost Nodes
        CogSolver.Compute_Bulk_Extended_Properties(nNode, nNode+nBNode);
        
        AddTime = TRUE;
        // Compute Residuals
        Compute_Residual(AddTime);
        
        // Compute Local Time Stepping
        Compute_DeltaT(iter);
        
        // Smooth the Residual: DeltaT Computation is required
        if (ResidualSmoothMethod != RESIDUAL_SMOOTH_METHOD_NONE)
            Residual_Smoothing();
        
        // Write RMS
        RMS_Res = RMS_Writer(iter+1);
        
        // Check for Residual NAN
        if (isnan(RMS_Res)) {
            info("Solve_Steady_Implicit: NAN Encountered ! - Abort");
            iter = SolverNIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }

        // Compute Jacobian and Fill CRS Matrix
        AddTime = TRUE;
        Compute_Jacobian(AddTime, iter);
        
        // Solve for Solution
//        lrms = MC_Iterative_Block_LU_Jacobi_CRS(LinearSolverNIteration, 2, SolverBlockMatrix);       
        lrms = MC_Iterative_Block_LUSGS_CRS(LinearSolverNIteration, SolverBlockMatrix);
        
        // Check for Linear Solver NAN
        if (isnan(lrms)) {
            info("Solve_Steady_Implicit: Liner Solver Iteration: NAN Encountered ! - Abort");
            iter = SolverNIteration + 1;
        }
        
        // Update Conservative Variables
        for (int i = 0; i < nNode; i++) {
            Q1[i] += Relaxation*SolverBlockMatrix.X[i][0];
            Q2[i] += Relaxation*SolverBlockMatrix.X[i][1];
            Q3[i] += Relaxation*SolverBlockMatrix.X[i][2];
            Q4[i] += Relaxation*SolverBlockMatrix.X[i][3];
            Q5[i] += Relaxation*SolverBlockMatrix.X[i][4];
        }
        
        // Limit the Solution in Valid Range
        Material_Limit_Solution();
        
        // Smooth the Solution if Required
        if (SolutionSmooth == TRUE)
            CogSolver.Smooth_Solution();
        
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

    // Finalize the Jacobian Data Structure
    Jacobian_Finalize();
    
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

