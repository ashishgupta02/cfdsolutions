/*******************************************************************************
 * File:        Solver_Unsteady_Implicit.cpp
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
//! Solver in Unsteady State Mode with Implicit Mode
//! Dual Time: Euler(tau)-Euler(t) Newton Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Implicit_Euler_Euler(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Unsteady State Mode with Implicit Mode
//! Dual Time: Euler(tau)-BDF(t) Newton Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Implicit_Euler_BDF(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Implicit Unsteady Mode
//------------------------------------------------------------------------------
int Solver_Unsteady_Implicit(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    double theta, eta, lrms;
    double *Qn01, *Qn02, *Qn03, *Qn04, *Qn05;
    double *Qn11, *Qn12, *Qn13, *Qn14, *Qn15;
    int giter;
    
    // Check if Unsteady Computations
    Qn01 = Qn02 = Qn03 = Qn04 = Qn05 = NULL;
    Qn11 = Qn12 = Qn13 = Qn14 = Qn15 = NULL;
    // Qn0
    Qn01 = new double[nNode];
    Qn02 = new double[nNode];
    Qn03 = new double[nNode];
    Qn04 = new double[nNode];
    Qn05 = new double[nNode];

    // Check if Second Order in Time Requested
    if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF) {
        // Qn1
        Qn11 = new double[nNode];
        Qn12 = new double[nNode];
        Qn13 = new double[nNode];
        Qn14 = new double[nNode];
        Qn15 = new double[nNode];
    }
    
     // Coefficients for First Order Physical Time Discretization 
    theta = 0.0;
    
    // Save Solver Parameters
    SaveOrder   = SolverOrder;
    SaveLimiter = LimiterMethod;
    
    Unsteady_Initialization();
    
    // Set Qn
    for (int i = 0; i < nNode; i++) {
        Qn01[i] = Q1[i];
        Qn02[i] = Q2[i];
        Qn03[i] = Q3[i];
        Qn04[i] = Q4[i];
        Qn05[i] = Q5[i];
    }
    
    giter = 0;
    // Physical Time Loop
    for (int t_iter = 0; t_iter < SolverNIteration; t_iter++) {
        
        // Check if Second Order in Time Requested
        if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF)
            if (t_iter > 0)
                theta = 0.5;
        
        // Set Qn and Qn-1
        for (int i = 0; i < nNode; i++) {
            // Check if Second Order in Time Requested
            if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF) {
                // Qn-1
                Qn11[i] = Qn01[i];
                Qn12[i] = Qn02[i];
                Qn13[i] = Qn03[i];
                Qn14[i] = Qn04[i];
                Qn15[i] = Qn05[i];
            }
            
            // Qn
            Qn01[i] = Q1[i];
            Qn02[i] = Q2[i];
            Qn03[i] = Q3[i];
            Qn04[i] = Q4[i];
            Qn05[i] = Q5[i];
        }
        
        // Inner Iterations: Newton Iterations
        for (int iter = 0; iter < InnerNIteration; iter++) {
            giter++;
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
                Compute_Least_Square_Gradient(0);
                if (LimiterMethod != LIMITER_METHOD_NONE)
                    Compute_Limiter();
            }

            // Apply boundary conditions
            Apply_Boundary_Condition(iter);

            AddTime = FALSE;
            // Compute Residuals
            Compute_Residual(AddTime);

//            // Compute Local Time Stepping
//            Compute_DeltaT(iter);
            
            // Compute Jacobian and Fill CRS Matrix
            AddTime = FALSE;
            Compute_Jacobian(AddTime, iter);

            // Add Physical Time Contribution to Residual and Jacobian
            // And Compute RMS with Time Term 
            int idgn;
            RMS[0] = RMS[1] = RMS[2] = RMS[3] = RMS[4] = 0.0;
            for (int i = 0; i < nNode; i++) {
                // Get the diagonal location
                idgn = SolverBlockMatrix.IAU[i];
                // Get the LHS
                eta = cVolume[i]/PhysicalDeltaTime;
                // Check if Second Order in Time Requested
                if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF) {
                    SolverBlockMatrix.B[i][0] -= eta*((1.0 + theta)*(Q1[i] - Qn01[i]) - theta*(Qn01[i] - Qn11[i]));
                    SolverBlockMatrix.B[i][1] -= eta*((1.0 + theta)*(Q2[i] - Qn02[i]) - theta*(Qn02[i] - Qn12[i]));
                    SolverBlockMatrix.B[i][2] -= eta*((1.0 + theta)*(Q3[i] - Qn03[i]) - theta*(Qn03[i] - Qn13[i]));
                    SolverBlockMatrix.B[i][3] -= eta*((1.0 + theta)*(Q4[i] - Qn04[i]) - theta*(Qn04[i] - Qn14[i]));
                    SolverBlockMatrix.B[i][4] -= eta*((1.0 + theta)*(Q5[i] - Qn05[i]) - theta*(Qn05[i] - Qn15[i]));
                } else if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_IMPLICIT_EULER_EULER) {
                    SolverBlockMatrix.B[i][0] -= eta*(Q1[i] - Qn01[i]);
                    SolverBlockMatrix.B[i][1] -= eta*(Q2[i] - Qn02[i]);
                    SolverBlockMatrix.B[i][2] -= eta*(Q3[i] - Qn03[i]);
                    SolverBlockMatrix.B[i][3] -= eta*(Q4[i] - Qn04[i]);
                    SolverBlockMatrix.B[i][4] -= eta*(Q5[i] - Qn05[i]);
                }
                
                for (int j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
                    for (int k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                        if (k == j)
                            SolverBlockMatrix.A[idgn][j][k] += (1.0 + theta)*eta;
                }
                
                // Compute Residual With Time Term
                for (int i = 0; i < nNode; i++) {
                    for (int j = 0; j < NEQUATIONS; j++)
                        RMS[j] += SolverBlockMatrix.B[i][j]*SolverBlockMatrix.B[i][j];
                }
            }
            
            // Compute RMS
            RMS_Res = RMS[0] + RMS[1] + RMS[2] + RMS[3] + RMS[4];
            RMS_Res = sqrt(RMS_Res/(5.0 * (double)nNode));
            RMS[0] = sqrt(RMS[0]/(double)nNode);
            RMS[1] = sqrt(RMS[1]/(double)nNode);
            RMS[2] = sqrt(RMS[2]/(double)nNode);
            RMS[3] = sqrt(RMS[3]/(double)nNode);
            RMS[4] = sqrt(RMS[4]/(double)nNode);

            // Write RMS
            RMS_Writer(giter, RMS);

            // Check for Residual NAN
            if (isnan(RMS_Res)) {
                info("Solve_Implicit_Unsteady: NAN Encountered ! - Abort");
                iter = InnerNIteration + 1;
                CheckNAN = 1;
                VTK_Writer("SolutionBeforeNAN.vtk", 1);
                break;
            }
            
            // Solve for Solution
            lrms = MC_Iterative_Block_LU_Jacobi_CRS(LinearSolverNIteration, 0, SolverBlockMatrix);
            printf("Titer: %d, GIter: %d, NIter: %d, LRMS: %10.5e\n", t_iter, giter, iter, lrms);
            
            // Check for Linear Solver NAN
            if (isnan(lrms)) {
                info("Solve_Implicit_Unsteady: Linear Solver Iteration: NAN Encountered ! - Abort");
                iter = InnerNIteration + 1;
                CheckNAN = 1;
                VTK_Writer("SolutionBeforeNAN.vtk", 1);
                break;
            }

            // Update Conservative/Primitive Variables
            for (int i = 0; i < nNode; i++) {
                Q1[i] += Relaxation*SolverBlockMatrix.X[i][0];
                Q2[i] += Relaxation*SolverBlockMatrix.X[i][1];
                Q3[i] += Relaxation*SolverBlockMatrix.X[i][2];
                Q4[i] += Relaxation*SolverBlockMatrix.X[i][3];
                Q5[i] += Relaxation*SolverBlockMatrix.X[i][4];
            }
        }
        
        // Check for NAN
        if (CheckNAN == 1) {
            t_iter = SolverNIteration + 1;
            break;
        }
        
        // Write the Unsteady Solution
        RestartIteration = t_iter+1;
        Check_Restart(t_iter);
    }
    
    // Check if Solution Restart is Requested
    if (RestartOutput && CheckNAN != 1)
        Restart_Writer(RestartOutputFilename, 1);

    // Debug the NAN
    if (CheckNAN)
        DebugNAN();
    
    // Free Memory
    // Qn0
    if (Qn01 != NULL)
        delete[] Qn01;
    if (Qn02 != NULL)
        delete[] Qn02;
    if (Qn03 != NULL)
        delete[] Qn03;
    if (Qn04 != NULL)
        delete[] Qn04;
    if (Qn05 != NULL)
        delete[] Qn05;
    
    // Check if Second Order in Time Requested
    if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_IMPLICIT_EULER_BDF) {
        // Qn1
        if (Qn11 != NULL)
            delete[] Qn11;
        if (Qn12 != NULL)
            delete[] Qn12;
        if (Qn13 != NULL)
            delete[] Qn13;
        if (Qn14 != NULL)
            delete[] Qn14;
        if (Qn15 != NULL)
            delete[] Qn15;
    }
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

