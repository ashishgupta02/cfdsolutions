/*******************************************************************************
 * File:        Solver_Unsteady_Explicit.cpp
 * Author:      Ashish Gupta
 * Revision:    4
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
//! Solver in Unsteady State Mode with Explicit Mode
//! Dual Time: Euler(tau)-Euler(t) Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Explicit_Euler_Euler(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Unsteady State Mode with Explicit Mode
//! Dual Time: Euler(tau)-BDF(t) Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Explicit_Euler_BDF(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Unsteady State Mode with Explicit Mode
//! Dual Time: RK4(tau)-Euler(t) Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Explicit_RK4_Euler(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Unsteady State Mode with Explicit Mode
//! Dual Time: RK4(tau)-BDF(t) Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Explicit_RK4_BDF(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Unsteady State Mode with Explicit Mode
//! Dual Time: RK5(tau)-Euler(t) Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Explicit_RK5_Euler(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Unsteady State Mode with Explicit Mode
//! Dual Time: RK5(tau)-BDF(t) Time Integration
//------------------------------------------------------------------------------
int Solver_Unsteady_Explicit_RK5_BDF(void) {
    return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! Solver in Explicit Unsteady Mode
//------------------------------------------------------------------------------
int Solver_Unsteady_Explicit(void) {
    int AddTime     = FALSE;
    int CheckNAN    = 0;
    int SaveOrder   = 0;
    int SaveLimiter = 0;
    double *W01, *W02, *W03, *W04, *W05;
    double *WDeltaT;
    double alpha, beta, theta;
    double phi1, phi2, phi3, phi4;
    double eta, reta, zeta;
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

    // Qn1
    Qn11 = new double[nNode];
    Qn12 = new double[nNode];
    Qn13 = new double[nNode];
    Qn14 = new double[nNode];
    Qn15 = new double[nNode];
    
    // Check if Euler or Runge-Kutta Scheme
    W01 = W02 = W03 = W04 = W05 = NULL;
    WDeltaT = NULL;
    if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_EXPLICIT_RK4) {
        // Wo
        W01 = new double[nNode];
        W02 = new double[nNode];
        W03 = new double[nNode];
        W04 = new double[nNode];
        W05 = new double[nNode];
        
        // WDeltaT : Needed because DeltaT computations is inside Residual Computation
        WDeltaT = new double[nNode];
    }
    
    // Coefficients for First Order Physical Time Discretization 
    alpha = 1.0;
    beta  = 1.0;
    theta = 0.0;
    
    // Helper Variables
    eta  = 0.0;
    reta = 0.0;
    zeta = 0.0;
    
    // RK4 Coefficients
    phi1 = 0.25;
    phi2 = 0.333333333333333333;
    phi3 = 0.5;
    phi4 = 1.0;
    
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
        
        // Second Order Physical Time Coefficients
        if (t_iter > 0) {
            alpha = 3.0;
            beta  = 4.0;
            theta = 1.0;
        }
        
        // Set Qn and Qn-1
        for (int i = 0; i < nNode; i++) {
            // Qn-1
            Qn11[i] = Qn01[i];
            Qn12[i] = Qn02[i];
            Qn13[i] = Qn03[i];
            Qn14[i] = Qn04[i];
            Qn15[i] = Qn05[i];
            
            // Qn
            Qn01[i] = Q1[i];
            Qn02[i] = Q2[i];
            Qn03[i] = Q3[i];
            Qn04[i] = Q4[i];
            Qn05[i] = Q5[i];
        }
        
        // Inner Iterations: Dual Time Iterations
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
                MinPrecondSigma = DBL_MAX;
                MaxPrecondSigma = DBL_MIN;
            }

            // Compute Far Field Conditions with Mach Ramping
            ComputeFarFieldCondition(iter);

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

            AddTime = TRUE;
            // Compute Residuals
            Compute_Residual(AddTime);
            
            // Compute Local Time Stepping
            Compute_DeltaT(iter);

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

            // Write RMS
            RMS_Writer(giter, RMS);

            // Check for Residual NAN
            if (isnan(RMS_Res)) {
                info("Solve_Explicit_Unsteady: NAN Encountered ! - Abort");
                iter = InnerNIteration + 1;
                CheckNAN = 1;
                VTK_Writer("SolutionBeforeNAN.vtk", 1);
                break;
            }

            // Euler Time Integration Method
            if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_EXPLICIT_EULER) {
                // Update Conservative or Primitive Variables
                for (int i = 0; i < nNode; i++) {
                    eta    = DeltaT[i]/cVolume[i];
                    reta   = DeltaT[i]/PhysicalDeltaTime;
                    if (t_iter > 0) { // BDF
                        zeta = 1.0 + 1.5*reta;
                        reta = 0.5*reta;
                    } else // Euler
                        zeta = 1.0 + reta;
                    Q1[i] -= (eta * Res1[i] + reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] -= (eta * Res2[i] + reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] -= (eta * Res3[i] + reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] -= (eta * Res4[i] + reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] -= (eta * Res5[i] + reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;
                }
            }

            // Runge-Kutta 4 Time Integration Method
            if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_EXPLICIT_RK4) {
                for (int i = 0; i < nNode; i++) {
                    // W0 = Wn
                    W01[i]     = Q1[i];
                    W02[i]     = Q2[i];
                    W03[i]     = Q3[i];
                    W04[i]     = Q4[i];
                    W05[i]     = Q5[i];
                    WDeltaT[i] = DeltaT[i];
                    eta        = WDeltaT[i]/cVolume[i];
                    reta       = WDeltaT[i]/PhysicalDeltaTime;
                    if (t_iter > 0) { // BDF
                        zeta = 1.0 + 1.5*reta;
                        reta = 0.5*reta;
                    } else // Euler
                        zeta = 1.0 + reta;
                    // W1 = W0 - (phi1*dtau/zeta)*(RES0/vol + (0.5/dT)*(alpha*W0 - beta*Qn0 + theta*Qn1))
                    Q1[i] = W01[i] - phi1*(eta*Res1[i] + reta*(alpha*W01[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi1*(eta*Res2[i] + reta*(alpha*W02[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi1*(eta*Res3[i] + reta*(alpha*W03[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi1*(eta*Res4[i] + reta*(alpha*W04[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi1*(eta*Res5[i] + reta*(alpha*W05[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;

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

                AddTime = FALSE;
                // Compute Residuals
                Compute_Residual(AddTime);

                for (int i = 0; i < nNode; i++) {
                    // W2 = W0 - (phi2*dtau/zeta)*(RES1/vol + (0.5/dT)*(alpha*W1 - beta*Qn0 + theta*Qn1))
                    eta   = WDeltaT[i]/cVolume[i];
                    reta  = WDeltaT[i]/PhysicalDeltaTime;
                    if (t_iter > 0) { // BDF
                        zeta = 1.0 + 1.5*reta;
                        reta = 0.5*reta;
                    } else // Euler
                        zeta = 1.0 + reta;
                    Q1[i] = W01[i] - phi2*(eta*Res1[i] + reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi2*(eta*Res2[i] + reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi2*(eta*Res3[i] + reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi2*(eta*Res4[i] + reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi2*(eta*Res5[i] + reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;

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

                AddTime = FALSE;
                // Compute Residuals
                Compute_Residual(AddTime);

                for (int i = 0; i < nNode; i++) {
                    // W3 = W0 - (phi3*dtau/zeta)*(RES2/vol + (0.5/dT)*(alpha*W2 - beta*Qn0 + theta*Qn1))
                    eta   = WDeltaT[i]/cVolume[i];
                    reta  = WDeltaT[i]/PhysicalDeltaTime;
                    if (t_iter > 0) { // BDF
                        zeta = 1.0 + 1.5*reta;
                        reta = 0.5*reta;
                    } else // Euler
                        zeta = 1.0 + reta;
                    Q1[i] = W01[i] - phi3*(eta*Res1[i] + reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi3*(eta*Res2[i] + reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi3*(eta*Res3[i] + reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi3*(eta*Res4[i] + reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi3*(eta*Res5[i] + reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;

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

                AddTime = FALSE;
                // Compute Residuals
                Compute_Residual(AddTime);

                for (int i = 0; i < nNode; i++) {
                    // W4 = W0 - (phi4*dtau/zeta)*(RES3/vol + (0.5/dT)*(alpha*W2 - beta*Qn0 + theta*Qn1))
                    eta   = WDeltaT[i]/cVolume[i];
                    reta  = WDeltaT[i]/PhysicalDeltaTime;
                    if (t_iter > 0) { // BDF
                        zeta = 1.0 + 1.5*reta;
                        reta = 0.5*reta;
                    } else // Euler
                        zeta = 1.0 + reta;
                    Q1[i] = W01[i] - phi4*(eta*Res1[i] + reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi4*(eta*Res2[i] + reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi4*(eta*Res3[i] + reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi4*(eta*Res4[i] + reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi4*(eta*Res5[i] + reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;
                }
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
    
    if (TimeIntegrationMethod == TIME_INTEGRATION_METHOD_EXPLICIT_RK4) {
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
    }
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}
