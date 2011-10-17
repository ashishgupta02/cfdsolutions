/*******************************************************************************
 * File:        Solver.cpp
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
#include "DebugSolver.h"

// Local Time
double *DeltaT;

// Conservative Variables
double *Q1;
double *Q2;
double *Q3;
double *Q4;
double *Q5;

// Conservative Variables Gradients
double *Q1x;
double *Q1y;
double *Q1z;

double *Q2x;
double *Q2y;
double *Q2z;

double *Q3x;
double *Q3y;
double *Q3z;

double *Q4x;
double *Q4y;
double *Q4z;

double *Q5x;
double *Q5y;
double *Q5z;

// Residuals (Convective and Dissipative)
double *Res1;
double *Res2;
double *Res3;
double *Res4;
double *Res5;

double *Res1_Diss;
double *Res2_Diss;
double *Res3_Diss;
double *Res4_Diss;
double *Res5_Diss;

// Limiters
double *Limiter_Phi1;
double *Limiter_Phi2;
double *Limiter_Phi3;
double *Limiter_Phi4;
double *Limiter_Phi5;

//RMS
double RMS[5];
double RMS_Res;

// MC CRS Matrix
MC_CRS SolverBlockMatrix;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Init(void) {
    // Local Time
    DeltaT          = NULL;

    // Conservative/Primitive Variables
    Q1              = NULL;
    Q2              = NULL;
    Q3              = NULL;
    Q4              = NULL;
    Q5              = NULL;

    // Conservative/Primitive Variables Gradients
    Q1x             = NULL;
    Q1y             = NULL;
    Q1z             = NULL;

    Q2x             = NULL;
    Q2y             = NULL;
    Q2z             = NULL;

    Q3x             = NULL;
    Q3y             = NULL;
    Q3z             = NULL;

    Q4x             = NULL;
    Q4y             = NULL;
    Q4z             = NULL;

    Q5x             = NULL;
    Q5y             = NULL;
    Q5z             = NULL;
    
    // Residuals (Convective and Dissipative)
    Res1            = NULL;
    Res2            = NULL;
    Res3            = NULL;
    Res4            = NULL;
    Res5            = NULL;

    Res1_Diss       = NULL;
    Res2_Diss       = NULL;
    Res3_Diss       = NULL;
    Res4_Diss       = NULL;
    Res5_Diss       = NULL;
    
    // Limiters
    Limiter_Phi1    = NULL;
    Limiter_Phi2    = NULL;
    Limiter_Phi3    = NULL;
    Limiter_Phi4    = NULL;
    Limiter_Phi5    = NULL;
    
    // RMS
    RMS[0]          = 0.0;
    RMS[1]          = 0.0;
    RMS[2]          = 0.0;
    RMS[3]          = 0.0;
    RMS[4]          = 0.0;
    RMS_Res         = 0.0;
    
    // MC CRS Matrix
    SolverBlockMatrix.Block_nRow = 0;
    SolverBlockMatrix.Block_nCol = 0;
    SolverBlockMatrix.nROW       = 0;
    SolverBlockMatrix.nCOL       = 0;
    SolverBlockMatrix.DIM        = 0;
    SolverBlockMatrix.A          = NULL;
    SolverBlockMatrix.B          = NULL;
    SolverBlockMatrix.X          = NULL;
    SolverBlockMatrix.IA         = NULL;
    SolverBlockMatrix.IAU        = NULL;
    SolverBlockMatrix.JA         = NULL;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Finalize(void) {
    // Conservative Variables
    if (Q1 != NULL)
        delete[] Q1;
    if (Q2 != NULL)
        delete[] Q2;
    if (Q3 != NULL)
        delete[] Q3;
    if (Q4 != NULL)
        delete[] Q4;
    if (Q5 != NULL)
        delete[] Q5;
    Q1 = NULL;
    Q2 = NULL;
    Q3 = NULL;
    Q4 = NULL;
    Q5 = NULL;

     // Conservative Variables Gradient
    if (Q1x != NULL)
        delete[] Q1x;
    if (Q1y != NULL)
        delete[] Q1y;
    if (Q1z != NULL)
        delete[] Q1z;
    Q1x = NULL;
    Q1y = NULL;
    Q1z = NULL;

    if (Q2x != NULL)
        delete[] Q2x;
    if (Q2y != NULL)
        delete[] Q2y;
    if (Q2z != NULL)
        delete[] Q2z;
    Q2x = NULL;
    Q2y = NULL;
    Q2z = NULL;

    if (Q3x != NULL)
        delete[] Q3x;
    if (Q3y != NULL)
        delete[] Q3y;
    if (Q3z != NULL)
        delete[] Q3z;
    Q3x = NULL;
    Q3y = NULL;
    Q3z = NULL;
    
    if (Q4x != NULL)
        delete[] Q4x;
    if (Q4y != NULL)
        delete[] Q4y;
    if (Q4z != NULL)
        delete[] Q4z;
    Q4x = NULL;
    Q4y = NULL;
    Q4z = NULL;

    if (Q5x != NULL)
        delete[] Q5x;
    if (Q5y != NULL)
        delete[] Q5y;
    if (Q5z != NULL)
        delete[] Q5z;
    Q5x = NULL;
    Q5y = NULL;
    Q5z = NULL;

    // Limiters
    if (Limiter_Phi1 != NULL)
        delete[] Limiter_Phi1;
    if (Limiter_Phi2 != NULL)
        delete[] Limiter_Phi2;
    if (Limiter_Phi3 != NULL)
        delete[] Limiter_Phi3;
    if (Limiter_Phi4 != NULL)
        delete[] Limiter_Phi4;
    if (Limiter_Phi5 != NULL)
        delete[] Limiter_Phi5;
    Limiter_Phi1 = NULL;
    Limiter_Phi2 = NULL;
    Limiter_Phi3 = NULL;
    Limiter_Phi4 = NULL;
    Limiter_Phi5 = NULL;
    
    // Residuals (Convective and Dissipative)
    if (Res1 != NULL)
        delete[] Res1;
    if (Res2 != NULL)
        delete[] Res2;
    if (Res3 != NULL)
        delete[] Res3;
    if (Res4 != NULL)
        delete[] Res4;
    if (Res5 != NULL)
        delete[] Res5;
    Res1 = NULL;
    Res2 = NULL;
    Res3 = NULL;
    Res4 = NULL;
    Res5 = NULL;

    if (Res1_Diss != NULL)
        delete[] Res1_Diss;
    if (Res2_Diss != NULL)
        delete[] Res2_Diss;
    if (Res3_Diss != NULL)
        delete[] Res3_Diss;
    if (Res4_Diss != NULL)
        delete[] Res4_Diss;
    if (Res5_Diss != NULL)
        delete[] Res5_Diss;
    Res1_Diss = NULL;
    Res2_Diss = NULL;
    Res3_Diss = NULL;
    Res4_Diss = NULL;
    Res5_Diss = NULL;
    
    // Local Time
    if (DeltaT != NULL)
        delete[] DeltaT;
    DeltaT = NULL;
    
    // Finalize Gradient Infrastructure
    if (Order == 2)
        Gradient_Finalize();
    
    // Check if Implicit Method
    if (SolverMethod == SOLVER_METHOD_IMPLICIT)
        Delete_CRS_SolverBlockMatrix();
    
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//! Set Initial Conditions
//------------------------------------------------------------------------------
void Solver_Set_Initial_Conditions(void) {

    // Allocate Memory to Store Conservative Variables
    Q1 = new double[nNode + nBNode];
    Q2 = new double[nNode + nBNode];
    Q3 = new double[nNode + nBNode];
    Q4 = new double[nNode + nBNode];
    Q5 = new double[nNode + nBNode];

    // Compute Free Stream Conditions
    ComputeFreeStreamCondition(0);
    
    // Initialize the variable with reference conditions
    for (int i = 0; i < (nNode + nBNode); i++) {
        Q1[i] = Inf_Rho;
        Q2[i] = Inf_Rho*Inf_U;
        Q3[i] = Inf_Rho*Inf_V;
        Q4[i] = Inf_Rho*Inf_W;
        Q5[i] = Inf_Rho*Inf_Et;
    }

    // Allocate Memory to Store Conservative Variables Gradients
    if (Order == 2) {
        // Initialize Gradient Infrastructure
        Gradient_Init();

        Q1x = new double[nNode];
        Q1y = new double[nNode];
        Q1z = new double[nNode];
        Gradient_Add_Function(Q1, Q1x, Q1y, Q1z, nNode);
        
        Q2x = new double[nNode];
        Q2y = new double[nNode];
        Q2z = new double[nNode];
        Gradient_Add_Function(Q2, Q2x, Q2y, Q2z, nNode);

        Q3x = new double[nNode];
        Q3y = new double[nNode];
        Q3z = new double[nNode];
        Gradient_Add_Function(Q3, Q3x, Q3y, Q3z, nNode);

        Q4x = new double[nNode];
        Q4y = new double[nNode];
        Q4z = new double[nNode];
        Gradient_Add_Function(Q4, Q4x, Q4y, Q4z, nNode);

        Q5x = new double[nNode];
        Q5y = new double[nNode];
        Q5z = new double[nNode];
        Gradient_Add_Function(Q5, Q5x, Q5y, Q5z, nNode);

        // Check if Limiter Memory is Required
        if (Limiter > 0) {
            Limiter_Phi1 = new double[nNode];
            Limiter_Phi2 = new double[nNode];
            Limiter_Phi3 = new double[nNode];
            Limiter_Phi4 = new double[nNode];
            Limiter_Phi5 = new double[nNode];
        }
    }
    
    // Allocate  Memory to Store Residuals (Convective and Dissipative)
    Res1 = new double[nNode];
    Res2 = new double[nNode];
    Res3 = new double[nNode];
    Res4 = new double[nNode];
    Res5 = new double[nNode];
    for (int i = 0; i < nNode; i++) {
        Res1[i] = 0.0;
        Res2[i] = 0.0;
        Res3[i] = 0.0;
        Res4[i] = 0.0;
        Res5[i] = 0.0;
    }
    
    Res1_Diss = new double[nNode];
    Res2_Diss = new double[nNode];
    Res3_Diss = new double[nNode];
    Res4_Diss = new double[nNode];
    Res5_Diss = new double[nNode];
    for (int i = 0; i < nNode; i++) {
        Res1_Diss[i] = 0.0;
        Res2_Diss[i] = 0.0;
        Res3_Diss[i] = 0.0;
        Res4_Diss[i] = 0.0;
        Res5_Diss[i] = 0.0;
    }
    
    // Allocate Memory to Store Time Step
    DeltaT = new double[nNode];
    // Initialize
    for (int i = 0; i < nNode; i++)
        DeltaT[i] = 0.0;
    
    // Set Boundary Conditions
    Initialize_Boundary_Condition();
    
    // Check if Implicit Method
    if (SolverMethod == SOLVER_METHOD_IMPLICIT)
        Create_CRS_SolverBlockMatrix();
}

//------------------------------------------------------------------------------
//! Solver in Explicit Mode
//------------------------------------------------------------------------------
int Solve_Explicit(void) {
    int CheckNAN  = 0;
    int SaveOrder = 0;
    int SaveLimiter = 0;
    double dtmp;
    double *W01, *W02, *W03, *W04, *W05;
    double *W41, *W42, *W43, *W44, *W45;
    
    // Check if Euler or Runge-Kutta Scheme
    W01 = W02 = W03 = W04 = W05 = NULL;
    W41 = W42 = W43 = W44 = W45 = NULL;
    if (TimeStepScheme == 2) {
        // Wo
        W01 = new double[nNode];
        W02 = new double[nNode];
        W03 = new double[nNode];
        W04 = new double[nNode];
        W05 = new double[nNode];

        // W4
        W41 = new double[nNode];
        W42 = new double[nNode];
        W43 = new double[nNode];
        W44 = new double[nNode];
        W45 = new double[nNode];
    }

    // Save Solver Parameters
    SaveOrder   = Order;
    SaveLimiter = Limiter;

    printf("=============================================================================\n");
    printf("-----------------------------------------------------------------------------\n");
    printf(" Iter        RMS_RHO    RMS_RHOU    RMS_RHOV    RMS_RHOW    RMS_E     RMS_RES\n");
    printf("-----------------------------------------------------------------------------\n");
    for (int iter = RestartIteration; iter < NIteration; iter++) {
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

        // Compute Free Stream Conditions with Mach Ramping
        ComputeFreeStreamCondition(iter);
        
        // Check if First Order Iterations are Required
        Order = SaveOrder;
        if (FirstOrderNIteration > iter)
            Order = 1;

        // Set the Start and End of Limiter Iterations
        Limiter = 0;
        if (StartLimiterNIteration < EndLimiterNIteration) {
            if ((StartLimiterNIteration <= iter+1) && (EndLimiterNIteration > iter))
                Limiter = SaveLimiter;
        }

        // Compute Least Square Gradient -- Unweighted
        if (Order == 2) {
            Compute_Least_Square_Gradient(0);
            if (Limiter > 0)
                Compute_Limiter();
        }
        
        // Apply boundary conditions
        Apply_Boundary_Condition(iter);

        // Compute Residuals
        Compute_Residual();
        
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

        printf("%5d %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",
                iter+1, RMS[0], RMS[1], RMS[2], RMS[3], RMS[4], RMS_Res);

        if (isnan(RMS_Res)) {
            info("Solve: NAN Encountered ! - Abort");
            iter = NIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }

        // Euler Time Stepping Scheme
        if (TimeStepScheme == 1) {
            // Update Conservative Variables
            for (int i = 0; i < nNode; i++) {
                // Compute Weiss Smith Precondition Variables
                if (PrecondMethod == SOLVER_PRECOND_ROE_WS) {
                    double rho, u, v, w, ht, rhoet, p, c, Ur, eps, mach, Cp, Theta, drhodt;
                    double lres1, lres2, lres3, lres4, lres5;
                    double Qc[5], Qp[5];
                    rho   = Q1[i];
                    u     = Q2[i]/rho;
                    v     = Q3[i]/rho;
                    w     = Q4[i]/rho;
                    rhoet = Q5[i];
                    p     = (Gamma - 1.0)*(rhoet - 0.5*rho*(u*u + v*v + w*w)) + Gauge_Pressure;
                    ht    = (rhoet + p)/rho;
                    c     = sqrt((Gamma * p) / rho);
                    
                    mach = sqrt(u*u + v*v + w*w)/c;
                    eps  = MIN(1.0, MAX(1e-10, mach*mach));
                    
                    // Compute Ur for Ideal Gas
                    dtmp = sqrt(u*u + v*v + w*w);
                    if (dtmp < eps*c)
                        Ur = eps*c;
                    else if ((eps*c < dtmp) && (dtmp < c))
                        Ur = dtmp;
                    else
                        Ur = c;
                    
                    Theta  = (1.0/(Ur*Ur) + (Gamma - 1.0)/(c*c));
                    Cp     = 1/(Gamma - 1.0);
                    drhodt = -(rho*rho)/(Gamma*p);
                    
                    dtmp = DeltaT[i]/cVolume[i];
                    lres1 = ht*Res1[i] - Res5[i] + u*Res2[i] - u*u*Res1[i] + v*Res3[i] - v*v*Res1[i] + w*Res4[i] - w*w*Res1[i];
                    lres1 = dtmp*(drhodt*lres1 + Cp*rho*Res1[i])/(drhodt + Cp*Theta*rho);
                    lres2 = dtmp*(Res2[i] - u*Res1[i])/rho;
                    lres3 = dtmp*(Res3[i] - v*Res1[i])/rho;
                    lres4 = dtmp*(Res4[i] - w*Res1[i])/rho;
                    lres5 = Theta*(Res5[i] - u*Res2[i] - v*Res3[i] - w*Res4[i]) + Res1[i]*(1.0 + Theta*(u*u + v*v + w*w - ht));
                    lres5 = dtmp*lres5/(drhodt + Cp*Theta*rho);
                    
                    // Convert Conservative to Primitive
                    Qc[0] = Q1[i];
                    Qc[1] = Q2[i];
                    Qc[2] = Q3[i];
                    Qc[3] = Q4[i];
                    Qc[4] = Q5[i];
                    ConservativeToPressureVelocityTemperature(Qc, Qp, Gamma, Gauge_Pressure);
                    Qp[0] -= lres1;
                    Qp[1] -= lres2;
                    Qp[2] -= lres3;
                    Qp[3] -= lres4;
                    Qp[4] -= lres5;
                    PressureVelocityTemperatureToConservative(Qp, Qc, Gamma, Gauge_Pressure);
                    Q1[i] = Qc[0];
                    Q2[i] = Qc[1];
                    Q3[i] = Qc[2];
                    Q4[i] = Qc[3];
                    Q5[i] = Qc[4];
                } else {
                    dtmp = DeltaT[i]/cVolume[i];
                    Q1[i] -= dtmp * Res1[i];
                    Q2[i] -= dtmp * Res2[i];
                    Q3[i] -= dtmp * Res3[i];
                    Q4[i] -= dtmp * Res4[i];
                    Q5[i] -= dtmp * Res5[i];
                }
            }
        }

        // Runge-Kutta Time Stepping Scheme
        if (TimeStepScheme == 2) {
            for (int i = 0; i < nNode; i++) {
                // W0 = Wn
                W01[i] = Q1[i];
                W02[i] = Q2[i];
                W03[i] = Q3[i];
                W04[i] = Q4[i];
                W05[i] = Q5[i];

                // W1 = W0 - (dT/2)*RES0
                dtmp = DeltaT[i]/(2.0*cVolume[i]);
                Q1[i] = W01[i] - dtmp*Res1[i];
                Q2[i] = W02[i] - dtmp*Res2[i];
                Q3[i] = W03[i] - dtmp*Res3[i];
                Q4[i] = W04[i] - dtmp*Res4[i];
                Q5[i] = W05[i] - dtmp*Res5[i];

                // W4 = W0 - (dT/6)*RES0
                dtmp = DeltaT[i]/(6.0*cVolume[i]);
                W41[i] = W01[i] - dtmp*Res1[i];
                W42[i] = W02[i] - dtmp*Res2[i];
                W43[i] = W03[i] - dtmp*Res3[i];
                W44[i] = W04[i] - dtmp*Res4[i];
                W45[i] = W05[i] - dtmp*Res5[i];

                // Reset the Residual
                Res1[i] = 0.0;
                Res2[i] = 0.0;
                Res3[i] = 0.0;
                Res4[i] = 0.0;
                Res5[i] = 0.0;
            }
            
            // Compute RES1
            // Compute Least Square Gradient -- Unweighted
            if (Order == 2) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
                    Compute_Limiter();
            }

            // Apply boundary conditions
            Apply_Boundary_Condition(iter);

            // Compute Residuals
            Compute_Residual();
            
            for (int i = 0; i < nNode; i++) {
                // W2 = W0 - (dT/2)*RES1
                dtmp = DeltaT[i]/(2.0*cVolume[i]);
                Q1[i] = W01[i] - dtmp*Res1[i];
                Q2[i] = W02[i] - dtmp*Res2[i];
                Q3[i] = W03[i] - dtmp*Res3[i];
                Q4[i] = W04[i] - dtmp*Res4[i];
                Q5[i] = W05[i] - dtmp*Res5[i];

                // W4 = W0 - (dT/6)*RES0 - (dT/3)*RES1
                dtmp = DeltaT[i]/(3.0*cVolume[i]);
                W41[i] = W41[i] - dtmp*Res1[i];
                W42[i] = W42[i] - dtmp*Res2[i];
                W43[i] = W43[i] - dtmp*Res3[i];
                W44[i] = W44[i] - dtmp*Res4[i];
                W45[i] = W45[i] - dtmp*Res5[i];

                // Reset the Residual
                Res1[i] = 0.0;
                Res2[i] = 0.0;
                Res3[i] = 0.0;
                Res4[i] = 0.0;
                Res5[i] = 0.0;
            }

            // Compute RES2
            // Compute Least Square Gradient -- Unweighted
            if (Order == 2) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
                    Compute_Limiter();
            }

            // Apply boundary conditions
            Apply_Boundary_Condition(iter);

            // Compute Residuals
            Compute_Residual();

            for (int i = 0; i < nNode; i++) {
                // W3 = W0 - dT*RES2
                dtmp = DeltaT[i]/cVolume[i];
                Q1[i] = W01[i] - dtmp*Res1[i];
                Q2[i] = W02[i] - dtmp*Res2[i];
                Q3[i] = W03[i] - dtmp*Res3[i];
                Q4[i] = W04[i] - dtmp*Res4[i];
                Q5[i] = W05[i] - dtmp*Res5[i];

                // W4 = W0 - (dT/6)*RES0 - (dT/3)*RES1 - (dT/3)*RES2
                dtmp = DeltaT[i]/(3.0*cVolume[i]);
                W41[i] = W41[i] - dtmp*Res1[i];
                W42[i] = W42[i] - dtmp*Res2[i];
                W43[i] = W43[i] - dtmp*Res3[i];
                W44[i] = W44[i] - dtmp*Res4[i];
                W45[i] = W45[i] - dtmp*Res5[i];

                // Reset the Residual
                Res1[i] = 0.0;
                Res2[i] = 0.0;
                Res3[i] = 0.0;
                Res4[i] = 0.0;
                Res5[i] = 0.0;
            }

            // Compute RES3
            // Compute Least Square Gradient -- Unweighted
            if (Order == 2) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
                    Compute_Limiter();
            }

            // Apply boundary conditions
            Apply_Boundary_Condition(iter);

            // Compute Residuals
            Compute_Residual();

            for (int i = 0; i < nNode; i++) {
                // W4 = W0 - (dT/6)*RES0 - (dT/3)*RES1 - (dT/3)*RES2 - (dT/6)*RES3
                dtmp = DeltaT[i]/(6.0*cVolume[i]);
                Q1[i] = W41[i] - dtmp*Res1[i];
                Q2[i] = W42[i] - dtmp*Res2[i];
                Q3[i] = W43[i] - dtmp*Res3[i];
                Q4[i] = W44[i] - dtmp*Res4[i];
                Q5[i] = W45[i] - dtmp*Res5[i];
            }
        }

        // Check Cyclic Restart is Requested
        RestartIteration = iter+1;
        Check_Restart(iter);
        
        if ((RMS_Res < (DBL_EPSILON*10.0))|| ((iter+1) == NIteration)) {
            iter = NIteration + 1;
        }
    }

    // Check if Solution Restart is Requested
    if (RestartOutput && CheckNAN != 1)
        Restart_Writer(RestartOutputFilename, 1);

    // Debug the NAN
    if (CheckNAN)
        DebugNAN();

    // Free Memory
    if (TimeStepScheme == 2) {
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

        // W4
        if (W41 != NULL)
            delete[] W41;
        if (W42 != NULL)
            delete[] W42;
        if (W43 != NULL)
            delete[] W43;
        if (W44 != NULL)
            delete[] W44;
        if (W45 != NULL)
            delete[] W45;
        
        W01 = W02 = W03 = W04 = W05 = NULL;
        W41 = W42 = W43 = W44 = W45 = NULL;
    }
    
    // Return Solver State
    if (CheckNAN)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
int Solve_Implicit(void) {
    int AddTime   = 0;
    int CheckNAN  = 0;
    int SaveOrder = 0;
    int SaveLimiter = 0;
    double lrms;
    
    // Save Solver Parameters
    SaveOrder   = Order;
    SaveLimiter = Limiter;

    printf("=============================================================================\n");
    printf("-----------------------------------------------------------------------------\n");
    printf(" Iter        RMS_RHO    RMS_RHOU    RMS_RHOV    RMS_RHOW    RMS_E     RMS_RES\n");
    printf("-----------------------------------------------------------------------------\n");
    for (int iter = RestartIteration; iter < NIteration; iter++) {
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

        // Compute Free Stream Conditions with Mach Ramping
        ComputeFreeStreamCondition(iter);
    
        // Check if First Order Iterations are Required
        Order = SaveOrder;
        if (FirstOrderNIteration > iter)
            Order = 1;

        // Set the Start and End of Limiter Iterations
        Limiter = 0;
        if (StartLimiterNIteration < EndLimiterNIteration) {
            if ((StartLimiterNIteration <= iter+1) && (EndLimiterNIteration > iter))
                Limiter = SaveLimiter;
        }

        // Compute Least Square Gradient -- Unweighted
        if (Order == 2) {
            Compute_Least_Square_Gradient(0);
            if (Limiter > 0)
                Compute_Limiter();
        }
        
        // Apply boundary conditions
        Apply_Boundary_Condition(iter);

        // Compute Residuals
        Compute_Residual();
        
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

        printf("%5d %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n",
                iter+1, RMS[0], RMS[1], RMS[2], RMS[3], RMS[4], RMS_Res);

        if (isnan(RMS_Res)) {
            info("Solve: NAN Encountered ! - Abort");
            iter = NIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }

        // Compute Jacobian and Fill CRS Matrix
        AddTime = 1;
        Compute_Jacobian(AddTime, iter);
        
        // Solve for Solution
        lrms = MC_Iterative_Block_LU_Jacobi_CRS(InnerNIteration, 0, SolverBlockMatrix);
        
        if (isnan(lrms)) {
            info("Solve: Liner Solver Iteration: NAN Encountered ! - Abort");
            iter = NIteration + 1;
        }
        
        // Update Conservative Variables
        for (int i = 0; i < nNode; i++) {
            Q1[i] += Relaxation*SolverBlockMatrix.X[i][0];
            Q2[i] += Relaxation*SolverBlockMatrix.X[i][1];
            Q3[i] += Relaxation*SolverBlockMatrix.X[i][2];
            Q4[i] += Relaxation*SolverBlockMatrix.X[i][3];
            Q5[i] += Relaxation*SolverBlockMatrix.X[i][4];
        }

        // Check Cyclic Restart is Requested
        RestartIteration = iter+1;
        Check_Restart(iter);
        
        if ((RMS_Res < (DBL_EPSILON*10.0))|| ((iter+1) == NIteration)) {
            iter = NIteration + 1;
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
//! 
//------------------------------------------------------------------------------
int Solve(void) {
    int rvalue = EXIT_FAILURE;
    
    // Initialize the Solver Scheme Data Structure
    switch (SolverScheme) {
        case SOLVER_SCHEME_ROE: // Roe
            Roe_Init();
            break;
        case SOLVER_SCHEME_LMROE: // LMRoe
            Roe_Init();
            break;
        case SOLVER_SCHEME_ROE_WS: // Roe Weiss Smith Precondition
            Roe_Init();
            break;
        case SOLVER_SCHEME_ROE_CV: // Roe Cecile Voizat Precondition
            Roe_Init();
            break;
        default:
            error("Solve: Invalid Solver Scheme - %d", SolverScheme);
            break;
    }
    
    // Check if Solution Restart is Requested
    if (RestartInput)
        Restart_Reader(RestartInputFilename);
    
    // Select the Solver Method Type
    switch (SolverMethod) {
        case SOLVER_METHOD_EXPLICIT:
            rvalue = Solve_Explicit();
            break;
        case SOLVER_METHOD_IMPLICIT:
            rvalue = Solve_Implicit();
            break;
        default:
            error("Solve: Invalid Solver Method - %d", SolverMethod);
            break;
    }
    
    // Finalize the Solver Scheme Data Structure
    switch (SolverScheme) {
        case SOLVER_SCHEME_ROE: // Roe
            Roe_Finalize();
            break;
        case SOLVER_SCHEME_LMROE: // LMRoe
            Roe_Finalize();
            break;
        case SOLVER_SCHEME_ROE_WS: // Roe Weiss Smith Precondition
            Roe_Finalize();
            break;
        case SOLVER_SCHEME_ROE_CV: // Roe Cecile Voizat Precondition
            Roe_Finalize();
            break;
        default:
            error("Solve: Invalid Solver Scheme - %d", SolverScheme);
            break;
    }
    
    return rvalue;
}

//------------------------------------------------------------------------------
//! Compute Free Stream Condition with Mach Ramping
//------------------------------------------------------------------------------
void ComputeFreeStreamCondition(int Iteration) {
    double tmpMach = 0.0;

    // Compute the Free Stream Condition Ramping
    if ((Mach_Ramp > 1) && (Mach_MAX > Mach_MIN)) {
        if (Iteration < Mach_Ramp)
            tmpMach = Mach_MIN + (Mach_MAX - Mach_MIN)*(((double)Iteration)/((double)(Mach_Ramp-1)));
        else
            tmpMach = Mach_MAX;
    } else
        tmpMach = Mach_MAX;

    // Set the Free Stream Conditions
    Inf_Mach     = tmpMach;
    Inf_Rho      = Ref_Rho;
    Inf_Pressure = Ref_Pressure - Gauge_Pressure;
    Inf_U        = Inf_Mach*cos(Ref_Alpha)*cos(Ref_Beta);
    Inf_V        = Inf_Mach*sin(Ref_Beta);
    Inf_W        = Inf_Mach*sin(Ref_Alpha)*cos(Ref_Beta);
    Inf_Et       = (Ref_Pressure - Gauge_Pressure)/((Gamma - 1.0)*Inf_Rho) + 0.5 *(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W);
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
void Create_CRS_SolverBlockMatrix(void) {
    int i, j, jstart, jend, k, ksave;
    int degree, index, min, minsave;
    
    SolverBlockMatrix.nROW       = nNode;
    SolverBlockMatrix.nCOL       = nNode;
    SolverBlockMatrix.Block_nRow = NEQUATIONS;
    SolverBlockMatrix.Block_nCol = NEQUATIONS;
    
    // Allocate Memory of IA Array to Store Start and End Location of Row
    SolverBlockMatrix.IA = new int[nNode+1];
    
    // Start Filling the Row Location and
    // Get the No of Non Zero Entries for SolverBlockMatrix
    SolverBlockMatrix.DIM = 0;
    SolverBlockMatrix.IA[0] = 0;
    for (i = 0; i < nNode; i++) {
        degree = crs_IA_Node2Node[i+1] - crs_IA_Node2Node[i];
        SolverBlockMatrix.DIM += degree + 1;
        SolverBlockMatrix.IA[i+1] = SolverBlockMatrix.IA[i] + degree + 1;
    }
    
    // Allocate Memory to IAU Array to store Diagonal Location
    SolverBlockMatrix.IAU = new int[nNode];
    
    // Allocate Memory to JA Array to store location of Non Zero Entries
    SolverBlockMatrix.JA = new int[SolverBlockMatrix.DIM];
    // Get the values for JA
    for (i = 0; i < nNode; i++) {
        // Make the first start entry as main node id
        index = SolverBlockMatrix.IA[i];
        SolverBlockMatrix.JA[index] = i;
        
        // Now add the nodes connected to main node
        for (j = crs_IA_Node2Node[i]; j < crs_IA_Node2Node[i+1]; j++) {
            index++;
            SolverBlockMatrix.JA[index] = crs_JA_Node2Node[j];
        }
    }
    
    /* Now Sort JA and Find IAU */
    // This step is necessary for computation of Transpose in Design Code - Adjoint
    for (i = 0; i < nNode; i++) {
        jstart = SolverBlockMatrix.IA[i];
        jend   = SolverBlockMatrix.IA[i + 1];
        for (j = jstart; j < jend; j++) {
            min = SolverBlockMatrix.JA[j];
            minsave = SolverBlockMatrix.JA[j];
            ksave = j;
            for (k = j + 1; k < jend; k++) {
                if (SolverBlockMatrix.JA[k] < min) {
                    min = SolverBlockMatrix.JA[k];
                    ksave = k;
                }
            }
            SolverBlockMatrix.JA[j] = min;
            SolverBlockMatrix.JA[ksave] = minsave;
            if (SolverBlockMatrix.JA[j] == i)
                SolverBlockMatrix.IAU[i] = j;
        }
    }
    
    // Allocate Memory for CRS Matrix
    SolverBlockMatrix.A = (double ***) malloc (SolverBlockMatrix.DIM*sizeof(double**));
    for (i = 0; i < SolverBlockMatrix.DIM; i++) {
        SolverBlockMatrix.A[i] = NULL;
        SolverBlockMatrix.A[i] = (double **) malloc (SolverBlockMatrix.Block_nRow*sizeof(double*));
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
            SolverBlockMatrix.A[i][j] = NULL;
            SolverBlockMatrix.A[i][j] = (double *) malloc (SolverBlockMatrix.Block_nCol*sizeof(double));
            for (k = 0; k < SolverBlockMatrix.Block_nCol; k++)
                SolverBlockMatrix.A[i][j][k] = 0.0;
        }
    }
    
    // Allocate Memory of RHS
    SolverBlockMatrix.B = (double **) malloc (SolverBlockMatrix.nROW*sizeof(double*));
    for (i = 0; i < SolverBlockMatrix.nROW; i++) {
        SolverBlockMatrix.B[i] = NULL;
        SolverBlockMatrix.B[i] = (double *) malloc (SolverBlockMatrix.Block_nRow*sizeof(double));
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++)
            SolverBlockMatrix.B[i][j] = 0.0;
    }

    // Allocate Memory for X
    SolverBlockMatrix.X = (double **) malloc (SolverBlockMatrix.nROW*sizeof(double*));
    for (i = 0; i < SolverBlockMatrix.nROW; i++) {
        SolverBlockMatrix.X[i] = NULL;
        SolverBlockMatrix.X[i] = (double *) malloc (SolverBlockMatrix.Block_nRow*sizeof(double));
        for (j = 0; j < SolverBlockMatrix.Block_nRow; j++)
            SolverBlockMatrix.X[i][j] = 0.0;
    }
}

//------------------------------------------------------------------------------
//! 
//------------------------------------------------------------------------------
void Delete_CRS_SolverBlockMatrix(void) {
    int i, j;
    
    if (SolverBlockMatrix.A != NULL) {
        for (i = 0; i < SolverBlockMatrix.DIM; i++) {
            if (SolverBlockMatrix.A[i] != NULL) {
                for (j = 0; j < SolverBlockMatrix.Block_nRow; j++) {
                    if (SolverBlockMatrix.A[i][j] != NULL)
                        free(SolverBlockMatrix.A[i][j]);
                }
                free(SolverBlockMatrix.A[i]);
            }
        }
        free(SolverBlockMatrix.A);
    }
    
    if (SolverBlockMatrix.B != NULL) {
        for (i = 0; i < SolverBlockMatrix.nROW; i++) {
            if (SolverBlockMatrix.B[i] != NULL)
                free(SolverBlockMatrix.B[i]);
        }
        free(SolverBlockMatrix.B);
    }
    
    if (SolverBlockMatrix.X != NULL) {
        for (i = 0; i < SolverBlockMatrix.nROW; i++) {
            if (SolverBlockMatrix.X[i] != NULL)
                free(SolverBlockMatrix.X[i]);
        }
        free(SolverBlockMatrix.X);
    }
    
    if (SolverBlockMatrix.IA != NULL)
        delete[] SolverBlockMatrix.IA;
    
    if (SolverBlockMatrix.IAU != NULL)
        delete[] SolverBlockMatrix.IAU;
    
    if (SolverBlockMatrix.JA != NULL)
        delete[] SolverBlockMatrix.JA;
}

