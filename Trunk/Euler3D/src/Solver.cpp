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
#include "Residual_Smoothing.h"
#include "Material.h"
#include "DebugSolver.h"

int SolverIteration;

// Local Time
double *DeltaT;
double MinDeltaT;
double MaxDeltaT;

// Conservative or Primitive Variables
double *Q1;
double *Q2;
double *Q3;
double *Q4;
double *Q5;

// Conservative or Primitive Variables Gradients
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

// Min and Max EigenValue Value
double MinEigenLamda1;
double MaxEigenLamda1;
double MinEigenLamda4;
double MaxEigenLamda4;
double MinEigenLamda5;
double MaxEigenLamda5;

// Min and Max Precondition Variable
double *PrecondSigma;
double MinPrecondSigma;
double MaxPrecondSigma;

// MC CRS Matrix
MC_CRS SolverBlockMatrix;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Init(void) {
    // Local Time
    DeltaT          = NULL;
    MinDeltaT       = DBL_MAX;
    MaxDeltaT       = DBL_MIN;

    // Conservative or Primitive Variables
    Q1              = NULL;
    Q2              = NULL;
    Q3              = NULL;
    Q4              = NULL;
    Q5              = NULL;

    // Conservative or Primitive Variables Gradients
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
    
    // Min and Max EigenValue Value
    MinEigenLamda1   = DBL_MAX;
    MaxEigenLamda1   = DBL_MIN;
    MinEigenLamda4   = DBL_MAX;
    MaxEigenLamda4   = DBL_MIN;
    MinEigenLamda5   = DBL_MAX;
    MaxEigenLamda5   = DBL_MIN;
    
    // Min and Max Precondition Variable
    PrecondSigma    = NULL;
    MinPrecondSigma = DBL_MAX;
    MaxPrecondSigma = DBL_MIN;
    
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
    
    // RMS Writer Initialization
    RMS_Writer_Init();
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Finalize(void) {
    // Conservative or Primitive Variables
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

     // Conservative or Primitive Variables Gradient
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
    
    // Precondition Variable
    if (PrecondSigma != NULL)
        delete[] PrecondSigma;
    PrecondSigma = NULL;
    
    // Finalize Gradient Infrastructure
    if (Order == SOLVER_ORDER_SECOND)
        Gradient_Finalize();
    
    // Check if Implicit Method
    if ((SolverMethod == SOLVER_METHOD_IMPLICIT) || (SolverMethod == SOLVER_METHOD_IMPLICIT_UNSTEADY))
        Delete_CRS_SolverBlockMatrix();
    
    // RMS Writer Finalize
    RMS_Writer_Finalize();
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//! Set Initial Conditions
//------------------------------------------------------------------------------
void Solver_Set_Initial_Conditions(void) {

    // Allocate Memory to Store Conservative or Primitive Variables
    Q1 = new double[nNode + nBNode];
    Q2 = new double[nNode + nBNode];
    Q3 = new double[nNode + nBNode];
    Q4 = new double[nNode + nBNode];
    Q5 = new double[nNode + nBNode];

    // Compute Free Stream Conditions
    ComputeFreeStreamCondition(0);
    
    // Initialize the variable with reference conditions
    switch (Variable_Type) {
        case VARIABLE_CONSERVATIVE:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Rho;
                Q2[i] = Inf_Rho*Inf_U;
                Q3[i] = Inf_Rho*Inf_V;
                Q4[i] = Inf_Rho*Inf_W;
                Q5[i] = Inf_Rho*Inf_Et;
            }
            break;
        case VARIABLE_PRIMITIVE_PUT:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Pressure;
                Q2[i] = Inf_U;
                Q3[i] = Inf_V;
                Q4[i] = Inf_W;
                Q5[i] = Inf_Temperature;
            }
            break;
        case VARIABLE_PRIMITIVE_RUP:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Rho;
                Q2[i] = Inf_U;
                Q3[i] = Inf_V;
                Q4[i] = Inf_W;
                Q5[i] = Inf_Pressure;
            }
            break;
         case VARIABLE_PRIMITIVE_RUT:
            for (int i = 0; i < (nNode + nBNode); i++) {
                Q1[i] = Inf_Rho;
                Q2[i] = Inf_U;
                Q3[i] = Inf_V;
                Q4[i] = Inf_W;
                Q5[i] = Inf_Temperature;
            }
            break;
        default:
            error("Solver_Set_Initial_Conditions: Undefined Variable Type - %d", Variable_Type);
            break;
    }
    

    // Allocate Memory to Store Conservative or Primitive Variables Gradients
    if (Order == SOLVER_ORDER_SECOND) {
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
    
    // Min and Max EigenValue Value
    MinEigenLamda1   = DBL_MAX;
    MaxEigenLamda1   = DBL_MIN;
    MinEigenLamda4   = DBL_MAX;
    MaxEigenLamda4   = DBL_MIN;
    MinEigenLamda5   = DBL_MAX;
    MaxEigenLamda5   = DBL_MIN;
    
    // Min and Max Precondition Variable
    if (PrecondMethod != SOLVER_PRECOND_NONE) {
        PrecondSigma    = new double[nNode];
        for (int i = 0; i < nNode; i++)
            PrecondSigma[i] = 0.0;
        MinPrecondSigma = DBL_MAX;
        MaxPrecondSigma = DBL_MIN;
    } else {
        MinPrecondSigma = 1.0;
        MaxPrecondSigma = 1.0;
    }
    
    // Allocate Memory to Store Time Step
    DeltaT = new double[nNode];
    // Initialize
    for (int i = 0; i < nNode; i++)
        DeltaT[i] = 0.0;
    MinDeltaT = DBL_MAX;
    MaxDeltaT = DBL_MIN;
    
    // Set Boundary Conditions
    Initialize_Boundary_Condition();
    
    // Check if Implicit Method
    if ((SolverMethod == SOLVER_METHOD_IMPLICIT) || (SolverMethod == SOLVER_METHOD_IMPLICIT_UNSTEADY))
        Create_CRS_SolverBlockMatrix();
}

//------------------------------------------------------------------------------
//! Solver in Explicit Steady State Mode
//------------------------------------------------------------------------------
int Solve_Explicit(void) {
    int AddTime   = FALSE;
    int CheckNAN  = 0;
    int SaveOrder = 0;
    int SaveLimiter = 0;
    double *W01, *W02, *W03, *W04, *W05;
    double *WDeltaT, *Dummy;
    double phi1, phi2, phi3, phi4;
    double eta, scale_mach;
    
    // Check if Euler or Runge-Kutta Scheme
    W01 = W02 = W03 = W04 = W05 = NULL;
    WDeltaT = NULL;
    Dummy   = NULL;
    if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_RK4) {
        // Wo
        W01 = new double[nNode];
        W02 = new double[nNode];
        W03 = new double[nNode];
        W04 = new double[nNode];
        W05 = new double[nNode];
        
        // WDeltaT : Needed because DeltaT computations is inside Residual Computation
        WDeltaT = new double[nNode];
    }
    
    // Helper Variables
    eta  = 0.0;
    scale_mach = 0.0;
    
    // RK4 Coefficients
    phi1 = 0.25;
    phi2 = 0.333333333333333333;
    phi3 = 0.5;
    phi4 = 1.0;
    
    // Save Solver Parameters
    SaveOrder   = Order;
    SaveLimiter = Limiter;
    
    for (int iter = RestartIteration; iter < NIteration; iter++) {
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
        if (PrecondMethod != SOLVER_PRECOND_NONE) {
            for (int i = 0; i < nNode; i++)
                PrecondSigma[i] = 0.0;
            MinPrecondSigma = DBL_MAX;
            MaxPrecondSigma = DBL_MIN;
        }
        
        // Compute Free Stream Conditions with Mach Ramping
        ComputeFreeStreamCondition(iter);
        
        // Check if First Order Iterations are Required
        Order = SaveOrder;
        if (FirstOrderNIteration > iter)
            Order = SOLVER_ORDER_FIRST;

        // Set the Start and End of Limiter Iterations
        Limiter = 0;
        if (StartLimiterNIteration < EndLimiterNIteration) {
            if ((StartLimiterNIteration <= iter+1) && (EndLimiterNIteration > iter))
                Limiter = SaveLimiter;
        }

        // Compute Least Square Gradient -- Unweighted
        if (Order == SOLVER_ORDER_SECOND) {
            Compute_Least_Square_Gradient(0);
            if (Limiter > 0)
                Compute_Limiter();
        }
        
        // Apply boundary conditions
        Apply_Boundary_Condition(iter);

        // Reset the Residual Smoothing Variables: Required before Residual Computation
        if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE)
            Residual_Smoothing_Reset();
        
        // Compute Residuals
        AddTime = TRUE;
        Compute_Residual(AddTime);
        
        // Compute Local Time Stepping
        Compute_DeltaT(iter);

        // Smooth the Residual: DeltaT Computation is required
        if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE)
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
        
        // Normalize the RMS with Mach Number
        if (Ref_Mach > 1.0)
            scale_mach = Ref_Mach;
        else
            scale_mach = 1.0/Ref_Mach;
        
        if (PrecondMethod != SOLVER_PRECOND_NONE && PrecondMethod != SOLVER_PRECOND_ROE_LMFIX) {
            scale_mach = scale_mach*scale_mach;
        }
        RMS_Res *= scale_mach;
        RMS[0]  *= scale_mach;
        RMS[1]  *= scale_mach;
        RMS[2]  *= scale_mach;
        RMS[3]  *= scale_mach;
        RMS[4]  *= scale_mach;
        
        // Write RMS
        RMS_Writer(iter+1, RMS);
        
        // Check for Residual NAN
        if (isnan(RMS_Res)) {
            info("Solve_Explicit: NAN Encountered ! - Abort");
            iter = NIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }

        // Euler Time Stepping Scheme
        if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_EULER) {
            // Update Conservative or Primitive Variables
            for (int i = 0; i < nNode; i++) {
                eta    = DeltaT[i]/cVolume[i];
                Q1[i] -= eta * Res1[i];
                Q2[i] -= eta * Res2[i];
                Q3[i] -= eta * Res3[i];
                Q4[i] -= eta * Res4[i];
                Q5[i] -= eta * Res5[i]; 
            }
        }

        // Runge-Kutta 4 Time Stepping Scheme
        if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_RK4) {
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
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
                    Compute_Limiter();
            }

            // Apply boundary conditions
            Apply_Boundary_Condition(iter);

            // Reset the Residual Smoothing Variables: Required before Residual Computation
            if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE)
                Residual_Smoothing_Reset();
            
            // Compute Residuals
            AddTime = FALSE;
            Compute_Residual(AddTime);
            
            // Smooth the Residual: DeltaT Computation is required
            if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE) {
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
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
                    Compute_Limiter();
            }

            // Apply boundary conditions
            Apply_Boundary_Condition(iter);
            
            // Reset the Residual Smoothing Variables: Required before Residual Computation
            if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE)
                Residual_Smoothing_Reset();
            
            // Compute Residuals
            AddTime = FALSE;
            Compute_Residual(AddTime);
            
            // Smooth the Residual: DeltaT Computation is required
            if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE) {
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
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
                    Compute_Limiter();
            }

            // Apply boundary conditions
            Apply_Boundary_Condition(iter);

            // Reset the Residual Smoothing Variables: Required before Residual Computation
            if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE)
                Residual_Smoothing_Reset();
            
            // Compute Residuals
            AddTime = FALSE;
            Compute_Residual(AddTime);
            
            // Smooth the Residual: DeltaT Computation is required
            if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE) {
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
        }

        // Precondition Variable Normalize
        if (PrecondMethod != SOLVER_PRECOND_NONE) {
            for (int i = 0; i < nNode; i++)
                PrecondSigma[i] /= (crs_IA_Node2Node[i+1] - crs_IA_Node2Node[i]);
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
    if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_RK4) {
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

//------------------------------------------------------------------------------
//! Solver in Explicit Unsteady Mode
//------------------------------------------------------------------------------
int Solve_Explicit_Unsteady(void) {
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
    if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_RK4) {
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
    SaveOrder   = Order;
    SaveLimiter = Limiter;
    
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
    for (int t_iter = 0; t_iter < NIteration; t_iter++) {
        
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
            if (PrecondMethod != SOLVER_PRECOND_NONE) {
                MinPrecondSigma = DBL_MAX;
                MaxPrecondSigma = DBL_MIN;
            }

            // Compute Free Stream Conditions with Mach Ramping
            ComputeFreeStreamCondition(iter);

            // Check if First Order Iterations are Required
            Order = SaveOrder;
            if (FirstOrderNIteration > iter)
                Order = SOLVER_ORDER_FIRST;

            // Set the Start and End of Limiter Iterations
            Limiter = 0;
            if (StartLimiterNIteration < EndLimiterNIteration) {
                if ((StartLimiterNIteration <= iter+1) && (EndLimiterNIteration > iter))
                    Limiter = SaveLimiter;
            }

            // Compute Least Square Gradient -- Unweighted
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
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

            // Euler Time Stepping Scheme
            if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_EULER) {
                // Update Conservative or Primitive Variables
                for (int i = 0; i < nNode; i++) {
                    eta    = DeltaT[i]/cVolume[i];
                    reta   = DeltaT[i]/PhysicalDeltaTime;
                    zeta   = (1.0 + 1.5*reta);
                    Q1[i] -= (eta * Res1[i] + 0.5*reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] -= (eta * Res2[i] + 0.5*reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] -= (eta * Res3[i] + 0.5*reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] -= (eta * Res4[i] + 0.5*reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] -= (eta * Res5[i] + 0.5*reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;
                }
            }

            // Runge-Kutta 4 Time Stepping Scheme
            if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_RK4) {
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
                    zeta       = (1.0 + 1.5*reta);
                    
                    // W1 = W0 - (phi1*dtau/zeta)*(RES0/vol + (0.5/dT)*(alpha*W0 - beta*Qn0 + theta*Qn1))
                    Q1[i] = W01[i] - phi1*(eta*Res1[i] + 0.5*reta*(alpha*W01[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi1*(eta*Res2[i] + 0.5*reta*(alpha*W02[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi1*(eta*Res3[i] + 0.5*reta*(alpha*W03[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi1*(eta*Res4[i] + 0.5*reta*(alpha*W04[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi1*(eta*Res5[i] + 0.5*reta*(alpha*W05[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;

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
                if (Order == SOLVER_ORDER_SECOND) {
                    Compute_Least_Square_Gradient(0);
                    if (Limiter > 0)
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
                    zeta  = (1.0 + 1.5*reta);
                    Q1[i] = W01[i] - phi2*(eta*Res1[i] + 0.5*reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi2*(eta*Res2[i] + 0.5*reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi2*(eta*Res3[i] + 0.5*reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi2*(eta*Res4[i] + 0.5*reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi2*(eta*Res5[i] + 0.5*reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;

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
                if (Order == SOLVER_ORDER_SECOND) {
                    Compute_Least_Square_Gradient(0);
                    if (Limiter > 0)
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
                    zeta  = (1.0 + 1.5*reta);
                    Q1[i] = W01[i] - phi3*(eta*Res1[i] + 0.5*reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi3*(eta*Res2[i] + 0.5*reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi3*(eta*Res3[i] + 0.5*reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi3*(eta*Res4[i] + 0.5*reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi3*(eta*Res5[i] + 0.5*reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;

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
                if (Order == SOLVER_ORDER_SECOND) {
                    Compute_Least_Square_Gradient(0);
                    if (Limiter > 0)
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
                    zeta  = (1.0 + 1.5*reta);
                    Q1[i] = W01[i] - phi4*(eta*Res1[i] + 0.5*reta*(alpha*Q1[i] - beta*Qn01[i] + theta*Qn11[i]))/zeta;
                    Q2[i] = W02[i] - phi4*(eta*Res2[i] + 0.5*reta*(alpha*Q2[i] - beta*Qn02[i] + theta*Qn12[i]))/zeta;
                    Q3[i] = W03[i] - phi4*(eta*Res3[i] + 0.5*reta*(alpha*Q3[i] - beta*Qn03[i] + theta*Qn13[i]))/zeta;
                    Q4[i] = W04[i] - phi4*(eta*Res4[i] + 0.5*reta*(alpha*Q4[i] - beta*Qn04[i] + theta*Qn14[i]))/zeta;
                    Q5[i] = W05[i] - phi4*(eta*Res5[i] + 0.5*reta*(alpha*Q5[i] - beta*Qn05[i] + theta*Qn15[i]))/zeta;
                }
            }
        }
        
        // Check for NAN
        if (CheckNAN == 1) {
            t_iter = NIteration + 1;
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
    
    if (TimeStepScheme == SOLVER_EXPLICIT_TIME_STEP_SCHEME_RK4) {
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

//------------------------------------------------------------------------------
//! Solver in Implicit Steady State Mode
//------------------------------------------------------------------------------
int Solve_Implicit(void) {
    int AddTime   = FALSE;
    int CheckNAN  = 0;
    int SaveOrder = 0;
    int SaveLimiter = 0;
    double lrms;
    
    // Save Solver Parameters
    SaveOrder   = Order;
    SaveLimiter = Limiter;

    // Pseduo-Time Loop
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
        if (PrecondMethod != SOLVER_PRECOND_NONE) {
            MinPrecondSigma = DBL_MAX;
            MaxPrecondSigma = DBL_MIN;
        }
        
        // Compute Free Stream Conditions with Mach Ramping
        ComputeFreeStreamCondition(iter);
    
        // Check if First Order Iterations are Required
        Order = SaveOrder;
        if (FirstOrderNIteration > iter)
            Order = SOLVER_ORDER_FIRST;

        // Set the Start and End of Limiter Iterations
        Limiter = 0;
        if (StartLimiterNIteration < EndLimiterNIteration) {
            if ((StartLimiterNIteration <= iter+1) && (EndLimiterNIteration > iter))
                Limiter = SaveLimiter;
        }

        // Compute Least Square Gradient -- Unweighted
        if (Order == SOLVER_ORDER_SECOND) {
            Compute_Least_Square_Gradient(0);
            if (Limiter > 0)
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
        RMS_Writer(iter+1, RMS);
        
        // Check for Residual NAN
        if (isnan(RMS_Res)) {
            info("Solve_Implicit: NAN Encountered ! - Abort");
            iter = NIteration + 1;
            CheckNAN = 1;
            VTK_Writer("SolutionBeforeNAN.vtk", 1);;
        }

        // Compute Jacobian and Fill CRS Matrix
        AddTime = TRUE;
        Compute_Jacobian(AddTime, iter);
        
        // Solve for Solution
        lrms = MC_Iterative_Block_LU_Jacobi_CRS(LinearSolverNIteration, 0, SolverBlockMatrix);
        
        // Check for Linear Solver NAN
        if (isnan(lrms)) {
            info("Solve_Implicit: Liner Solver Iteration: NAN Encountered ! - Abort");
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
//! Solver in Implicit Unsteady Mode
//------------------------------------------------------------------------------
int Solve_Implicit_Unsteady(void) {
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
    if (TimeStepScheme == SOLVER_IMPLICIT_TIME_STEP_SCHEME_BDF) {
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
    SaveOrder   = Order;
    SaveLimiter = Limiter;
    
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
    for (int t_iter = 0; t_iter < NIteration; t_iter++) {
        
        // Check if Second Order in Time Requested
        if (TimeStepScheme == SOLVER_IMPLICIT_TIME_STEP_SCHEME_BDF)
            if (t_iter > 0)
                theta = 0.5;
        
        // Set Qn and Qn-1
        for (int i = 0; i < nNode; i++) {
            // Check if Second Order in Time Requested
            if (TimeStepScheme == SOLVER_IMPLICIT_TIME_STEP_SCHEME_BDF) {
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
            if (PrecondMethod != SOLVER_PRECOND_NONE) {
                MinPrecondSigma = DBL_MAX;
                MaxPrecondSigma = DBL_MIN;
            }

            // Compute Free Stream Conditions with Mach Ramping
            ComputeFreeStreamCondition(iter);

            // Check if First Order Iterations are Required
            Order = SaveOrder;
            if (FirstOrderNIteration > iter)
                Order = SOLVER_ORDER_FIRST;

            // Set the Start and End of Limiter Iterations
            Limiter = 0;
            if (StartLimiterNIteration < EndLimiterNIteration) {
                if ((StartLimiterNIteration <= iter+1) && (EndLimiterNIteration > iter))
                    Limiter = SaveLimiter;
            }

            // Compute Least Square Gradient -- Unweighted
            if (Order == SOLVER_ORDER_SECOND) {
                Compute_Least_Square_Gradient(0);
                if (Limiter > 0)
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
                if (TimeStepScheme == SOLVER_IMPLICIT_TIME_STEP_SCHEME_BDF) {
                    SolverBlockMatrix.B[i][0] -= eta*((1.0 + theta)*(Q1[i] - Qn01[i]) - theta*(Qn01[i] - Qn11[i]));
                    SolverBlockMatrix.B[i][1] -= eta*((1.0 + theta)*(Q2[i] - Qn02[i]) - theta*(Qn02[i] - Qn12[i]));
                    SolverBlockMatrix.B[i][2] -= eta*((1.0 + theta)*(Q3[i] - Qn03[i]) - theta*(Qn03[i] - Qn13[i]));
                    SolverBlockMatrix.B[i][3] -= eta*((1.0 + theta)*(Q4[i] - Qn04[i]) - theta*(Qn04[i] - Qn14[i]));
                    SolverBlockMatrix.B[i][4] -= eta*((1.0 + theta)*(Q5[i] - Qn05[i]) - theta*(Qn05[i] - Qn15[i]));
                } else if (TimeStepScheme == SOLVER_IMPLICIT_TIME_STEP_SCHEME_EULER) {
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
            t_iter = NIteration + 1;
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
    if (TimeStepScheme == SOLVER_IMPLICIT_TIME_STEP_SCHEME_BDF) {
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
        case SOLVER_SCHEME_HLLC: // HLLC
            HLLC_Init();
            break;
        case SOLVER_SCHEME_AUSM: // AUSM
            AUSM_Init();
            break;
        case SOLVER_SCHEME_VANLEER: // Van Leer
            VanLeer_Init();
            break;
        case SOLVER_SCHEME_LDFSS: // LDFSS
            LDFSS_Init();
            break;
        case SOLVER_SCHEME_OSHER: // Osher
            Osher_Init();
            break;
        case SOLVER_SCHEME_STEGERWARMING: // Steger Warming
            StegerWarming_Init();
            break;
        default:
            error("Solve: Invalid Solver Scheme - %d", SolverScheme);
            break;
    }
    
    // Check if Solution Restart is Requested
    if (RestartInput)
        Restart_Reader(RestartInputFilename);
    
    // Check if Residual Smoothing is Requested
    if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE)
        Residual_Smoothing_Init();
    
    // Select the Solver Method Type
    switch (SolverMethod) {
        case SOLVER_METHOD_EXPLICIT:
            rvalue = Solve_Explicit();
            break;
        case SOLVER_METHOD_IMPLICIT:
            rvalue = Solve_Implicit();
            break;
        case SOLVER_METHOD_EXPLICIT_UNSTEADY:
            rvalue = Solve_Explicit_Unsteady();
            break;
        case SOLVER_METHOD_IMPLICIT_UNSTEADY:
            rvalue = Solve_Implicit_Unsteady();
            break;
        default:
            error("Solve: Invalid Solver Method - %d", SolverMethod);
            break;
    }
    
    // Finalize Residual Smoothing Data Structure
    if (ResidualSmoothType != RESIDUAL_SMOOTH_NONE)
        Residual_Smoothing_Finalize();
    
    // Finalize the Solver Scheme Data Structure
    switch (SolverScheme) {
        case SOLVER_SCHEME_ROE: // Roe
            Roe_Finalize();
            break;
        case SOLVER_SCHEME_HLLC: // HLLC
            HLLC_Finalize();
            break;
        case SOLVER_SCHEME_AUSM: // AUSM
            AUSM_Finalize();
            break;
        case SOLVER_SCHEME_VANLEER: // Van Leer
            VanLeer_Finalize();
            break;
        case SOLVER_SCHEME_LDFSS: // LDFSS
            LDFSS_Finalize();
            break;
        case SOLVER_SCHEME_OSHER: // Osher
            Osher_Finalize();
            break;
        case SOLVER_SCHEME_STEGERWARMING: // Steger Warming
            StegerWarming_Finalize();
            break;
        default:
            error("Solve: Invalid Solver Scheme - %d", SolverScheme);
            break;
    }
    
    return rvalue;
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

