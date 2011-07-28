/*******************************************************************************
 * File:        Solver.h
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
#include "DebugSolver.h"

// Linear Solver Parameters
int    SolverMethod;
int    SolverScheme;
int    TimeAccuracy;
int    TimeStepScheme;
int    SolverBCScheme;
int    JacobianMethod;
int    JacobianUpdate;
int    Order;
int    NIteration;
int    InnerNIteration;
int    FirstOrderNIteration;
double Relaxation;

// Flux Limiter
int    Limiter;
int    LimiterSmooth;
int    LimiterOrder;
int    StartLimiterNIteration;
int    EndLimiterNIteration;
double Venkat_KThreshold;

// Entropy Fix
int EntropyFix;

// Low Mach Fix
int LMRoeFix;

// Restart Parameters
int  RestartInput;
int  RestartOutput;
int  RestartIteration;
int  RestartCycle;
char RestartInputFilename[256];
char RestartOutputFilename[256];

// Reference Conditions
double Ref_Rho;
double Ref_Mach;
double Ref_Alpha;
double Ref_Pressure;
double Ref_Temperature;

// Free Stream Conditions
double Inf_Rho;
double Inf_U;
double Inf_V;
double Inf_W;
double Inf_Et;
double Inf_Pressure;
double Inf_Mach;

// Constants
double Gamma;

// Solver Tuning Parameters
// CFL Conditions
int    CFL_Ramp;
double CFL_MAX;
double CFL_MIN;
double CFL;

// Mach Ramping
int    Mach_Ramp;
double Mach_MAX;
double Mach_MIN;

// Zero Pressure Gradient No of Iterations
int    ZPGIteration;

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

// Residuals
double *Res1;
double *Res2;
double *Res3;
double *Res4;
double *Res5;

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
    // Linear Solver Parameters
    SolverMethod    = SOLVER_METHOD_NONE;
    SolverScheme    = SOLVER_SCHEME_NONE;
    TimeAccuracy    = 0;
    TimeStepScheme  = 0;
    SolverBCScheme  = SOLVER_BC_SCHEME_NONE;
    JacobianMethod  = SOLVER_JACOBIAN_NONE;
    JacobianUpdate  = 0;
    Order           = 0;
    NIteration      = 0;
    InnerNIteration = 0;
    FirstOrderNIteration    = 0;
    Relaxation      = 0.0;

    // Flux Limiter
    Limiter         = 0;
    LimiterSmooth   = 0;
    LimiterOrder    = 0;
    StartLimiterNIteration  = 0;
    EndLimiterNIteration    = 0;
    Venkat_KThreshold       = 0.0;

    // Entropy Fix
    EntropyFix      = 0;
    
    // Low Mach Fix
    LMRoeFix        = 0;
    
    // Restart Parameters
    RestartInput    = 0;
    RestartOutput   = 0;
    RestartIteration= 0;
    RestartCycle    = 0;
    str_blank(RestartInputFilename);
    str_blank(RestartOutputFilename);

    // Reference Conditions
    Ref_Rho         = 0.0;
    Ref_Mach        = 0.0;
    Ref_Alpha       = 0.0;
    Ref_Pressure    = 0.0;
    Ref_Temperature = 0.0;

    // Free Stream Conditions
    Inf_Rho         = 0.0;
    Inf_U           = 0.0;
    Inf_V           = 0.0;
    Inf_W           = 0.0;
    Inf_Et          = 0.0;
    Inf_Pressure    = 0.0;
    Inf_Mach        = 0.0;
    
    // Constants
    Gamma           = 0.0;

    // CFL Conditions
    CFL_Ramp        = 0;
    CFL_MAX         = 0.0;
    CFL_MIN         = 0.0;
    CFL             = 0.0;

    // Mach Ramping
    Mach_Ramp       = 0;
    Mach_MAX        = 0.0;
    Mach_MIN        = 0.0;

    // Zero Pressure Gradient No of Iterations
    ZPGIteration    = 0;
    
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
    
    // Residuals
    Res1            = NULL;
    Res2            = NULL;
    Res3            = NULL;
    Res4            = NULL;
    Res5            = NULL;

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
    
    // Residuals
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
//!
//------------------------------------------------------------------------------
void Solver_Read_Params(const char *filename) {
    FILE *fp;
    int bdim = 256;
    char buff[256];
    char *dummy;

    if ((fp = fopen(filename, "r")) == NULL)
        error("Solver_Read_Params: Unable to Read Parameter File - %s", filename);
    
    // Number of Boundaries
    int nb;
    // Read Number of Boundaries
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &nb);

    // Allocate memory to store boundary type and read
    bndType = new int[nBC];
    for (int i = 0; i < nb; i++) {
        dummy = fgets(buff, bdim, fp);
        dummy = fgets(buff, bdim, fp);
        sscanf(buff, "%d", &bndType[i]);
    }

    // 1) Get the Solver Scheme
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &SolverMethod);
    
    // 2) Get the Solver Scheme
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &SolverScheme);

    // 3) Get the Time Accuracy
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &TimeAccuracy);

    // 4) Get the Time Stepping Scheme
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &TimeStepScheme);
    
    // 5) Get the Solver Scheme
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &SolverBCScheme);
    
    // 6) Get the Jacobian Method
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &JacobianMethod);
    
    // 7) Get the Jacobian Update Frequency
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &JacobianUpdate);
    
    // 8) Get the Order
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &Order);
    
    // 9) Get Number of Outer Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &NIteration);

    // 10) Get Number of Inner Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &InnerNIteration);

    // 11) Get Number of First Order Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &FirstOrderNIteration);
    
    // 12) Get the Relaxation Factor
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Relaxation);

    // 13) Get the Flux Limiter
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &Limiter);

    // 14) Get the Flux Limiter Smooth
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &LimiterSmooth);

    // 15) Get the Flux Limiter Order
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &LimiterOrder);
    
    // 16) Get the Start Limiter Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &StartLimiterNIteration);

    // 17) Get the End Limiter Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &EndLimiterNIteration);

    // 18) Read Ventakakrishanan K Threshold
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Venkat_KThreshold);

    // 19) Read Entropy Fix
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &EntropyFix);
    
    // 20) Read Gamma
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Gamma);

    // 21) Read Reference Density
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Rho);

    // 22) Read Reference Mach
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Mach);

    // 23) Read Reference Pressure
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Pressure);

    // 24) Read Alpha
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Alpha);

    // 25) Read CFL Ramp
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &CFL_Ramp);

    // 26) Read CFL MIN
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &CFL_MIN);

    // 27) Read CFL MAX
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &CFL_MAX);

    // 28) Read Mach Ramp
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &Mach_Ramp);

    // 29) Read Mach MIN
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Mach_MIN);

    // 30) Read Mach MAX
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Mach_MAX);

    // 31) Read No of Zero Pressure Gradient ZPG Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &ZPGIteration);

    // 32) Read Restart Solution
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &RestartInput);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &RestartOutput);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &RestartCycle);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%s", RestartInputFilename);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%s", RestartOutputFilename);
    
    CFL             = CFL_MIN;
    Ref_Alpha       = Ref_Alpha * M_PI / 180.0;
    Ref_Temperature = Gamma*Ref_Pressure/Ref_Rho;
    
    // Avoid Invalid Mach
    if (Ref_Mach < Mach_MAX)
        Mach_MAX = Ref_Mach;
    if (Mach_MIN == 0.0)
        Mach_MIN = Ref_Mach;
    
    // Close file
    fclose(fp);

    printf("=============================================================================\n");
    info("Input Solver Parameters");
    printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n");
    info("Solver Method ------------------------: %d",  SolverMethod);
    info("Solver Scheme ------------------------: %d",  SolverScheme);
    info("Time Accuracy ------------------------: %d",  TimeAccuracy);
    info("Time Stepping Scheme -----------------: %d",  TimeStepScheme);
    info("Solver Boundary Condition Scheme -----: %d",  SolverBCScheme);
    info("Solver Jacobian Method ---------------: %d",  JacobianMethod);
    info("Solver Jacobian Update Frequency -----: %d",  JacobianUpdate);
    info("Order --------------------------------: %d",  Order);
    info("No of Iterations ---------------------: %d",  NIteration);
    info("No of Inner Iterations ---------------: %d",  InnerNIteration);
    info("No of First Order Iterations ---------: %d",  FirstOrderNIteration);
    info("Relaxation Factor --------------------: %lf", Relaxation);
    info("Limiter Type -------------------------: %d",  Limiter);
    info("Limiter Smooth -----------------------: %d",  LimiterSmooth);
    info("Limiter Order ------------------------: %d",  LimiterOrder);
    info("Limiter Start Iteration --------------: %d",  StartLimiterNIteration);
    info("Limiter End Iteration ----------------: %d",  EndLimiterNIteration);
    info("Venkatakrishan Limiter Threshold -----: %lf", Venkat_KThreshold);
    info("Entropy Fix --------------------------: %d",  EntropyFix);
    info("Gamma --------------------------------: %lf", Gamma);
    info("Reference Density --------------------: %lf", Ref_Rho);
    info("Reference Mach Number ----------------: %lf", Ref_Mach);
    info("Reference Pressure -------------------: %lf", Ref_Pressure);
    info("Reference Alpha ----------------------: %lf", Ref_Alpha);
    info("CFL Ramp -----------------------------: %d",  CFL_Ramp);
    info("CFL_Min ------------------------------: %lf", CFL_MIN);
    info("CFL_Max ------------------------------: %lf", CFL_MAX);
    info("Mach Ramp ----------------------------: %d",  Mach_Ramp);
    info("Mach_Min -----------------------------: %lf", Mach_MIN);
    info("Mach_Max -----------------------------: %lf", Mach_MAX);
    info("No of ZPG Iterations -----------------: %d",  ZPGIteration);
    info("Restart Input ------------------------: %d",  RestartInput);
    info("Restart Input Filename ---------------: %s",  RestartInputFilename);
    info("Restart Output -----------------------: %d",  RestartOutput);
    info("Restart Output Filename --------------: %s",  RestartOutputFilename);
    info("Restart Cycle ------------------------: %d",  RestartCycle);
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
    
    // Allocate  Memory to Store Residuals
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
            Res1[i]   = 0.0;
            Res2[i]   = 0.0;
            Res3[i]   = 0.0;
            Res4[i]   = 0.0;
            Res5[i]   = 0.0;
            DeltaT[i] = 0.0;
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
                dtmp = DeltaT[i]/cVolume[i];
                Q1[i] -= dtmp * Res1[i];
                Q2[i] -= dtmp * Res2[i];
                Q3[i] -= dtmp * Res3[i];
                Q4[i] -= dtmp * Res4[i];
                Q5[i] -= dtmp * Res5[i];
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
            Res1[i]   = 0.0;
            Res2[i]   = 0.0;
            Res3[i]   = 0.0;
            Res4[i]   = 0.0;
            Res5[i]   = 0.0;
            DeltaT[i] = 0.0;
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
    Inf_Pressure = Ref_Pressure;
    Inf_U        = Inf_Mach*cos(Ref_Alpha);
    Inf_V        = Inf_Mach*sin(Ref_Alpha);
    Inf_W        = 0.0;
    Inf_Et       = Inf_Pressure/((Gamma - 1.0)*Inf_Rho) + 0.5 *(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W);
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

