/*******************************************************************************
 * File:        Solver.h
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

#ifndef SOLVER_H
#define	SOLVER_H

// Boundary Conditions Type
#define BC_NONE                -1
#define BC_SOLID_WALL           0
#define BC_FREE_STREAM          1
#define BC_SUBSONIC_INFLOW      2
#define BC_SUBSONIC_OUTFLOW     3
#define BC_SUPERSONIC_INFLOW    4
#define BC_SUPERSONIC_OUTFLOW   5

// Solver Type
#define SOLVER_NONE             0
#define SOLVER_ROE              1
#define SOLVER_LMROE            2

// Variable Type
#define VARIABLE_NONE           0
#define VARIABLE_CONSERVATIVE   1
#define VARIABLE_RUP            2
#define VARIABLE_RUT            3
#define VARIABLE_PUT            4

// Linear Solver Parameters
extern int    SolverScheme;
extern int    TimeAccuracy;
extern int    TimeStepScheme;
extern int    Order;
extern int    NIteration;
extern int    InnerNIteration;
extern int    FirstOrderNIteration;
extern double Relaxation;

// Flux Limiters
extern int    Limiter;
extern int    LimiterSmooth;
extern int    LimiterOrder;
extern int    StartLimiterNIteration;
extern int    EndLimiterNIteration;
extern double Venkat_KThreshold;

// Entropy Fix
extern int EntropyFix;

// Low Mach Fix
extern int LMRoeFix;

// Restart Parameters
extern int  RestartInput;
extern int  RestartOutput;
extern int  RestartIteration;
extern int  RestartCycle;
extern char RestartInputFilename[256];
extern char RestartOutputFilename[256];

// Reference Conditions
extern double Ref_Rho;
extern double Ref_Mach;
extern double Ref_Alpha;
extern double Ref_Pressure;
extern double Ref_Temperature;

// Free Stream Conditions
extern double Inf_Rho;
extern double Inf_U;
extern double Inf_V;
extern double Inf_W;
extern double Inf_Et;
extern double Inf_Pressure;
extern double Inf_Mach;

// Constants
extern double Gamma;

// Solver Tuning Parameters
// CFL Conditions
extern int    CFL_Ramp;
extern double CFL_MAX;
extern double CFL_MIN;
extern double CFL;

// Mach Ramping
extern int    Mach_Ramp;
extern double Mach_MAX;
extern double Mach_MIN;

// Zero Pressure Gradient No of Iterations
extern int    ZPGIteration;

// Conservative Variables
extern double *Q1;
extern double *Q2;
extern double *Q3;
extern double *Q4;
extern double *Q5;

// Conservative Variables Gradients
extern double *Q1x;
extern double *Q1y;
extern double *Q1z;

extern double *Q2x;
extern double *Q2y;
extern double *Q2z;

extern double *Q3x;
extern double *Q3y;
extern double *Q3z;

extern double *Q4x;
extern double *Q4y;
extern double *Q4z;

extern double *Q5x;
extern double *Q5y;
extern double *Q5z;

// Residuals
extern double *Res1;
extern double *Res2;
extern double *Res3;
extern double *Res4;
extern double *Res5;

// Local Time
extern double *DeltaT;

// Limiters
extern double *Limiter_Phi1;
extern double *Limiter_Phi2;
extern double *Limiter_Phi3;
extern double *Limiter_Phi4;
extern double *Limiter_Phi5;

// RMS
extern double RMS[5];
extern double RMS_Res;

// Compute Free Stream Condition with Mach Ramping
void ComputeFreeStreamCondition(int Iteration);

// Initialize the Solver Data Structure
void Solver_Init(void);
// Finalize the Solver Data Structure
void Solver_Finalize(void);
// Get Solver Parameters
void Solver_Read_Params(const char *paramfile);
void Solver_Set_Initial_Conditions(void);
int  Solve(void);

void Initialize_Boundary_Condition(void);
void Apply_Boundary_Condition(int Iteration);
void Compute_DeltaT(int Iteration);

// Roe Scheme Functions
void Compute_RoeVariables(double *Q_L, double *Q_R, double *Q_Roe);
void Compute_Residual_Roe(void);

// Limiter
void Compute_Limiter(void);
void Compute_Limiter_Barth_Jespersen(int node_C, double *Phi_C);
void Compute_Limiter_Venkatakrishnan(int node_C, double *Phi_C);
void Compute_Limiter_PressureCorrection(int node_C, double *Phi_C);
void Compute_Limiter_RoePressureCorrection(int node_L, int node_R, double *Phi_L, double *Phi_R);

// Higher Order Reconstruction of Q's
void Compute_SecondOrderReconstructQ(int node_L, int node_R, double *Q_L, double *Q_R);

// Entropy Fix
void Roe_EntropyFix(double ubar_L, double c_L, double ubar_R, double c_R, double ubar, double c, double **Eigen);

#endif	/* SOLVER_H */

