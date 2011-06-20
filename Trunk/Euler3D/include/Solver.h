/*******************************************************************************
 * File:        Solver.h
 * Author:      Ashish Gupta
 * Revision:    2
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

// Linear Solver Params
extern int    TimeStepScheme;
extern int    Order;
extern int    NIteration;
extern int    InnerNIteration;
extern int    FirstOrderNIteration;
extern double Relaxation;

// Flux Limiters
extern int    Limiter;
extern int    StartLimiterNIteration;
extern int    EndLimiterNIteration;
extern double Venkat_KThreshold;

// Entropy Fix
extern int EntropyFix;

// Restart Params
extern int  RestartInput;
extern int  RestartOutput;
extern int  RestartIteration;
extern char RestartInputFilename[256];
extern char RestartOutputFilename[256];

// Reference Conditions
extern double Ref_Mach;
extern double Ref_Alpha;
extern double Ref_Pressure;

// Free Stream Conditions
extern double Inf_Rho;
extern double Inf_U;
extern double Inf_V;
extern double Inf_W;
extern double Inf_Et;
extern double Inf_Pressure;

// Constants
extern double Gamma;

// CFL Conditons
extern int    CFL_Ramp;
extern double CFL_MAX;
extern double CFL_MIN;
extern double CFL;

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

// RMS
extern double RMS[5];
extern double RMS_Res;

// Initialize the Solver Data Structure
void Solver_Init(void);
// Finalize the Solver Data Structure
void Solver_Finalize(void);
// Get Solver Parameters
void Solver_Read_Params(const char *paramfile);
void Solver_Set_Initial_Conditions(void);
int  Solve(void);

void Initialize_Boundary_Condition(void);
void Apply_Boundary_Condition(void);
void Compute_Residual(void);
void Compute_Residual_LMRoe(void);
void Compute_DeltaT(int Iteration);

// Flux Limiter
void Compute_Flux_Limiter_Barth_Jespersen(int node_L, int node_R, double *Phi_L, double *Phi_R);
void Compute_Flux_Limiter_Venkatakrishnan(int node_L, int node_R, double *Phi_L, double *Phi_R);

// Entropy Fix
void Roe_EntropyFix(double ubar_L, double c_L, double ubar_R, double c_R, double ubar, double c, double **Eigen);

#endif	/* SOLVER_H */

