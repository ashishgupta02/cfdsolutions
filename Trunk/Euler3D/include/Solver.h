/*******************************************************************************
 * File:        Solver.h
 * Author:      Ashish Gupta
 * Revision:    1
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
extern int    Order;
extern int    NIteration;
extern int    InnerNIteration;
extern double Relaxation;

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
void Solver_Set_Initial_Conditions();
void Solve();

void Initialize_Boundary_Condition();
void Apply_Boundary_Condition();
void Compute_Residual();
void Compute_DeltaT(int Iteration);

#endif	/* SOLVER_H */

