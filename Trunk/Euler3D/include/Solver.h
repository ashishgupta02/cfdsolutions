/*******************************************************************************
 * File:        Solver.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _SOLVER_H
#define	_SOLVER_H

// Custom header files
#include "SolverParameters.h"
#include "MC.h"

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

// Residuals (Convective and Dissipative)
extern double *Res1;
extern double *Res2;
extern double *Res3;
extern double *Res4;
extern double *Res5;

extern double *Res1_Diss;
extern double *Res2_Diss;
extern double *Res3_Diss;
extern double *Res4_Diss;
extern double *Res5_Diss;

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

// MC CRS Matrix
extern MC_CRS SolverBlockMatrix;

// Compute Free Stream Condition with Mach Ramping
void ComputeFreeStreamCondition(int Iteration);

// Initialize the Solver Data Structure
void Solver_Init(void);
// Finalize the Solver Data Structure
void Solver_Finalize(void);
void Solver_Set_Initial_Conditions(void);
int  Solve(void);
int  Solve_Explicit(void);
void Compute_Residual(void);

// Implicit Method Routines
void Create_CRS_SolverBlockMatrix(void);
void Delete_CRS_SolverBlockMatrix(void);
int  Solve_Implicit(void);
void Compute_Jacobian(int AddTime, int Iteration);

// Boundary Conditions
void Initialize_Boundary_Condition(void);
void Compute_Characteristic_BoundaryCondition(int BEdgeID, int Iteration);
int  Compute_Characteristic_BoundaryCondition_New(double Q_L[NEQUATIONS], double Q_R[NEQUATIONS], Vector3D AreaVec, int BCType, int Iteration, double Q_B[NEQUATIONS]);
void Compute_Pressure_BoundaryCondition(int BEdgeID, int Iteration);
void Apply_Characteristic_Boundary_Condition(int Iteration);
void Apply_Pressure_Boundary_Condition(int Iteration);
void Apply_Boundary_Condition(int Iteration);

// Time Stepping
void Compute_DeltaT(int Iteration);

// Roe Scheme Functions
void Roe_Init(void);
void Roe_Finalize(void);
void Roe_Reset(void);
void Compute_RoeVariables(double *Q_L, double *Q_R, double *Q_Roe);
void Compute_RoeAJacobian(double *Q_L, double *Q_R, Vector3D areavec, double **AJacobian_Roe);
void Compute_RoeFlux(int node_L, int node_R, Vector3D areavec, double *Flux_Roe);
void Compute_Residual_Roe(void);
void Compute_Jacobian_Exact_Roe(int AddTime, int Iteration);
void Compute_Jacobian_Approximate_Roe(int AddTime, int Iteration);
void Compute_Jacobian_FiniteDifference_Roe(int AddTime, int Iteration);
void Compute_Jacobian_Roe(int AddTime, int Iteration);

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

#endif	/* _SOLVER_H */

