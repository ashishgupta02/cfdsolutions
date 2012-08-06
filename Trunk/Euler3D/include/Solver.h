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

extern int SolverIteration;

// Conservative or Primitive Variables
extern double *Q1; /* rho,    P   */
extern double *Q2; /* rho*u,  u   */
extern double *Q3; /* rho*v,  v   */
extern double *Q4; /* rho*w,  w   */
extern double *Q5; /* rho*et, T   */

// Conservative or Primitive Variables Gradients
extern double *Q1x; /* rho,     P */
extern double *Q1y;
extern double *Q1z;

extern double *Q2x; /* rho*u,   u */
extern double *Q2y;
extern double *Q2z;

extern double *Q3x; /* rho*v,   v */
extern double *Q3y;
extern double *Q3z;

extern double *Q4x; /* rho*w,    w */
extern double *Q4y;
extern double *Q4z;

extern double *Q5x; /* rho*et,   T */
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
extern double MinDeltaT;
extern double MaxDeltaT;

// Limiters
extern double *Limiter_Phi1;
extern double *Limiter_Phi2;
extern double *Limiter_Phi3;
extern double *Limiter_Phi4;
extern double *Limiter_Phi5;

// RMS
extern double RMS[5];
extern double RMS_Res;

// Min and Max EigenValue Value
extern double MinEigenLamda1;
extern double MaxEigenLamda1;
extern double MinEigenLamda4;
extern double MaxEigenLamda4;
extern double MinEigenLamda5;
extern double MaxEigenLamda5;

// Min and Max Precondition Variable
extern double *PrecondSigma;
extern double MinPrecondSigma;
extern double MaxPrecondSigma;

// MC CRS Matrix
extern MC_CRS SolverBlockMatrix;

// Initialize the Solver Data Structure
void Solver_Init(void);
// Finalize the Solver Data Structure
void Solver_Finalize(void);
void Solver_Set_Initial_Conditions(void);
int  Solve(void);
int  Solve_Explicit(void);
int  Solve_Explicit_Unsteady(void);
int  Solve_Implicit(void);
int  Solve_Implicit_Unsteady(void);
void Compute_Residual(int AddTime);
void Unsteady_Initialization(void);

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
void Compute_RoeFlux(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_RoeFlux_Optimized(int node_L, int node_R, Vector3D areavec, double *Flux_Roe_Conv, double *Flux_Roe_Diss, int AddTime);
void Compute_Residual_Roe(int AddTime);
void Compute_Jacobian_Exact_Roe(int AddTime, int Iteration);
void Compute_Jacobian_Approximate_Roe(int AddTime, int Iteration);
void Compute_Jacobian_FiniteDifference_Roe(int AddTime, int Iteration);
void Compute_Jacobian_Roe(int AddTime, int Iteration);

// HLLC Scheme Functions
void HLLC_Init(void);
void HLLC_Finalize(void);
void HLLC_Reset(void);
void Compute_Residual_HLLC(int AddTime);

// AUSM Scheme Functions
void AUSM_Init(void);
void AUSM_Finalize(void);
void AUSM_Reset(void);
void Compute_Residual_AUSM(int AddTime);

// Van Leer Scheme Functions
void VanLeer_Init(void);
void VanLeer_Finalize(void);
void VanLeer_Reset(void);
void Compute_Residual_VanLeer(int AddTime);

// LDFSS Scheme Functions
void LDFSS_Init(void);
void LDFSS_Finalize(void);
void LDFSS_Reset(void);
void Compute_Residual_LDFSS(int AddTime);

// Osher Scheme Functions
void Osher_Init(void);
void Osher_Finalize(void);
void Osher_Reset(void);
void Compute_Residual_Osher(int AddTime);

// Steger Warming Scheme Functions
void StegerWarming_Init(void);
void StegerWarming_Finalize(void);
void StegerWarming_Reset(void);
void Compute_Residual_StegerWarming(int AddTime);

// JST Scheme Functions
void JST_Init(void);
void JST_Finalize(void);
void JST_Reset(void);
void Compute_Residual_JST(int AddTime);

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

