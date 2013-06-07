/*******************************************************************************
 * File:        Solver.h
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

#ifndef _SOLVER_H
#define	_SOLVER_H

// Custom header files
#include "SolverDB.h"
#include "SolverParameters.h"
#include "MC.h"
#include "BC.h"
#include "Roe_Fluxes.h"
#include "HLLC_Fluxes.h"
#include "VanLeer_Fluxes.h"
#include "AUSM_Fluxes.h"
#include "JST_Fluxes.h"
#include "LDFSS_Fluxes.h"
#include "Osher_Fluxes.h"
#include "StegerWarming_Fluxes.h"

// Solver Object
extern CSolver CogSolver;

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
extern double *Res1_Conv;
extern double *Res2_Conv;
extern double *Res3_Conv;
extern double *Res4_Conv;
extern double *Res5_Conv;

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
int  Solver(void);
int  Solver_Steady_Explicit(void);
int  Solver_Steady_Implicit(void);
int  Solver_Unsteady_Explicit(void);
int  Solver_Unsteady_Implicit(void);
void Compute_Flux(int node_L, int node_R, Vector3D areavec, double *Flux_Conv, double *Flux_Diss, int AddTime);
void Compute_Transformed_Preconditioner_Matrix(int nodeID, int reqInv, double **PrecondMatrix);
void Compute_Residual(int AddTime);
void Unsteady_Initialization(void);

// Implicit Method Routines
void Create_CRS_SolverBlockMatrix(void);
void Delete_CRS_SolverBlockMatrix(void);
int  Solve_Implicit(void);

// Time Stepping
void Compute_DeltaT(int Iteration);

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

