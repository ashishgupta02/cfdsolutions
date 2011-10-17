/*******************************************************************************
 * File:        SolverParameters.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _SOLVERPARAMETERS_H
#define	_SOLVERPARAMETERS_H

// Boundary Conditions Type
#define BC_NONE                -1
#define BC_SOLID_WALL           0
#define BC_FREE_STREAM          1
#define BC_SUBSONIC_INFLOW      2
#define BC_SUBSONIC_OUTFLOW     3
#define BC_SUPERSONIC_INFLOW    4
#define BC_SUPERSONIC_OUTFLOW   5

// Solver Method
#define SOLVER_METHOD_NONE      0
#define SOLVER_METHOD_EXPLICIT  1
#define SOLVER_METHOD_IMPLICIT  2

// Solver Scheme
#define SOLVER_SCHEME_NONE      0
#define SOLVER_SCHEME_ROE       1
#define SOLVER_SCHEME_LMROE     2
#define SOLVER_SCHEME_ROE_WS    3
#define SOLVER_SCHEME_ROE_CV    4
#define SOLVER_SCHEME_ROE_BTW   5

// Precondition Method
#define SOLVER_PRECOND_NONE             0
#define SOLVER_PRECOND_ROE_LMFIX        1
#define SOLVER_PRECOND_ROE_WS           2
#define SOLVER_PRECOND_ROE_CV           3
#define SOLVER_PRECOND_ROE_BTW          4

// Solver Boundary Condition Scheme
#define SOLVER_BC_SCHEME_NONE           0
#define SOLVER_BC_SCHEME_CHARCTERISTIC  1
#define SOLVER_BC_SCHEME_PRESSURE       2

// Solver Time Accuracy
#define SOLVER_TIME_ACCURACY_NONE               0
#define SOLVER_TIME_ACCURACY_GLOBAL             1
#define SOLVER_TIME_ACCURACY_LOCAL              2
#define SOLVER_TIME_ACCURACY_CHARCTERISTIC      3

// Solver Explict Time Step Scheme 
#define SOLVER_EXPLICT_TIME_STEP_SCHEME_NONE    0
#define SOLVER_EXPLICT_TIME_STEP_SCHEME_EULER   1
#define SOLVER_EXPLICT_TIME_STEP_SCHEME_RK4     2

// Solver Jacobian Computation Method
#define SOLVER_JACOBIAN_NONE           -1
#define SOLVER_JACOBIAN_CENTRAL         0
#define SOLVER_JACOBIAN_FORWARD         1
#define SOLVER_JACOBIAN_BACKWARD        2
#define SOLVER_JACOBIAN_ALTERNATE       3
#define SOLVER_JACOBIAN_APPROX          4
#define SOLVER_JACOBIAN_EXACT           5

// Variable Type
#define VARIABLE_NONE           0
#define VARIABLE_CONSERVATIVE   1
#define VARIABLE_RUP            2
#define VARIABLE_RUT            3
#define VARIABLE_PUT            4

// No of Equations
#define NEQUATIONS              5

// Mesh Parameters
extern char   MeshInputFilename[256];
extern int    MeshReorder;      /* Cuthill Mckee */

// Boundary Condition Parameters
extern char   BCInputFilename[256];

// Solution Parameters
extern char   SolutionOutputFilename[256];

// Linear Solver Parameters   
extern int    SolverMethod;     /* Explicit, Implicit */
extern int    SolverScheme;     /* Roe, LMRoe */
extern int    PrecondMethod;    /* Roe_LMFix, Roe_WS, Roe_CV, Roe_BTW */
extern int    TimeAccuracy;     /* Global Time, Local Time, Characteristic */
extern int    TimeStepScheme;   /* Euler, Runge-Kutta 4 */
extern int    SolverBCScheme;   /* Characteristic, Pressure */
extern int    JacobianMethod;   /* Central, Forward, Backward, Approx, Exact */
extern int    JacobianUpdate;   
extern int    Order;            /* First, Second */
extern int    NIteration;       /* Number of Main Iteration */ 
extern int    InnerNIteration;  /* Linear Solver Iterations */ 
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
extern double Ref_Beta;
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

// Gauge Pressure
extern double Gauge_Pressure;

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


// Function Definitions
void Solver_Parameters_Init(void);
void Solver_Parameters_Read(const char *filename);
void Solver_BC_Parameters_Read(const char *filename);

#endif	/* _SOLVERPARAMETERS_H */

