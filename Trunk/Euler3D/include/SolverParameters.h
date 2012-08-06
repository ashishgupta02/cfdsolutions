/*******************************************************************************
 * File:        SolverParameters.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _SOLVERPARAMETERS_H
#define	_SOLVERPARAMETERS_H

// Boundary Conditions Type
#define BC_NONE                                 -1
#define BC_SOLID_WALL                           0
#define BC_FREE_STREAM                          1
#define BC_INFLOW                               2
#define BC_OUTFLOW                              3
#define BC_SUBSONIC_INFLOW                      4
#define BC_SUBSONIC_OUTFLOW                     5
#define BC_SUPERSONIC_INFLOW                    6
#define BC_SUPERSONIC_OUTFLOW                   7
#define BC_FREE_STREAM_SUBSONIC_INFLOW          8
#define BC_FREE_STREAM_SUBSONIC_OUTFLOW         9
#define BC_FREE_STREAM_SUPERSONIC_INFLOW        10
#define BC_FREE_STREAM_SUPERSONIC_OUTFLOW       11

// Solver Method
#define SOLVER_METHOD_NONE                      0
#define SOLVER_METHOD_EXPLICIT                  1
#define SOLVER_METHOD_IMPLICIT                  2
#define SOLVER_METHOD_EXPLICIT_UNSTEADY         3
#define SOLVER_METHOD_IMPLICIT_UNSTEADY         4

// Solver Scheme
#define SOLVER_SCHEME_NONE                      0
#define SOLVER_SCHEME_ROE                       1
#define SOLVER_SCHEME_HLLC                      2
#define SOLVER_SCHEME_AUSM                      3
#define SOLVER_SCHEME_VANLEER                   4
#define SOLVER_SCHEME_LDFSS                     5
#define SOLVER_SCHEME_OSHER                     6
#define SOLVER_SCHEME_STEGERWARMING             7
#define SOLVER_SCHEME_JST                       8

// Precondition Method
#define SOLVER_PRECOND_NONE                     0
#define SOLVER_PRECOND_ROE_LMFIX                1
#define SOLVER_PRECOND_ROE_WS                   2
#define SOLVER_PRECOND_ROE_CV                   3
#define SOLVER_PRECOND_ROE_CV_ORIGINAL          4
#define SOLVER_PRECOND_ROE_BTW                  5
#define SOLVER_PRECOND_ROE_BTW_ORIGINAL         6
#define SOLVER_PRECOND_ROE_ERIKSSON             7
#define SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL    8
#define SOLVER_PRECOND_HLLC_LMFIX               9

// Precondition Type
#define PRECOND_TYPE_NONE                       0
#define PRECOND_TYPE_LOCAL                      1
#define PRECOND_TYPE_GLOBAL                     2

// Residual Smoother Type
#define RESIDUAL_SMOOTH_NONE                    0
#define RESIDUAL_SMOOTH_EXPLICIT                1
#define RESIDUAL_SMOOTH_EXPLICIT_WEIGHTED       2
#define RESIDUAL_SMOOTH_IMPLICIT                3
#define RESIDUAL_SMOOTH_IMPLICIT_CONSTANT       4

// Residual Smoother Type
#define RESIDUAL_SMOOTH_REGION_NONE             0
#define RESIDUAL_SMOOTH_REGION_LOCAL            1
#define RESIDUAL_SMOOTH_REGION_GLOBAL           2

// Solver Boundary Condition Scheme
#define SOLVER_BC_SCHEME_NONE                   0
#define SOLVER_BC_SCHEME_CHARCTERISTIC          1
#define SOLVER_BC_SCHEME_CHARCTERISTIC_WEAK     2
#define SOLVER_BC_SCHEME_PRESSURE               3

// Solver Time Accuracy
#define SOLVER_TIME_ACCURACY_NONE               0
#define SOLVER_TIME_ACCURACY_GLOBAL             1
#define SOLVER_TIME_ACCURACY_LOCAL              2
#define SOLVER_TIME_ACCURACY_CHARCTERISTIC      3

// Solver Time Step Scheme 
#define SOLVER_TIME_STEP_SCHEME_NONE            0
#define SOLVER_EXPLICIT_TIME_STEP_SCHEME_EULER  1
#define SOLVER_EXPLICIT_TIME_STEP_SCHEME_RK4    2
#define SOLVER_IMPLICIT_TIME_STEP_SCHEME_EULER  3
#define SOLVER_IMPLICIT_TIME_STEP_SCHEME_BDF    4

// Solver Jacobian Computation Method
#define SOLVER_JACOBIAN_NONE                    -1
#define SOLVER_JACOBIAN_CENTRAL                 0
#define SOLVER_JACOBIAN_FORWARD                 1
#define SOLVER_JACOBIAN_BACKWARD                2
#define SOLVER_JACOBIAN_ALTERNATE               3
#define SOLVER_JACOBIAN_APPROX                  4
#define SOLVER_JACOBIAN_EXACT                   5

// Variable Type
#define VARIABLE_NONE                           0
#define VARIABLE_CONSERVATIVE                   1
#define VARIABLE_PRIMITIVE_PUT                  2
#define VARIABLE_PRIMITIVE_RUP                  3
#define VARIABLE_PRIMITIVE_RUT                  4

// Solver Order
#define SOLVER_ORDER_NONE                       0
#define SOLVER_ORDER_FIRST                      1
#define SOLVER_ORDER_SECOND                     2

// No of Equations
#define NEQUATIONS                              5

// Mesh Parameters
extern char   MeshInputFilename[256];
extern int    MeshReorder;      /* Cuthill Mckee */

// Boundary Condition Parameters
extern char   BCInputFilename[256];

// Solution Parameters
extern char   SolutionOutputFilename[256];
extern char   RMSOutputFilename[256];

// Variable Type
extern int    Variable_Type;

// Linear Solver Parameters   
extern int    SolverMethod;     /* Explicit, Implicit */
extern int    SolverScheme;     /* Roe, LMRoe */
extern int    PrecondMethod;    /* Roe_LMFix, Roe_WS, Roe_CV, Roe_CV_Org, Roe_BTW */
extern int    TimeAccuracy;     /* Global Time, Local Time, Characteristic */
extern int    TimeStepScheme;   /* Euler, Runge-Kutta 4 */
extern int    SolverBCScheme;   /* Characteristic, Pressure */
extern int    JacobianMethod;   /* Central, Forward, Backward, Approx, Exact */
extern int    JacobianUpdate;   
extern int    Order;            /* First, Second */
extern int    NIteration;       /* Number of Main Iteration */ 
extern int    InnerNIteration;  /* Inner Iterations: Newton or Dual Time */
extern int    LinearSolverNIteration; /* Linear Solver Iterations */
extern int    FirstOrderNIteration;
extern double PhysicalDeltaTime; /* Physical Time Step For Unsteady Calculations */
extern double Relaxation;

// Residual Smoother Variables
extern int    ResidualSmoothType;
extern int    ResidualSmoothRegion;
extern int    ResidualSmoothNIteration;
extern double ResidualSmoothRelaxation;

// Flux Limiters
extern int    Limiter;
extern int    LimiterSmooth;
extern int    LimiterOrder;
extern int    StartLimiterNIteration;
extern int    EndLimiterNIteration;
extern double Venkat_KThreshold;

// Entropy Fix
extern int EntropyFix;

// Precondition Variables
extern int PrecondSmooth;
extern int PrecondType;

// Restart Parameters
extern int  RestartInput;
extern int  RestartOutput;
extern int  RestartIteration;
extern int  RestartCycle;
extern char RestartInputFilename[256];
extern char RestartOutputFilename[256];

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

void RMS_Writer_Init(void);
void RMS_Writer_Finalize(void);
void RMS_Writer(int Iteration, double *RMS);

#endif	/* _SOLVERPARAMETERS_H */

