/*******************************************************************************
 * File:        SolverParameters.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Trim_Utils.h"
#include "Commons.h"
#include "SolverParameters.h"
#include "Material.h"

// Mesh Parameters
char   MeshInputFilename[256];
int    MeshReorder;

// Boundary Condition Parameters
char   BCInputFilename[256];

// Solution Parameters
char   SolutionOutputFilename[256];
char   RMSOutputFilename[256];

// Variable Type
int    Variable_Type;

// Linear Solver Parameters
int    SolverMethod;
int    SolverScheme;
int    PrecondMethod;
int    TimeAccuracy;
int    TimeStepScheme;
int    SolverBCScheme;
int    JacobianMethod;
int    JacobianUpdate;
int    Order;
int    NIteration;
int    InnerNIteration;
int    LinearSolverNIteration;
int    FirstOrderNIteration;
double PhysicalDeltaTime;
double Relaxation;

// Residual Smoother Variables
int    ResidualSmoothType;
int    ResidualSmoothRegion;
int    ResidualSmoothNIteration;
double ResidualSmoothRelaxation;

// Flux Limiter
int    Limiter;
int    LimiterSmooth;
int    LimiterOrder;
int    StartLimiterNIteration;
int    EndLimiterNIteration;
double Venkat_KThreshold;

// Entropy Fix
int    EntropyFix;

// Precondition Variables
int    PrecondSmooth;
int    PrecondType;

// Restart Parameters
int    RestartInput;
int    RestartOutput;
int    RestartIteration;
int    RestartCycle;
char   RestartInputFilename[256];
char   RestartOutputFilename[256];

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


//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Parameters_Init(void) {
    // Mesh Parameters
    str_blank(MeshInputFilename);
    MeshReorder             = 0;
    
    // Boundary Condition Parameters
    str_blank(BCInputFilename);
    
    // Solution Parameters
    str_blank(SolutionOutputFilename);
    str_blank(RMSOutputFilename);
    
    // Variable Type
    Variable_Type            = VARIABLE_NONE;
    
    // Linear Solver Parameters
    SolverMethod             = SOLVER_METHOD_NONE;
    SolverScheme             = SOLVER_SCHEME_NONE;
    PrecondMethod            = SOLVER_PRECOND_NONE;
    TimeAccuracy             = SOLVER_TIME_ACCURACY_NONE;
    TimeStepScheme           = SOLVER_TIME_STEP_SCHEME_NONE;
    SolverBCScheme           = SOLVER_BC_SCHEME_NONE;
    JacobianMethod           = SOLVER_JACOBIAN_NONE;
    JacobianUpdate           = 0;
    Order                    = SOLVER_ORDER_NONE;
    NIteration               = 0;
    InnerNIteration          = 0;
    LinearSolverNIteration   = 0;
    FirstOrderNIteration     = 0;
    PhysicalDeltaTime        = 0.0;
    Relaxation               = 0.0;

    // Residual Smoother Variables
    ResidualSmoothType       = RESIDUAL_SMOOTH_NONE;
    ResidualSmoothRegion     = RESIDUAL_SMOOTH_REGION_NONE;
    ResidualSmoothNIteration = 0;
    ResidualSmoothRelaxation = 0.0;
    
    // Flux Limiter
    Limiter                  = 0;
    LimiterSmooth            = 0;
    LimiterOrder             = 0;
    StartLimiterNIteration   = 0;
    EndLimiterNIteration     = 0;
    Venkat_KThreshold        = 0.0;

    // Entropy Fix
    EntropyFix               = 0;
    
    // Precondition Variables
    PrecondSmooth            = 0;
    PrecondType              = PRECOND_TYPE_NONE;
    
    // Restart Parameters
    RestartInput             = 0;
    RestartOutput            = 0;
    RestartIteration         = 0;
    RestartCycle             = 0;
    str_blank(RestartInputFilename);
    str_blank(RestartOutputFilename);

    // Reference and Material Properties
    Material_Init();
    
    // CFL Conditions
    CFL_Ramp                 = 0;
    CFL_MAX                  = 0.0;
    CFL_MIN                  = 0.0;
    CFL                      = 0.0;

    // Mach Ramping
    Mach_Ramp                = 0;
    Mach_MAX                 = 0.0;
    Mach_MIN                 = 0.0;

    // Zero Pressure Gradient No of Iterations
    ZPGIteration             = 0;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Parameters_Read(const char *filename) {
    FILE *fp;
    int bdim = 256;
    char buff[256];
    char *dummy;

    if ((fp = fopen(filename, "r")) == NULL)
        error("Solver_Parameters_Read: Unable to Read Parameter File - %s", filename);
    
    // Mesh Input Filename
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%s", MeshInputFilename);
    
    // Get Mesh Reorder 
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &MeshReorder);
    
    // BC Input Filename
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%s", BCInputFilename);
    
    // Solution Output Filename
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%s", SolutionOutputFilename);
    
    // Solver Parameters
    // 1) Get the Solver Method
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &SolverMethod);
    
    // 2) Get the Solver Scheme
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &SolverScheme);

    // 3) Get the Solver Precondition Method
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &PrecondMethod);
    
    // 4) Get the Time Accuracy
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &TimeAccuracy);

    // 5) Get the Time Stepping Scheme
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &TimeStepScheme);
    
    // 6) Get the Boundary Condition Scheme
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &SolverBCScheme);
    
    // 7) Get the Jacobian Method
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &JacobianMethod);
    
    // 8) Get the Jacobian Update Frequency
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &JacobianUpdate);
    
    // 9) Get the Order
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &Order);
    
    // 10) Get Number of Main Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &NIteration);

    // 11) Get Number of Inner Iterations: Newton or Dual Time
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &InnerNIteration);

    // 12) Get Number of Linear Solver Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &LinearSolverNIteration);
    
    // 13) Get Number of First Order Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &FirstOrderNIteration);
    
    // 14) Get the Physical Delta Time
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &PhysicalDeltaTime);
    
    // 15) Get the Relaxation Factor
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Relaxation);

    // 16) Residual Smoother Type
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &ResidualSmoothType);
    
    // 17) Residual Smoother Region Type
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &ResidualSmoothRegion);
    
    // 18) Residual Smoother Number of Iteration
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &ResidualSmoothNIteration);
    
    // 19) Get the Residual Smoother Relaxation Factor
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &ResidualSmoothRelaxation);
    
    // 20) Get the Flux Limiter
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &Limiter);

    // 21) Get the Flux Limiter Smooth
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &LimiterSmooth);

    // 22) Get the Flux Limiter Order
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &LimiterOrder);
    
    // 23) Get the Start Limiter Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &StartLimiterNIteration);

    // 24) Get the End Limiter Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &EndLimiterNIteration);

    // 25) Read Ventakakrishanan K Threshold
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Venkat_KThreshold);

    // 26) Read Entropy Fix
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &EntropyFix);
    
    // 27) Read Precondition Smooth
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &PrecondSmooth);
    
    // 28) Read Precondition Type
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &PrecondType);
    
    // Reference and Material Properties
    // 29) Read Gamma
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Gamma);

    // 30) Read Reference Density
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Rho);

    // 31) Read Reference Mach
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Mach);

    // 32) Read Reference Pressure
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Pressure);

    // 33) Read Alpha
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Alpha);

    // 34) Read Beta
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Ref_Beta);
    
    // 35) Read CFL Ramp
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &CFL_Ramp);

    // 36) Read CFL MIN
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &CFL_MIN);

    // 37) Read CFL MAX
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &CFL_MAX);

    // 38) Read Mach Ramp
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &Mach_Ramp);

    // 39) Read Mach MIN
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Mach_MIN);

    // 40) Read Mach MAX
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%lf", &Mach_MAX);

    // 41) Read No of Zero Pressure Gradient ZPG Iterations
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &ZPGIteration);

    // 42) Read Restart Solution Input
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &RestartInput);
    
    // 43) Read Restart Solution Input Filename
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%s", RestartInputFilename);
    
    // 44) Read Restart Solution Output
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &RestartOutput);
    
    // 45) Read Restart Solution Output Filename
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%s", RestartOutputFilename);
    
    // 46) Read Restart Solution Cycle
    dummy = fgets(buff, bdim, fp);
    dummy = fgets(buff, bdim, fp);
    sscanf(buff, "%d", &RestartCycle);
    
    // Set the Remaining Quantities
    CFL             = CFL_MIN;
    
    // Reference and Material Properties
    if (PrecondMethod != SOLVER_PRECOND_NONE) {
        switch (PrecondMethod) {
            case SOLVER_PRECOND_ROE_LMFIX:
                Variable_Type = VARIABLE_CONSERVATIVE;
                break;
            case SOLVER_PRECOND_ROE_WS:
                Variable_Type = VARIABLE_PRIMITIVE_PUT;
                break;
            case SOLVER_PRECOND_ROE_CV:
                Variable_Type = VARIABLE_CONSERVATIVE;
                break;
            case SOLVER_PRECOND_ROE_CV_ORIGINAL:
                Variable_Type = VARIABLE_CONSERVATIVE;
                break;
            case SOLVER_PRECOND_ROE_BTW:
                Variable_Type = VARIABLE_PRIMITIVE_RUP;
                break;
            case SOLVER_PRECOND_ROE_BTW_ORIGINAL:
                Variable_Type = VARIABLE_PRIMITIVE_RUP;
                break;
            case SOLVER_PRECOND_ROE_ERIKSSON:
                Variable_Type = VARIABLE_PRIMITIVE_RUP;
                break;
            case SOLVER_PRECOND_ROE_ERIKSSON_ORIGINAL:
                Variable_Type = VARIABLE_PRIMITIVE_RUP;
                break;
            default:
                Variable_Type = VARIABLE_CONSERVATIVE;
                break;
        }
    } else
        Variable_Type = VARIABLE_CONSERVATIVE;
    
    // Avoid Invalid Mach
    if (Ref_Mach < Mach_MAX)
        Mach_MAX = Ref_Mach;
    if (Mach_MIN == 0.0)
        Mach_MIN = Ref_Mach;
    
    // Set the Reference, Non-Dimensional and Material Properties
    Material_Set_Properties();
    
    // Close file
    fclose(fp);

    printf("=============================================================================\n");
    info("Mesh Parameters");
    info("Mesh Input Filename ----------------------: %s",  MeshInputFilename);
    info("Mesh Reorder -----------------------------: %d",  MeshReorder);
    printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n");
    info("Boundary Condition Parameters");
    info("Boundary Condition input Filename --------: %s",  BCInputFilename);
    printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n");
    info("Solution Parameters");
    info("Solution Output Filename -----------------: %s",  SolutionOutputFilename);
    printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n");
    int count = 0;
    info("Solver Parameters");
    info("%2d) Solver Method ------------------------: %d",  ++count, SolverMethod);
    info("%2d) Solver Scheme ------------------------: %d",  ++count, SolverScheme);
    info("%2d) Precondition Method ------------------: %d",  ++count, PrecondMethod);
    info("%2d) Time Accuracy ------------------------: %d",  ++count, TimeAccuracy);
    info("%2d) Time Stepping Scheme -----------------: %d",  ++count, TimeStepScheme);
    info("%2d) Solver Boundary Condition Scheme -----: %d",  ++count, SolverBCScheme);
    info("%2d) Solver Jacobian Method ---------------: %d",  ++count, JacobianMethod);
    info("%2d) Solver Jacobian Update Frequency -----: %d",  ++count, JacobianUpdate);
    info("%2d) Order --------------------------------: %d",  ++count, Order);
    info("%2d) No of Iterations ---------------------: %d",  ++count, NIteration);
    info("%2d) No of Inner Iterations ---------------: %d",  ++count, InnerNIteration);
    info("%2d) No of Linear Solver Iterations -------: %d",  ++count, LinearSolverNIteration);
    info("%2d) No of First Order Iterations ---------: %d",  ++count, FirstOrderNIteration);
    info("%2d) Physical Delta Time ------------------: %lf", ++count, PhysicalDeltaTime);
    info("%2d) Relaxation Factor --------------------: %lf", ++count, Relaxation);
    info("%2d) Residual Smoother Type ---------------: %d",  ++count, ResidualSmoothType);
    info("%2d) Residual Smoother Region -------------: %d",  ++count, ResidualSmoothRegion);
    info("%2d) No of Residual Smoother Iterations ---: %d",  ++count, ResidualSmoothNIteration);
    info("%2d) Residual Smoother Relaxation Factor --: %lf", ++count, ResidualSmoothRelaxation);
    info("%2d) Limiter Type -------------------------: %d",  ++count, Limiter);
    info("%2d) Limiter Smooth -----------------------: %d",  ++count, LimiterSmooth);
    info("%2d) Limiter Order ------------------------: %d",  ++count, LimiterOrder);
    info("%2d) Limiter Start Iteration --------------: %d",  ++count, StartLimiterNIteration);
    info("%2d) Limiter End Iteration ----------------: %d",  ++count, EndLimiterNIteration);
    info("%2d) Venkatakrishan Limiter Threshold -----: %lf", ++count, Venkat_KThreshold);
    info("%2d) Entropy Fix --------------------------: %d",  ++count, EntropyFix);
    info("%2d) Preconditioner Smooth  ---------------: %d",  ++count, PrecondSmooth);
    info("%2d) Preconditioner Type ------------------: %d",  ++count, PrecondType);
    info("%2d) Gamma --------------------------------: %lf", ++count, Gamma);
    info("%2d) Reference Density --------------------: %lf", ++count, Ref_Rho);
    info("%2d) Reference Mach Number ----------------: %lf", ++count, Ref_Mach);
    info("%2d) Reference Pressure -------------------: %lf", ++count, Ref_Pressure);
    info("%2d) Reference Alpha ----------------------: %lf", ++count, Ref_Alpha);
    info("%2d) Reference Beta -----------------------: %lf", ++count, Ref_Beta);
    info("%2d) CFL Ramp -----------------------------: %d",  ++count, CFL_Ramp);
    info("%2d) CFL_Min ------------------------------: %lf", ++count, CFL_MIN);
    info("%2d) CFL_Max ------------------------------: %lf", ++count, CFL_MAX);
    info("%2d) Mach Ramp ----------------------------: %d",  ++count, Mach_Ramp);
    info("%2d) Mach_Min -----------------------------: %lf", ++count, Mach_MIN);
    info("%2d) Mach_Max -----------------------------: %lf", ++count, Mach_MAX);
    info("%2d) No of ZPG Iterations -----------------: %d",  ++count, ZPGIteration);
    info("%2d) Restart Input ------------------------: %d",  ++count, RestartInput);
    info("%2d) Restart Input Filename ---------------: %s",  ++count, RestartInputFilename);
    info("%2d) Restart Output -----------------------: %d",  ++count, RestartOutput);
    info("%2d) Restart Output Filename --------------: %s",  ++count, RestartOutputFilename);
    info("%2d) Restart Cycle ------------------------: %d",  ++count, RestartCycle);
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_BC_Parameters_Read(const char *filename) {
    FILE *fp;
    int bdim = 256;
    char buff[256];
    char *dummy;

    if ((fp = fopen(filename, "r")) == NULL)
        error("Solver_BC_Parameters_Read: Unable to Read Boundary Condition Parameter File - %s", filename);
    
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
    
    // Close file
    fclose(fp);
}