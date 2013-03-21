/*******************************************************************************
 * File:        SolverParameters.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#include <fstream>

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
int    VariableType;

// Linear Solver Parameters
int    SolverMethod;
int    SolverScheme;
int    FluxScheme;
int    TimeIntegrationMethod;
int    TimeIntegrationType;
int    BCMethod;
int    JacobianMethod;
int    JacobianUpdate;
int    SolverOrder;
int    SolverNIteration;
int    InnerNIteration;
int    LinearSolverNIteration;
int    FirstOrderNIteration;
double PhysicalDeltaTime;
double Relaxation;

// Residual Smoother Variables
int    ResidualSmoothMethod;
int    ResidualSmoothType;
int    ResidualSmoothScheme;
int    ResidualSmoothNIteration;
double ResidualSmoothRelaxation;

// Flux Limiter
int    LimiterMethod;
int    LimiterSmooth;
int    LimiterOrder;
int    LimiterStartSolverIteration;
int    LimiterEndSolverIteration;
double Venkat_KThreshold;

// Entropy Fix
int    EntropyFix;

// Precondition Variables
int    PrecondMethod;
int    PrecondType;
int    PrecondSmooth;
int    PrecondVariableType;

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
    VariableType                = VARIABLE_NONE;
    
    // Linear Solver Parameters
    SolverMethod                = SOLVER_METHOD_NONE;
    SolverScheme                = SOLVER_SCHEME_NONE;
    FluxScheme                  = FLUX_SCHEME_NONE;
    PrecondMethod               = PRECOND_METHOD_NONE;
    TimeIntegrationType         = TIME_INTEGRATION_TYPE_NONE;
    TimeIntegrationMethod       = TIME_INTEGRATION_METHOD_NONE;
    BCMethod                    = BC_METHOD_NONE;
    JacobianMethod              = JACOBIAN_METHOD_NONE;
    JacobianUpdate              = 0;
    SolverOrder                 = SOLVER_ORDER_NONE;
    SolverNIteration            = 0;
    InnerNIteration             = 0;
    LinearSolverNIteration      = 0;
    FirstOrderNIteration        = 0;
    PhysicalDeltaTime           = 0.0;
    Relaxation                  = 0.0;

    // Residual Smoother Variables
    ResidualSmoothMethod        = RESIDUAL_SMOOTH_METHOD_NONE;
    ResidualSmoothType          = RESIDUAL_SMOOTH_TYPE_NONE;
    ResidualSmoothScheme        = RESIDUAL_SMOOTH_SCHEME_NONE;
    ResidualSmoothNIteration    = 0;
    ResidualSmoothRelaxation    = 0.0;
    
    // Flux Limiter
    LimiterMethod               = LIMITER_METHOD_NONE;
    LimiterSmooth               = 0;
    LimiterOrder                = LIMITER_ORDER_NONE;
    LimiterStartSolverIteration = 0;
    LimiterEndSolverIteration   = 0;
    Venkat_KThreshold           = 0.0;

    // Entropy Fix
    EntropyFix                  = 0;
    
    // Precondition Variables
    PrecondSmooth               = 0;
    PrecondType                 = PRECOND_TYPE_NONE;
    PrecondVariableType         = VARIABLE_NONE;
    
    // Restart Parameters
    RestartInput                = 0;
    RestartOutput               = 0;
    RestartIteration            = 0;
    RestartCycle                = 0;
    str_blank(RestartInputFilename);
    str_blank(RestartOutputFilename);

    // Reference and Material Properties
    Material_Init();
    
    // CFL Conditions
    CFL_Ramp                    = 0;
    CFL_MAX                     = 0.0;
    CFL_MIN                     = 0.0;
    CFL                         = 0.0;

    // Mach Ramping
    Mach_Ramp                   = 0;
    Mach_MAX                    = 0.0;
    Mach_MIN                    = 0.0;

    // Zero Pressure Gradient No of Iterations
    ZPGIteration                = 0;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_Parameters_Read(const char *filename) {
    ifstream param_file;
    string text_line;
    string option_name;
    string::size_type position;
    size_t strlength;
    
    // Open the Parameter File
    param_file.open(filename, ios::in);
    if (param_file.fail())
        error("Solver_Parameters_Read: Unable to Read Parameter File - %s", filename);
    
    // Start Reading the Parameter File
    while (getline(param_file, text_line)) {
        // Check if line is comment
        position = text_line.find("#", 0);
        if (position != string::npos)
            continue;
        
        // Remove Any Space, Return or End Characters
        for (size_t i = 0 ; i < text_line.size(); i++) {
            position = text_line.find(" ", 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find("\r", 0);
            if (position != string::npos) text_line.erase(position, 1);
            position = text_line.find("\n", 0);
            if (position != string::npos) text_line.erase(position, 1);
        }
        
        // Get the Mesh Input File Name
        position = text_line.find("MESH_FILENAME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            strlength = text_line.copy(MeshInputFilename, text_line.size(), 0);
            MeshInputFilename[strlength] = '\0';
        }
        
        // Get the Boundary Condition Input File Name
        position = text_line.find("BC_FILENAME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 12);
            strlength = text_line.copy(BCInputFilename, text_line.size(), 0);
            BCInputFilename[strlength] = '\0';
        }
        
        // Get the Solution Output File Name
        position = text_line.find("SOLUTION_FILENAME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 18);
            strlength = text_line.copy(SolutionOutputFilename, text_line.size(), 0);
            SolutionOutputFilename[strlength] = '\0';
        }
        
        // Get the Restart Input File Name
        position = text_line.find("RESTART_INPUT_FILENAME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 23);
            strlength = text_line.copy(RestartInputFilename, text_line.size(), 0);
            RestartInputFilename[strlength] = '\0';
        }
        
        // Get the Restart Output File Name
        position = text_line.find("RESTART_OUTPUT_FILENAME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 24);
            strlength = text_line.copy(RestartOutputFilename, text_line.size(), 0);
            RestartOutputFilename[strlength] = '\0';
        }
        
        // Get the Mesh Reorder
        position = text_line.find("MESH_REORDER=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, MeshReorder, SolverBooleanMap);
        }
        
        // Get the Solver Method
        position = text_line.find("SOLVER_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, SolverMethod, SolverMethodMap);
        }
        
        // Get the Solver Scheme
        position = text_line.find("SOLVER_SCHEME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, SolverScheme, SolverSchemeMap);
        }
        
        // Get the Flux Scheme
        position = text_line.find("FLUX_SCHEME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 12);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, FluxScheme, FluxSchemeMap);
        }
        
        // Get the Solver Main Iterations 
        position = text_line.find("SOLVER_NITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 18);
            option_name.assign(text_line);
            SolverNIteration = atoi(option_name.c_str());
        }
        
        // Get the Solver Inner Iterations: Newton or Dual Time 
        position = text_line.find("SOLVER_INNER_NITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 24);
            option_name.assign(text_line);
            InnerNIteration = atoi(option_name.c_str());
        }
        
        // Get the Linear Solver Iterations 
        position = text_line.find("LINEAR_SOLVER_NITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 25);
            option_name.assign(text_line);
            LinearSolverNIteration = atoi(option_name.c_str());
        }
        
        // Get the Solver Order
        position = text_line.find("SOLVER_ORDER=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, SolverOrder, SolverOrderMap);
        }
        
        // Get the Time Integration Method
        position = text_line.find("TIME_INTEGRATION_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 24);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, TimeIntegrationMethod, TimeIntegrationMethodMap);
        }
        
        // Get the Time Integration Type
        position = text_line.find("TIME_INTEGRATION_TYPE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 22);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, TimeIntegrationType, TimeIntegrationTypeMap);
        }
        
        // Get the Precondition Method
        position = text_line.find("PRECOND_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 15);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, PrecondMethod, PrecondMethodMap);
        }
        
        // Get the Precondition Type
        position = text_line.find("PRECOND_TYPE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, PrecondType, PrecondTypeMap);
        }
        
        // Get the Precondition Smooth
        position = text_line.find("PRECOND_SMOOTH=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 15);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, PrecondSmooth, SolverBooleanMap);
        }
        
        // Get the Boundary Condition Method
        position = text_line.find("BC_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 10);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, BCMethod, BCMethodMap);
        }
        
        // Read No of Zero Pressure Gradient ZPG Iterations
        position = text_line.find("ZPG_ITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            ZPGIteration = atoi(option_name.c_str());
        }
        
        // Get the Variable Type for Computation
        position = text_line.find("VARIABLE_TYPE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, VariableType, VariableTypeMap);
        }
        
        // Get the Entropy Fix
        position = text_line.find("ENTROPY_FIX=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 12);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, EntropyFix, SolverBooleanMap);
        }
        
        // Get the Jacobian Method
        position = text_line.find("JACOBIAN_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 16);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, JacobianMethod, JacobianMethodMap);
        }
        
        // Get the Jacobian Update Frequency
        position = text_line.find("JACOBIAN_UPDATE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 16);
            option_name.assign(text_line);
            JacobianUpdate = atoi(option_name.c_str());
        }
        
        // Get the Number of First Order Iterations for Higher Order Solver
        position = text_line.find("FIRST_ORDER_NITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 23);
            option_name.assign(text_line);
            FirstOrderNIteration = atoi(option_name.c_str());
        }
        
        // Get the Physical Delta Time
        position = text_line.find("PHYSICAL_DELTA_TIME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 20);
            option_name.assign(text_line);
            PhysicalDeltaTime = atof(option_name.c_str());
        }
        
        // Get the Relaxation Factor
        position = text_line.find("SOLVER_RELAXATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 18);
            option_name.assign(text_line);
            Relaxation = atof(option_name.c_str());
        }
        
        // Residual Smoother Method
        position = text_line.find("RESIDUAL_SMOOTH_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 23);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, ResidualSmoothMethod, ResidualSmoothMethodMap);
        }
        
        // Residual Smoother Type
        position = text_line.find("RESIDUAL_SMOOTH_TYPE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 21);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, ResidualSmoothType, ResidualSmoothTypeMap);
        }
        
        // Residual Smoother Scheme
        position = text_line.find("RESIDUAL_SMOOTH_SCHEME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 23);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, ResidualSmoothScheme, ResidualSmoothSchemeMap);
        }
        
        // Get the Number of Residual Smoother Iterations
        position = text_line.find("RESIDUAL_SMOOTH_NITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 27);
            option_name.assign(text_line);
            ResidualSmoothNIteration = atoi(option_name.c_str());
        }
        
        // Get the Residual Smoother Relaxation Factor
        position = text_line.find("RESIDUAL_SMOOTH_RELAXATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 27);
            option_name.assign(text_line);
            ResidualSmoothRelaxation = atof(option_name.c_str());
        }
        
        // Get the Flux Limiter Method
        position = text_line.find("LIMITER_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 15);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, LimiterMethod, LimiterMethodMap);
        }
        
        // Get the Flux Limiter Smooth Boolean
        position = text_line.find("LIMITER_SMOOTH=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 15);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, LimiterSmooth, SolverBooleanMap);
        }
        
        // Get the Flux Limiter Order
        position = text_line.find("LIMITER_ORDER=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, LimiterOrder, LimiterOrderMap);
        }
        
        // Get the Limiter Start Solver Iteration
        position = text_line.find("LIMITER_START_SOLVER_ITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 31);
            option_name.assign(text_line);
            LimiterStartSolverIteration = atoi(option_name.c_str());
        }
        
        // Get the Limiter End Solver Iteration
        position = text_line.find("LIMITER_END_SOLVER_ITERATION=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 29);
            option_name.assign(text_line);
            LimiterEndSolverIteration = atoi(option_name.c_str());
        }
        
        // Read Venkatakrishanan Limiter K Threshold
        position = text_line.find("VENKAT_K_THRESHOLD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 19);
            option_name.assign(text_line);
            Venkat_KThreshold = atof(option_name.c_str());
        }
        
        // Get the Non Dimensionalization Method for Computation
        position = text_line.find("NONDIMENSIONAL_METHOD=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 22);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, NonDimensionalMethod, NonDimensionalMethodMap);
        }
        
        // Read Gamma
        position = text_line.find("GAMMA=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 6);
            option_name.assign(text_line);
            Gamma = atof(option_name.c_str());
        }
        
        // Read Reference Density
        position = text_line.find("REF_DENSITY=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 12);
            option_name.assign(text_line);
            Ref_Rho = atof(option_name.c_str());
        }
        
        // Read Reference Mach
        position = text_line.find("REF_MACH=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 9);
            option_name.assign(text_line);
            Ref_Mach = atof(option_name.c_str());
        }
        
        // Read Reference Pressure
        position = text_line.find("REF_PRESSURE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            Ref_Pressure = atof(option_name.c_str());
        }
        
        // Read Reference Alpha
        position = text_line.find("REF_ALPHA=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 10);
            option_name.assign(text_line);
            Ref_Alpha = atof(option_name.c_str());
        }
        
        // Read Reference Beta
        position = text_line.find("REF_BETA=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 9);
            option_name.assign(text_line);
            Ref_Beta = atof(option_name.c_str());
        }
        
        // Read CFL Ramp
        position = text_line.find("CFL_RAMP=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 9);
            option_name.assign(text_line);
            CFL_Ramp = atoi(option_name.c_str());
        }
        
        // Read CFL Minimum
        position = text_line.find("CFL_MIN=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 8);
            option_name.assign(text_line);
            CFL_MIN = atof(option_name.c_str());
        }
        
        // Read CFL Maximum
        position = text_line.find("CFL_MAX=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 8);
            option_name.assign(text_line);
            CFL_MAX = atof(option_name.c_str());
        }
        
        // Read Mach Ramp
        position = text_line.find("MACH_RAMP=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 10);
            option_name.assign(text_line);
            Mach_Ramp = atoi(option_name.c_str());
        }
        
        // Read Mach Minimum
        position = text_line.find("MACH_MIN=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 9);
            option_name.assign(text_line);
            Mach_MIN = atof(option_name.c_str());
        }
        
        // Read Mach Maximum
        position = text_line.find("MACH_MAX=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 9);
            option_name.assign(text_line);
            Mach_MAX = atof(option_name.c_str());
        }
        
        // Read Restart Solution Input Boolean
        position = text_line.find("RESTART_INPUT=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, RestartInput, SolverBooleanMap);
        }
        
        // Read Restart Solution Output Boolean
        position = text_line.find("RESTART_OUTPUT=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 15);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, RestartOutput, SolverBooleanMap);
        }
        
        // Read Restart Solution Cycle
        position = text_line.find("RESTART_DUMP_CYCLE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 19);
            option_name.assign(text_line);
            RestartCycle = atoi(option_name.c_str());
        }
    }
    
    // Close the Parameter File
    param_file.close();
    
    // Set the Remaining Quantities
    CFL             = CFL_MIN;
    
    // Reference and Material Properties
    if (PrecondMethod != PRECOND_METHOD_NONE) {
        switch (PrecondMethod) {
            case PRECOND_METHOD_ROE_LMFIX:
                PrecondVariableType = VARIABLE_PRIMITIVE_RUP;
                break;
            case PRECOND_METHOD_ROE_WS:
                PrecondVariableType = VARIABLE_PRIMITIVE_PUT;
                break;
            case PRECOND_METHOD_ROE_CV:
                PrecondVariableType = VARIABLE_PRIMITIVE_PUS;
                break;
            case PRECOND_METHOD_ROE_CV_ORIGINAL:
                PrecondVariableType = VARIABLE_PRIMITIVE_PUS;
                break;
            case PRECOND_METHOD_ROE_BTW:
                PrecondVariableType = VARIABLE_PRIMITIVE_RUP;
                break;
            case PRECOND_METHOD_ROE_BTW_ORIGINAL:
                PrecondVariableType = VARIABLE_PRIMITIVE_RUP;
                break;
            case PRECOND_METHOD_ROE_ERIKSSON:
                PrecondVariableType = VARIABLE_PRIMITIVE_RUP;
                break;
            case PRECOND_METHOD_ROE_ERIKSSON_ORIGINAL:
                PrecondVariableType = VARIABLE_PRIMITIVE_RUP;
                break;
            default:
                PrecondVariableType = VARIABLE_CONSERVATIVE;
                break;
        }
    } else {
        PrecondVariableType = VARIABLE_PRIMITIVE_RUP;
    }
    
    // Avoid Invalid Mach
    if (Ref_Mach < Mach_MAX)
        Mach_MAX = Ref_Mach;
    if (Mach_MIN == 0.0)
        Mach_MIN = Ref_Mach;
    
    // Set the Reference, Non-Dimensional and Material Properties
    Material_Set_Properties();

    int count = 0;
    printf("=============================================================================\n");
    info("%2d) Mesh Input Filename ------------------: %s", ++count, MeshInputFilename);
    info("%2d) Boundary Condition Input Filename ----: %s", ++count, BCInputFilename);
    info("%2d) Solution Output Filename -------------: %s", ++count, SolutionOutputFilename);
    info("%2d) Restart Input Filename ---------------: %s", ++count, RestartInputFilename);
    info("%2d) Restart Output Filename --------------: %s", ++count, RestartOutputFilename);
    GetOptionValueName(option_name, MeshReorder, SolverBooleanMap);
    info("%2d) Mesh Reorder -------------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, SolverMethod, SolverMethodMap);
    info("%2d) Solver Method ------------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, SolverScheme, SolverSchemeMap);
    info("%2d) Solver Scheme ------------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, FluxScheme, FluxSchemeMap);
    info("%2d) Flux Scheme --------------------------: %s",  ++count, option_name.c_str());
    info("%2d) No of Solver Iterations --------------: %d",  ++count, SolverNIteration);
    info("%2d) No of Inner (Newton/Dual) Iterations -: %d",  ++count, InnerNIteration);
    info("%2d) No of Linear Solver Iterations -------: %d",  ++count, LinearSolverNIteration);
    info("%2d) Solver Order -------------------------: %d",  ++count, SolverOrder);
    GetOptionValueName(option_name, TimeIntegrationMethod, TimeIntegrationMethodMap);
    info("%2d) Time Integration Method --------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, TimeIntegrationType, TimeIntegrationTypeMap);
    info("%2d) Time Integration Type ----------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, PrecondMethod, PrecondMethodMap);
    info("%2d) Precondition Method ------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, PrecondType, PrecondTypeMap);
    info("%2d) Preconditioner Type ------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, PrecondSmooth, SolverBooleanMap);
    info("%2d) Preconditioner Smooth  ---------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, BCMethod, BCMethodMap);
    info("%2d) Boundary Condition Method ------------: %s",  ++count, option_name.c_str());
    info("%2d) No of ZPG Iterations -----------------: %d",  ++count, ZPGIteration);
    GetOptionValueName(option_name, VariableType, VariableTypeMap);
    info("%2d) Variable Type for Computation ---- ---: %s",  ++count, option_name.c_str());
    info("%2d) Entropy Fix --------------------------: %d",  ++count, EntropyFix);
    GetOptionValueName(option_name, JacobianMethod, JacobianMethodMap);
    info("%2d) Jacobian Method ----------------------: %s",  ++count, option_name.c_str());
    info("%2d) Jacobian Update Frequency ------------: %d",  ++count, JacobianUpdate);
    info("%2d) No of First Order Iterations ---------: %d",  ++count, FirstOrderNIteration);
    info("%2d) Physical Delta Time ------------------: %lf", ++count, PhysicalDeltaTime);
    info("%2d) Relaxation Factor --------------------: %lf", ++count, Relaxation);
    GetOptionValueName(option_name, ResidualSmoothMethod, ResidualSmoothMethodMap);
    info("%2d) Residual Smoother Method -------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, ResidualSmoothType, ResidualSmoothTypeMap);
    info("%2d) Residual Smoother Type ---------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, ResidualSmoothScheme, ResidualSmoothSchemeMap);
    info("%2d) Residual Smoother Scheme -------------: %s",  ++count, option_name.c_str());
    info("%2d) No of Residual Smoother Iterations ---: %d",  ++count, ResidualSmoothNIteration);
    info("%2d) Residual Smoother Relaxation Factor --: %lf", ++count, ResidualSmoothRelaxation);
    GetOptionValueName(option_name, LimiterMethod, LimiterMethodMap);
    info("%2d) Limiter Method  ----------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, LimiterSmooth, SolverBooleanMap);
    info("%2d) Limiter Smooth -----------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, LimiterOrder, LimiterOrderMap);
    info("%2d) Limiter Order ------------------------: %s",  ++count, option_name.c_str());
    info("%2d) Limiter Start Solver Iteration -------: %d",  ++count, LimiterStartSolverIteration);
    info("%2d) Limiter End Solver Iteration ---------: %d",  ++count, LimiterEndSolverIteration);
    GetOptionValueName(option_name, NonDimensionalMethod, NonDimensionalMethodMap);
    info("%2d) Non Dimensional Method ---------------: %s",  ++count, option_name.c_str());
    info("%2d) Venkatakrishan Limiter Threshold -----: %lf", ++count, Venkat_KThreshold);
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
    info("%2d) Restart Input ------------------------: %d",  ++count, RestartInput);
    info("%2d) Restart Output -----------------------: %d",  ++count, RestartOutput);
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