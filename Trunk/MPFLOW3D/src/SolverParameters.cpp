/*******************************************************************************
 * File:        SolverParameters.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
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
    
    // CFL Conditions
    CFL_Ramp                    = 0;
    CFL_MAX                     = 0.0;
    CFL_MIN                     = 0.0;
    CFL                         = 0.0;

    // Zero Pressure Gradient No of Iterations
    ZPGIteration                = 0;
    
    // Material Properties
    Material_Init();
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
            // Check the input
            if (VariableType == VARIABLE_PUS || VariableType == VARIABLE_NONE)
                error("Solver_Parameters_Read: VARIABLE_TYPE=%s is invalid", option_name.c_str());
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
        
        //-------------------Residual Smoother Properties----------------------------
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
        
        //-------------------Flux Limiter Properties----------------------------
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
        
        //-------------------Material Properties--------------------------------
        // Read Material Type
        position = text_line.find("MATERIAL_TYPE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, MaterialType, MaterialTypeMap);
        }
        
        // Get the Material Name
        position = text_line.find("MATERIAL_NAME=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            strlength = text_line.copy(MaterialName, text_line.size(), 0);
            MaterialName[strlength] = '\0';
        }
        
        //-------------------Reference Values-----------------------------------
        // Read Reference Pressure
        position = text_line.find("REF_PRESSURE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            Ref_Pressure = atof(option_name.c_str());
        }
        
        // Read Reference Temperature
        position = text_line.find("REF_TEMPERATURE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 16);
            option_name.assign(text_line);
            Ref_Temperature = atof(option_name.c_str());
        }
        
        // Read Reference Temperature
        position = text_line.find("REF_LENGTH=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 11);
            option_name.assign(text_line);
            Ref_Length = atof(option_name.c_str());
        }
        
        //-------------------Far Field (Infinity) Values------------------------
        // Read Far Field Flow Pressure
        position = text_line.find("INF_PRESSURE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            Inf_Pressure = atof(option_name.c_str());
        }
        
        // Read Far Field Flow Temperature
        position = text_line.find("INF_TEMPERATURE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 16);
            option_name.assign(text_line);
            Inf_Temperature = atof(option_name.c_str());
        }
        
        // Read Far Field Flow Mach Number
        position = text_line.find("INF_MACH=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 9);
            option_name.assign(text_line);
            Inf_Mach = atof(option_name.c_str());
        }
        
        // Read Far Field Mach Ramp
        position = text_line.find("INF_MACH_RAMP=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 14);
            option_name.assign(text_line);
            Inf_Mach_Ramp = atoi(option_name.c_str());
        }
        
        // Read Far Field Mach Minimum
        position = text_line.find("INF_MACH_MIN=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            Inf_Mach_MIN = atof(option_name.c_str());
        }
        
        // Read Far Field Mach Maximum
        position = text_line.find("INF_MACH_MAX=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 13);
            option_name.assign(text_line);
            Inf_Mach_MAX = atof(option_name.c_str());
        }
        
        // Read Far Field Flow Angle Alpha
        position = text_line.find("INF_ALPHA=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 10);
            option_name.assign(text_line);
            Inf_Alpha = atof(option_name.c_str());
        }
        
        // Read Far Field Flow Angle Beta
        position = text_line.find("INF_BETA=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 9);
            option_name.assign(text_line);
            Inf_Beta = atof(option_name.c_str());
        }
        
        //-------------------Pressure Conditions------------------------
        // Read Outflow Flow Pressure
        position = text_line.find("OUTFLOW_PRESSURE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 17);
            option_name.assign(text_line);
            Outflow_Pressure = atof(option_name.c_str());
        }
        
        // Read Gauge Pressure
        position = text_line.find("GUAGE_PRESSURE=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 15);
            option_name.assign(text_line);
            Gauge_Pressure = atof(option_name.c_str());
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
    CFL = CFL_MIN;
    PrecondVariableType = VARIABLE_RUP;
    
    // Avoid Invalid Infinity Mach Ramping Conditions
    if (Inf_Mach_MAX != Inf_Mach)
        Inf_Mach_MAX = Inf_Mach;
    if (Inf_Mach_Ramp <= 1) {
        Inf_Mach_MIN  = Inf_Mach;
        Inf_Mach_Ramp = 1;
    } else if ((Inf_Mach_MIN <= 0.0) || (Inf_Mach_MIN >= Inf_Mach_MAX)) {
        Inf_Mach_MIN = Inf_Mach_MAX/((double)(Inf_Mach_Ramp - 1));
    }

    int count = 0;
    printf("=============================================================================\n");
    info("%2d) Mesh Input Filename ------------------: %s",  ++count, MeshInputFilename);
    info("%2d) Boundary Condition Input Filename ----: %s",  ++count, BCInputFilename);
    info("%2d) Solution Output Filename -------------: %s",  ++count, SolutionOutputFilename);
    info("%2d) Restart Input Filename ---------------: %s",  ++count, RestartInputFilename);
    info("%2d) Restart Output Filename --------------: %s",  ++count, RestartOutputFilename);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Solver Properties--------------------------------------------------\n");
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
    GetOptionValueName(option_name, SolverOrder, SolverOrderMap);
    info("%2d) Solver Order -------------------------: %s",  ++count, option_name.c_str());
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
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Residual Smoother Properties---------------------------------------\n");
    GetOptionValueName(option_name, ResidualSmoothMethod, ResidualSmoothMethodMap);
    info("%2d) Residual Smoother Method -------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, ResidualSmoothType, ResidualSmoothTypeMap);
    info("%2d) Residual Smoother Type ---------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, ResidualSmoothScheme, ResidualSmoothSchemeMap);
    info("%2d) Residual Smoother Scheme -------------: %s",  ++count, option_name.c_str());
    info("%2d) No of Residual Smoother Iterations ---: %d",  ++count, ResidualSmoothNIteration);
    info("%2d) Residual Smoother Relaxation Factor --: %lf", ++count, ResidualSmoothRelaxation);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Flux Limiter Properties--------------------------------------------\n");
    GetOptionValueName(option_name, LimiterMethod, LimiterMethodMap);
    info("%2d) Limiter Method  ----------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, LimiterSmooth, SolverBooleanMap);
    info("%2d) Limiter Smooth -----------------------: %s",  ++count, option_name.c_str());
    GetOptionValueName(option_name, LimiterOrder, LimiterOrderMap);
    info("%2d) Limiter Order ------------------------: %s",  ++count, option_name.c_str());
    info("%2d) Limiter Start Solver Iteration -------: %d",  ++count, LimiterStartSolverIteration);
    info("%2d) Limiter End Solver Iteration ---------: %d",  ++count, LimiterEndSolverIteration);
    info("%2d) Venkatakrishan Limiter Threshold -----: %lf", ++count, Venkat_KThreshold);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Material Type------------------------------------------------------\n");
    GetOptionValueName(option_name, MaterialType, MaterialTypeMap);
    info("%2d) Material Type ------------------------: %s",  ++count, option_name.c_str());
    info("%2d) Material Name ------------------------: %s",  ++count, MaterialName);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Reference Conditions-----------------------------------------------\n");
    info("%2d) Reference Pressure -------------------: %lf", ++count, Ref_Pressure);
    info("%2d) Reference Temperature ----------------: %lf", ++count, Ref_Temperature);
    info("%2d) Reference Length ---------------------: %lf", ++count, Ref_Length);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Far Field (Infinity) Conditions------------------------------------\n");
    info("%2d) Infinity Pressure --------------------: %lf", ++count, Inf_Pressure);
    info("%2d) Infinity Temperature -----------------: %lf", ++count, Inf_Temperature);
    info("%2d) Infinity Mach Number -----------------: %lf", ++count, Inf_Mach);
    info("%2d) Infinity Mach Ramp -------------------: %d",  ++count, Inf_Mach_Ramp);
    info("%2d) Infinity Mach_Min --------------------: %lf", ++count, Inf_Mach_MIN);
    info("%2d) Infinity Mach_Max --------------------: %lf", ++count, Inf_Mach_MAX);
    info("%2d) Infinity Flow Angle Alpha ------------: %lf", ++count, Inf_Alpha);
    info("%2d) Infinity Flow Angle Beta -------------: %lf", ++count, Inf_Beta);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Pressure Conditions------------------------------------------------\n");
    info("%2d) Outflow Pressure ---------------------: %lf", ++count, Outflow_Pressure);
    info("%2d) Gauge Pressure -----------------------: %lf", ++count, Gauge_Pressure);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: CFL Conditions-----------------------------------------------------\n");
    info("%2d) CFL Ramp -----------------------------: %d",  ++count, CFL_Ramp);
    info("%2d) CFL_Min ------------------------------: %lf", ++count, CFL_MIN);
    info("%2d) CFL_Max ------------------------------: %lf", ++count, CFL_MAX);
    printf("-----------------------------------------------------------------------------\n");
    printf("---INPUT: Restart Conditions-------------------------------------------------\n");
    info("%2d) Restart Input ------------------------: %d",  ++count, RestartInput);
    info("%2d) Restart Output -----------------------: %d",  ++count, RestartOutput);
    info("%2d) Restart Cycle ------------------------: %d",  ++count, RestartCycle);
    
    // Set the Material Properties
    Material_Set_Properties();
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Solver_BC_Parameters_Read(const char *filename) {
    ifstream bc_file;
    string text_line;
    string option_name;
    string::size_type position;
    int nb  = 0;
    int bid = -1;
    
    // Open the Boundary Condition File
    bc_file.open(filename, ios::in);
    if (bc_file.fail())
        error("Solver_BC_Parameters_Read: Unable to Read Boundary Condition File - %s", filename);
    
    // Start Reading the Boundary Condition File
    while (getline(bc_file, text_line)) {
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
        
        // Get the Number of Boundaries from BC File
        position = text_line.find("NUMBER_OF_BOUNDARY=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 19);
            option_name.assign(text_line);
            nb = atoi(option_name.c_str());
            break;
        }
    }
    // Close the Boundary Condition File
    bc_file.close();
    
    // Validate with number of boundaries in the mesh file
    if (nb != nBC)
        error("Solver_BC_Parameters_Read: Insufficient/Excess Boundary Conditions in  - %s", filename);
    
    // Reopen the Boundary Condition File
    bc_file.open(filename, ios::in);
    if (bc_file.fail())
        error("Solver_BC_Parameters_Read: Unable to Reopen Boundary Condition File - %s", filename);
    
    // Allocate memory to store boundary type
    bndType = new int[nBC];
    for (int i = 0; i < nBC; i++)
        bndType[i] = BC_TYPE_NONE;
    
    // Start Reading the Boundary Condition File
    while (getline(bc_file, text_line)) {
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
        
        // Get the Number of Boundaries from BC File
        position = text_line.find("BND_ID=", 0);
        if ((position != string::npos) && (position == 0)) {
            text_line.erase(0, 7);
            
            // Get the Boundary ID from the string
            bid = -1;
            position = text_line.find("=", 0);
            // Check if Boundary ID is provided
            if ((position == 0) || (position == string::npos))
                error("Solver_BC_Parameters_Read: No Boundary ID Found - %s", text_line.c_str());
            option_name.assign(text_line, 0, position);
            bid = atoi(option_name.c_str());
            
            // Validate the Boundary ID
            if ((bid <= 0) || (bid > nBC))
                error("Solver_BC_Parameters_Read: Invalid Boundary ID Found - %s", text_line.c_str());
            
            // Check if Boundary Condition is Already Updated
            if (bndType[bid-1] != BC_TYPE_NONE) {
                error("Solver_BC_Parameters_Read: Multiple Boundary ID Found - %s", text_line.c_str());
                break;
            }
            
            // Erase to get Boundary Type
            text_line.erase(0, position+1);
            
            // Finally Assign the Boundary Type
            option_name.assign(text_line);
            StringToUpperCase(option_name);
            GetOptionNameValue(option_name, bndType[bid-1], BCTypeMap);
        }
    }
    // Close the Boundary Condition File
    bc_file.close();
}

