/*******************************************************************************
 * File:        Material.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Material.h"
#include "SolverParameters.h"
#include "Trim_Utils.h"

// Material Variables
int    MaterialType;

// Non Dimensionalization Method
int    NonDimensionalMethod;

// Dimensional Properties
double Ref_R;
double Ref_MolecularWeight;
double Ref_Length;
double Ref_Area;
double Ref_Time;
double Ref_Rho;
double Ref_Pressure;
double Ref_Temperature;
double Ref_Velocity;
double Ref_Mach;

// Other Reference Properties
double Gamma;
double Ref_Alpha;
double Ref_Beta;

// Non Dimensional Properties
double NonDim_R;
double NonDim_Time;
double NonDim_Rho;
double NonDim_Pressure;
double NonDim_Temperature;

// Far Field Conditions
double Inf_Rho;
double Inf_U;
double Inf_V;
double Inf_W;
double Inf_Et;
double Inf_Pressure;
double Inf_Temperature;
double Inf_SpeedSound;
double Inf_Mach;

// Gauge Variables
double Gauge_Pressure;
double Gauge_Temperature;
double Outflow_Pressure;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Material_Init(void) {
    // Initialize
    NonDimensionalMethod = NONDIMENSIONAL_METHOD_NONE;
    MaterialType         = MATERIAL_TYPE_NONE;
    
    // Dimensional Reference Properties
    Ref_R               = 0.0;
    Ref_MolecularWeight = 0.0;
    Ref_Length          = 0.0;
    Ref_Area            = 0.0;
    Ref_Time            = 0.0;
    Ref_Rho             = 0.0;
    Ref_Pressure        = 0.0;
    Ref_Temperature     = 0.0;
    Ref_Velocity        = 0.0;
    Ref_Mach            = 0.0;
    
    // Other Reference Properties
    Gamma               = 0.0;
    Ref_Alpha           = 0.0;
    Ref_Beta            = 0.0;
    
    // Non Dimensional Properties
    NonDim_R            = 0.0;
    NonDim_Time         = 0.0;
    NonDim_Rho          = 0.0;
    NonDim_Pressure     = 0.0;
    NonDim_Temperature  = 0.0;
    
    // Far Field Conditions
    Inf_Rho             = 0.0;
    Inf_U               = 0.0;
    Inf_V               = 0.0;
    Inf_W               = 0.0;
    Inf_Et              = 0.0;
    Inf_Pressure        = 0.0;
    Inf_Temperature     = 0.0;
    Inf_SpeedSound      = 0.0;
    Inf_Mach            = 0.0;
    
    // Gauge Variables
    Gauge_Pressure      = 0.0;
    Gauge_Temperature   = 0.0;
    Outflow_Pressure    = 0.0;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Material_Set_Properties(void) {
    MaterialType       = MATERIAL_TYPE_IDEAL_GAS; // Only Support for now
    
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            Ref_MolecularWeight = 28.97;
            Ref_R               = UNIV_GAS_CONST/Ref_MolecularWeight;
            
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    NonDim_R           = 1.0/Gamma;
                    NonDim_Rho         = 1.0;
                    NonDim_Temperature = Gamma;
                    NonDim_Pressure    = NonDim_R*NonDim_Rho*NonDim_Temperature;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    NonDim_R           = 1.0/(Gamma*Ref_Mach*Ref_Mach);
                    NonDim_Rho         = 1.0;
                    NonDim_Temperature = Gamma*Ref_Mach*Ref_Mach;
                    NonDim_Pressure    = NonDim_R*NonDim_Rho*NonDim_Temperature;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    NonDim_R           = 1.0;
                    NonDim_Rho         = 1.0;
                    NonDim_Temperature = 1.0;
                    NonDim_Pressure    = NonDim_R*NonDim_Rho*NonDim_Temperature;
                    break;
                default:
                    error("Material_Set_Properties: Undefined Non-Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Material_Set_Properties: Undefined Material Type - %d", MaterialType);
            break;
    }
    
    // Gauge Variables
    switch (VariableType) {
        case VARIABLE_CONSERVATIVE:
            Gauge_Pressure    = 0.0;
            Gauge_Temperature = 0.0;
            break;
        case VARIABLE_PRIMITIVE_PUT:
            Gauge_Pressure    = NonDim_Pressure;
            Gauge_Temperature = NonDim_Temperature;
            break;
        case VARIABLE_PRIMITIVE_RUP:
            Gauge_Pressure    = NonDim_Pressure;
            Gauge_Temperature = 0.0;
            break;
        case VARIABLE_PRIMITIVE_RUT:
            Gauge_Pressure    = 0.0;
            Gauge_Temperature = NonDim_Temperature;
            break;
        default:
            error("Material_Set_Properties: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    Outflow_Pressure = 0.935*NonDim_Pressure;
    
    // Set Flow Angles in Radians
    Ref_Alpha  = Ref_Alpha * M_PI / 180.0;
    Ref_Beta   = Ref_Beta * M_PI / 180.0;
    Ref_Length = 1.0;
    Ref_Area   = 1.0;
    Ref_Time   = 1.0;
}

//------------------------------------------------------------------------------
//! Compute Far Field Condition with Mach Ramping
//------------------------------------------------------------------------------
void ComputeFarFieldCondition(int Iteration) {
    double tmpMach = 0.0;

    // Compute the Far Field Condition Ramping
    if ((Mach_Ramp > 1) && (Mach_MAX > Mach_MIN)) {
        if (Iteration < Mach_Ramp)
            tmpMach = Mach_MIN + (Mach_MAX - Mach_MIN)*(((double)Iteration)/((double)(Mach_Ramp-1)));
        else
            tmpMach = Mach_MAX;
    } else
        tmpMach = Mach_MAX;

    // Set the Far Field Conditions
    Inf_Mach        = tmpMach;
    Inf_Rho         = NonDim_Rho;
    Inf_Pressure    = NonDim_Pressure - Gauge_Pressure;
    Inf_Temperature = NonDim_Temperature - Gauge_Temperature;
    Inf_SpeedSound  = Get_SpeedSound(Inf_Temperature);
    Inf_U           = Inf_Mach*Inf_SpeedSound*cos(Ref_Alpha)*cos(Ref_Beta);
    Inf_V           = Inf_Mach*Inf_SpeedSound*sin(Ref_Alpha);
    Inf_W           = Inf_Mach*Inf_SpeedSound*sin(Ref_Beta)*cos(Ref_Alpha);
    Inf_Et          = Get_TotalEnergy(Inf_Rho, Inf_Pressure, Inf_U, Inf_V, Inf_W);
    
    // Set Other Conditions
    Outflow_Pressure -= Gauge_Pressure;
}


//------------------------------------------------------------------------------
//! Compute Equation of State Variables at Face/Edge
//------------------------------------------------------------------------------
void Compute_EOS_Variables_Face(double *Q, double nx, double ny, double nz,
        double &rho, double &pressure, double &temperature,
        double &velocity_u, double &velocity_v, double &velocity_w, double &q2,
        double &speed_sound, double &mach, double &ubar, double &total_energy, 
        double &total_enthalpy) {
    
    // Compute the EOS Based on Non Dimensionalization and Variable Type of Q
    switch (NonDimensionalMethod) {
        // Default Non Dimensionalization
        case NONDIMENSIONAL_METHOD_GENERIC:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    rho            = Q[0];
                    velocity_u     = Q[1]/rho;
                    velocity_v     = Q[2]/rho;
                    velocity_w     = Q[3]/rho;
                    total_energy   = Q[4]/rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    pressure       = (Gamma - 1.0)*(rho*total_energy - 0.5*rho*q2);
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    temperature    = pressure/(NonDim_R*rho);
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    pressure       = Q[0] + Gauge_Pressure;
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    rho            = pressure/(NonDim_R*temperature);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    pressure       = Q[0]; // Only Perturbation
                    temperature    = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    pressure       = Q[4] + Gauge_Pressure;
                    temperature    = pressure/(NonDim_R*rho);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    pressure       = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    pressure       = temperature*NonDim_R*rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    temperature    = Q[4]; // Only Perturbation
                    break;
                default:
                    error("Compute_EOS_Variables_Face: Undefined Variable Type - %d - Error-1", VariableType);
                    break;
            }
            break;
        // Briley Taylor Whitfield Non Dimensionalization
        case NONDIMENSIONAL_METHOD_BTW:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    rho            = Q[0];
                    velocity_u     = Q[1]/rho;
                    velocity_v     = Q[2]/rho;
                    velocity_w     = Q[3]/rho;
                    total_energy   = Q[4]/rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    pressure       = rho*(total_energy/(Ref_Mach*Ref_Mach) - 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    temperature    = pressure/(NonDim_R*rho);
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    pressure       = Q[0] + Gauge_Pressure;
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    rho            = pressure/(NonDim_R*temperature);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = (Ref_Mach*Ref_Mach)*(pressure/rho + 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    pressure       = Q[0]; // Only Perturbation
                    temperature    = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    pressure       = Q[4] + Gauge_Pressure;
                    temperature    = pressure/(NonDim_R*rho);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = (Ref_Mach*Ref_Mach)*(pressure/rho + 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    pressure       = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    pressure       = temperature*NonDim_R*rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = (Ref_Mach*Ref_Mach)*(pressure/rho + 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    temperature    = Q[4]; // Only Perturbation
                    break;
                default:
                    error("Compute_EOS_Variables_Face: Undefined Variable Type - %d - Error-2", VariableType);
                    break;
            }
            break;
        // LMRoe Non Dimensionalization
        case NONDIMENSIONAL_METHOD_LMROE:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    rho            = Q[0];
                    velocity_u     = Q[1]/rho;
                    velocity_v     = Q[2]/rho;
                    velocity_w     = Q[3]/rho;
                    total_energy   = Q[4]/rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    pressure       = (Gamma - 1.0)*rho*(total_energy - 0.5*q2);
                    total_enthalpy = total_energy + pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    temperature    = pressure/(NonDim_R*rho);
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    pressure       = Q[0] + Gauge_Pressure;
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    rho            = pressure/(NonDim_R*temperature);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    pressure       = Q[0]; // Only Perturbation
                    temperature    = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    pressure       = Q[4] + Gauge_Pressure;
                    temperature    = pressure/(NonDim_R*rho);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    pressure       = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    pressure       = temperature*NonDim_R*rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    ubar           = velocity_u*nx + velocity_v*ny + velocity_w*nz;
                    temperature    = Q[4]; // Only Perturbation
                    break;
                default:
                    error("Compute_EOS_Variables_Face: Undefined Variable Type - %d - Error-3", VariableType);
                    break;
            }
            break;
        default:
            error("Compute_EOS_Variables_Face: Undefined Non Dimensional Method - %d - Error-4", NonDimensionalMethod);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute Equation of State Variables at Control Volume
//------------------------------------------------------------------------------
void Compute_EOS_Variables_ControlVolume(double *Q,
        double &rho, double &pressure, double &temperature,
        double &velocity_u, double &velocity_v, double &velocity_w, double &q2,
        double &speed_sound, double &mach, double &total_energy, double &total_enthalpy) {
    
    // Compute the EOS Based on Non Dimensionalization and Variable Type of Q
    switch (NonDimensionalMethod) {
        // Default Non Dimensionalization
        case NONDIMENSIONAL_METHOD_GENERIC:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    rho            = Q[0];
                    velocity_u     = Q[1]/rho;
                    velocity_v     = Q[2]/rho;
                    velocity_w     = Q[3]/rho;
                    total_energy   = Q[4]/rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    pressure       = (Gamma - 1.0)*(rho*total_energy - 0.5*rho*q2);
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    temperature    = pressure/(NonDim_R*rho);
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    pressure       = Q[0] + Gauge_Pressure;
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    rho            = pressure/(NonDim_R*temperature);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    pressure       = Q[0]; // Only Perturbation
                    temperature    = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    pressure       = Q[4] + Gauge_Pressure;
                    temperature    = pressure/(NonDim_R*rho);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    pressure       = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    pressure       = temperature*NonDim_R*rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    temperature    = Q[4]; // Only Perturbation
                    break;
                default:
                    error("Compute_EOS_Variables_ControlVolume: Undefined Variable Type - %d - Error-1", VariableType);
                    break;
            }
            break;
        // Briley Taylor Whitfield Non Dimensionalization
        case NONDIMENSIONAL_METHOD_BTW:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    rho            = Q[0];
                    velocity_u     = Q[1]/rho;
                    velocity_v     = Q[2]/rho;
                    velocity_w     = Q[3]/rho;
                    total_energy   = Q[4]/rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    pressure       = rho*(total_energy/(Ref_Mach*Ref_Mach) - 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    temperature    = pressure/(NonDim_R*rho);
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    pressure       = Q[0] + Gauge_Pressure;
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    rho            = pressure/(NonDim_R*temperature);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = (Ref_Mach*Ref_Mach)*(pressure/rho + 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    pressure       = Q[0]; // Only Perturbation
                    temperature    = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    pressure       = Q[4] + Gauge_Pressure;
                    temperature    = pressure/(NonDim_R*rho);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = (Ref_Mach*Ref_Mach)*(pressure/rho + 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    pressure       = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    pressure       = temperature*NonDim_R*rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = (Ref_Mach*Ref_Mach)*(pressure/rho + 0.5*(Gamma - 1.0)*q2);
                    total_enthalpy = total_energy + (Gamma - 1.0)*Ref_Mach*Ref_Mach*pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    temperature    = Q[4]; // Only Perturbation
                    break;
                default:
                    error("Compute_EOS_Variables_ControlVolume: Undefined Variable Type - %d - Error-2", VariableType);
                    break;
            }
            break;
        // LMRoe Non Dimensionalization
        case NONDIMENSIONAL_METHOD_LMROE:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    rho            = Q[0];
                    velocity_u     = Q[1]/rho;
                    velocity_v     = Q[2]/rho;
                    velocity_w     = Q[3]/rho;
                    total_energy   = Q[4]/rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    pressure       = (Gamma - 1.0)*rho*(total_energy - 0.5*q2);
                    total_enthalpy = total_energy + pressure/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    temperature    = pressure/(NonDim_R*rho);
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    pressure       = Q[0] + Gauge_Pressure;
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    rho            = pressure/(NonDim_R*temperature);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    pressure       = Q[0]; // Only Perturbation
                    temperature    = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    pressure       = Q[4] + Gauge_Pressure;
                    temperature    = pressure/(NonDim_R*rho);
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    pressure       = Q[4]; // Only Perturbation
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    rho            = Q[0];
                    velocity_u     = Q[1];
                    velocity_v     = Q[2];
                    velocity_w     = Q[3];
                    temperature    = Q[4] + Gauge_Temperature;
                    pressure       = temperature*NonDim_R*rho;
                    q2             = velocity_u*velocity_u + velocity_v*velocity_v + velocity_w*velocity_w;
                    total_energy   = pressure/((Gamma - 1.0)*rho) + 0.5*q2;
                    total_enthalpy = (rho*total_energy + pressure)/rho;
                    speed_sound    = sqrt((Gamma * pressure) / rho);
                    mach           = sqrt(q2)/speed_sound;
                    temperature    = Q[4]; // Only Perturbation
                    break;
                default:
                    error("Compute_EOS_Variables_ControlVolume: Undefined Variable Type - %d - Error-3", VariableType);
                    break;
            }
            break;
        default:
            error("Compute_EOS_Variables_ControlVolume: Undefined Non Dimensional Method - %d - Error-4", NonDimensionalMethod);
            break;
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_Rho(double Pressure, double Temperature) {
    double Rho = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            Rho = (Pressure + Gauge_Pressure)/(NonDim_R*(Temperature + Gauge_Temperature));
            break;
        default:
            error("Get_Rho: Undefined Material Type - %d", MaterialType);
            break;
    }
    return Rho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_Pressure(double Rho, double Temperature) {
    double Pressure = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            Pressure = NonDim_R*Rho*(Temperature + Gauge_Temperature);
            break;
        default:
            error("Get_Pressure: Undefined Material Type - %d", MaterialType);
            break;
    }
    return Pressure;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_Pressure(double *Q) {
    double Pressure = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            switch (NonDimensionalMethod) {
                // Default Non Dimensionalization
                case NONDIMENSIONAL_METHOD_GENERIC:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            Pressure = (Gamma - 1.0)*(Q[4] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0]);
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Pressure = Q[0];
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure = Q[4];
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            Pressure = NonDim_R*Q[0]*(Q[4] + Gauge_Temperature);
                            break;
                        default:
                            error("Get_Pressure: Undefined Variable Type - %d - Error-1", VariableType);
                            break;
                    }
                    break;
                // Briley Taylor Whitfield Non Dimensionalization
                case NONDIMENSIONAL_METHOD_BTW:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            Pressure = Q[4]/(Ref_Mach*Ref_Mach) - 0.5*(Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0];
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Pressure = Q[0];
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure = Q[4];
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            Pressure = NonDim_R*Q[0]*(Q[4] + Gauge_Temperature);
                            break;
                        default:
                            error("Get_Pressure: Undefined Variable Type - %d - Error-2", VariableType);
                            break;
                    }
                    break;
                // LMRoe Non Dimensionalization
                case NONDIMENSIONAL_METHOD_LMROE:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            Pressure = (Gamma - 1.0)*(Q[4] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0]);
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Pressure = Q[0];
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure = Q[4];
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            Pressure = NonDim_R*Q[0]*(Q[4] + Gauge_Temperature);
                            break;
                        default:
                            error("Get_Pressure: Undefined Variable Type - %d - Error-3", VariableType);
                            break;
                    }
                    break;
                default:
                    error("Get_Pressure: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Get_Pressure: Undefined Material Type - %d", MaterialType);
            break;
    }
    return Pressure;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_Temperature(double Rho, double Pressure) {
    double Temperature = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            Temperature = (Pressure + Gauge_Pressure)/(NonDim_R*Rho) - Gauge_Temperature;
            break;
        default:
            error("Get_Temperature: Undefined Material Type - %d", MaterialType);
            break;
    }
    return Temperature;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_Temperature(double *Q) {
    double Temperature = 0.0, Pressure = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            switch (NonDimensionalMethod) {
                // Default Non Dimensionalization
                case NONDIMENSIONAL_METHOD_GENERIC:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            Pressure    = (Gamma - 1.0)*(Q[4] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0]);
                            Temperature = Pressure/(NonDim_R*Q[0]);
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Temperature = Q[4];
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure    = Q[4] + Gauge_Pressure;
                            Temperature = Pressure/(NonDim_R*Q[0]) - Gauge_Temperature;
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            Temperature = Q[4];
                            break;
                        default:
                            error("Get_Temperature: Undefined Variable Type - %d - Error-1", VariableType);
                            break;
                    }
                    break;
                // Briley Taylor Whitfield Non Dimensionalization
                case NONDIMENSIONAL_METHOD_BTW:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            Pressure    = Q[4]/(Ref_Mach*Ref_Mach) - 0.5*(Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0];
                            Temperature = Pressure/(NonDim_R*Q[0]);
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Temperature = Q[4];
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure    = Q[4] + Gauge_Pressure;
                            Temperature = Pressure/(NonDim_R*Q[0]) - Gauge_Temperature;
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            Temperature = Q[4];
                            break;
                        default:
                            error("Get_Temperature: Undefined Variable Type - %d - Error-2", VariableType);
                            break;
                    }
                    break;
                // LMRoe Non Dimensionalization
                case NONDIMENSIONAL_METHOD_LMROE:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            Pressure    = (Gamma - 1.0)*(Q[4] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0]);
                            Temperature = Pressure/(NonDim_R*Q[0]);
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Temperature = Q[4];
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure   = Q[4] + Gauge_Pressure;
                            Temperature = Pressure/(NonDim_R*Q[0]) - Gauge_Temperature;
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            Temperature = Q[4];
                            break;
                        default:
                            error("Get_Temperature: Undefined Variable Type - %d - Error-3", VariableType);
                            break;
                    }
                    break;
                default:
                    error("Get_Temperature: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Get_Temperature: Undefined Material Type - %d", MaterialType);
            break;
    }
    return Temperature;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_SpeedSound(double Rho, double Pressure) {
    double SpeedSound = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            SpeedSound = sqrt(Gamma*(Pressure + Gauge_Pressure)/Rho);
            break;
        default:
            error("Get_SpeedSound: Undefined Material Type - %d", MaterialType);
            break;
    }
    return SpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_SpeedSound(double Temperature) {
    double SpeedSound = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            SpeedSound = sqrt(Gamma*NonDim_R*(Temperature + Gauge_Temperature));
            break;
        default:
            error("Get_SpeedSound: Undefined Material Type - %d", MaterialType);
            break;
    }
    return SpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_SpeedSound(double Velocity_U, double Velocity_V, double Velocity_W, double Total_Enthalpy) {
    double SpeedSound = 0.0;
    double phi = 0.0;
    
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            phi = 0.5*(Gamma - 1.0)*(Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            switch (NonDimensionalMethod) {
                // Default Non Dimensionalization
                case NONDIMENSIONAL_METHOD_GENERIC:
                    SpeedSound = (Gamma - 1.0) * Total_Enthalpy - phi;
                    SpeedSound = sqrt(SpeedSound);
                    break;
                // Briley Taylor Whitfield Non Dimensionalization
                case NONDIMENSIONAL_METHOD_BTW:
                    SpeedSound = Total_Enthalpy / (Ref_Mach * Ref_Mach) - phi;
                    SpeedSound = sqrt(SpeedSound);
                    break;
                // LMRoe Non Dimensionalization
                case NONDIMENSIONAL_METHOD_LMROE:
                    SpeedSound = (Gamma - 1.0) * Total_Enthalpy - phi;
                    SpeedSound = sqrt(SpeedSound);
                    break;
                default:
                    error("Get_SpeedSound: Undefined Non-Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Get_SpeedSound: Undefined Material Type - %d", MaterialType);
            break;
    }
    return SpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_SpeedSoundSquare(double Velocity_U, double Velocity_V, double Velocity_W, double Total_Enthalpy) {
    double SpeedSoundSquare = 0.0;
    double phi = 0.0;
    
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            phi = 0.5*(Gamma - 1.0)*(Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            switch (NonDimensionalMethod) {
                // Default Non Dimensionalization
                case NONDIMENSIONAL_METHOD_GENERIC:
                    SpeedSoundSquare = (Gamma - 1.0) * Total_Enthalpy - phi;
                    break;
                // Briley Taylor Whitfield Non Dimensionalization
                case NONDIMENSIONAL_METHOD_BTW:
                    SpeedSoundSquare = Total_Enthalpy / (Ref_Mach * Ref_Mach) - phi;
                    break;
                // LMRoe Non Dimensionalization
                case NONDIMENSIONAL_METHOD_LMROE:
                    SpeedSoundSquare = (Gamma - 1.0) * Total_Enthalpy - phi;
                    break;
                default:
                    error("Get_SpeedSoundSquare: Undefined Non-Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Get_SpeedSoundSquare: Undefined Material Type - %d", MaterialType);
            break;
    }
    return SpeedSoundSquare;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_TotalEnergy(double Rho, double Pressure, double Velocity_U, double Velocity_V, double Velocity_W) {
    double TotalEnergy = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            switch (NonDimensionalMethod) {
                // Default Non Dimensionalization
                case NONDIMENSIONAL_METHOD_GENERIC:
                    TotalEnergy = Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W;
                    TotalEnergy = (Pressure + Gauge_Pressure)/((Gamma - 1.0)*Rho) + 0.5*TotalEnergy;
                    break;
                // Briley Taylor Whitfield Non Dimensionalization
                case NONDIMENSIONAL_METHOD_BTW:
                    TotalEnergy = Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W;
                    TotalEnergy = (Ref_Mach*Ref_Mach)*((Pressure + Gauge_Pressure)/Rho + 0.5*(Gamma - 1.0)*TotalEnergy);
                    break;
                // LMRoe Non Dimensionalization
                case NONDIMENSIONAL_METHOD_LMROE:
                    TotalEnergy = Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W;
                    TotalEnergy = (Pressure + Gauge_Pressure)/((Gamma - 1.0)*Rho) + 0.5*TotalEnergy;
                    break;
                default:
                    error("Get_TotalEnergy: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Get_TotalEnergy: Undefined Material Type - %d", MaterialType);
            break;
    }
    return TotalEnergy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Get_Mach(double *Q) {
    double q2 = 0.0, var = 0.0, Pressure = 0.0, Mach = 0.0;
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            switch (NonDimensionalMethod) {
                // Default Non Dimensionalization
                case NONDIMENSIONAL_METHOD_GENERIC:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            q2       = (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/(Q[0]*Q[0]);
                            Pressure = (Gamma - 1.0)*(Q[4] - 0.5*Q[0]*q2);
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Pressure = Q[0] + Gauge_Pressure;
                            var      = Q[4] + Gauge_Temperature;
                            var      = Pressure/(NonDim_R*var);
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*var/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure = Q[4] + Gauge_Pressure;
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            var      = Q[4] + Gauge_Temperature;
                            Pressure = NonDim_R*Q[0]*var;
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        default:
                            error("Get_Mach: Undefined Variable Type - %d - Error-1", VariableType);
                            break;
                    }
                    break;
                // Briley Taylor Whitfield Non Dimensionalization
                case NONDIMENSIONAL_METHOD_BTW:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            q2       = (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/(Q[0]*Q[0]);
                            Pressure = Q[4]/(Ref_Mach*Ref_Mach) - 0.5*(Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0];
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Pressure = Q[0] + Gauge_Pressure;
                            var      = Q[4] + Gauge_Temperature;
                            var      = Pressure/(NonDim_R*var);
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*var/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure = Q[4] + Gauge_Pressure;
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            var      = Q[4] + Gauge_Temperature;
                            Pressure = NonDim_R*Q[0]*var;
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        default:
                            error("Get_Mach: Undefined Variable Type - %d - Error-2", VariableType);
                            break;
                    }
                    break;
                // LMRoe Non Dimensionalization
                case NONDIMENSIONAL_METHOD_LMROE:
                    switch (VariableType) {
                        // Conservative Variable Formulation
                        case VARIABLE_CONSERVATIVE:
                            q2       = (Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/(Q[0]*Q[0]);
                            Pressure = (Gamma - 1.0)*(Q[4] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/Q[0]);
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Pressure Velocity Temperature
                        case VARIABLE_PRIMITIVE_PUT:
                            Pressure = Q[0] + Gauge_Pressure;
                            var      = Q[4] + Gauge_Temperature;
                            var      = Pressure/(NonDim_R*var);
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*var/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Density Velocity Pressure
                        case VARIABLE_PRIMITIVE_RUP:
                            Pressure = Q[4] + Gauge_Pressure;
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        // Primitive Variable Formulation Density Velocity Temperature
                        case VARIABLE_PRIMITIVE_RUT:
                            var      = Q[4] + Gauge_Temperature;
                            Pressure = NonDim_R*Q[0]*var;
                            q2       = Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3];
                            Mach     = sqrt(q2*Q[0]/(Gamma * Pressure));
                            break;
                        default:
                            error("Get_Mach: Undefined Variable Type - %d - Error-3", VariableType);
                            break;
                    }
                    break;
                default:
                    error("Get_Mach: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Get_Mach: Undefined Material Type - %d", MaterialType);
            break;
    }
    return Mach;
}

//------------------------------------------------------------------------------
//! Transformation: M01 = dQ/dq1
//------------------------------------------------------------------------------
void Compute_Transformation_Conservative_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            phi = 0.5*(Gamma - 1.0)*(Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            
            Matrix[0][0] = 1.0;
            Matrix[0][1] = 0.0;
            Matrix[0][2] = 0.0;
            Matrix[0][3] = 0.0;
            Matrix[0][4] = 0.0;

            Matrix[1][0] = Velocity_U;
            Matrix[1][1] = Rho;
            Matrix[1][2] = 0.0;
            Matrix[1][3] = 0.0;
            Matrix[1][4] = 0.0;

            Matrix[2][0] = Velocity_V;
            Matrix[2][1] = 0.0;
            Matrix[2][2] = Rho;
            Matrix[2][3] = 0.0;
            Matrix[2][4] = 0.0;

            Matrix[3][0] = Velocity_W;
            Matrix[3][1] = 0.0;
            Matrix[3][2] = 0.0;
            Matrix[3][3] = Rho;
            Matrix[3][4] = 0.0;

            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[4][0] = phi/(Gamma - 1.0);
                    Matrix[4][1] = Rho * Velocity_U;
                    Matrix[4][2] = Rho * Velocity_V;
                    Matrix[4][3] = Rho * Velocity_W;
                    Matrix[4][4] = 1.0/(Gamma - 1.0);
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[4][0] = phi * Ref_Mach * Ref_Mach;
                    Matrix[4][1] = (Gamma - 1.0) * Rho * Velocity_U * Ref_Mach * Ref_Mach;
                    Matrix[4][2] = (Gamma - 1.0) * Rho * Velocity_V * Ref_Mach * Ref_Mach;
                    Matrix[4][3] = (Gamma - 1.0) * Rho * Velocity_W * Ref_Mach * Ref_Mach;
                    Matrix[4][4] = Ref_Mach * Ref_Mach;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[4][0] = phi/(Gamma - 1.0);
                    Matrix[4][1] = Rho * Velocity_U;
                    Matrix[4][2] = Rho * Velocity_V;
                    Matrix[4][3] = Rho * Velocity_W;
                    Matrix[4][4] = 1.0/(Gamma - 1.0);
                    break;
                default:
                    error("Compute_Transformation_Conservative_To_RUP: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_Conservative_To_RUP: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M10 = dq1/dQ
//------------------------------------------------------------------------------
void Compute_Transformation_RUP_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            phi = 0.5*(Gamma - 1.0)*(Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            
            Matrix[0][0] = 1.0;
            Matrix[0][1] = 0.0;
            Matrix[0][2] = 0.0;
            Matrix[0][3] = 0.0;
            Matrix[0][4] = 0.0;

            Matrix[1][0] = -Velocity_U/Rho;
            Matrix[1][1] = 1.0/Rho;
            Matrix[1][2] = 0.0;
            Matrix[1][3] = 0.0;
            Matrix[1][4] = 0.0;

            Matrix[2][0] = -Velocity_V/Rho;
            Matrix[2][1] = 0.0;
            Matrix[2][2] = 1.0/Rho;
            Matrix[2][3] = 0.0;
            Matrix[2][4] = 0.0;

            Matrix[3][0] = -Velocity_W/Rho;
            Matrix[3][1] = 0.0;
            Matrix[3][2] = 0.0;
            Matrix[3][3] = 1.0/Rho;
            Matrix[3][4] = 0.0;

            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[4][0] = phi;
                    Matrix[4][1] = -Velocity_U * (Gamma - 1.0);
                    Matrix[4][2] = -Velocity_V * (Gamma - 1.0);
                    Matrix[4][3] = -Velocity_W * (Gamma - 1.0);
                    Matrix[4][4] = (Gamma - 1.0);
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[4][0] = phi;
                    Matrix[4][1] = -Velocity_U * (Gamma - 1.0);
                    Matrix[4][2] = -Velocity_V * (Gamma - 1.0);
                    Matrix[4][3] = -Velocity_W * (Gamma - 1.0);
                    Matrix[4][4] = 1.0/(Ref_Mach * Ref_Mach);
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[4][0] = phi;
                    Matrix[4][1] = -Velocity_U * (Gamma - 1.0);
                    Matrix[4][2] = -Velocity_V * (Gamma - 1.0);
                    Matrix[4][3] = -Velocity_W * (Gamma - 1.0);
                    Matrix[4][4] = (Gamma - 1.0);
                    break;
                default:
                    error("Compute_Transformation_RUP_To_Conservative: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUP_To_Conservative: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M02 = dQ/dq2
//------------------------------------------------------------------------------
void Compute_Transformation_Conservative_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    double c2  = 0.0;
    double q2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            q2  = (Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            phi = 0.5*(Gamma - 1.0)*q2;
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/c2;

                    Matrix[1][0] = (Gamma*Velocity_U)/c2;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = -(Rho*Velocity_U)/c2;

                    Matrix[2][0] = (Gamma*Velocity_V)/c2;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = -(Rho*Velocity_V)/c2;

                    Matrix[3][0] = (Gamma*Velocity_W)/c2;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = -(Rho*Velocity_W)/c2;
                    
                    Matrix[4][0] = 1.0/(Gamma - 1.0) + (Gamma*q2)/(2*c2);
                    Matrix[4][1] = Rho*Velocity_U;
                    Matrix[4][2] = Rho*Velocity_V;
                    Matrix[4][3] = Rho*Velocity_W;
                    Matrix[4][4] = -(Rho*q2)/(2*c2);
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/(c2*Ref_Mach*Ref_Mach);

                    Matrix[1][0] = (Gamma*Velocity_U)/c2;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = -(Rho*Velocity_U)/(c2*Ref_Mach*Ref_Mach);

                    Matrix[2][0] = (Gamma*Velocity_V)/c2;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = -(Rho*Velocity_V)/(c2*Ref_Mach*Ref_Mach);

                    Matrix[3][0] = (Gamma*Velocity_W)/c2;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = -(Rho*Velocity_W)/(c2*Ref_Mach*Ref_Mach);
                    
                    Matrix[4][0] = (1.0 + (Gamma*phi)/c2)*Ref_Mach*Ref_Mach;
                    Matrix[4][1] = (Gamma - 1.0)*Rho*Velocity_U*Ref_Mach*Ref_Mach;
                    Matrix[4][2] = (Gamma - 1.0)*Rho*Velocity_V*Ref_Mach*Ref_Mach;
                    Matrix[4][3] = (Gamma - 1.0)*Rho*Velocity_W*Ref_Mach*Ref_Mach;
                    Matrix[4][4] = -(Rho*phi)/c2;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -(Gamma*Rho)/c2;

                    Matrix[1][0] = (Gamma*Velocity_U)/c2;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = -(Gamma*Rho*Velocity_U)/c2;

                    Matrix[2][0] = (Gamma*Velocity_V)/c2;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = -(Gamma*Rho*Velocity_V)/c2;

                    Matrix[3][0] = (Gamma*Velocity_W)/c2;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = -(Gamma*Rho*Velocity_W)/c2;
                    
                    Matrix[4][0] = 1.0/(Gamma - 1.0) + (Gamma*q2)/(2*c2);
                    Matrix[4][1] = Rho*Velocity_U;
                    Matrix[4][2] = Rho*Velocity_V;
                    Matrix[4][3] = Rho*Velocity_W;
                    Matrix[4][4] = -(Gamma*Rho*q2)/(2*c2);
                    break;
                default:
                    error("Compute_Transformation_Conservative_To_PUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_Conservative_To_PUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M20 = dq2/dQ
//------------------------------------------------------------------------------
void Compute_Transformation_PUT_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    double c2  = 0.0;
    double q2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            q2  = (Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            phi = 0.5*(Gamma - 1.0)*q2;
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = phi;
                    Matrix[0][1] = -(Gamma - 1.0)*Velocity_U;
                    Matrix[0][2] = -(Gamma - 1.0)*Velocity_V;
                    Matrix[0][3] = -(Gamma - 1.0)*Velocity_W;
                    Matrix[0][4] = Gamma - 1.0;

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma*phi - c2)/Rho;
                    Matrix[4][1] = -Gamma*(Gamma - 1.0)*Velocity_U/Rho;
                    Matrix[4][2] = -Gamma*(Gamma - 1.0)*Velocity_V/Rho;
                    Matrix[4][3] = -Gamma*(Gamma - 1.0)*Velocity_W/Rho;
                    Matrix[4][4] = (Gamma*(Gamma - 1.0))/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = phi;
                    Matrix[0][1] = -(Gamma - 1.0)*Velocity_U;
                    Matrix[0][2] = -(Gamma - 1.0)*Velocity_V;
                    Matrix[0][3] = -(Gamma - 1.0)*Velocity_W;
                    Matrix[0][4] = 1.0/(Ref_Mach*Ref_Mach);

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = ((Gamma*phi - c2)/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][1] = -(Gamma*(Gamma - 1.0)*Velocity_U/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][2] = -(Gamma*(Gamma - 1.0)*Velocity_V/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][3] = -(Gamma*(Gamma - 1.0)*Velocity_W/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][4] = Gamma/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = phi;
                    Matrix[0][1] = -(Gamma - 1.0)*Velocity_U;
                    Matrix[0][2] = -(Gamma - 1.0)*Velocity_V;
                    Matrix[0][3] = -(Gamma - 1.0)*Velocity_W;
                    Matrix[0][4] = Gamma - 1.0;

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma*phi - c2)/(Gamma*Rho);
                    Matrix[4][1] = -(Gamma - 1.0)*Velocity_U/Rho;
                    Matrix[4][2] = -(Gamma - 1.0)*Velocity_V/Rho;
                    Matrix[4][3] = -(Gamma - 1.0)*Velocity_W/Rho;
                    Matrix[4][4] = (Gamma - 1.0)/Rho;
                    break;
                default:
                    error("Compute_Transformation_PUT_To_Conservative: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUT_To_Conservative: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M03 = dQ/dq3
//------------------------------------------------------------------------------
void Compute_Transformation_Conservative_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    double c2  = 0.0;
    double q2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            q2  = (Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            phi = 0.5*(Gamma - 1.0)*q2;
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = Velocity_U/c2;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = -(Rho*Velocity_U)/Gamma;

                    Matrix[2][0] = Velocity_V/c2;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = -(Rho*Velocity_V)/Gamma;

                    Matrix[3][0] = Velocity_W/c2;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = -(Rho*Velocity_W)/Gamma;
                    
                    Matrix[4][0] = phi/((Gamma - 1.0)*c2) + 1.0/(Gamma - 1.0);
                    Matrix[4][1] = Rho*Velocity_U;
                    Matrix[4][2] = Rho*Velocity_V;
                    Matrix[4][3] = Rho*Velocity_W;
                    Matrix[4][4] = -(Rho*phi)/(Gamma*(Gamma - 1.0));
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = Velocity_U/c2;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = -(Rho*Velocity_U)/Gamma;

                    Matrix[2][0] = Velocity_V/c2;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = -(Rho*Velocity_V)/Gamma;

                    Matrix[3][0] = Velocity_W/c2;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = -(Rho*Velocity_W)/Gamma;
                    
                    Matrix[4][0] = (1.0 + phi/c2)*Ref_Mach*Ref_Mach;
                    Matrix[4][1] = (Gamma - 1.0)*Rho*Velocity_U*Ref_Mach*Ref_Mach;
                    Matrix[4][2] = (Gamma - 1.0)*Rho*Velocity_V*Ref_Mach*Ref_Mach;
                    Matrix[4][3] = (Gamma - 1.0)*Rho*Velocity_W*Ref_Mach*Ref_Mach;
                    Matrix[4][4] = -((Rho*phi)/Gamma)*Ref_Mach*Ref_Mach;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = Velocity_U/c2;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = -(Rho*Velocity_U)/Gamma;

                    Matrix[2][0] = Velocity_V/c2;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = -(Rho*Velocity_V)/Gamma;

                    Matrix[3][0] = Velocity_W/c2;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = -(Rho*Velocity_W)/Gamma;
                    
                    Matrix[4][0] = phi/((Gamma - 1.0)*c2) + 1.0/(Gamma - 1.0);
                    Matrix[4][1] = Rho*Velocity_U;
                    Matrix[4][2] = Rho*Velocity_V;
                    Matrix[4][3] = Rho*Velocity_W;
                    Matrix[4][4] = -(Rho*phi)/(Gamma*(Gamma - 1.0));
                    break;
                default:
                    error("Compute_Transformation_Conservative_To_PUS: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_Conservative_To_PUS: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M30 = dq3/dQ
//------------------------------------------------------------------------------
void Compute_Transformation_PUS_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    double c2  = 0.0;
    double q2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            q2  = (Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            phi = 0.5*(Gamma - 1.0)*q2;
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = phi;
                    Matrix[0][1] = -(Gamma - 1.0)*Velocity_U;
                    Matrix[0][2] = -(Gamma - 1.0)*Velocity_V;
                    Matrix[0][3] = -(Gamma - 1.0)*Velocity_W;
                    Matrix[0][4] = Gamma - 1.0;

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = Gamma*(phi - c2)/(Rho*c2);
                    Matrix[4][1] = -Gamma*(Gamma - 1.0)*Velocity_U/(Rho*c2);
                    Matrix[4][2] = -Gamma*(Gamma - 1.0)*Velocity_V/(Rho*c2);
                    Matrix[4][3] = -Gamma*(Gamma - 1.0)*Velocity_W/(Rho*c2);
                    Matrix[4][4] = (Gamma*(Gamma - 1.0))/(Rho*c2);
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = phi;
                    Matrix[0][1] = -(Gamma - 1.0)*Velocity_U;
                    Matrix[0][2] = -(Gamma - 1.0)*Velocity_V;
                    Matrix[0][3] = -(Gamma - 1.0)*Velocity_W;
                    Matrix[0][4] = 1.0/(Ref_Mach*Ref_Mach);

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = Gamma*(phi - c2)/(Rho*c2);
                    Matrix[4][1] = -Gamma*(Gamma - 1.0)*Velocity_U/(Rho*c2);
                    Matrix[4][2] = -Gamma*(Gamma - 1.0)*Velocity_V/(Rho*c2);
                    Matrix[4][3] = -Gamma*(Gamma - 1.0)*Velocity_W/(Rho*c2);
                    Matrix[4][4] = Gamma/(Rho*c2*Ref_Mach*Ref_Mach);
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = phi;
                    Matrix[0][1] = -(Gamma - 1.0)*Velocity_U;
                    Matrix[0][2] = -(Gamma - 1.0)*Velocity_V;
                    Matrix[0][3] = -(Gamma - 1.0)*Velocity_W;
                    Matrix[0][4] = Gamma - 1.0;

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = Gamma*(phi - c2)/(Rho*c2);
                    Matrix[4][1] = -Gamma*(Gamma - 1.0)*Velocity_U/(Rho*c2);
                    Matrix[4][2] = -Gamma*(Gamma - 1.0)*Velocity_V/(Rho*c2);
                    Matrix[4][3] = -Gamma*(Gamma - 1.0)*Velocity_W/(Rho*c2);
                    Matrix[4][4] = (Gamma*(Gamma - 1.0))/(Rho*c2);
                    break;
                default:
                    error("Compute_Transformation_PUS_To_Conservative: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUS_To_Conservative: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M04 = dQ/dq4
//------------------------------------------------------------------------------
void Compute_Transformation_Conservative_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    double c2  = 0.0;
    double q2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            q2  = (Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            phi = 0.5*(Gamma - 1.0)*q2;
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = Velocity_U;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = Velocity_V;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = Velocity_W;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (c2  + Gamma*phi)/(Gamma*(Gamma - 1.0));
                    Matrix[4][1] = Rho*Velocity_U;
                    Matrix[4][2] = Rho*Velocity_V;
                    Matrix[4][3] = Rho*Velocity_W;
                    Matrix[4][4] = Rho/(Gamma*(Gamma - 1.0));
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = Velocity_U;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = Velocity_V;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = Velocity_W;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (c2/Gamma + phi)*Ref_Mach*Ref_Mach;
                    Matrix[4][1] = (Gamma - 1.0)*Rho*Velocity_U*Ref_Mach*Ref_Mach;
                    Matrix[4][2] = (Gamma - 1.0)*Rho*Velocity_V*Ref_Mach*Ref_Mach;
                    Matrix[4][3] = (Gamma - 1.0)*Rho*Velocity_W*Ref_Mach*Ref_Mach;
                    Matrix[4][4] = Rho/Gamma;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = Velocity_U;
                    Matrix[1][1] = Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = Velocity_V;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = Velocity_W;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (c2  + Gamma*phi)/(Gamma*(Gamma - 1.0));
                    Matrix[4][1] = Rho*Velocity_U;
                    Matrix[4][2] = Rho*Velocity_V;
                    Matrix[4][3] = Rho*Velocity_W;
                    Matrix[4][4] = Rho/(Gamma - 1.0);
                    break;
                default:
                    error("Compute_Transformation_Conservative_To_RUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_Conservative_To_RUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M40 = dq4/dQ
//------------------------------------------------------------------------------
void Compute_Transformation_RUT_To_Conservative(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double phi = 0.0;
    double c2  = 0.0;
    double q2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            q2  = (Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W);
            phi = 0.5*(Gamma - 1.0)*q2;
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma*phi - c2)/Rho;
                    Matrix[4][1] = -Gamma*(Gamma - 1.0)*Velocity_U/Rho;
                    Matrix[4][2] = -Gamma*(Gamma - 1.0)*Velocity_V/Rho;
                    Matrix[4][3] = -Gamma*(Gamma - 1.0)*Velocity_W/Rho;
                    Matrix[4][4] = (Gamma*(Gamma - 1.0))/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = ((Gamma*phi - c2)/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][1] = -(Gamma*(Gamma - 1.0)*Velocity_U/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][2] = -(Gamma*(Gamma - 1.0)*Velocity_V/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][3] = -(Gamma*(Gamma - 1.0)*Velocity_W/Rho)*Ref_Mach*Ref_Mach;
                    Matrix[4][4] = Gamma/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = -Velocity_U/Rho;
                    Matrix[1][1] = 1.0/Rho;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = -Velocity_V/Rho;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0/Rho;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = -Velocity_W/Rho;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0/Rho;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma*phi - c2)/(Gamma*Rho);
                    Matrix[4][1] = -(Gamma - 1.0)*Velocity_U/Rho;
                    Matrix[4][2] = -(Gamma - 1.0)*Velocity_V/Rho;
                    Matrix[4][3] = -(Gamma - 1.0)*Velocity_W/Rho;
                    Matrix[4][4] = (Gamma - 1.0)/Rho;
                    break;
                default:
                    error("Compute_Transformation_RUT_To_Conservative: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUT_To_Conservative: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M12 = dq1/dq2
//------------------------------------------------------------------------------
void Compute_Transformation_RUP_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) { 
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/c2;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 1.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 0.0;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/(c2*Ref_Mach*Ref_Mach);

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 1.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 0.0;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Gamma*Rho/c2;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 1.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 0.0;
                    break;
                default:
                    error("Compute_Transformation_RUP_To_PUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUP_To_PUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M21 = dq2/dq1
//------------------------------------------------------------------------------
void Compute_Transformation_PUT_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 0.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 1.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -c2/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 0.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 1.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -(c2*Ref_Mach*Ref_Mach)/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma*Ref_Mach*Ref_Mach/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 0.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 1.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -c2/(Gamma*Rho);
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0/Rho;
                    break;
                default:
                    error("Compute_Transformation_PUT_To_RUP: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUT_To_RUP: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M13 = dq1/dq3
//------------------------------------------------------------------------------
void Compute_Transformation_RUP_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 1.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 0.0;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 1.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 0.0;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 1.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 0.0;
                    break;
                default:
                    error("Compute_Transformation_RUP_To_PUS: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUP_To_PUS: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M31 = dq3/dq1
//------------------------------------------------------------------------------
void Compute_Transformation_PUS_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 0.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 1.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -Gamma/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/(c2*Rho);
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 0.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 1.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -Gamma/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/(c2*Rho);
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 0.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 1.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -Gamma/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/(c2*Rho);
                    break;
                default:
                    error("Compute_Transformation_PUS_To_RUP: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUS_To_RUP: Undefined Material Type - %d", MaterialType);
            break;
    }
}


//------------------------------------------------------------------------------
//! Transformation: M14 = dq1/dq4
//------------------------------------------------------------------------------
void Compute_Transformation_RUP_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = c2/Gamma;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Rho/Gamma;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = c2/Gamma;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Rho/(Gamma*Ref_Mach*Ref_Mach);
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = c2/Gamma;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Rho;
                    break;
                default:
                    error("Compute_Transformation_RUP_To_RUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUP_To_RUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M41 = dq4/dq1
//------------------------------------------------------------------------------
void Compute_Transformation_RUT_To_RUP(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -c2/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -(c2*Ref_Mach*Ref_Mach)/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = (Gamma*Ref_Mach*Ref_Mach)/Rho;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -c2/(Gamma*Rho);
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0/Rho;
                    break;
                default:
                    error("Compute_Transformation_RUT_To_RUP: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUT_To_RUP: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M23 = dq2/dq3
//------------------------------------------------------------------------------
void Compute_Transformation_PUT_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma - 1.0)/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = c2/Gamma;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma - 1.0)*Ref_Mach*Ref_Mach/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = c2*Ref_Mach*Ref_Mach/Gamma;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma - 1.0)/(Gamma*Rho);
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = c2/(Gamma*Gamma);
                    break;
                default:
                    error("Compute_Transformation_PUT_To_PUS: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUT_To_PUS: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M32 = dq3/dq2
//------------------------------------------------------------------------------
void Compute_Transformation_PUS_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -(Gamma - 1.0)*Gamma/(c2*Rho);
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/c2;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -(Gamma - 1.0)*Gamma/(c2*Rho);
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/(c2*Ref_Mach*Ref_Mach);
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = 0.0;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = -(Gamma - 1.0)*Gamma/(c2*Rho);
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma*Gamma/c2;
                    break;
                default:
                    error("Compute_Transformation_PUS_To_PUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUS_To_PUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M24 = dq2/dq4
//------------------------------------------------------------------------------
void Compute_Transformation_PUT_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = c2/Gamma;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 0.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = c2/Gamma;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = Rho/(Gamma*Ref_Mach*Ref_Mach);

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 0.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = c2/Gamma;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = Rho;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 0.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0;
                    break;
                default:
                    error("Compute_Transformation_PUT_To_RUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUT_To_RUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M42 = dq4/dq2
//------------------------------------------------------------------------------
void Compute_Transformation_RUT_To_PUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/c2;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 0.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/(c2*Ref_Mach*Ref_Mach);

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 0.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = Gamma/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Gamma*Rho/c2;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = 0.0;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0;
                    break;
                default:
                    error("Compute_Transformation_RUT_To_PUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUT_To_PUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M34 = dq3/dq4
//------------------------------------------------------------------------------
void Compute_Transformation_PUS_To_RUT(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = c2/Gamma;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (1.0 - Gamma)/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0/c2;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = c2/Gamma;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = Rho/(Gamma*Ref_Mach*Ref_Mach);

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (1.0 - Gamma)/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = 1.0/(c2*Ref_Mach*Ref_Mach);
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = c2/Gamma;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = Rho;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (1.0 - Gamma)/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = Gamma/c2;
                    break;
                default:
                    error("Compute_Transformation_PUS_To_RUT: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_PUS_To_RUT: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation: M43 = dq4/dq3
//------------------------------------------------------------------------------
void Compute_Transformation_RUT_To_PUS(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    double c2  = 0.0;
    
    // Matrix
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEAL_GAS:
            c2  = SpeedSound*SpeedSound;
            switch (NonDimensionalMethod) {
                case NONDIMENSIONAL_METHOD_GENERIC:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma - 1.0)/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = c2/Gamma;
                    break;
                case NONDIMENSIONAL_METHOD_BTW:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma - 1.0)*Ref_Mach*Ref_Mach/Rho;
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = c2*Ref_Mach*Ref_Mach/Gamma;
                    break;
                case NONDIMENSIONAL_METHOD_LMROE:
                    Matrix[0][0] = 1.0/c2;
                    Matrix[0][1] = 0.0;
                    Matrix[0][2] = 0.0;
                    Matrix[0][3] = 0.0;
                    Matrix[0][4] = -Rho/Gamma;

                    Matrix[1][0] = 0.0;
                    Matrix[1][1] = 1.0;
                    Matrix[1][2] = 0.0;
                    Matrix[1][3] = 0.0;
                    Matrix[1][4] = 0.0;

                    Matrix[2][0] = 0.0;
                    Matrix[2][1] = 0.0;
                    Matrix[2][2] = 1.0;
                    Matrix[2][3] = 0.0;
                    Matrix[2][4] = 0.0;

                    Matrix[3][0] = 0.0;
                    Matrix[3][1] = 0.0;
                    Matrix[3][2] = 0.0;
                    Matrix[3][3] = 1.0;
                    Matrix[3][4] = 0.0;
                    
                    Matrix[4][0] = (Gamma - 1.0)/(Gamma*Rho);
                    Matrix[4][1] = 0.0;
                    Matrix[4][2] = 0.0;
                    Matrix[4][3] = 0.0;
                    Matrix[4][4] = c2/(Gamma*Gamma);
                    break;
                default:
                    error("Compute_Transformation_RUT_To_PUS: Undefined Non Dimensional Method - %d", NonDimensionalMethod);
                    break;
            }
            break;
        default:
            error("Compute_Transformation_RUT_To_PUS: Undefined Material Type - %d", MaterialType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation Matrix: Mpr = dqp/dqr
//------------------------------------------------------------------------------
void Compute_Transformation_Matrix(int VarType1, int VarType2, double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double SpeedSound, double **Matrix) {
    // Check the Variable Type and Return the proper Transformation Matrix
    // Return Identity
    if (VarType1 == VarType2) {
        for (int i = 0; i < NEQUATIONS; i++) {
            for (int j = 0; j < NEQUATIONS; j++) {
                if (i == j)
                    Matrix[i][j] = 1.0;
                else
                    Matrix[i][j] = 0.0;
            }
        }
    } else {
        switch (VarType1) {
            case VARIABLE_CONSERVATIVE:
                switch (VarType2) {
                    case VARIABLE_PRIMITIVE_PUT:
                        Compute_Transformation_Conservative_To_PUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUP:
                        Compute_Transformation_Conservative_To_RUP(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUT:
                        Compute_Transformation_Conservative_To_RUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_PUS:
                        Compute_Transformation_Conservative_To_PUS(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    default:
                        error("Compute_Transformation_Matrix: Undefined Variable Type - %d - Error-1", VarType2);
                        break;
                }
                break;
            case VARIABLE_PRIMITIVE_PUT:
                switch (VarType2) {
                    case VARIABLE_CONSERVATIVE:
                        Compute_Transformation_PUT_To_Conservative(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUP:
                        Compute_Transformation_PUT_To_RUP(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUT:
                        Compute_Transformation_PUT_To_RUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_PUS:
                        Compute_Transformation_PUT_To_PUS(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    default:
                        error("Compute_Transformation_Matrix: Undefined Variable Type - %d - Error-2", VarType2);
                        break;
                }
                break;
            case VARIABLE_PRIMITIVE_RUP:
                switch (VarType2) {
                    case VARIABLE_CONSERVATIVE:
                        Compute_Transformation_RUP_To_Conservative(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_PUT:
                        Compute_Transformation_RUP_To_PUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUT:
                        Compute_Transformation_RUP_To_RUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_PUS:
                        Compute_Transformation_RUP_To_PUS(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    default:
                        error("Compute_Transformation_Matrix: Undefined Variable Type - %d - Error-3", VarType2);
                        break;
                }
                break;
            case VARIABLE_PRIMITIVE_RUT:
                switch (VarType2) {
                    case VARIABLE_CONSERVATIVE:
                        Compute_Transformation_RUT_To_Conservative(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_PUT:
                        Compute_Transformation_RUT_To_PUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUP:
                        Compute_Transformation_RUT_To_RUP(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_PUS:
                        Compute_Transformation_RUT_To_PUS(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    default:
                        error("Compute_Transformation_Matrix: Undefined Variable Type - %d - Error-4", VarType2);
                        break;
                }
                break;
            case VARIABLE_PRIMITIVE_PUS:
                switch (VarType2) {
                    case VARIABLE_CONSERVATIVE:
                        Compute_Transformation_PUS_To_Conservative(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_PUT:
                        Compute_Transformation_PUS_To_PUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUP:
                        Compute_Transformation_PUS_To_RUP(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    case VARIABLE_PRIMITIVE_RUT:
                        Compute_Transformation_PUS_To_RUT(Rho, Velocity_U, Velocity_V, Velocity_W, SpeedSound, Matrix);
                        break;
                    default:
                        error("Compute_Transformation_Matrix: Undefined Variable Type - %d - Error-5", VarType2);
                        break;
                }
                break;
            default:
                error("Compute_Transformation_Matrix: Undefined Variable Type - %d - Error-6", VarType1);
                break;
        }
    }
}

