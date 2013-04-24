/*******************************************************************************
 * File:        Material.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Material.h"
#include "SolverParameters.h"
#include "Trim_Utils.h"
#include "EOS.h"

// Material Variables
int    MaterialType;
int    MaterialCompType;
int    MaterialEOS_IO_Type;
char   MaterialName[256];

// Dimensional Properties (SI Units)
double Ref_Rho;
double Ref_Pressure;
double Ref_Temperature;
double Ref_Velocity;
double Ref_TotalEnergy;
double Ref_Length;
double Ref_Time;
double Ref_Mach;

// Other Reference Properties
double Gamma;

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
double Inf_Alpha;
double Inf_Beta;
// Infinity Mach Ramping
int    Inf_Mach_Ramp;
double Inf_Mach_MAX;
double Inf_Mach_MIN;

// Gauge Variables
double Gauge_Pressure;
double Outflow_Pressure;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Material_Init(void) {
    // Initialize
    MaterialType        = MATERIAL_TYPE_NONE;
    MaterialCompType    = MATERIAL_COMP_TYPE_NONE;
    MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_NONE;
    str_blank(MaterialName);
    
    // Dimensional Reference Properties
    Ref_Rho             = 0.0;
    Ref_Pressure        = 0.0;
    Ref_Temperature     = 0.0;
    Ref_Velocity        = 0.0;
    Ref_TotalEnergy     = 0.0;
    Ref_Length          = 0.0;
    Ref_Time            = 0.0;
    Ref_Mach            = 0.0;
    
    // Other Reference Properties
    Gamma               = 0.0;
    
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
    Inf_Alpha           = 0.0;
    Inf_Beta            = 0.0;
    // Mach Ramping
    Inf_Mach_Ramp       = 0;
    Inf_Mach_MAX        = 0.0;
    Inf_Mach_MIN        = 0.0;
    
    // Gauge Variables
    Gauge_Pressure      = 0.0;
    Outflow_Pressure    = 0.0;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void Material_Set_Properties(void) {
    double Q[NEQUATIONS];
    double Density[3];
    double Prop[EOS_THERM_DIM];
    
    Gamma = 1.4;
    // Initialize and Set the EOS based on Material Type
    switch (MaterialType) {
        case MATERIAL_TYPE_IDEALGAS:
            EOS_Init(EOS_MODEL_IDEALGAS);
            break;
        case MATERIAL_TYPE_NIST:
            EOS_Init(EOS_MODEL_NIST);
            break;
        default:
            error("Material_Set_Properties:1: Undefined Material Type - %d", MaterialType);
            break;
    }
    EOS_Set(MaterialName);
    
    // Set the EOS Reference Properties
    EOS_Set_Reference_Properties(Ref_Pressure, Ref_Temperature, Ref_Length);
    EOS_Print_Reference_Properties();
    // Get the Reference Properties
    EOS_Get_Reference_Properties(Prop);
    
    // Update the Final Reference Properties (SI Units)
    Ref_Rho             = Prop[ 0];
    Ref_Pressure        = Prop[ 1];
    Ref_Temperature     = Prop[ 2];
    Ref_Velocity        = Prop[ 3];
    Ref_Length          = Prop[ 4];
    Ref_Mach            = Prop[ 6];
    Ref_Time            = Prop[ 7];
    Ref_TotalEnergy     = Prop[12];
    
    // Compute the Infinity Conditions (SI Units)
    EOS_Get_PT_Density(EOS_DIMENSIONAL_IO_D_D, Inf_Pressure, Inf_Temperature, Density);
    Inf_Rho = Density[0];
    Q[0]    = Inf_Rho;
    Q[1]    = Q[2] = Q[3] = 0.0;
    Q[4]    = Inf_Temperature;
    EOS_Get_Properties(EOS_DIMENSIONAL_IO_D_D, EOS_VARIABLE_RUT, Q, Prop);
    Inf_Alpha       = Inf_Alpha * M_PI / 180.0;
    Inf_Beta        = Inf_Beta * M_PI / 180.0;
    Inf_SpeedSound  = Prop[9];
    Inf_U           = Inf_Mach*Inf_SpeedSound*cos(Inf_Alpha)*cos(Inf_Beta);
    Inf_V           = Inf_Mach*Inf_SpeedSound*sin(Inf_Alpha);
    Inf_W           = Inf_Mach*Inf_SpeedSound*sin(Inf_Beta)*cos(Inf_Alpha);
    Inf_Et          = Prop[13] + 0.5*(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W);
    
    // Print the Infinity Conditions
    info("Infinity Conditions:");
    printf("-----------------------------------------------------------------------------\n");
    info("Density_Inf -------------------------------: %15.6f kg/m^3", Inf_Rho);
    info("Pressure_Inf ------------------------------: %15.6f kPa",    Inf_Pressure/1000.0);
    info("Temperature_Inf ---------------------------: %15.6f K",      Inf_Temperature);
    info("Velocity_U_Inf ----------------------------: %15.6f m/s",    Inf_U);
    info("Velocity_V_Inf ----------------------------: %15.6f m/s",    Inf_V);
    info("Velocity_W_Inf ----------------------------: %15.6f m/s",    Inf_W);
    info("SpeedSound_Inf ----------------------------: %15.6f m/s",    Inf_SpeedSound);
    info("Mach_Inf ----------------------------------: %15.6f ",       Inf_Mach);
    info("TotalEnergy_Inf ---------------------------: %15.6f J/kg",   Inf_Et);
    printf("=============================================================================\n");
    
    // Check the Material Computation Type
    switch (MaterialCompType) {
        case MATERIAL_COMP_TYPE_NDIM:
            //--START Non-Dimensionalization 
            // Infinity Conditions
            Inf_Rho         /= Ref_Rho;
            Inf_Pressure    /= Ref_Pressure;
            Inf_Temperature /= Ref_Temperature;
            Inf_U           /= Ref_Velocity;
            Inf_V           /= Ref_Velocity;
            Inf_W           /= Ref_Velocity;
            Inf_SpeedSound  /= Ref_Velocity;
            Inf_Et          /= Ref_TotalEnergy;
            // Pressure Conditions
            Gauge_Pressure   /= Ref_Pressure;
            Outflow_Pressure /= Ref_Pressure;
            // Set Computations I/O Mode
            MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_ND_ND;
            //--END Non-Dimensionalization
            break;
        case MATERIAL_COMP_TYPE_DIM:
            // Set Computations I/O Mode
            MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D;
            break;
        default:
            error("Material_Set_Properties:2: Undefined Material Computation Type - %d", MaterialCompType);
            break;
    }
    
    // Set All Pressure be Perturbations
    Inf_Pressure     -= Gauge_Pressure;
    Outflow_Pressure -= Gauge_Pressure;
    
    PrecondGlobalMach = Inf_Mach_MAX;
}

//------------------------------------------------------------------------------
//! Set Far Field Condition with Mach Ramping
//------------------------------------------------------------------------------
void Material_Set_InfinityCondition(int Iteration) {
    double tmpMach = 0.0;

    // Compute the Far Field Condition Ramping
    if ((Inf_Mach_Ramp > 1) && (Inf_Mach_MAX > Inf_Mach_MIN)) {
        if (Iteration < Inf_Mach_Ramp)
            tmpMach = Inf_Mach_MIN + (Inf_Mach_MAX - Inf_Mach_MIN)*(((double)Iteration)/((double)(Inf_Mach_Ramp-1)));
        else
            tmpMach = Inf_Mach_MAX;
    } else
        tmpMach = Inf_Mach_MAX;

    // Set the Far Field Conditions
    Inf_Mach = tmpMach;
    Inf_Et   = Inf_Et - 0.5*(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W); 
    Inf_U    = Inf_Mach*Inf_SpeedSound*cos(Inf_Alpha)*cos(Inf_Beta);
    Inf_V    = Inf_Mach*Inf_SpeedSound*sin(Inf_Alpha);
    Inf_W    = Inf_Mach*Inf_SpeedSound*sin(Inf_Beta)*cos(Inf_Alpha);
    Inf_Et   = Inf_Et + 0.5*(Inf_U*Inf_U + Inf_V*Inf_V + Inf_W*Inf_W);
}

//------------------------------------------------------------------------------
//! Compute Equation of State Variables at Control Volume
//! Note: MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D or EOS_DIMENSIONAL_IO_ND_ND
//       is only supported. Done to improve the accuracy
//------------------------------------------------------------------------------
void Material_Get_ControlVolume_Properties(double *dpVariableIn,
        double &Rho, double &Pressure, double &Temperature,
        double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
        double &SpeedSound, double &Mach, double &TotalEnergy, double &TotalEnthalpy) {
    double Prop[EOS_THERM_DIM];
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    // Compute the EOS Based on Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable, Prop);
            // Update the variables
            Rho         = dpVariableIn[0];
            Velocity_U  = dpVariableIn[1]/Rho;
            Velocity_V  = dpVariableIn[2]/Rho;
            Velocity_W  = dpVariableIn[3]/Rho;
            TotalEnergy = dpVariableIn[4]/Rho;
            Pressure    = Prop[ 3] - Gauge_Pressure;
            Temperature = Prop[ 4];
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable, Prop);
            Rho         = dpVariableIn[0];
            Velocity_U  = dpVariableIn[1];
            Velocity_V  = dpVariableIn[2];
            Velocity_W  = dpVariableIn[3];
            Pressure    = dpVariableIn[4]; // Perturbation
            Temperature = Prop[ 4];
            TotalEnergy = Prop[15];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable, Prop);
            Pressure    = dpVariableIn[0]; // Perturbation
            Velocity_U  = dpVariableIn[1];
            Velocity_V  = dpVariableIn[2];
            Velocity_W  = dpVariableIn[3];
            Temperature = dpVariableIn[4];
            Rho         = Prop[ 0];
            TotalEnergy = Prop[15];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            EOS_Get_Properties(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable, Prop);
            Rho         = dpVariableIn[0];
            Velocity_U  = dpVariableIn[1];
            Velocity_V  = dpVariableIn[2];
            Velocity_W  = dpVariableIn[3];
            Temperature = dpVariableIn[4];
            Pressure    = Prop[ 3] - Gauge_Pressure;
            TotalEnergy = Prop[15];
            break;
        default:
            error("Material_Get_ControlVolume_Properties:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    // Update the Variables
    Q2            = Prop[ 8];
    SpeedSound    = Prop[ 9];
    Mach          = Prop[10];
    TotalEnthalpy = Prop[14];
}

//------------------------------------------------------------------------------
//! Compute Equation of State Variables at Face/Edge
//! Note: MaterialEOS_IO_Type = EOS_DIMENSIONAL_IO_D_D or EOS_DIMENSIONAL_IO_ND_ND
//       is only supported. Done to improve the accuracy
//------------------------------------------------------------------------------
void Material_Get_Face_Properties(double *dpVariableIn, double nx, double ny, double nz,
        double &Rho, double &Pressure, double &Temperature,
        double &Velocity_U, double &Velocity_V, double &Velocity_W, double &Q2,
        double &SpeedSound, double &Mach, double &Ubar, double &TotalEnergy, 
        double &TotalEnthalpy) {
    
    Material_Get_ControlVolume_Properties(dpVariableIn, Rho, Pressure, Temperature,
            Velocity_U, Velocity_V, Velocity_W, Q2, SpeedSound, Mach, 
            TotalEnergy, TotalEnthalpy);
    
    // Compute the Ubar
    Ubar = Velocity_U*nx + Velocity_V*ny + Velocity_W*nz;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Density(double *dpVariableIn) {
    double Rho = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Rho = EOS_Get_Density(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Density:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return Rho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Pressure(double *dpVariableIn) {
    double Pressure     = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Pressure = EOS_Get_Pressure(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Pressure:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return (Pressure - Gauge_Pressure);
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Temperature(double *dpVariableIn) {
    double Temperature = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Temperature = EOS_Get_Temperature(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Temperature:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return Temperature;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_Mach(double *dpVariableIn) {
    double Mach = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Mach = EOS_Get_Mach(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_Mach:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return Mach;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_TotalEnergy(double *dpVariableIn) {
    double TotalEnergy = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            TotalEnergy = EOS_Get_TotalEnergy(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_TotalEnergy:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return TotalEnergy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_SpeedSound(double *dpVariableIn) {
    double SpeedSound = 0.0;
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_CON, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_RUP, daEOSVariable);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_PUT, daEOSVariable);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            SpeedSound = EOS_Get_SpeedSound(MaterialEOS_IO_Type, EOS_VARIABLE_RUT, daEOSVariable);
            break;
        default:
            error("Material_Get_SpeedSound:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    return SpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_DH_SpeedSound(double dvDensity, double dvEnthalpy) {
    double SpeedSound = 0.0;
    
    SpeedSound = EOS_Get_DH_SpeedSound(MaterialEOS_IO_Type, dvDensity, dvEnthalpy);
    
    return SpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_DH_Pressure(double dvDensity, double dvEnthalpy) {
    double Pressure = 0.0;
    
    Pressure = EOS_Get_DH_Pressure(MaterialEOS_IO_Type, dvDensity, dvEnthalpy);
    
    return Pressure;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double Material_Get_DH_Temperature(double dvDensity, double dvEnthalpy) {
    double Temperature = 0.0;
    
    Temperature = EOS_Get_DH_Temperature(MaterialEOS_IO_Type, dvDensity, dvEnthalpy);
    
    return Temperature;
}

//------------------------------------------------------------------------------
//! Compute Q's from EOS Variables
//------------------------------------------------------------------------------
void Material_Get_RUH_To_Q(double Rho, double Velocity_U, double Velocity_V, double Velocity_W, double Enthalpy, double *Q) {
    double Pressure, Q2;
    
    // Compute the EOS Based on Non Dimensionalization and Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            Q2       = Velocity_U*Velocity_U + Velocity_V*Velocity_V + Velocity_W*Velocity_W;
            Pressure = Material_Get_DH_Pressure(Rho, Enthalpy);
            Q[0]     = Rho;
            Q[1]     = Rho*Velocity_U;
            Q[2]     = Rho*Velocity_V;
            Q[3]     = Rho*Velocity_W;
            Q[4]     = Rho*(Enthalpy + 0.5*Q2) - Pressure;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            Pressure = Material_Get_DH_Pressure(Rho, Enthalpy);
            Q[0]     = Pressure - Gauge_Pressure;
            Q[1]     = Velocity_U;
            Q[2]     = Velocity_V;
            Q[3]     = Velocity_W;
            Q[4]     = Material_Get_DH_Temperature(Rho, Enthalpy);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            Pressure = Material_Get_DH_Pressure(Rho, Enthalpy);
            Q[0]     = Rho;
            Q[1]     = Velocity_U;
            Q[2]     = Velocity_V;
            Q[3]     = Velocity_W;
            Q[4]     = Pressure - Gauge_Pressure;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            Q[0]     = Rho;
            Q[1]     = Velocity_U;
            Q[2]     = Velocity_V;
            Q[3]     = Velocity_W;
            Q[4]     = Material_Get_DH_Temperature(Rho, Enthalpy);
            break;
        default:
            error("Compute_Q_From_EOS_Variables:1: Undefined Variable Type - %d", VariableType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Transformation Matrix: Mpr = dqp/dqr
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix(double *dpVariableIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix) {
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    // Compute the EOS Based on Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            // Nothing to Do
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            // Nothing to Do
            break;
        default:
            error("Material_Get_Transformation_Matrix:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    // Compute the Transformation Matrix
    EOS_Get_Transformation_Matrix(MaterialEOS_IO_Type, VariableType, daEOSVariable, ivVarTypeFrom, ivVarTypeTo, Matrix);
}

