/*******************************************************************************
 * File:        EOS_Internal.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "Trim_Utils.h"
#include "EOS_Internal.h"

/* Common Shared Data */
/* NIST Thermodynamics Input-Output Helpers */
EOS_S_NISTHelper SogNISTHelper;
/* Fluid Information */
EOS_S_FluidInfo SogFluidInfo;
/* Reference Properties */
int    ivgSet_Ref;
double dvgTemperature_Ref;
double dvgPressure_Ref;
double dvgDensity_Ref;
double dvgVelocity_Ref;
double dvgLength_Ref;
double dvgMach_Ref;
double dvgSpeedSound_Ref;
double dvgTime_Ref;
double dvgEnthalpy_Ref;
double dvgTotalEnthalpy_Ref;
double dvgInternalEnergy_Ref;
double dvgTotalEnergy_Ref;
double dvgEntropy_Ref;
double dvgEntropyConst_Ref;
double dvgGasConstant_Ref;
double dvgHeatCapacityCv_Ref;
double dvgHeatCapacityCp_Ref;
double dvgRatioSpecificHeat_Ref;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void EOS_Internal_Init() {
    int i;
    // NIST Thermodynamics Input-Output Helpers
    SogNISTHelper.ivNumberComponent = 0;
    SogNISTHelper.ivError           = 0;
    str_blank(SogNISTHelper.csFiles);
    str_blank(SogNISTHelper.csFmix);
    str_blank(SogNISTHelper.csError);
    str_blank(SogNISTHelper.csPath);
    str_blank(SogNISTHelper.csRef);
    for (i = 0; i < NISTTHERMO_MAX_COMPONENT; i++) {
        SogNISTHelper.daX[i]    = 0.0;
        SogNISTHelper.daXLiq[i] = 0.0;
        SogNISTHelper.daXVap[i] = 0.0;
    }
    
    // Initialize the Fluid Information
    EOS_Internal_Init_Fluid_Information();
    
    // Initialize the Fluid Reference Properties
    EOS_Internal_Init_Reference_Properties();
}

//------------------------------------------------------------------------------
//! Initialize Fluid Information
//------------------------------------------------------------------------------
void EOS_Internal_Init_Fluid_Information() {
    str_blank(SogFluidInfo.csFluidName);
    SogFluidInfo.dvTemperature_Crit = 0.0;
    SogFluidInfo.dvPressure_Crit    = 0.0;
    SogFluidInfo.dvDensity_Crit     = 0.0;
    SogFluidInfo.dvMolecularWeight  = 0.0;
    SogFluidInfo.dvTPTemperature    = 0.0;
    SogFluidInfo.dvNBPTemperature   = 0.0;
    SogFluidInfo.dvZ_Crit           = 0.0;
    SogFluidInfo.dvAcentric_Fact    = 0.0;
    SogFluidInfo.dvDipole_Mom       = 0.0;
    SogFluidInfo.dvGasConstant      = 0.0;
    SogFluidInfo.dvTemperature_Min  = 0.0;
    SogFluidInfo.dvTemperature_Max  = 0.0;
    SogFluidInfo.dvPressure_Max     = 0.0;
    SogFluidInfo.dvDensity_Max      = 0.0;
}

//------------------------------------------------------------------------------
//! Initialize Fluid Information from NIST Thermodynamic Library
//  All Thermodynamic variables are in SI Units
//------------------------------------------------------------------------------
void EOS_Internal_Set_NIST_Fluid_Information(int ivComp) {
    // Get the Information about the fluid
    INFO(&ivComp, &SogFluidInfo.dvMolecularWeight, &SogFluidInfo.dvTPTemperature, 
            &SogFluidInfo.dvNBPTemperature, &SogFluidInfo.dvTemperature_Crit, 
            &SogFluidInfo.dvPressure_Crit, &SogFluidInfo.dvDensity_Crit, 
            &SogFluidInfo.dvZ_Crit, &SogFluidInfo.dvAcentric_Fact, &SogFluidInfo.dvDipole_Mom, 
            &SogFluidInfo.dvGasConstant);
    // Get the fluid operating range
    str_blank(SogNISTHelper.csRef);
    strncpy(SogNISTHelper.csRef, "EOS", 3);
#if defined(WIN32) && !defined(IFORT)
    LIMITK(SogNISTHelper.csRef, NISTTHERMO_REF_STR_LEN, &ivComp, &SogFluidInfo.dvTemperature_Crit, &SogFluidInfo.dvDensity_Crit,
            &SogFluidInfo.dvPressure_Crit, &SogFluidInfo.dvTemperature_Min, &SogFluidInfo.dvTemperature_Max,
            &SogFluidInfo.dvDensity_Max, &SogFluidInfo.dvPressure_Max, &SogNISTHelper.ivError, SogNISTHelper.csError, 
            NISTTHERMO_GEN_STR_LEN);
#else
    LIMITK(SogNISTHelper.csRef, &ivComp, &SogFluidInfo.dvTemperature_Crit, &SogFluidInfo.dvDensity_Crit,
            &SogFluidInfo.dvPressure_Crit, &SogFluidInfo.dvTemperature_Min, &SogFluidInfo.dvTemperature_Max,
            &SogFluidInfo.dvDensity_Max, &SogFluidInfo.dvPressure_Max, &SogNISTHelper.ivError, SogNISTHelper.csError,
            NISTTHERMO_REF_STR_LEN, NISTTHERMO_GEN_STR_LEN);
#endif
    if (SogNISTHelper.ivError != 0)
        error("EOS_Internal_Set_NIST_Fluid_Information: %s", SogNISTHelper.csError);
    
    // Convert into International Standard Units (SI)
    SogFluidInfo.dvPressure_Crit = EOS_Internal_Get_NIST_To_SI_Pressure(SogFluidInfo.dvPressure_Crit);
    SogFluidInfo.dvPressure_Max  = EOS_Internal_Get_NIST_To_SI_Pressure(SogFluidInfo.dvPressure_Max);
    SogFluidInfo.dvDensity_Crit  = EOS_Internal_Get_NIST_To_SI_Density(SogFluidInfo.dvDensity_Crit);
    SogFluidInfo.dvDensity_Max   = EOS_Internal_Get_NIST_To_SI_Density(SogFluidInfo.dvDensity_Max);
    SogFluidInfo.dvGasConstant  /= 0.001*SogFluidInfo.dvMolecularWeight;
}

//------------------------------------------------------------------------------
//! Print Fluid Information from NIST Thermodynamic Library
//------------------------------------------------------------------------------
void EOS_Internal_Print_NIST_Fluid_Information(int ivComp) {
    if (ivComp != 1)
        error("Only Single fluid implemented !");
    printf("=============================================================================\n");
    info("FLUID INFORMATION : %s", SogFluidInfo.csFluidName);
    printf("=============================================================================\n");
    info("Molecular Weight                  = %15.6f kg/kmol", SogFluidInfo.dvMolecularWeight);
    info("Gas Constant                      = %15.6f J/kg-K",  SogFluidInfo.dvGasConstant);
    info("Temperature Critical              = %15.6f K",       SogFluidInfo.dvTemperature_Crit);
    info("Pressure Critical                 = %15.6f kPa",     EOS_Internal_Get_SI_To_NIST_Pressure(SogFluidInfo.dvPressure_Crit));
    info("Compressibility @ Critical Point  = %15.6f ",        SogFluidInfo.dvZ_Crit);
    info("Density Critical                  = %15.6f kg/m^3",  SogFluidInfo.dvDensity_Crit);
    info("Triple Point Temperature          = %15.6f K",       SogFluidInfo.dvTPTemperature);
    info("Normal Boiling Point Temperature  = %15.6f K",       SogFluidInfo.dvNBPTemperature);
    info("Acentric Factor                   = %15.6f ",        SogFluidInfo.dvAcentric_Fact);
    info("Dipole Moment                     = %15.6f debye",   SogFluidInfo.dvDipole_Mom);
    info("Temperature Minimum               = %15.6f K",       SogFluidInfo.dvTemperature_Min);
    info("Temperature Maximum               = %15.6f K",       SogFluidInfo.dvTemperature_Max);
    info("Pressure Maximum                  = %15.6f kPa",     EOS_Internal_Get_SI_To_NIST_Pressure(SogFluidInfo.dvPressure_Max));
    info("Density Maximum                   = %15.6f kg/m^3",  SogFluidInfo.dvDensity_Max);
}

//------------------------------------------------------------------------------
//! Initialize Reference Property for Fluid
//------------------------------------------------------------------------------
void EOS_Internal_Init_Reference_Properties() {
    /* Initialize the Reference Properties */
    ivgSet_Ref               = 0;
    dvgTemperature_Ref       = 0.0;
    dvgPressure_Ref          = 0.0;
    dvgGasConstant_Ref       = 0.0;
    dvgRatioSpecificHeat_Ref = 0.0;
    dvgLength_Ref            = 0.0;
    dvgDensity_Ref           = 0.0;
    dvgSpeedSound_Ref        = 0.0;
    dvgVelocity_Ref          = 0.0;
    dvgMach_Ref              = 0.0;
    dvgTime_Ref              = 0.0;
    dvgEnthalpy_Ref          = 0.0;
    dvgTotalEnthalpy_Ref     = 0.0;
    dvgInternalEnergy_Ref    = 0.0;
    dvgTotalEnergy_Ref       = 0.0;
    dvgEntropy_Ref           = 0.0;
    dvgHeatCapacityCv_Ref    = 0.0;
    dvgHeatCapacityCp_Ref    = 0.0;
    dvgEntropyConst_Ref      = 0.0;
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Density
// Density in kg/m^3 to mol/L
//------------------------------------------------------------------------------
double EOS_Internal_Get_SI_To_NIST_Density(double dvDensity) {
    return dvDensity/SogFluidInfo.dvMolecularWeight;
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Density
// Density in mol/L to kg/m^3 
//------------------------------------------------------------------------------
double EOS_Internal_Get_NIST_To_SI_Density(double dvDensity) {
    return dvDensity*SogFluidInfo.dvMolecularWeight;
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Pressure
// Pressure in Pa to kPa
//------------------------------------------------------------------------------
double EOS_Internal_Get_SI_To_NIST_Pressure(double dvPressure) {
    return 0.001*dvPressure;
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Pressure
// Pressure in kPa to Pa
//------------------------------------------------------------------------------
double EOS_Internal_Get_NIST_To_SI_Pressure(double dvPressure) {
    return 1000.0*dvPressure;
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Variables
//  Note: Reference values and dpVariableOut are in SI units
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_Variables(int ivVariableType, double *dpVariableIn, double *dpVariableOut) {
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            dpVariableOut[0] = dpVariableIn[0]*dvgDensity_Ref;
            dpVariableOut[1] = dpVariableIn[1]*dvgDensity_Ref*dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]*dvgDensity_Ref*dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]*dvgDensity_Ref*dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]*dvgDensity_Ref*dvgTotalEnergy_Ref;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dpVariableOut[0] = dpVariableIn[0]*dvgDensity_Ref;
            dpVariableOut[1] = dpVariableIn[1]*dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]*dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]*dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]*dvgPressure_Ref;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dpVariableOut[0] = dpVariableIn[0]*dvgPressure_Ref;
            dpVariableOut[1] = dpVariableIn[1]*dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]*dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]*dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]*dvgTemperature_Ref;
            break;
        // Primitive Variable Formulation Pressure Velocity Entropy
        case EOS_VARIABLE_PUS:
            dpVariableOut[0] = dpVariableIn[0]*dvgPressure_Ref;
            dpVariableOut[1] = dpVariableIn[1]*dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]*dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]*dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]*dvgEntropy_Ref;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dpVariableOut[0] = dpVariableIn[0]*dvgDensity_Ref;
            dpVariableOut[1] = dpVariableIn[1]*dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]*dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]*dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]*dvgTemperature_Ref;
            break;
        default:
            error("EOS_Internal_Dimensionalize_Variables: Undefined Variable Type - %d", ivVariableType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute the Non-Dimensionalize the Properties
//  Note: Reference values and dpVariableIn are in SI units
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Variables(int ivVariableType, double *dpVariableIn, double *dpVariableOut) {
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            dpVariableOut[0] = dpVariableIn[0]/dvgDensity_Ref;
            dpVariableOut[1] = dpVariableIn[1]/(dvgDensity_Ref*dvgVelocity_Ref);
            dpVariableOut[2] = dpVariableIn[2]/(dvgDensity_Ref*dvgVelocity_Ref);
            dpVariableOut[3] = dpVariableIn[3]/(dvgDensity_Ref*dvgVelocity_Ref);
            dpVariableOut[4] = dpVariableIn[4]/(dvgDensity_Ref*dvgTotalEnergy_Ref);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dpVariableOut[0] = dpVariableIn[0]/dvgDensity_Ref;
            dpVariableOut[1] = dpVariableIn[1]/dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]/dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]/dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]/dvgPressure_Ref;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dpVariableOut[0] = dpVariableIn[0]/dvgPressure_Ref;
            dpVariableOut[1] = dpVariableIn[1]/dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]/dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]/dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]/dvgTemperature_Ref;
            break;
        // Primitive Variable Formulation Pressure Velocity Entropy
        case EOS_VARIABLE_PUS:
            dpVariableOut[0] = dpVariableIn[0]/dvgPressure_Ref;
            dpVariableOut[1] = dpVariableIn[1]/dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]/dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]/dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]/dvgEntropy_Ref;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dpVariableOut[0] = dpVariableIn[0]/dvgDensity_Ref;
            dpVariableOut[1] = dpVariableIn[1]/dvgVelocity_Ref;
            dpVariableOut[2] = dpVariableIn[2]/dvgVelocity_Ref;
            dpVariableOut[3] = dpVariableIn[3]/dvgVelocity_Ref;
            dpVariableOut[4] = dpVariableIn[4]/dvgTemperature_Ref;
            break;
        default:
            error("EOS_Internal_NonDimensionalize_Variables: Undefined Variable Type - %d", ivVariableType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Properties
//  Note: Reference values and dpPropertyOut are in SI units
//--------------------------------------------------
void EOS_Internal_Dimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut) {
    dpPropertyOut[ 0] = dpPropertyIn[ 0] * dvgDensity_Ref;
    dpPropertyOut[ 1] = dpPropertyIn[ 1] * dvgPressure_Ref;
    dpPropertyOut[ 2] = dpPropertyIn[ 2] * dvgTemperature_Ref;
    dpPropertyOut[ 3] = dpPropertyIn[ 3] * dvgVelocity_Ref;
    dpPropertyOut[ 4] = dpPropertyIn[ 4] * dvgVelocity_Ref;
    dpPropertyOut[ 5] = dpPropertyIn[ 5] * dvgVelocity_Ref;
    dpPropertyOut[ 6] = dpPropertyIn[ 6] * (dvgVelocity_Ref*dvgVelocity_Ref);
    dpPropertyOut[ 7] = dpPropertyIn[ 7] * dvgSpeedSound_Ref;
    dpPropertyOut[ 8] = dpPropertyIn[ 8] * dvgMach_Ref;
    dpPropertyOut[ 9] = dpPropertyIn[ 9] * dvgEntropy_Ref;
    dpPropertyOut[10] = dpPropertyIn[10] * dvgEnthalpy_Ref;
    dpPropertyOut[11] = dpPropertyIn[11] * dvgInternalEnergy_Ref;
    dpPropertyOut[12] = dpPropertyIn[12] * dvgTotalEnthalpy_Ref;
    dpPropertyOut[13] = dpPropertyIn[13] * dvgTotalEnergy_Ref;
    dpPropertyOut[14] = dpPropertyIn[14] * dvgHeatCapacityCv_Ref;
    dpPropertyOut[15] = dpPropertyIn[15] * dvgHeatCapacityCp_Ref;
}

//------------------------------------------------------------------------------
//! Compute the Non-Dimensionalize the Properties
//  Note: Reference values and dpPropertyIn are in SI units
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut) {
    dpPropertyOut[ 0] = dpPropertyIn[ 0] / dvgDensity_Ref;
    dpPropertyOut[ 1] = dpPropertyIn[ 1] / dvgPressure_Ref;
    dpPropertyOut[ 2] = dpPropertyIn[ 2] / dvgTemperature_Ref;
    dpPropertyOut[ 3] = dpPropertyIn[ 3] / dvgVelocity_Ref;
    dpPropertyOut[ 4] = dpPropertyIn[ 4] / dvgVelocity_Ref;
    dpPropertyOut[ 5] = dpPropertyIn[ 5] / dvgVelocity_Ref;
    dpPropertyOut[ 6] = dpPropertyIn[ 6] / (dvgVelocity_Ref*dvgVelocity_Ref);
    dpPropertyOut[ 7] = dpPropertyIn[ 7] / dvgSpeedSound_Ref;
    dpPropertyOut[ 8] = dpPropertyIn[ 8] / dvgMach_Ref;
    dpPropertyOut[ 9] = dpPropertyIn[ 9] / dvgEntropy_Ref;
    dpPropertyOut[10] = dpPropertyIn[10] / dvgEnthalpy_Ref;
    dpPropertyOut[11] = dpPropertyIn[11] / dvgInternalEnergy_Ref;
    dpPropertyOut[12] = dpPropertyIn[12] / dvgTotalEnthalpy_Ref;
    dpPropertyOut[13] = dpPropertyIn[13] / dvgTotalEnergy_Ref;
    dpPropertyOut[14] = dpPropertyIn[14] / dvgHeatCapacityCv_Ref;
    dpPropertyOut[15] = dpPropertyIn[15] / dvgHeatCapacityCp_Ref;
}

//------------------------------------------------------------------------------
//! Dimensionalize Density and Temperature
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_DT(double *dpDensity, double *dpTemperature) {
    *dpDensity     = *dpDensity*dvgDensity_Ref;
    *dpTemperature = *dpTemperature*dvgTemperature_Ref;
}

//------------------------------------------------------------------------------
//! Dimensionalize Density and Pressure
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_DP(double *dpDensity, double *dpPressure) {
    *dpDensity  = *dpDensity*dvgDensity_Ref;
    *dpPressure = *dpPressure*dvgPressure_Ref;
}

//------------------------------------------------------------------------------
//! Dimensionalize Pressure and Temperature
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_PT(double *dpPressure, double *dpTemperature) {
    *dpPressure    = *dpPressure*dvgPressure_Ref;
    *dpTemperature = *dpTemperature*dvgTemperature_Ref;
}

//------------------------------------------------------------------------------
//! Non-Dimensionalize Pressure
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Pressure(double *dpPressure) {
    *dpPressure = *dpPressure/dvgPressure_Ref;
}

//------------------------------------------------------------------------------
//! Non-Dimensionalize Density
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Density(double *dpDensity) {
    *dpDensity = *dpDensity/dvgDensity_Ref;
}

//------------------------------------------------------------------------------
//! Non-Dimensionalize Temperature
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Temperature(double *dpTemperature) {
    *dpTemperature = *dpTemperature/dvgTemperature_Ref;
}

