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
double dvTemperature_Ref;
double dvPressure_Ref;
double dvDensity_Ref;
double dvVelocity_Ref;
double dvLength_Ref;
double dvMach_Ref;
double dvSpeedSound_Ref;
double dvTime_Ref;
double dvEnthalpy_Ref;
double dvTotalEnthalpy_Ref;
double dvInternalEnergy_Ref;
double dvTotalEnergy_Ref;
double dvEntropy_Ref;
double dvEntropyConst_Ref;
double dvGasConstant_Ref;
double dvRatioSpecificHeat_Ref;

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
    info("Molecular Weight                  = %10.4e", SogFluidInfo.dvMolecularWeight);
    info("Gas Constant                      = %10.4e", SogFluidInfo.dvGasConstant);
    info("Temperature Critical              = %10.4e", SogFluidInfo.dvTemperature_Crit);
    info("Pressure Critical                 = %10.4e", SogFluidInfo.dvPressure_Crit);
    info("Compressibility @ Critical Point  = %10.4e", SogFluidInfo.dvZ_Crit);
    info("Density Critical                  = %10.4e", SogFluidInfo.dvDensity_Crit);
    info("Triple Point Temperature          = %10.4e", SogFluidInfo.dvTPTemperature);
    info("Normal Boiling Point Temperature  = %10.4e", SogFluidInfo.dvNBPTemperature);
    info("Acentric Factor                   = %10.4e", SogFluidInfo.dvAcentric_Fact);
    info("Dipole Moment                     = %10.4e", SogFluidInfo.dvDipole_Mom);
    info("Temperature Minimum               = %10.4e", SogFluidInfo.dvTemperature_Min);
    info("Temperature Maximum               = %10.4e", SogFluidInfo.dvTemperature_Max);
    info("Pressure Maximum                  = %10.4e", SogFluidInfo.dvPressure_Max);
    info("Density Maximum                   = %10.4e", SogFluidInfo.dvDensity_Max);
}

//------------------------------------------------------------------------------
//! Initialize Reference Property for Fluid
//! Initialization is based for Ideal Gas
//------------------------------------------------------------------------------
void EOS_Internal_Init_Reference_Properties() {
    /* Properties of Ideal Gas */
    dvTemperature_Ref       = 273.15;
    dvPressure_Ref          = 101325;
    dvGasConstant_Ref       = 287.04;
    dvRatioSpecificHeat_Ref = 1.4;
    dvLength_Ref            = 1.0;
    dvDensity_Ref           = dvPressure_Ref/(dvGasConstant_Ref*dvTemperature_Ref);
    dvSpeedSound_Ref        = sqrt(dvRatioSpecificHeat_Ref*dvGasConstant_Ref*dvTemperature_Ref);
    dvVelocity_Ref          = dvSpeedSound_Ref;
    dvMach_Ref              = 1.0;
    dvTime_Ref              = dvLength_Ref/dvSpeedSound_Ref;
    dvEnthalpy_Ref          = dvSpeedSound_Ref*dvSpeedSound_Ref;
    dvTotalEnthalpy_Ref     = dvEnthalpy_Ref;
    dvInternalEnergy_Ref    = dvEnthalpy_Ref;
    dvTotalEnergy_Ref       = dvEnthalpy_Ref;
    dvEntropy_Ref           = dvGasConstant_Ref/(dvRatioSpecificHeat_Ref - 1.0);
    dvEntropyConst_Ref      = 1.0;
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Properties
//
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_Properties(int ivVariableType, double *dpPropertyIn, double *dpPropertyOut) {
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            dpPropertyOut[0] = dpPropertyIn[0]*dvDensity_Ref;
            dpPropertyOut[1] = dpPropertyIn[1]*dvDensity_Ref*dvVelocity_Ref;
            dpPropertyOut[2] = dpPropertyIn[2]*dvDensity_Ref*dvVelocity_Ref;
            dpPropertyOut[3] = dpPropertyIn[3]*dvDensity_Ref*dvVelocity_Ref;
            dpPropertyOut[4] = dpPropertyIn[4]*dvDensity_Ref*dvTotalEnergy_Ref;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dpPropertyOut[0] = dpPropertyIn[0]*dvPressure_Ref;
            dpPropertyOut[1] = dpPropertyIn[1]*dvVelocity_Ref;
            dpPropertyOut[2] = dpPropertyIn[2]*dvVelocity_Ref;
            dpPropertyOut[3] = dpPropertyIn[3]*dvVelocity_Ref;
            dpPropertyOut[4] = dpPropertyIn[4]*dvTemperature_Ref;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dpPropertyOut[0] = dpPropertyIn[0]*dvDensity_Ref;
            dpPropertyOut[1] = dpPropertyIn[1]*dvVelocity_Ref;
            dpPropertyOut[2] = dpPropertyIn[2]*dvVelocity_Ref;
            dpPropertyOut[3] = dpPropertyIn[3]*dvVelocity_Ref;
            dpPropertyOut[4] = dpPropertyIn[4]*dvPressure_Ref;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dpPropertyOut[0] = dpPropertyIn[0]*dvDensity_Ref;
            dpPropertyOut[1] = dpPropertyIn[1]*dvVelocity_Ref;
            dpPropertyOut[2] = dpPropertyIn[2]*dvVelocity_Ref;
            dpPropertyOut[3] = dpPropertyIn[3]*dvVelocity_Ref;
            dpPropertyOut[4] = dpPropertyIn[4]*dvTemperature_Ref;
            break;
        default:
            error("EOS_Internal_Dimensionalize_Properties: Undefined Variable Type - %d", ivVariableType);
            break;
    }
}


//------------------------------------------------------------------------------
//! Dimensionalize Density and Temperature
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_DT(double *dpDensity, double *dpTemperature) {
    *dpDensity     = *dpDensity*dvDensity_Ref;
    *dpTemperature = *dpTemperature*dvTemperature_Ref;
}

//------------------------------------------------------------------------------
//! Dimensionalize Density and Pressure
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_DP(double *dpDensity, double *dpPressure) {
    *dpDensity  = *dpDensity*dvDensity_Ref;
    *dpPressure = *dpPressure*dvPressure_Ref;
}

//------------------------------------------------------------------------------
//! Dimensionalize Pressure and Temperature
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_PT(double *dpPressure, double *dpTemperature) {
    *dpPressure    = *dpPressure*dvPressure_Ref;
    *dpTemperature = *dpTemperature*dvTemperature_Ref;
}

//------------------------------------------------------------------------------
//! Non-Dimensionalize Pressure
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Pressure(double *dpPressure) {
    *dpPressure = *dpPressure/dvPressure_Ref;
}

//------------------------------------------------------------------------------
//! Non-Dimensionalize Density
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Density(double *dpDensity) {
    *dpDensity = *dpDensity/dvDensity_Ref;
}

//------------------------------------------------------------------------------
//! Non-Dimensionalize Temperature
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Temperature(double *dpTemperature) {
    *dpTemperature = *dpTemperature/dvTemperature_Ref;
}

