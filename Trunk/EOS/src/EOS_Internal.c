/*******************************************************************************
 * File:        EOS_Internal.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "Trim_Utils.h"
#include "EOS.h"
#include "EOS_Internal.h"

/* Common Shared Data */
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
    // Initialize the Fluid Information
    EOS_Internal_Init_Fluid_Information();
    
    // Initialize the Fluid Reference Properties
    EOS_Internal_Init_Reference_Properties();
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void EOS_Internal_Finalize() {
    
}

//------------------------------------------------------------------------------
//! Initialize Fluid Information
//------------------------------------------------------------------------------
void EOS_Internal_Init_Fluid_Information() {
    SogFluidInfo.ivEOSModelType     = EOS_MODEL_NONE;
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
//! Compute the Dimensionalize the Properties
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut) {
    dpPropertyOut[ 0] = dpPropertyIn[ 0] * dvgDensity_Ref;
    dpPropertyOut[ 1] = dpPropertyIn[ 1] * dvgDensity_Ref;
    dpPropertyOut[ 2] = dpPropertyIn[ 2] * dvgDensity_Ref;
    dpPropertyOut[ 3] = dpPropertyIn[ 3] * dvgPressure_Ref;
    dpPropertyOut[ 4] = dpPropertyIn[ 4] * dvgTemperature_Ref;
    dpPropertyOut[ 5] = dpPropertyIn[ 5] * dvgVelocity_Ref;
    dpPropertyOut[ 6] = dpPropertyIn[ 6] * dvgVelocity_Ref;
    dpPropertyOut[ 7] = dpPropertyIn[ 7] * dvgVelocity_Ref;
    dpPropertyOut[ 8] = dpPropertyIn[ 8] * (dvgVelocity_Ref*dvgVelocity_Ref);
    dpPropertyOut[ 9] = dpPropertyIn[ 9] * dvgSpeedSound_Ref;
    dpPropertyOut[10] = dpPropertyIn[10] * dvgMach_Ref;
    dpPropertyOut[11] = dpPropertyIn[11] * dvgEntropy_Ref;
    dpPropertyOut[12] = dpPropertyIn[12] * dvgEnthalpy_Ref;
    dpPropertyOut[13] = dpPropertyIn[13] * dvgInternalEnergy_Ref;
    dpPropertyOut[14] = dpPropertyIn[14] * dvgTotalEnthalpy_Ref;
    dpPropertyOut[15] = dpPropertyIn[15] * dvgTotalEnergy_Ref;
    dpPropertyOut[16] = dpPropertyIn[16] * dvgHeatCapacityCv_Ref;
    dpPropertyOut[17] = dpPropertyIn[17] * dvgHeatCapacityCp_Ref;
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Extended Properties
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
void EOS_Internal_Dimensionalize_Extended_Properties(double *dpPropertyIn, double *dpPropertyOut) {
    dpPropertyOut[ 0] = dpPropertyIn[ 0] * dvgDensity_Ref;
    dpPropertyOut[ 1] = dpPropertyIn[ 1] * dvgDensity_Ref;
    dpPropertyOut[ 2] = dpPropertyIn[ 2] * dvgDensity_Ref;
    dpPropertyOut[ 3] = dpPropertyIn[ 3] * dvgPressure_Ref;
    dpPropertyOut[ 4] = dpPropertyIn[ 4] * dvgTemperature_Ref;
    dpPropertyOut[ 5] = dpPropertyIn[ 5] * dvgVelocity_Ref;
    dpPropertyOut[ 6] = dpPropertyIn[ 6] * dvgVelocity_Ref;
    dpPropertyOut[ 7] = dpPropertyIn[ 7] * dvgVelocity_Ref;
    dpPropertyOut[ 8] = dpPropertyIn[ 8] * (dvgVelocity_Ref*dvgVelocity_Ref);
    dpPropertyOut[ 9] = dpPropertyIn[ 9] * dvgSpeedSound_Ref;
    dpPropertyOut[10] = dpPropertyIn[10] * dvgMach_Ref;
    dpPropertyOut[11] = dpPropertyIn[11] * dvgEntropy_Ref;
    dpPropertyOut[12] = dpPropertyIn[12] * dvgEnthalpy_Ref;
    dpPropertyOut[13] = dpPropertyIn[13] * dvgInternalEnergy_Ref;
    dpPropertyOut[14] = dpPropertyIn[14] * dvgTotalEnthalpy_Ref;
    dpPropertyOut[15] = dpPropertyIn[15] * dvgTotalEnergy_Ref;
    dpPropertyOut[16] = dpPropertyIn[16] * dvgHeatCapacityCv_Ref;
    dpPropertyOut[17] = dpPropertyIn[17] * dvgHeatCapacityCp_Ref;
    dpPropertyOut[18] = dpPropertyIn[18] * (dvgPressure_Ref/dvgDensity_Ref);
    dpPropertyOut[19] = dpPropertyIn[19] * (dvgPressure_Ref/dvgTemperature_Ref);
    dpPropertyOut[20] = dpPropertyIn[20] * (dvgDensity_Ref/dvgTemperature_Ref);
    dpPropertyOut[21] = dpPropertyIn[21] * (dvgDensity_Ref/dvgPressure_Ref);
    dpPropertyOut[22] = dpPropertyIn[22] * (dvgPressure_Ref/(dvgDensity_Ref*dvgDensity_Ref));
    dpPropertyOut[23] = dpPropertyIn[23] * (dvgPressure_Ref/(dvgTemperature_Ref*dvgTemperature_Ref));
    dpPropertyOut[24] = dpPropertyIn[24] * (dvgPressure_Ref/(dvgTemperature_Ref*dvgDensity_Ref));
    dpPropertyOut[25] = dpPropertyIn[25] * (dvgEnthalpy_Ref/dvgTemperature_Ref);
    dpPropertyOut[26] = dpPropertyIn[26] * (dvgEnthalpy_Ref/dvgTemperature_Ref);
    dpPropertyOut[27] = dpPropertyIn[27] * (dvgEnthalpy_Ref/dvgDensity_Ref);
    dpPropertyOut[28] = dpPropertyIn[28] * (dvgEnthalpy_Ref/dvgDensity_Ref);
    dpPropertyOut[29] = dpPropertyIn[29] * (dvgEnthalpy_Ref/dvgPressure_Ref);
    dpPropertyOut[30] = dpPropertyIn[30] * (dvgEnthalpy_Ref/dvgPressure_Ref);
    dpPropertyOut[31] = dpPropertyIn[31] * (dvgInternalEnergy_Ref/dvgTemperature_Ref);
    dpPropertyOut[32] = dpPropertyIn[32] * (dvgInternalEnergy_Ref/dvgTemperature_Ref);
    dpPropertyOut[33] = dpPropertyIn[33] * (dvgInternalEnergy_Ref/dvgDensity_Ref);
    dpPropertyOut[34] = dpPropertyIn[34] * (dvgInternalEnergy_Ref/dvgDensity_Ref);
    dpPropertyOut[35] = dpPropertyIn[35] * (dvgInternalEnergy_Ref/dvgPressure_Ref);
    dpPropertyOut[36] = dpPropertyIn[36] * (dvgInternalEnergy_Ref/dvgPressure_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_Density(double dvDensity) {
    return (dvDensity * dvgDensity_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_Pressure(double dvPressure) {
    return (dvPressure * dvgPressure_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_Temperature(double dvTemperature) {
    return (dvTemperature * dvgTemperature_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_SpeedSound(double dvSpeedSound) {
   return (dvSpeedSound * dvgSpeedSound_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_Entropy(double dvEntropy) {
    return (dvEntropy * dvgEntropy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_Enthalpy(double dvEnthalpy) {
    return (dvEnthalpy * dvgEnthalpy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_InternalEnergy(double dvInternalEnergy) {
    return (dvInternalEnergy * dvgInternalEnergy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_TotalEnthalpy(double dvTotalEnthalpy) {
    return (dvTotalEnthalpy * dvgTotalEnthalpy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_TotalEnergy(double dvTotalEnergy) {
    return (dvTotalEnergy * dvgTotalEnergy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_HeatCapacityCv(double dvHeatCapacityCv) {
    return (dvHeatCapacityCv * dvgHeatCapacityCv_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_HeatCapacityCp(double dvHeatCapacityCp) {
    return (dvHeatCapacityCp * dvgHeatCapacityCp_Ref);
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DPressureDDensity(double dvDPressureDDensity) {
    return (dvDPressureDDensity * (dvgPressure_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DPressureDTemperature(double dvDPressureDTemperature) {
    return (dvDPressureDTemperature * (dvgPressure_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DDensityDTemperature(double dvDDensityDTemperature) {
    return (dvDDensityDTemperature * (dvgDensity_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_D2PressureDDensity2(double dvD2PressureDDensity2) {
    return (dvD2PressureDDensity2 * (dvgPressure_Ref/(dvgDensity_Ref*dvgDensity_Ref)));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_D2PressureDTemperature2(double dvD2PressureDTemperature2) {
    return (dvD2PressureDTemperature2 * (dvgPressure_Ref/(dvgTemperature_Ref*dvgTemperature_Ref)));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_D2PressureDTemperatureDensity(double dvD2PressureDTemperatureDensity) {
    return (dvD2PressureDTemperatureDensity * (dvgPressure_Ref/(dvgTemperature_Ref*dvgDensity_Ref)));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dvPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DEnthalpyDTemperature_CDensity(double dvDEnthalpyDTemperature_CDensity) {
    return (dvDEnthalpyDTemperature_CDensity * (dvgEnthalpy_Ref/dvgTemperature_Ref));   
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DEnthalpyDTemperature_CPressure(double dvDEnthalpyDTemperature_CPressure) {
    return (dvDEnthalpyDTemperature_CPressure * (dvgEnthalpy_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DEnthalpyDDensity_CTemperature(double dvDEnthalpyDDensity_CTemperature) {
    return (dvDEnthalpyDDensity_CTemperature * (dvgEnthalpy_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DEnthalpyDDensity_CPressure(double dvDEnthalpyDDensity_CPressure) {
    return (dvDEnthalpyDDensity_CPressure * (dvgEnthalpy_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DEnthalpyDPressure_CTemperature(double dvDEnthalpyDPressure_CTemperature) {
    return (dvDEnthalpyDPressure_CTemperature * (dvgEnthalpy_Ref/dvgPressure_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DEnthalpyDPressure_CDensity(double dvDEnthalpyDPressure_CDensity) {
    return (dvDEnthalpyDPressure_CDensity * (dvgEnthalpy_Ref/dvgPressure_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DInternalEnergyDTemperature_CDensity(double dvDInternalEnergyDTemperature_CDensity) {
    return (dvDInternalEnergyDTemperature_CDensity * (dvgInternalEnergy_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DInternalEnergyDTemperature_CPressure(double dvDInternalEnergyDTemperature_CPressure) {
    return (dvDInternalEnergyDTemperature_CPressure * (dvgInternalEnergy_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DInternalEnergyDDensity_CTemperature(double dvDInternalEnergyDDensity_CTemperature) {
    return (dvDInternalEnergyDDensity_CTemperature * (dvgInternalEnergy_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DInternalEnergyDDensity_CPressure(double dvDInternalEnergyDDensity_CPressure) {
    return (dvDInternalEnergyDDensity_CPressure * (dvgInternalEnergy_Ref/dvgDensity_Ref));
}    

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DInternalEnergyDPressure_CTemperature(double dvDInternalEnergyDPressure_CTemperature) {
    return (dvDInternalEnergyDPressure_CTemperature * (dvgInternalEnergy_Ref/dvgPressure_Ref));
}

//------------------------------------------------------------------------------
//! Compute the Dimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_Dimensionalize_DInternalEnergyDPressure_CDensity(double dvDInternalEnergyDPressure_CDensity) {
    return (dvDInternalEnergyDPressure_CDensity * (dvgInternalEnergy_Ref/dvgPressure_Ref));
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
//! Compute the Non-Dimensionalize the Properties
//  Note: Reference values and dpPropertyIn are in SI units
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut) {
    dpPropertyOut[ 0] = dpPropertyIn[ 0] / dvgDensity_Ref;
    dpPropertyOut[ 1] = dpPropertyIn[ 1] / dvgDensity_Ref;
    dpPropertyOut[ 2] = dpPropertyIn[ 2] / dvgDensity_Ref;
    dpPropertyOut[ 3] = dpPropertyIn[ 3] / dvgPressure_Ref;
    dpPropertyOut[ 4] = dpPropertyIn[ 4] / dvgTemperature_Ref;
    dpPropertyOut[ 5] = dpPropertyIn[ 5] / dvgVelocity_Ref;
    dpPropertyOut[ 6] = dpPropertyIn[ 6] / dvgVelocity_Ref;
    dpPropertyOut[ 7] = dpPropertyIn[ 7] / dvgVelocity_Ref;
    dpPropertyOut[ 8] = dpPropertyIn[ 8] / (dvgVelocity_Ref*dvgVelocity_Ref);
    dpPropertyOut[ 9] = dpPropertyIn[ 9] / dvgSpeedSound_Ref;
    dpPropertyOut[10] = dpPropertyIn[10] / dvgMach_Ref;
    dpPropertyOut[11] = dpPropertyIn[11] / dvgEntropy_Ref;
    dpPropertyOut[12] = dpPropertyIn[12] / dvgEnthalpy_Ref;
    dpPropertyOut[13] = dpPropertyIn[13] / dvgInternalEnergy_Ref;
    dpPropertyOut[14] = dpPropertyIn[14] / dvgTotalEnthalpy_Ref;
    dpPropertyOut[15] = dpPropertyIn[15] / dvgTotalEnergy_Ref;
    dpPropertyOut[16] = dpPropertyIn[16] / dvgHeatCapacityCv_Ref;
    dpPropertyOut[17] = dpPropertyIn[17] / dvgHeatCapacityCp_Ref;
}

//------------------------------------------------------------------------------
//! Compute the Non-Dimensionalize the Extended Properties
//  Note: Reference values and dpPropertyIn are in SI units
//------------------------------------------------------------------------------
void EOS_Internal_NonDimensionalize_Extended_Properties(double *dpPropertyIn, double *dpPropertyOut) {
    dpPropertyOut[ 0] = dpPropertyIn[ 0] / dvgDensity_Ref;
    dpPropertyOut[ 1] = dpPropertyIn[ 1] / dvgDensity_Ref;
    dpPropertyOut[ 2] = dpPropertyIn[ 2] / dvgDensity_Ref;
    dpPropertyOut[ 3] = dpPropertyIn[ 3] / dvgPressure_Ref;
    dpPropertyOut[ 4] = dpPropertyIn[ 4] / dvgTemperature_Ref;
    dpPropertyOut[ 5] = dpPropertyIn[ 5] / dvgVelocity_Ref;
    dpPropertyOut[ 6] = dpPropertyIn[ 6] / dvgVelocity_Ref;
    dpPropertyOut[ 7] = dpPropertyIn[ 7] / dvgVelocity_Ref;
    dpPropertyOut[ 8] = dpPropertyIn[ 8] / (dvgVelocity_Ref*dvgVelocity_Ref);
    dpPropertyOut[ 9] = dpPropertyIn[ 9] / dvgSpeedSound_Ref;
    dpPropertyOut[10] = dpPropertyIn[10] / dvgMach_Ref;
    dpPropertyOut[11] = dpPropertyIn[11] / dvgEntropy_Ref;
    dpPropertyOut[12] = dpPropertyIn[12] / dvgEnthalpy_Ref;
    dpPropertyOut[13] = dpPropertyIn[13] / dvgInternalEnergy_Ref;
    dpPropertyOut[14] = dpPropertyIn[14] / dvgTotalEnthalpy_Ref;
    dpPropertyOut[15] = dpPropertyIn[15] / dvgTotalEnergy_Ref;
    dpPropertyOut[16] = dpPropertyIn[16] / dvgHeatCapacityCv_Ref;
    dpPropertyOut[17] = dpPropertyIn[17] / dvgHeatCapacityCp_Ref;
    dpPropertyOut[18] = dpPropertyIn[18] / (dvgPressure_Ref/dvgDensity_Ref);
    dpPropertyOut[19] = dpPropertyIn[19] / (dvgPressure_Ref/dvgTemperature_Ref);
    dpPropertyOut[20] = dpPropertyIn[20] / (dvgDensity_Ref/dvgTemperature_Ref);
    dpPropertyOut[21] = dpPropertyIn[21] / (dvgDensity_Ref/dvgPressure_Ref);
    dpPropertyOut[22] = dpPropertyIn[22] / (dvgPressure_Ref/(dvgDensity_Ref*dvgDensity_Ref));
    dpPropertyOut[23] = dpPropertyIn[23] / (dvgPressure_Ref/(dvgTemperature_Ref*dvgTemperature_Ref));
    dpPropertyOut[24] = dpPropertyIn[24] / (dvgPressure_Ref/(dvgTemperature_Ref*dvgDensity_Ref));
    dpPropertyOut[25] = dpPropertyIn[25] / (dvgEnthalpy_Ref/dvgTemperature_Ref);
    dpPropertyOut[26] = dpPropertyIn[26] / (dvgEnthalpy_Ref/dvgTemperature_Ref);
    dpPropertyOut[27] = dpPropertyIn[27] / (dvgEnthalpy_Ref/dvgDensity_Ref);
    dpPropertyOut[28] = dpPropertyIn[28] / (dvgEnthalpy_Ref/dvgDensity_Ref);
    dpPropertyOut[29] = dpPropertyIn[29] / (dvgEnthalpy_Ref/dvgPressure_Ref);
    dpPropertyOut[30] = dpPropertyIn[30] / (dvgEnthalpy_Ref/dvgPressure_Ref);
    dpPropertyOut[31] = dpPropertyIn[31] / (dvgInternalEnergy_Ref/dvgTemperature_Ref);
    dpPropertyOut[32] = dpPropertyIn[32] / (dvgInternalEnergy_Ref/dvgTemperature_Ref);
    dpPropertyOut[33] = dpPropertyIn[33] / (dvgInternalEnergy_Ref/dvgDensity_Ref);
    dpPropertyOut[34] = dpPropertyIn[34] / (dvgInternalEnergy_Ref/dvgDensity_Ref);
    dpPropertyOut[35] = dpPropertyIn[35] / (dvgInternalEnergy_Ref/dvgPressure_Ref);
    dpPropertyOut[36] = dpPropertyIn[36] / (dvgInternalEnergy_Ref/dvgPressure_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_Density(double dvDensity) {
    return (dvDensity / dvgDensity_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_Pressure(double dvPressure) {
    return (dvPressure / dvgPressure_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_Temperature(double dvTemperature) {
    return (dvTemperature / dvgTemperature_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_SpeedSound(double dvSpeedSound) {
   return (dvSpeedSound / dvgSpeedSound_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_Entropy(double dvEntropy) {
    return (dvEntropy / dvgEntropy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_Enthalpy(double dvEnthalpy) {
    return (dvEnthalpy / dvgEnthalpy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_InternalEnergy(double dvInternalEnergy) {
    return (dvInternalEnergy / dvgInternalEnergy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_TotalEnthalpy(double dvTotalEnthalpy) {
    return (dvTotalEnthalpy / dvgTotalEnthalpy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_TotalEnergy(double dvTotalEnergy) {
    return (dvTotalEnergy / dvgTotalEnergy_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_HeatCapacityCv(double dvHeatCapacityCv) {
    return (dvHeatCapacityCv / dvgHeatCapacityCv_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_HeatCapacityCp(double dvHeatCapacityCp) {
    return (dvHeatCapacityCp / dvgHeatCapacityCp_Ref);
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DPressureDDensity(double dvDPressureDDensity) {
    return (dvDPressureDDensity / (dvgPressure_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DPressureDTemperature(double dvDPressureDTemperature) {
    return (dvDPressureDTemperature / (dvgPressure_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DDensityDTemperature(double dvDDensityDTemperature) {
    return (dvDDensityDTemperature / (dvgDensity_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_D2PressureDDensity2(double dvD2PressureDDensity2) {
    return (dvD2PressureDDensity2 / (dvgPressure_Ref/(dvgDensity_Ref*dvgDensity_Ref)));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_D2PressureDTemperature2(double dvD2PressureDTemperature2) {
    return (dvD2PressureDTemperature2 / (dvgPressure_Ref/(dvgTemperature_Ref*dvgTemperature_Ref)));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_D2PressureDTemperatureDensity(double dvD2PressureDTemperatureDensity) {
    return (dvD2PressureDTemperatureDensity / (dvgPressure_Ref/(dvgTemperature_Ref*dvgDensity_Ref)));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dvPropertyIn are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CDensity(double dvDEnthalpyDTemperature_CDensity) {
    return (dvDEnthalpyDTemperature_CDensity / (dvgEnthalpy_Ref/dvgTemperature_Ref));   
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CPressure(double dvDEnthalpyDTemperature_CPressure) {
    return (dvDEnthalpyDTemperature_CPressure / (dvgEnthalpy_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DEnthalpyDDensity_CTemperature(double dvDEnthalpyDDensity_CTemperature) {
    return (dvDEnthalpyDDensity_CTemperature / (dvgEnthalpy_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DEnthalpyDDensity_CPressure(double dvDEnthalpyDDensity_CPressure) {
    return (dvDEnthalpyDDensity_CPressure / (dvgEnthalpy_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DEnthalpyDPressure_CTemperature(double dvDEnthalpyDPressure_CTemperature) {
    return (dvDEnthalpyDPressure_CTemperature / (dvgEnthalpy_Ref/dvgPressure_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DEnthalpyDPressure_CDensity(double dvDEnthalpyDPressure_CDensity) {
    return (dvDEnthalpyDPressure_CDensity / (dvgEnthalpy_Ref/dvgPressure_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CDensity(double dvDInternalEnergyDTemperature_CDensity) {
    return (dvDInternalEnergyDTemperature_CDensity / (dvgInternalEnergy_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CPressure(double dvDInternalEnergyDTemperature_CPressure) {
    return (dvDInternalEnergyDTemperature_CPressure / (dvgInternalEnergy_Ref/dvgTemperature_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CTemperature(double dvDInternalEnergyDDensity_CTemperature) {
    return (dvDInternalEnergyDDensity_CTemperature / (dvgInternalEnergy_Ref/dvgDensity_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CPressure(double dvDInternalEnergyDDensity_CPressure) {
    return (dvDInternalEnergyDDensity_CPressure / (dvgInternalEnergy_Ref/dvgDensity_Ref));
}    

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CTemperature(double dvDInternalEnergyDPressure_CTemperature) {
    return (dvDInternalEnergyDPressure_CTemperature / (dvgInternalEnergy_Ref/dvgPressure_Ref));
}

//------------------------------------------------------------------------------
//! Compute the NonDimensionalize the Property
//  Note: Reference values and dpPropertyOut are in SI units
//------------------------------------------------------------------------------
double EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CDensity(double dvDInternalEnergyDPressure_CDensity) {
    return (dvDInternalEnergyDPressure_CDensity / (dvgInternalEnergy_Ref/dvgPressure_Ref));
}

