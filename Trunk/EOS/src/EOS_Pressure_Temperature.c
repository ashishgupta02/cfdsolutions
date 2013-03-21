/*******************************************************************************
 * File:        EOS_Pressure_Temperature.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "Trim_Utils.h"
#include "EOS_Internal.h"

// ------------------ Pressure and Temperature Formulation ---------------------
// This functions are valid Only in this Regions:
// 1) Super Critical State
// 2) Subcooled Compressed Liquid
// 3) Super Heated Vapor
// -----------------------------------------------------------------------------

 // Compute EOS Properties
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_Pressure(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_SpeedSound(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_Entropy(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_Enthalpy(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_InternalEnergy(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_HeatCapacityCv(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_HeatCapacityCp(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_TotalEnergy(double dvPressure, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_Mach(double dvPressure, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W){
    double result = 0.0;
    
    return result;
}


// Compute First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DPressureDDensity(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DPressureDTemperature(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DDensityDPressure(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DDensityDTemperature(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DTemperatureDDensity(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DTemperatureDPressure(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

// Enthalpy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DEnthalpyDTemperature_CDensity(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DEnthalpyDTemperature_CPressure(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DEnthalpyDDensity_CTemperature(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DEnthalpyDDensity_CPressure(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DEnthalpyDPressure_CTemperature(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DEnthalpyDPressure_CDensity(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

// Internal Energy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DInternalEnergyDTemperature_CDensity(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DInternalEnergyDTemperature_CPressure(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DInternalEnergyDDensity_CTemperature(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DInternalEnergyDDensity_CPressure(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DInternalEnergyDPressure_CTemperature(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_DInternalEnergyDPressure_CDensity(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

// Compute Second Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_D2PressureDDensity2(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_D2PressureDTemperature2(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_PT_Get_D2PressureDTemperatureDensity(double dvPressure, double dvTemperature){
    double result = 0.0;
    
    return result;
}

