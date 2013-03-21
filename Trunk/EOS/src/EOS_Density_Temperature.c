/*******************************************************************************
 * File:        EOS_Density_Temperature.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "Trim_Utils.h"
#include "EOS_Internal.h"

// ------------------ Density and Temperature Formulation ----------------------
// This functions are valid for All Regions:
// 1) Super Critical State
// 2) Subcooled Compressed Liquid
// 3) Super Heated Vapor
// 4) Multi-Phase Liquid-Vapor
// -----------------------------------------------------------------------------

// Compute EOS Properties
void EOS_DT_Get_Properties(double dvDensity, double dvTemperature, double *dpProperties) {
    // Initialize
    dpProperties[0] = 0.0; // Pressure
    dpProperties[1] = 0.0; // Internal Energy
    dpProperties[2] = 0.0; // Enthalpy
    dpProperties[3] = 0.0; // Entropy
    dpProperties[4] = 0.0; // Specific Heat at Constant Volume Cv
    dpProperties[5] = 0.0; // Specific Heat at Constant Pressure Cp
    dpProperties[6] = 0.0; // Speed of Sound
    
    // Dimensionalize Density and Temperature
    EOS_Internal_Dimensionalize_DT(&dvDensity, &dvTemperature);
    
}


//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_Pressure(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    // Dimensionalize Density and Temperature
    EOS_Internal_Dimensionalize_DT(&dvDensity, &dvTemperature);
    
    // Determine the Thermodynamic Region
    if ((dvTemperature <= SogFluidInfo.dvTemperature_Min) || (dvTemperature >= SogFluidInfo.dvTemperature_Max))
        error("EOS_DT_Get_Pressure: Temperature Out of Applicability Limit!");
    
    // Region: Super Critical State fluid and Super Heated Vapor
    if (dvTemperature > SogFluidInfo.dvTemperature_Crit) {
        PRESS(&dvTemperature, &dvDensity, SogNISTHelper.daX, &result);
    } else {
        int ivKph = 1; // For Liquid
        double dvDLiq, dvDVap, dvP;
        // Compute the Saturation Values for Liquid
        dvDLiq = dvDVap = dvP = 0.0;
        SATT(&dvTemperature, SogNISTHelper.daX, &ivKph, &dvP, &dvDLiq, &dvDVap,
                SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &SogNISTHelper.ivError, 
                SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            error("EOS_DT_Get_Pressure: %s", SogNISTHelper.csError);
        
        // Determine Thermodynamic Region based on Density
        // Region: Subcooled Compressed Liquid or Heated Vapor
        if ((dvDensity < dvDVap) || (dvDensity > dvDLiq)) {
            dvP = 0.0;
            PRESS(&dvTemperature, &dvDensity, SogNISTHelper.daX, &dvP);
            result = dvP;
        } else {
            // Region: Multiphase
            // Pressure is Saturated Value.
            result = dvP;
        }
    }
    
    // Non-Dimensionalize Pressure
    EOS_Internal_NonDimensionalize_Pressure(&result);
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_SpeedSound(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    // Dimensionalize Density and Temperature
    EOS_Internal_Dimensionalize_DT(&dvDensity, &dvTemperature);
    
    // Determine the Thermodynamic Region
    if ((dvTemperature <= SogFluidInfo.dvTemperature_Min) || (dvTemperature >= SogFluidInfo.dvTemperature_Max))
        error("EOS_DT_Get_SpeedSound: Temperature Out of Applicability Limit!");
    
    // Region: Super Critical State fluid and Super Heated Vapor
    if (dvTemperature > SogFluidInfo.dvTemperature_Crit) {
        PRESS(&dvTemperature, &dvDensity, SogNISTHelper.daX, &result);
    } else {
        int ivKph = 1; // For Liquid
        double dvDLiq, dvDVap, dvP;
        // Compute the Saturation Values for Liquid
        dvDLiq = dvDVap = dvP = 0.0;
        SATT(&dvTemperature, SogNISTHelper.daX, &ivKph, &dvP, &dvDLiq, &dvDVap,
                SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &SogNISTHelper.ivError, 
                SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            error("EOS_DT_Get_SpeedSound: %s", SogNISTHelper.csError);
        
        // Determine Thermodynamic Region based on Density
        // Region: Subcooled Compressed Liquid or Heated Vapor
        if ((dvDensity < dvDVap) || (dvDensity > dvDLiq)) {
            dvP = 0.0;
            PRESS(&dvTemperature, &dvDensity, SogNISTHelper.daX, &dvP);
            result = dvP;
        } else {
            // Region: Multiphase
            // Pressure is Saturated Value.
            result = dvP;
        }
    }
    
    // Non-Dimensionalize Pressure
    EOS_Internal_NonDimensionalize_Pressure(&result);
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_Entropy(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_Enthalpy(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_InternalEnergy(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_HeatCapacityCv(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_HeatCapacityCp(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_TotalEnergy(double dvDensity, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_Mach(double dvDensity, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W) {
    double result = 0.0;
    
    return result;
}

// Compute First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DPressureDDensity(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DPressureDTemperature(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DDensityDPressure(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DDensityDTemperature(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DTemperatureDDensity(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DTemperatureDPressure(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

// Enthalpy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DEnthalpyDTemperature_CDensity(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DEnthalpyDTemperature_CPressure(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DEnthalpyDDensity_CTemperature(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DEnthalpyDDensity_CPressure(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DEnthalpyDPressure_CTemperature(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DEnthalpyDPressure_CDensity(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

// Internal Energy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DInternalEnergyDTemperature_CDensity(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DInternalEnergyDTemperature_CPressure(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DInternalEnergyDDensity_CTemperature(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DInternalEnergyDDensity_CPressure(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DInternalEnergyDPressure_CTemperature(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_DInternalEnergyDPressure_CDensity(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

// Compute Second Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_D2PressureDDensity2(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_D2PressureDTemperature2(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_DT_Get_D2PressureDTemperatureDensity(double dvDensity, double dvTemperature) {
    double result = 0.0;
    
    return result;
}

