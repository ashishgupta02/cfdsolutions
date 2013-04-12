/*******************************************************************************
 * File:        EOS.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>
#include <math.h>

#include "EOS.h"
#include "EOS_NIST.h"
#include "EOS_IdealGas.h"
#include "EOS_Internal.h"
#include "EOS_Density_Pressure.h"
#include "EOS_Density_Temperature.h"
#include "EOS_Pressure_Temperature.h"
#include "Trim_Utils.h"

//------------------------------------------------------------------------------
//! Initialize EOS
//------------------------------------------------------------------------------
void EOS_Init(int ivEOSModelType) {
    
    // Initialize the Internals of EOS library
    EOS_Internal_Init();
    
    // Set the EOS Model Type
    SogFluidInfo.ivEOSModelType = ivEOSModelType;
    
    // Setup According to the EOS Model Type
    switch (ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            EOS_IdealGas_Init();
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            EOS_NIST_Init();
            break;
        default:
            error("EOS_Init:1: Undefined Equation of State Model");
            break;
    }
}

//------------------------------------------------------------------------------
//! Finalize EOS
//------------------------------------------------------------------------------
void EOS_Finalize() {
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            EOS_IdealGas_Finalize();
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            EOS_NIST_Finalize();
            break;
        default:
            error("EOS_Finalize:1: Undefined Equation of State Model");
            break;
    }
    
    // Finalize the Internal EOS Library
    EOS_Internal_Finalize();
}

//------------------------------------------------------------------------------
//! Setup only for pure fluid or ideal gas
//  Improve to handle mixture fluids
//------------------------------------------------------------------------------
void EOS_Set(const char* csFluidName) {
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            EOS_IdealGas_Set();
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            EOS_NIST_Set(csFluidName);
            break;
        default:
            error("EOS_Set:1: Undefined Equation of State Model");
            break;
    }
}

//------------------------------------------------------------------------------
//! Get the NIST Fluid Information
//------------------------------------------------------------------------------
void EOS_Get_Fluid_Information(double *dpProperty) {
    dpProperty[ 0] = SogFluidInfo.dvMolecularWeight;
    dpProperty[ 1] = SogFluidInfo.dvGasConstant;
    dpProperty[ 2] = SogFluidInfo.dvTemperature_Crit;
    dpProperty[ 3] = SogFluidInfo.dvPressure_Crit;
    dpProperty[ 4] = SogFluidInfo.dvZ_Crit;
    dpProperty[ 5] = SogFluidInfo.dvDensity_Crit;
    dpProperty[ 6] = SogFluidInfo.dvTPTemperature;
    dpProperty[ 7] = SogFluidInfo.dvNBPTemperature;
    dpProperty[ 8] = SogFluidInfo.dvAcentric_Fact;
    dpProperty[ 9] = SogFluidInfo.dvDipole_Mom;
    dpProperty[10] = SogFluidInfo.dvTemperature_Min;
    dpProperty[11] = SogFluidInfo.dvTemperature_Max;
    dpProperty[12] = SogFluidInfo.dvPressure_Max;
    dpProperty[13] = SogFluidInfo.dvDensity_Max;
}

//------------------------------------------------------------------------------
//! Set the Reference Property for Fluid : Generic (Reference Density is Computed)
//  Pressure (Pa), Temperature (K) and Length (m)
//------------------------------------------------------------------------------
void EOS_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength) {
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            EOS_IdealGas_Set_Reference_Properties(dvPressure, dvTemperature, dvLength);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            EOS_NIST_Set_Reference_Properties(dvPressure, dvTemperature, dvLength);
            break;
        default:
            error("EOS_Set_Reference_Properties:1: Undefined Equation of State Model");
            break;
    }
}

//------------------------------------------------------------------------------
//! Get the Reference Property for Fluid
//  Note: Reference values are in SI units
//------------------------------------------------------------------------------
void EOS_Get_Reference_Properties(double *dpProperty) {
    if (ivgSet_Ref == 0)
        warn("EOS_Get_Reference_Properties:1: Reference Properties are Not Set");
    
    dpProperty[ 0] = dvgDensity_Ref;                 // Density
    dpProperty[ 1] = dvgPressure_Ref;                // Pressure in Pa
    dpProperty[ 2] = dvgTemperature_Ref;             // Temperature
    dpProperty[ 3] = dvgVelocity_Ref;                // Velocity
    dpProperty[ 4] = dvgLength_Ref;                  // Length
    dpProperty[ 5] = dvgSpeedSound_Ref;              // Speed of Sound
    dpProperty[ 6] = dvgMach_Ref;                    // Mach
    dpProperty[ 7] = dvgTime_Ref;                    // Time
    dpProperty[ 8] = dvgEntropy_Ref;                 // Entropy
    dpProperty[ 9] = dvgEnthalpy_Ref;                // Enthalpy
    dpProperty[10] = dvgInternalEnergy_Ref;          // Internal Energy
    dpProperty[11] = dvgTotalEnthalpy_Ref;           // Total Enthalpy
    dpProperty[12] = dvgTotalEnergy_Ref;             // Total Energy
    dpProperty[13] = dvgHeatCapacityCv_Ref;          // Heat Capacity Cv
    dpProperty[14] = dvgHeatCapacityCp_Ref;          // Heat Capacity Cp
    dpProperty[15] = dvgGasConstant_Ref;             // Gas Constant
    dpProperty[16] = dvgRatioSpecificHeat_Ref;       // Specific Heat Ratio
    dpProperty[17] = dvgEntropyConst_Ref;            // Entropy Constant
}

//------------------------------------------------------------------------------
//! Print the Reference Property for Fluid
//
//------------------------------------------------------------------------------
void EOS_Print_Reference_Properties() {
    if (ivgSet_Ref == 0)
        warn("EOS_Print_Reference_Properties:1: Reference Properties are Not Set");
    
    printf("=============================================================================\n");
    info("EOS: Reference Properties:");
    printf("-----------------------------------------------------------------------------\n");
    info("Density_Ref -------------------------------: %15.6f kg/m^3", dvgDensity_Ref);
    info("Pressure_Ref ------------------------------: %15.6f kPa",    dvgPressure_Ref/1000.0); // Pressure in kPa
    info("Temperature_Ref ---------------------------: %15.6f K",      dvgTemperature_Ref);
    info("Velocity_Ref ------------------------------: %15.6f m/s",    dvgVelocity_Ref);
    info("Length_Ref --------------------------------: %15.6f m",      dvgLength_Ref);
    info("SpeedSound_Ref ----------------------------: %15.6f m/s",    dvgSpeedSound_Ref);
    info("Mach_Ref ----------------------------------: %15.6f",        dvgMach_Ref);
    info("Time_Ref ----------------------------------: %15.6f s",      dvgTime_Ref);
    info("Entropy_Ref -------------------------------: %15.6f J/kg-K", dvgEntropy_Ref);
    info("Enthalpy_Ref ------------------------------: %15.6f J/kg",   dvgEnthalpy_Ref);
    info("InternalEnergy_Ref ------------------------: %15.6f J/kg",   dvgInternalEnergy_Ref);
    info("TotalEnthalpy_Ref -------------------------: %15.6f J/kg",   dvgTotalEnthalpy_Ref);
    info("TotalEnergy_Ref ---------------------------: %15.6f J/kg",   dvgTotalEnergy_Ref);
    info("Heat Capacity Cv_Ref ----------------------: %15.6f J/kg-K", dvgHeatCapacityCv_Ref);
    info("Heat Capacity Cp_Ref ----------------------: %15.6f J/kg-K", dvgHeatCapacityCp_Ref);
    info("GasConstant_Ref ---------------------------: %15.6f J/kg-K", dvgGasConstant_Ref);
    info("RatioSpecificHeat_Ref ---------------------: %15.6f",        dvgRatioSpecificHeat_Ref);
    info("EntropyConst_Ref --------------------------: %15.6f",        dvgEntropyConst_Ref);
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//! Compute the EOS Properties based on Variable Type
//! Input and Output property are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_Get_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            EOS_IdealGas_Get_Properties(ivDimIOType, ivVariableType, dpVariableIn, dpPropertyOut);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            EOS_NIST_Get_Properties(ivDimIOType, ivVariableType, dpVariableIn, dpPropertyOut);
            break;
        default:
            error("EOS_Get_Properties:1: Undefined Equation of State Model");
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute the EOS Extended Properties based on Variable Type
//! Input and Output property are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_Get_Extended_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            EOS_IdealGas_Get_Extended_Properties(ivDimIOType, ivVariableType, dpVariableIn, dpPropertyOut);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            EOS_NIST_Get_Extended_Properties(ivDimIOType, ivVariableType, dpVariableIn, dpPropertyOut);
            break;
        default:
            error("EOS_Get_Extended_Properties:1: Undefined Equation of State Model");
            break;
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_Density(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_Density(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_Density(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_Density:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DensityLiquid(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DensityLiquid(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DensityLiquid(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DensityLiquid:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DensityVapor(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DensityVapor(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DensityVapor(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DensityVapor:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_Pressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_Pressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_Pressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_Pressure:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_Temperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_Temperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_Temperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_Temperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_SpeedSound(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_SpeedSound(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_SpeedSound(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_SpeedSound:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_Mach(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_Mach(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_Mach(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_Mach:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_Entropy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_Entropy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_Entropy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_Entropy:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_Enthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_Enthalpy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_Enthalpy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_Enthalpy:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_InternalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_InternalEnergy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_InternalEnergy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_InternalEnergy:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_TotalEnthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_TotalEnthalpy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_TotalEnthalpy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_TotalEnthalpy:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_TotalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_TotalEnergy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_TotalEnergy(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_TotalEnergy:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_HeatCapacityCv(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_HeatCapacityCv(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_HeatCapacityCv(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_HeatCapacityCv:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_HeatCapacityCp(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_HeatCapacityCp(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_HeatCapacityCp(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_HeatCapacityCp:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

// Compute First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DPressureDDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DPressureDDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DPressureDDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DPressureDDensity:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DPressureDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DPressureDTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DPressureDTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DPressureDTemperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DDensityDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DDensityDTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DDensityDTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DDensityDTemperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

// Compute Second Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_D2PressureDDensity2(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_D2PressureDDensity2(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_D2PressureDDensity2(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_D2PressureDDensity2:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_D2PressureDTemperature2(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_D2PressureDTemperature2(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_D2PressureDTemperature2(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_D2PressureDTemperature2:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_D2PressureDTemperatureDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_D2PressureDTemperatureDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_D2PressureDTemperatureDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_D2PressureDTemperatureDensity:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

// Enthalpy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DEnthalpyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DEnthalpyDTemperature_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DEnthalpyDTemperature_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DEnthalpyDTemperature_CDensity:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DEnthalpyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DEnthalpyDTemperature_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DEnthalpyDTemperature_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DEnthalpyDTemperature_CPressure:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DEnthalpyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DEnthalpyDDensity_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DEnthalpyDDensity_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DEnthalpyDDensity_CTemperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DEnthalpyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DEnthalpyDDensity_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DEnthalpyDDensity_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DEnthalpyDDensity_CPressure:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DEnthalpyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DEnthalpyDPressure_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DEnthalpyDPressure_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DEnthalpyDPressure_CTemperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DEnthalpyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DEnthalpyDPressure_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DEnthalpyDPressure_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DEnthalpyDPressure_CDensity:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

// Internal Energy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DInternalEnergyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DInternalEnergyDTemperature_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DInternalEnergyDTemperature_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DInternalEnergyDTemperature_CDensity:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DInternalEnergyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DInternalEnergyDTemperature_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DInternalEnergyDTemperature_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DInternalEnergyDTemperature_CPressure:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DInternalEnergyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DInternalEnergyDDensity_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DInternalEnergyDDensity_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DInternalEnergyDDensity_CTemperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DInternalEnergyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DInternalEnergyDDensity_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DInternalEnergyDDensity_CPressure(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DInternalEnergyDDensity_CPressure:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DInternalEnergyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DInternalEnergyDPressure_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DInternalEnergyDPressure_CTemperature(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DInternalEnergyDPressure_CTemperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DInternalEnergyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DInternalEnergyDPressure_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DInternalEnergyDPressure_CDensity(ivDimIOType, ivVariableType, dpVariableIn);
            break;
        default:
            error("EOS_Get_DInternalEnergyDPressure_CDensity:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void EOS_Get_PT_Density(int ivDimIOType, double dvPressure, double dvTemperature, double *dpDensity) {
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            EOS_IdealGas_Get_PT_Density(ivDimIOType, dvPressure, dvTemperature, dpDensity);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            EOS_NIST_Get_PT_Density(ivDimIOType, dvPressure, dvTemperature, dpDensity);
            break;
        default:
            error("EOS_Get_PT_Density:1: Undefined Equation of State Model");
            break;
    }
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DH_SpeedSound(int ivDimIOType, double dvDensity, double dvEnthalpy) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DH_SpeedSound(ivDimIOType, dvDensity, dvEnthalpy);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DH_SpeedSound(ivDimIOType, dvDensity, dvEnthalpy);
            break;
        default:
            error("EOS_Get_DH_SpeedSound:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DH_Pressure(int ivDimIOType, double dvDensity, double dvEnthalpy) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DH_Pressure(ivDimIOType, dvDensity, dvEnthalpy);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DH_Pressure(ivDimIOType, dvDensity, dvEnthalpy);
            break;
        default:
            error("EOS_Get_DH_Pressure:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_Get_DH_Temperature(int ivDimIOType, double dvDensity, double dvEnthalpy) {
    double result = 0.0;
    
    // Setup According to the EOS Model Type
    switch (SogFluidInfo.ivEOSModelType) {
        // Thermally Ideal Gas Model
        case EOS_MODEL_IDEALGAS:
            result = EOS_IdealGas_Get_DH_Temperature(ivDimIOType, dvDensity, dvEnthalpy);
            break;
        // Fluid Models provided by NIST
        case EOS_MODEL_NIST:
            result = EOS_NIST_Get_DH_Temperature(ivDimIOType, dvDensity, dvEnthalpy);
            break;
        default:
            error("EOS_Get_DH_Temperature:1: Undefined Equation of State Model");
            break;
    }
    
    return result;
}

//------------------------------------------------------------------------------
//! Compute the Transformation Matrix Based on Variable Type Input
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_Get_Transformation_Matrix(int ivDimIOType, int ivVarTypeIn, double *dpVariableIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix) {
    int i, j;
    
    // Check the Variable Type and Return the proper Transformation Matrix
    // Return Identity
    if (ivVarTypeFrom == ivVarTypeTo) {
        for (i = 0; i < EOS_NEQUATIONS; i++) {
            for (j = 0; j < EOS_NEQUATIONS; j++) {
                if (i == j)
                    Matrix[i][j] = 1.0;
                else
                    Matrix[i][j] = 0.0;
            }
        }
    } else {
        switch (ivVarTypeIn) {
            case EOS_VARIABLE_CON:
                error("EOS_Get_Transformation_Matrix:1: Not Implemented For EOS_VARIABLE_CON");
                break;
            case EOS_VARIABLE_RUP:
                switch (ivVarTypeFrom) {
                    case EOS_VARIABLE_CON:
                        switch (ivVarTypeTo) {
                            case EOS_VARIABLE_RUP:
                                EOS_DP_Get_Transformation_Matrix_CON_To_RUP(ivDimIOType, dpVariableIn, Matrix);
                                break;
                            case EOS_VARIABLE_PUT:
                                error("EOS_Get_Transformation_Matrix:2: Not Implemented For EOS_VARIABLE_PUT");
                                break;
                            case EOS_VARIABLE_PUS:
                                error("EOS_Get_Transformation_Matrix:3: Not Implemented For EOS_VARIABLE_PUS");
                                break;
                            case EOS_VARIABLE_RUT:
                                error("EOS_Get_Transformation_Matrix:4: Not Implemented For EOS_VARIABLE_RUT");
                                break;
                            default:
                                error("EOS_Get_Transformation_Matrix:5: Undefined Variable Type - %d", ivVarTypeTo);
                        }
                        break;
                    case EOS_VARIABLE_RUP:
                        switch (ivVarTypeTo) {
                            case EOS_VARIABLE_CON:
                                EOS_DP_Get_Transformation_Matrix_RUP_To_CON(ivDimIOType, dpVariableIn, Matrix);
                                break;
                            case EOS_VARIABLE_PUT:
                                error("EOS_Get_Transformation_Matrix:6: Not Implemented For EOS_VARIABLE_PUT");
                                break;
                            case EOS_VARIABLE_PUS:
                                error("EOS_Get_Transformation_Matrix:7: Not Implemented For EOS_VARIABLE_PUS");
                                break;
                            case EOS_VARIABLE_RUT:
                                error("EOS_Get_Transformation_Matrix:8: Not Implemented For EOS_VARIABLE_RUT");
                                break;
                            default:
                                error("EOS_Get_Transformation_Matrix:9: Undefined Variable Type - %d", ivVarTypeTo);
                        }
                        break;
                    case EOS_VARIABLE_PUT:
                        error("EOS_Get_Transformation_Matrix:10: Not Implemented For EOS_VARIABLE_PUT");
                        break;
                    case EOS_VARIABLE_PUS:
                        error("EOS_Get_Transformation_Matrix:11: Not Implemented For EOS_VARIABLE_PUS");
                        break;
                    case EOS_VARIABLE_RUT:
                        error("EOS_Get_Transformation_Matrix:12: Not Implemented For EOS_VARIABLE_RUT");
                        break;
                    default:
                        error("EOS_Get_Transformation_Matrix:13: Undefined Variable Type - %d", ivVarTypeFrom);
                        break;
                }
                break;
            case EOS_VARIABLE_PUT:
                error("EOS_Get_Transformation_Matrix:14: Not Implemented For EOS_VARIABLE_PUT");
                break;
            case EOS_VARIABLE_RUT:
                switch (ivVarTypeFrom) {
                    case EOS_VARIABLE_CON:
                        switch (ivVarTypeTo) {
                            case EOS_VARIABLE_RUP:
                                EOS_DT_Get_Transformation_Matrix_CON_To_RUP(ivDimIOType, dpVariableIn, Matrix);
                                break;
                            case EOS_VARIABLE_PUT:
                                error("EOS_Get_Transformation_Matrix:15: Not Implemented For EOS_VARIABLE_PUT");
                                break;
                            case EOS_VARIABLE_PUS:
                                error("EOS_Get_Transformation_Matrix:16: Not Implemented For EOS_VARIABLE_PUS");
                                break;
                            case EOS_VARIABLE_RUT:
                                EOS_DT_Get_Transformation_Matrix_CON_To_RUT(ivDimIOType, dpVariableIn, Matrix);
                                break;
                            default:
                                error("EOS_Get_Transformation_Matrix:17: Undefined Variable Type - %d", ivVarTypeTo);
                        }
                        break;
                    case EOS_VARIABLE_RUP:
                        switch (ivVarTypeTo) {
                            case EOS_VARIABLE_CON:
                                error("EOS_Get_Transformation_Matrix:18: Not Implemented For EOS_VARIABLE_CON");
                                break;
                            case EOS_VARIABLE_PUT:
                                error("EOS_Get_Transformation_Matrix:19: Not Implemented For EOS_VARIABLE_PUT");
                                break;
                            case EOS_VARIABLE_PUS:
                                error("EOS_Get_Transformation_Matrix:20: Not Implemented For EOS_VARIABLE_PUS");
                                break;
                            case EOS_VARIABLE_RUT:
                                EOS_DT_Get_Transformation_Matrix_RUP_To_RUT(ivDimIOType, dpVariableIn, Matrix);
                                break;
                            default:
                                error("EOS_Get_Transformation_Matrix:21: Undefined Variable Type - %d", ivVarTypeTo);
                        }
                        break;
                    case EOS_VARIABLE_PUT:
                        error("EOS_Get_Transformation_Matrix:22: Not Implemented For EOS_VARIABLE_PUT");
                        break;
                    case EOS_VARIABLE_PUS:
                        error("EOS_Get_Transformation_Matrix:23: Not Implemented For EOS_VARIABLE_PUS");
                        break;
                    case EOS_VARIABLE_RUT:
                        switch (ivVarTypeTo) {
                            case EOS_VARIABLE_CON:
                                EOS_DT_Get_Transformation_Matrix_RUT_To_CON(ivDimIOType, dpVariableIn, Matrix);
                                break;
                            case EOS_VARIABLE_RUP:
                                EOS_DT_Get_Transformation_Matrix_RUT_To_RUP(ivDimIOType, dpVariableIn, Matrix);
                                break;
                            case EOS_VARIABLE_PUT:
                                error("EOS_Get_Transformation_Matrix:24: Not Implemented For EOS_VARIABLE_PUT");
                                break;
                            case EOS_VARIABLE_PUS:
                                error("EOS_Get_Transformation_Matrix:25: Not Implemented For EOS_VARIABLE_PUS");
                                break;
                            default:
                                error("EOS_Get_Transformation_Matrix:26: Undefined Variable Type - %d", ivVarTypeTo);
                        }
                        break;
                    default:
                        error("EOS_Get_Transformation_Matrix:27: Undefined Variable Type - %d", ivVarTypeFrom);
                        break;
                }
                break;
            case EOS_VARIABLE_PUS:
                error("EOS_Get_Transformation_Matrix:28: Not Implemented For EOS_VARIABLE_PUS");
                break;
            default:
                error("EOS_Get_Transformation_Matrix:29: Undefined Variable Type - %d - Error-6", ivVarTypeIn);
                break;
        }
    }
}

