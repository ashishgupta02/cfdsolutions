/*******************************************************************************
 * File:        EOS.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>
#include <math.h>

#include "NISTThermo.h"
#include "NISTThermo_Extension.h"
#include "EOS.h"
#include "EOS_Internal.h"
#include "EOS_Density_Pressure.h"
#include "EOS_Density_Temperature.h"
#include "EOS_Pressure_Temperature.h"
#include "Trim_Utils.h"

//------------------------------------------------------------------------------
//! Initialize EOS
//------------------------------------------------------------------------------
void EOS_Init() {
    int i;
    char *NIST_PATH;
    
    // Initialize the Internals of EOS library
    EOS_Internal_Init();
    
    // Get the path of NIST Fluids
    NIST_PATH = getenv("NIST_PATH");
    if (NIST_PATH != NULL) { 
        i = strlen(NIST_PATH);
        // Make sure that length does not exceed the storage
        if (i > NISTTHERMO_GEN_STR_LEN)
            i = NISTTHERMO_GEN_STR_LEN;
        strncpy(SogNISTHelper.csPath, NIST_PATH, i);
        SETPATH(SogNISTHelper.csPath, NISTTHERMO_GEN_STR_LEN);
    } else
        error("EOS_Init: Environment Variable NIST_PATH Not Found");
}

//------------------------------------------------------------------------------
//! Finalize EOS
//------------------------------------------------------------------------------
void EOS_Finalize() {
    
}

//------------------------------------------------------------------------------
//! Setup only for pure fluid
//  Improve to handle mixture fluids
//------------------------------------------------------------------------------
void EOS_Set() {
    // Setup for nitrogen for now
    SogNISTHelper.ivNumberComponent = 1;
    SogNISTHelper.ivError           = 0;
    strncpy(SogNISTHelper.csFiles, "nitrogen.fld", 12);
    strncpy(SogNISTHelper.csFmix, "hmx.bnc", 7);
    strncpy(SogNISTHelper.csRef, "DEF", 3);
    strncpy(SogNISTHelper.csError, "Ok", 2);
    
    // Call NISTThermo Setup
#if defined(WIN32) && !defined(IFORT)
    SETUP(&SogNISTHelper.ivNumberComponent, SogNISTHelper.csFiles, NISTTHERMO_FILE_STR_LEN, SogNISTHelper.csFmix, NISTTHERMO_GEN_STR_LEN, 
            SogNISTHelper.csRef, NISTTHERMO_REF_STR_LEN, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
#else
    SETUP(&SogNISTHelper.ivNumberComponent, SogNISTHelper.csFiles, SogNISTHelper.csFmix, SogNISTHelper.csRef, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_FILE_STR_LEN, 
            NISTTHERMO_GEN_STR_LEN, NISTTHERMO_REF_STR_LEN, NISTTHERMO_GEN_STR_LEN);
#endif
    if (SogNISTHelper.ivError != 0)
        error("EOS_Set: %s", SogNISTHelper.csError);
    
    // Set the Fluid Info Name
    strncpy(SogFluidInfo.csFluidName, "nitrogen.fld", 12);
    
    // Get the NIST Fluid Information
    EOS_Internal_Set_NIST_Fluid_Information(SogNISTHelper.ivNumberComponent);
    
    // Print the NIST Fluid Information
    EOS_Internal_Print_NIST_Fluid_Information(SogNISTHelper.ivNumberComponent);
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
    int ivKph, ivRefState;
    double dvDensity, dvDensityLiq, dvDensityVap, dvPressureNIST;
    double dvQuality, dvInternalEnergy, dvEnthalpy, dvEntropy, dvCv, dvCp, dvSpeedSound;
    double daX[NISTTHERMO_MAX_COMPONENT], daXLiq[NISTTHERMO_MAX_COMPONENT], daXVap[NISTTHERMO_MAX_COMPONENT];
    
    printf("=============================================================================\n");
    // Check for Equation of State Validity Bounds
    if ((dvTemperature <= SogFluidInfo.dvTemperature_Min) || (dvTemperature >= SogFluidInfo.dvTemperature_Max))
        error("EOS_Set_Reference_Properties:1: Reference Temperature Out of Bounds of Applicability");
    if (dvPressure >= SogFluidInfo.dvPressure_Max)
        error("EOS_Set_Reference_Properties:2: Reference Pressure Out of Bounds of Applicability");
    
    // Identify the Fluid Thermodynamic Region @ Reference Conditions.
    ivRefState = EOS_REG_UNKNOWN;
    // Region: Gas, Liquid or Multi-Phase
    if ((dvTemperature >= SogFluidInfo.dvTemperature_Crit) && (dvPressure >= SogFluidInfo.dvPressure_Crit)) {
        // Region: Super Critical State
        ivRefState = EOS_REG_SUPER_CRITICAL_STATE;
        info("Reference Properties in Region: Super Critical State");
    } else if ((dvTemperature < SogFluidInfo.dvTemperature_Crit) && (dvPressure >= SogFluidInfo.dvPressure_Crit)) {
        // Region: Subcooled Compressed Liquid
        ivRefState = EOS_REG_SUBCOOOLED_COMPRESSED_LIQUID;
        info("Reference Properties in Region: Subcooled Compressed Liquid");
    } else if((dvTemperature > SogFluidInfo.dvTemperature_Crit) && (dvPressure <= SogFluidInfo.dvPressure_Crit)) {
        // Region: Super Heated Vapor
        ivRefState = EOS_REG_SUPER_HEATED_VAPOR;
        info("Reference Properties in Region: Super Heated Vapor");
    } else if((dvTemperature <= SogFluidInfo.dvTemperature_Crit) && (dvPressure <= SogFluidInfo.dvPressure_Crit)) {
        // Region: Multiphase Region
        ivRefState = EOS_REG_MULTIPHASE;
        info("Reference Properties in Region: Multiphase Region");
    } else {
        // Region: Unknown
        ivRefState = EOS_REG_UNKNOWN;
        error("EOS_Set_Reference_Properties:3: Reference Temperature and Pressure Belong to UnKnown Fluid State");
    }
    
     // Convert Pressure from SI to NIST
    dvPressureNIST = EOS_Internal_Get_SI_To_NIST_Pressure(dvPressure);
    
    // Compute reference properties based on Fluid State Region
    if (ivRefState == EOS_REG_MULTIPHASE) {
        int icount;
        double tmp, dP, SSp, SSm;
        
        // Start the perturbation test on pressure: 0.1% fluctuation
        dP = 0.001*dvPressureNIST;
        
        // Verify if Temperature or Pressure are Saturated Values
        TPFLSH(&dvTemperature, &dvPressureNIST, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
                daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
                &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            warn("EOS_Set_Reference_Properties:1: %s", SogNISTHelper.csError);
        // Compute +dP property values
        dvPressureNIST += dP;
        TPFLSH(&dvTemperature, &dvPressureNIST, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
                daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
                &dvCv, &dvCp, &SSp, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            warn("EOS_Set_Reference_Properties:2: %s", SogNISTHelper.csError);
        // Compute -dP property values
        dvPressureNIST -= 2.0*dP;
        TPFLSH(&dvTemperature, &dvPressureNIST, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
                daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
                &dvCv, &dvCp, &SSm, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            warn("EOS_Set_Reference_Properties:3: %s", SogNISTHelper.csError);
        // Compute the change in speed of sound
        tmp = fabs(SSp - SSm)/dvSpeedSound;
        dvPressureNIST += dP; /* Original Input Pressure */
        // If Speed of Sound Changes more then 1%: Adjust the Reference Pressure to Saturated Pressure
        if (tmp > 0.01) {
            info("Saturated Reference Condition Detected !");
            info("Adjusting the Reference Pressure to Saturated Pressure");
            // Note: This process is not valid for Pseudo Fluids
            icount = 0;
            // Get the Saturated Pressure from given Temperature
            tmp = 0.0;
            ivKph = 1; /* For Non-Pseudo Fluids ivKph = 1 or 2 have same answer */
            SATT(&dvTemperature, daX, &ivKph, &dP, &dvDensityLiq, &dvDensityVap,
                    daXLiq, daXVap, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
            // Parse the error
            if (SogNISTHelper.ivError != 0) {
                warn("EOS_Set_Reference_Properties:4: Unable to Compute Saturated Pressure");
                warn(":4: %s", SogNISTHelper.csError);
            }
            // Compute the change in pressure
            // If change in pressure is more then 1% register error
            tmp = fabs(dP - dvPressureNIST)/dvPressureNIST;
            if (tmp > 0.01)
                icount++;
            
            // Get the Saturated Temperature from given Pressure
            tmp = 0.0;
            ivKph = 1; /* For Non-Pseudo Fluids ivKph = 1 or 2 have same answer */
            SATP(&dvPressureNIST, daX, &ivKph, &tmp, &dvDensityLiq, &dvDensityVap,
                    daXLiq, daXVap, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
            // Parse the error
            if (SogNISTHelper.ivError != 0) {
                warn("EOS_Set_Reference_Properties:5: Unable to Compute Saturated Temperature");
                warn(":5: %s", SogNISTHelper.csError);
            }
            // Compute the change in Temperature
            // If change in Temperature is more then 1% register error
            tmp = fabs(tmp - dvTemperature)/dvTemperature;
            if (tmp > 0.01)
                icount++;
            
            if (icount >= 1)
                error("EOS_Set_Reference_Properties:3: Unable to Adjust the Reference Pressure to Saturated Pressure");

            // Finally Adjust the Reference Pressure to Saturated Pressure
            dvPressureNIST = dP;
            info("Adjusted Reference Pressure: %10.4f kPa", dvPressureNIST);
        }
    }
    
    // Get the various properties
    // All Region: Multi-Phase, Super Heated Compressed Gas, Super Compressed Liquid, Super Heated Gas
    TPFLSH(&dvTemperature, &dvPressureNIST, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
            daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
            &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
    if (SogNISTHelper.ivError != 0)
        error("EOS_Set_Reference_Properties:4: %s", SogNISTHelper.csError);
    
    // Set the Reference Properties in SI Units
    dvgTemperature_Ref       = dvTemperature;
    dvgLength_Ref            = dvLength;
    dvgDensity_Ref           = EOS_Internal_Get_NIST_To_SI_Density(dvDensity);
    dvgSpeedSound_Ref        = dvSpeedSound;
    dvgRatioSpecificHeat_Ref = dvCp/dvCv;
    // Derived Reference Properties
    dvgPressure_Ref          = dvgDensity_Ref*dvgSpeedSound_Ref*dvgSpeedSound_Ref; // Pressure in Pa
    dvgVelocity_Ref          = dvgSpeedSound_Ref;
    dvgMach_Ref              = 1.0;
    dvgTime_Ref              = dvgLength_Ref/dvgSpeedSound_Ref;
    dvgEnthalpy_Ref          = dvgSpeedSound_Ref*dvgSpeedSound_Ref;
    dvgTotalEnthalpy_Ref     = dvgEnthalpy_Ref;
    dvgInternalEnergy_Ref    = dvgEnthalpy_Ref;
    dvgTotalEnergy_Ref       = dvgEnthalpy_Ref;
    dvgGasConstant_Ref       = dvgEnthalpy_Ref/dvgTemperature_Ref;
    dvgEntropy_Ref           = dvgGasConstant_Ref;
    dvgHeatCapacityCv_Ref    = dvgGasConstant_Ref;
    dvgHeatCapacityCp_Ref    = dvgGasConstant_Ref;
    dvgEntropyConst_Ref      = 1.0;
    ivgSet_Ref               = 1;
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
    info("Density_Ref                       = %15.6f kg/m^3", dvgDensity_Ref);
    info("Pressure_Ref                      = %15.6f kPa",    EOS_Internal_Get_SI_To_NIST_Pressure(dvgPressure_Ref)); // Pressure in kPa
    info("Temperature_Ref                   = %15.6f K",      dvgTemperature_Ref);
    info("Velocity_Ref                      = %15.6f m/s",    dvgVelocity_Ref);
    info("Length_Ref                        = %15.6f m",      dvgLength_Ref);
    info("SpeedSound_Ref                    = %15.6f m/s",    dvgSpeedSound_Ref);
    info("Mach_Ref                          = %15.6f",        dvgMach_Ref);
    info("Time_Ref                          = %15.6f s",      dvgTime_Ref);
    info("Entropy_Ref                       = %15.6f J/kg-K", dvgEntropy_Ref);
    info("Enthalpy_Ref                      = %15.6f J/kg",   dvgEnthalpy_Ref);
    info("InternalEnergy_Ref                = %15.6f J/kg",   dvgInternalEnergy_Ref);
    info("TotalEnthalpy_Ref                 = %15.6f J/kg",   dvgTotalEnthalpy_Ref);
    info("TotalEnergy_Ref                   = %15.6f J/kg",   dvgTotalEnergy_Ref);
    info("Heat Capacity Cv                  = %15.6f J/kg-K", dvgHeatCapacityCv_Ref);
    info("Heat Capacity Cp                  = %15.6f J/kg-K", dvgHeatCapacityCp_Ref);
    info("GasConstant_Ref                   = %15.6f J/kg-K", dvgGasConstant_Ref);
    info("RatioSpecificHeat_Ref             = %15.6f",        dvgRatioSpecificHeat_Ref);
    info("EntropyConst_Ref                  = %15.6f",        dvgEntropyConst_Ref);
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//! Compute the EOS Properties based on Variable Type
//! Input and Output property are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_Get_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    int i;
    double dvRho, dvRhoL, dvRhoV, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2, dvSpeedSound, dvMach;
    double dvEntropy, dvEnthalpy, dvInternalEnergy, dvTotalEnergy, dvTotalEnthalpy;
    double dvCv, dvCp;
    double daVariableDimensional[EOS_NEQUATIONS];
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST, dvQualityNIST;
    
    // Dimensionalize the Input Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            EOS_Internal_Dimensionalize_Variables(ivVariableType, dpVariableIn, daVariableDimensional);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            for (i = 0; i < EOS_NEQUATIONS; i++)
                daVariableDimensional[i] = dpVariableIn[i];
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            EOS_Internal_Dimensionalize_Variables(ivVariableType, dpVariableIn, daVariableDimensional);
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            for (i = 0; i < EOS_NEQUATIONS; i++)
                daVariableDimensional[i] = dpVariableIn[i];
            break;
        default:
            error("EOS_Get_Properties:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_Get_Properties:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvVelocityU = daVariableDimensional[1];
            dvVelocityV = daVariableDimensional[2];
            dvVelocityW = daVariableDimensional[3];
            dvPressure  = daVariableDimensional[4];
            dvQ2        = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            // Compute NIST compatible values
            dvRhoNIST      = EOS_Internal_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_Internal_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_Get_Properties:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_Get_Properties:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvVelocityU   = daVariableDimensional[1];
            dvVelocityV   = daVariableDimensional[2];
            dvVelocityW   = daVariableDimensional[3];
            dvTemperature = daVariableDimensional[4];
            dvQ2          = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            // Compute NIST compatible values
            dvRhoNIST     = EOS_Internal_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_Get_Properties:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Properties in NIST Units
    dpPropertyOut[ 0] = dvRhoNIST;
    dpPropertyOut[ 1] = dvRhoLNIST;
    dpPropertyOut[ 2] = dvRhoVNIST;
    dpPropertyOut[ 3] = dvPressureNIST;
    dpPropertyOut[ 4] = dvTemperature;
    dpPropertyOut[ 5] = dvSpeedSound;
    dpPropertyOut[ 6] = dvEntropy;
    dpPropertyOut[ 7] = dvEnthalpy;
    dpPropertyOut[ 8] = dvInternalEnergy;
    dpPropertyOut[ 9] = dvCv;
    dpPropertyOut[10] = dvCp;
    
    // Convert Properties From NIST to SI units
    EOS_Internal_Get_NIST_To_SI_Property(dpPropertyOut);
    
    // Get the Properties in SI units
    dvRho            = dpPropertyOut[ 0];
    dvRhoL           = dpPropertyOut[ 1];
    dvRhoV           = dpPropertyOut[ 2];
    dvPressure       = dpPropertyOut[ 3];
    dvTemperature    = dpPropertyOut[ 4];
    dvSpeedSound     = dpPropertyOut[ 5];
    dvEntropy        = dpPropertyOut[ 6];
    dvEnthalpy       = dpPropertyOut[ 7];
    dvInternalEnergy = dpPropertyOut[ 8];
    dvCv             = dpPropertyOut[ 9];
    dvCp             = dpPropertyOut[10];
    
    // Compute the Remaining Properties in SI units
    dvMach          = sqrt(dvQ2)/dvSpeedSound;
    dvTotalEnthalpy = dvEnthalpy + 0.5*dvQ2;
    dvTotalEnergy   = dvInternalEnergy + 0.5*dvQ2;
    
    // Update the Output (is Dimensional)
    dpPropertyOut[ 0] = dvRho;
    dpPropertyOut[ 1] = dvRhoL;
    dpPropertyOut[ 2] = dvRhoV;
    dpPropertyOut[ 3] = dvPressure;
    dpPropertyOut[ 4] = dvTemperature;
    dpPropertyOut[ 5] = dvVelocityU;
    dpPropertyOut[ 6] = dvVelocityV;
    dpPropertyOut[ 7] = dvVelocityW;
    dpPropertyOut[ 8] = dvQ2;
    dpPropertyOut[ 9] = dvSpeedSound;
    dpPropertyOut[10] = dvMach;
    dpPropertyOut[11] = dvEntropy;
    dpPropertyOut[12] = dvEnthalpy;
    dpPropertyOut[13] = dvInternalEnergy;
    dpPropertyOut[14] = dvTotalEnthalpy;
    dpPropertyOut[15] = dvTotalEnergy;
    dpPropertyOut[16] = dvCv;
    dpPropertyOut[17] = dvCp;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            EOS_Internal_NonDimensionalize_Properties(dpPropertyOut, dpPropertyOut);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            // Do Nothing Output is Dimensional (SI Units)
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            // Do Nothing Output is Dimensional (SI Units)
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            EOS_Internal_NonDimensionalize_Properties(dpPropertyOut, dpPropertyOut);
            break;
        default:
            error("EOS_Get_Properties:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
}

//------------------------------------------------------------------------------
//! Compute the EOS Extended Properties based on Variable Type
//! Input and Output property are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_Get_Extended_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    int i;
    double dvRho, dvRhoL, dvRhoV, dvPressure, dvTemperature;
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2, dvSpeedSound, dvMach;
    double dvEntropy, dvEnthalpy, dvInternalEnergy, dvTotalEnergy, dvTotalEnthalpy;
    double dvCv, dvCp, dvDPDRho, dvDPDT, dvDRhoDT, dvDRhoDP, dvD2PDRho2, dvD2PDT2, dvD2PDTRho;
    double dvDHdT_Rho, dvDHdT_P, dvDHDRho_T, dvDHDRho_P, dvDHDP_T, dvDHDP_Rho;
    double dvDEdT_Rho, dvDEdT_P, dvDEDRho_T, dvDEDRho_P, dvDEDP_T, dvDEDP_Rho;
    double daVariableDimensional[EOS_NEQUATIONS];
    
    // Dimensionalize the Input Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            EOS_Internal_Dimensionalize_Variables(ivVariableType, dpVariableIn, daVariableDimensional);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            for (i = 0; i < EOS_NEQUATIONS; i++)
                daVariableDimensional[i] = dpVariableIn[i];
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            EOS_Internal_Dimensionalize_Variables(ivVariableType, dpVariableIn, daVariableDimensional);
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            for (i = 0; i < EOS_NEQUATIONS; i++)
                daVariableDimensional[i] = dpVariableIn[i];
            break;
        default:
            error("EOS_Get_Properties:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_Get_Extended_Properties:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvVelocityU = daVariableDimensional[1];
            dvVelocityV = daVariableDimensional[2];
            dvVelocityW = daVariableDimensional[3];
            dvPressure  = daVariableDimensional[4];
            dvQ2        = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            // Compute NIST compatible values
            dvRhoNIST      = EOS_Internal_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_Internal_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Extended Quantities
            PDEPDFLSH2(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, SogNISTHelper.daXLiq, SogNISTHelper.daXVap,
                    dpPropertyOut, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_Get_Extended_Properties:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_Get_Extended_Properties:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvVelocityU   = daVariableDimensional[1];
            dvVelocityV   = daVariableDimensional[2];
            dvVelocityW   = daVariableDimensional[3];
            dvTemperature = daVariableDimensional[4];
            dvQ2          = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            // Compute NIST compatible values
            dvRhoNIST     = EOS_Internal_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH2(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, SogNISTHelper.daXLiq, SogNISTHelper.daXVap, 
                    dpPropertyOut, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_Get_Extended_Properties:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Properties in NIST units
    switch (ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_Get_Extended_Properties:6: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvTemperature = dpPropertyOut[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_Get_Extended_Properties:7: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_Get_Extended_Properties:8: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvPressureNIST = dpPropertyOut[0];
            break;
        default:
            error("EOS_Get_Extended_Properties:9: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    dvRhoLNIST       = dpPropertyOut[1];
    dvRhoVNIST       = dpPropertyOut[2];
    dvInternalEnergy = dpPropertyOut[4];
    dvEnthalpy       = dpPropertyOut[5];
    dvEntropy        = dpPropertyOut[6];
    dvCv             = dpPropertyOut[7];
    dvCp             = dpPropertyOut[8];
    dvSpeedSound     = dpPropertyOut[9];
    dvDPDRho         = dpPropertyOut[16];
    dvDPDT           = dpPropertyOut[17];
    dvDRhoDT         = dpPropertyOut[18];
    dvDRhoDP         = dpPropertyOut[19];
    dvD2PDRho2       = dpPropertyOut[20];
    dvD2PDT2         = dpPropertyOut[21];
    dvD2PDTRho       = dpPropertyOut[22];
    dvDHdT_Rho       = dpPropertyOut[23];
    dvDHdT_P         = dpPropertyOut[24];
    dvDHDRho_T       = dpPropertyOut[25];
    dvDHDRho_P       = dpPropertyOut[26];
    dvDHDP_T         = dpPropertyOut[27];
    dvDHDP_Rho       = dpPropertyOut[28];
    
    // Convert the Properties From NIST to SI units
    dpPropertyOut[ 0] = dvRhoNIST;
    dpPropertyOut[ 1] = dvRhoLNIST;
    dpPropertyOut[ 2] = dvRhoVNIST;
    dpPropertyOut[ 3] = dvPressureNIST;
    dpPropertyOut[ 4] = dvTemperature;
    dpPropertyOut[ 5] = dvSpeedSound;
    dpPropertyOut[ 6] = dvEntropy;
    dpPropertyOut[ 7] = dvEnthalpy;
    dpPropertyOut[ 8] = dvInternalEnergy;
    dpPropertyOut[ 9] = dvCv;
    dpPropertyOut[10] = dvCp;
    dpPropertyOut[11] = dvDPDRho;
    dpPropertyOut[12] = dvDPDT;
    dpPropertyOut[13] = dvDRhoDT;
    dpPropertyOut[14] = dvDRhoDP;
    dpPropertyOut[15] = dvD2PDRho2;
    dpPropertyOut[16] = dvD2PDT2;
    dpPropertyOut[17] = dvD2PDTRho;
    dpPropertyOut[18] = dvDHdT_Rho;
    dpPropertyOut[19] = dvDHdT_P;
    dpPropertyOut[20] = dvDHDRho_T;
    dpPropertyOut[21] = dvDHDRho_P;
    dpPropertyOut[22] = dvDHDP_T;
    dpPropertyOut[23] = dvDHDP_Rho;
    EOS_Internal_Get_NIST_To_SI_Extended_Property(dpPropertyOut);
    
    // Get the Properties in SI units
    dvRho            = dpPropertyOut[ 0];
    dvRhoL           = dpPropertyOut[ 1];
    dvRhoV           = dpPropertyOut[ 2];
    dvPressure       = dpPropertyOut[ 3];
    dvTemperature    = dpPropertyOut[ 4];
    dvSpeedSound     = dpPropertyOut[ 5];
    dvEntropy        = dpPropertyOut[ 6];
    dvEnthalpy       = dpPropertyOut[ 7];
    dvInternalEnergy = dpPropertyOut[ 8];
    dvCv             = dpPropertyOut[ 9];
    dvCp             = dpPropertyOut[10];
    dvDPDRho         = dpPropertyOut[11];
    dvDPDT           = dpPropertyOut[12];
    dvDRhoDT         = dpPropertyOut[13];
    dvDRhoDP         = dpPropertyOut[14];
    dvD2PDRho2       = dpPropertyOut[15];
    dvD2PDT2         = dpPropertyOut[16];
    dvD2PDTRho       = dpPropertyOut[17];
    dvDHdT_Rho       = dpPropertyOut[18];
    dvDHdT_P         = dpPropertyOut[19];
    dvDHDRho_T       = dpPropertyOut[20];
    dvDHDRho_P       = dpPropertyOut[21];
    dvDHDP_T         = dpPropertyOut[22];
    dvDHDP_Rho       = dpPropertyOut[23];
    
    // Compute the Remaining Properties in SI units
    dvMach          = sqrt(dvQ2)/dvSpeedSound;
    dvTotalEnthalpy = dvEnthalpy + 0.5*dvQ2;
    dvTotalEnergy   = dvInternalEnergy + 0.5*dvQ2;
    dvDEdT_Rho      = dvDHdT_Rho - dvDPDT/dvRho;
    dvDEdT_P        = dvDHdT_P + dvPressure*dvDRhoDT/(dvRho*dvRho);
    dvDEDRho_T      = dvDHDRho_T + dvPressure/(dvRho*dvRho) - dvDPDRho/dvRho;
    dvDEDRho_P      = dvDHDRho_P + dvPressure/(dvRho*dvRho);
    dvDEDP_T        = dvDHDP_T - 1.0/dvRho + dvPressure*dvDRhoDP/(dvRho*dvRho);
    dvDEDP_Rho      = dvDHDP_Rho - 1.0/dvRho;
    
    // Update the Output (is Dimensional)
    dpPropertyOut[ 0] = dvRho;
    dpPropertyOut[ 1] = dvRhoL;
    dpPropertyOut[ 2] = dvRhoV;
    dpPropertyOut[ 3] = dvPressure;
    dpPropertyOut[ 4] = dvTemperature;
    dpPropertyOut[ 5] = dvVelocityU;
    dpPropertyOut[ 6] = dvVelocityV;
    dpPropertyOut[ 7] = dvVelocityW;
    dpPropertyOut[ 8] = dvQ2;
    dpPropertyOut[ 9] = dvSpeedSound;
    dpPropertyOut[10] = dvMach;
    dpPropertyOut[11] = dvEntropy;
    dpPropertyOut[12] = dvEnthalpy;
    dpPropertyOut[13] = dvInternalEnergy;
    dpPropertyOut[14] = dvTotalEnthalpy;
    dpPropertyOut[15] = dvTotalEnergy;
    dpPropertyOut[16] = dvCv;
    dpPropertyOut[17] = dvCp;
    dpPropertyOut[18] = dvDPDRho;
    dpPropertyOut[19] = dvDPDT;
    dpPropertyOut[20] = dvDRhoDT;
    dpPropertyOut[21] = dvDRhoDP;
    dpPropertyOut[22] = dvD2PDRho2;
    dpPropertyOut[23] = dvD2PDT2;
    dpPropertyOut[24] = dvD2PDTRho;
    dpPropertyOut[25] = dvDHdT_Rho;
    dpPropertyOut[26] = dvDHdT_P;
    dpPropertyOut[27] = dvDHDRho_T;
    dpPropertyOut[28] = dvDHDRho_P;
    dpPropertyOut[29] = dvDHDP_T;
    dpPropertyOut[30] = dvDHDP_Rho;
    dpPropertyOut[31] = dvDEdT_Rho;
    dpPropertyOut[32] = dvDEdT_P;
    dpPropertyOut[33] = dvDEDRho_T;
    dpPropertyOut[34] = dvDEDRho_P;
    dpPropertyOut[35] = dvDEDP_T;
    dpPropertyOut[36] = dvDEDP_Rho;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            EOS_Internal_NonDimensionalize_Extended_Properties(dpPropertyOut, dpPropertyOut);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            // Do Nothing Output is Dimensional (SI Units)
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            // Do Nothing Output is Dimensional (SI Units)
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            EOS_Internal_NonDimensionalize_Extended_Properties(dpPropertyOut, dpPropertyOut);
            break;
        default:
            error("EOS_Get_Properties:10: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
}

//------------------------------------------------------------------------------
//! Compute the Transformation Matrix Based on Variable Type Input
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_Get_Transformation_Matrix(int ivDimIOType, int ivVarTypeIn, double *dpPropertyIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix) {
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
                                EOS_DP_Get_Transformation_Matrix_CON_To_RUP(ivDimIOType, dpPropertyIn, Matrix);
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
                                EOS_DP_Get_Transformation_Matrix_RUP_To_CON(ivDimIOType, dpPropertyIn, Matrix);
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
                                EOS_DT_Get_Transformation_Matrix_CON_To_RUP(ivDimIOType, dpPropertyIn, Matrix);
                                break;
                            case EOS_VARIABLE_PUT:
                                error("EOS_Get_Transformation_Matrix:15: Not Implemented For EOS_VARIABLE_PUT");
                                break;
                            case EOS_VARIABLE_PUS:
                                error("EOS_Get_Transformation_Matrix:16: Not Implemented For EOS_VARIABLE_PUS");
                                break;
                            case EOS_VARIABLE_RUT:
                                EOS_DT_Get_Transformation_Matrix_CON_To_RUT(ivDimIOType, dpPropertyIn, Matrix);
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
                                EOS_DT_Get_Transformation_Matrix_RUP_To_RUT(ivDimIOType, dpPropertyIn, Matrix);
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
                                EOS_DT_Get_Transformation_Matrix_RUT_To_CON(ivDimIOType, dpPropertyIn, Matrix);
                                break;
                            case EOS_VARIABLE_RUP:
                                EOS_DT_Get_Transformation_Matrix_RUT_To_RUP(ivDimIOType, dpPropertyIn, Matrix);
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

