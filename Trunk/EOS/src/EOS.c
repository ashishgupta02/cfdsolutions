/*******************************************************************************
 * File:        EOS.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>
#include <math.h>

#include "NISTThermo.h"
#include "EOS_Internal.h"
#include "EOS.h"
#include "Trim_Utils.h"

// Static Shared Data
const int ncmax = 20; // Note: ncmax is the max number of components
const int numparams = 72;
const int maxcoefs = 50;



//------------------------------------------------------------------------------
//!
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
//!
//------------------------------------------------------------------------------
void EOS_Finalize() {
    
}

//------------------------------------------------------------------------------
//! Setup only for pure fluid
// Improve to handle mixture fluids
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
//! Set the Reference Property for Fluid
//
//------------------------------------------------------------------------------
void EOS_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength) {
    int ivKph, ivRefState;
    double dvDensity, dvDensityLiq, dvDensityVap;
    double dvQuality, dvInternalEnergy, dvEnthalpy, dvEntropy, dvCv, dvCp, dvSpeedSound;
    double daX[NISTTHERMO_MAX_COMPONENT], daXLiq[NISTTHERMO_MAX_COMPONENT], daXVap[NISTTHERMO_MAX_COMPONENT];
    
    printf("=============================================================================\n");
    // Check for Equation of State Validity Bounds
    if ((dvTemperature <= SogFluidInfo.dvTemperature_Min) || (dvTemperature >= SogFluidInfo.dvTemperature_Max))
        error("EOS_Set_Reference_Properties:1: Reference Temperature Out of Bounds of Applicability");
    if (dvPressure >= SogFluidInfo.dvPressure_Max)
        error("EOS_Set_Reference_Properties:2: Reference Pressure Out of Bounds of Applicability");
    
    // Identify the Fluid Thermodynamic Region @ Reference Conditions.
    ivRefState = EOS_REG_UKN;
    // Region: Gas, Liquid or Multi-Phase
    if ((dvTemperature >= SogFluidInfo.dvTemperature_Crit) && (dvPressure >= SogFluidInfo.dvPressure_Crit)) {
        // Region: Super Critical State
        ivRefState = EOS_REG_SCS;
        info("Reference Properties in Region: Super Critical State");
    } else if ((dvTemperature < SogFluidInfo.dvTemperature_Crit) && (dvPressure >= SogFluidInfo.dvPressure_Crit)) {
        // Region: Subcooled Compressed Liquid
        ivRefState = EOS_REG_SCL;
        info("Reference Properties in Region: Subcooled Compressed Liquid");
    } else if((dvTemperature > SogFluidInfo.dvTemperature_Crit) && (dvPressure <= SogFluidInfo.dvPressure_Crit)) {
        // Region: Super Heated Vapor
        ivRefState = EOS_REG_SHV;
        info("Reference Properties in Region: Super Heated Vapor");
    } else if((dvTemperature <= SogFluidInfo.dvTemperature_Crit) && (dvPressure <= SogFluidInfo.dvPressure_Crit)) {
        // Region: Multiphase Region
        ivRefState = EOS_REG_MP;
        info("Reference Properties in Region: Multiphase Region");
    } else {
        // Region: Unknown
        ivRefState = EOS_REG_UKN;
        error("EOS_Set_Reference_Properties:3: Reference Temperature and Pressure Belong to UnKnown Fluid State");
    }
    
    // Compute reference properties based on Fluid State Region
    if (ivRefState == EOS_REG_MP) {
        int icount;
        double tmp, dP, SSp, SSm;
        
        // Start the perturbation test on pressure: 0.1% fluctuation
        dP = 0.001*dvPressure;
        
        // Verify if Temperature or Pressure are Saturated Values
        TPFLSH(&dvTemperature, &dvPressure, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
                daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
                &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            warn("EOS_Set_Reference_Properties:1: %s", SogNISTHelper.csError);
        // Compute +dP property values
        dvPressure += dP;
        TPFLSH(&dvTemperature, &dvPressure, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
                daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
                &dvCv, &dvCp, &SSp, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            warn("EOS_Set_Reference_Properties:2: %s", SogNISTHelper.csError);
        // Compute -dP property values
        dvPressure -= 2.0*dP;
        TPFLSH(&dvTemperature, &dvPressure, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
                daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
                &dvCv, &dvCp, &SSm, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
        if (SogNISTHelper.ivError != 0)
            warn("EOS_Set_Reference_Properties:3: %s", SogNISTHelper.csError);
        // Compute the change in speed of sound
        tmp = fabs(SSp - SSm)/dvSpeedSound;
        dvPressure += dP; /* Original Input Pressure */
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
            tmp = fabs(dP - dvPressure)/dvPressure;
            if (tmp > 0.01)
                icount++;
            
            // Get the Saturated Temperature from given Pressure
            tmp = 0.0;
            ivKph = 1; /* For Non-Pseudo Fluids ivKph = 1 or 2 have same answer */
            SATP(&dvPressure, daX, &ivKph, &tmp, &dvDensityLiq, &dvDensityVap,
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
            dvPressure = dP;
            info("Adjusted Reference Pressure: %10.4f", dvPressure);
        }
    }
    
    // Get the various properties
    // All Region: Multi-Phase, Super Heated Compressed Gas, Super Compressed Liquid, Super Heated Gas
    TPFLSH(&dvTemperature, &dvPressure, daX, &dvDensity, &dvDensityLiq, &dvDensityVap,
            daXLiq, daXVap, &dvQuality, &dvInternalEnergy, &dvEnthalpy, &dvEntropy, 
            &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
    if (SogNISTHelper.ivError != 0)
        error("EOS_Set_Reference_Properties:4: %s", SogNISTHelper.csError);
    
    // Set the Reference Properties
    dvTemperature_Ref       = dvTemperature;
    dvPressure_Ref          = dvPressure;
    dvLength_Ref            = dvLength;
    dvDensity_Ref           = dvDensity;
    dvRatioSpecificHeat_Ref = dvCp/dvCv;
    dvGasConstant_Ref       = SogFluidInfo.dvGasConstant;
    dvSpeedSound_Ref        = dvSpeedSound;
    dvVelocity_Ref          = dvSpeedSound;
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
//! Get the Reference Property for Fluid
//
//------------------------------------------------------------------------------
void EOS_Get_Reference_Properties() {
    printf("=============================================================================\n");
    info("Temperature_Ref                   = %10.4e", dvTemperature_Ref);
    info("Pressure_Ref                      = %10.4e", dvPressure_Ref);
    info("Length_Ref                        = %10.4e", dvLength_Ref);
    info("Density_Ref                       = %10.4e", dvDensity_Ref);
    info("RatioSpecificHeat_Ref             = %10.4e", dvRatioSpecificHeat_Ref);
    info("GasConstant_Ref                   = %10.4e", dvGasConstant_Ref);
    info("SpeedSound_Ref                    = %10.4e", dvSpeedSound_Ref);
    info("Velocity_Ref                      = %10.4e", dvVelocity_Ref);
    info("Mach_Ref                          = %10.4e", dvMach_Ref);
    info("Time_Ref                          = %10.4e", dvTime_Ref);
    info("Enthalpy_Ref                      = %10.4e", dvEnthalpy_Ref);
    info("TotalEnthalpy_Ref                 = %10.4e", dvTotalEnthalpy_Ref);
    info("InternalEnergy_Ref                = %10.4e", dvInternalEnergy_Ref);
    info("TotalEnergy_Ref                   = %10.4e", dvTotalEnergy_Ref);
    info("Entropy_Ref                       = %10.4e", dvEntropy_Ref);
    info("EntropyConst_Ref                  = %10.4e", dvEntropyConst_Ref);
    printf("=============================================================================\n");
}

//------------------------------------------------------------------------------
//! Compute the EOS Properties based on Variable Type
//
//------------------------------------------------------------------------------
void EOS_Compute_Properties(int ivVariableType, double *dpPropertyIn, double *dpPropertyOut) {
    double dvRho, dvTemperature, dvPressure, dvVelocityU, dvVelocityV, dvVelocityW ;
    double dvQ2, dvTotalEnergy, dvTotalEnthalpy, dvSpeedSound, dvMach;
    double dvInternalEnergy, dvEnthalpy;
    double daPropertyDimensional[5];
    
    // Dimensionalize the Input Properties
    EOS_Internal_Dimensionalize_Properties(ivVariableType, dpPropertyIn, daPropertyDimensional);
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            dvRho            = daPropertyDimensional[0];
            dvVelocityU      = daPropertyDimensional[1]/dvRho;
            dvVelocityV      = daPropertyDimensional[2]/dvRho;
            dvVelocityW      = daPropertyDimensional[3]/dvRho;
            dvTotalEnergy    = daPropertyDimensional[4]/dvRho;
            dvQ2             = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvInternalEnergy = dvTotalEnergy - 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure    = daPropertyDimensional[0];
            dvVelocityU   = daPropertyDimensional[1];
            dvVelocityV   = daPropertyDimensional[2];
            dvVelocityW   = daPropertyDimensional[3];
            dvTemperature = daPropertyDimensional[4];
            dvQ2          = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daPropertyDimensional[0];
            dvVelocityU = daPropertyDimensional[1];
            dvVelocityV = daPropertyDimensional[2];
            dvVelocityW = daPropertyDimensional[3];
            dvPressure  = daPropertyDimensional[4];
            dvQ2        = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daPropertyDimensional[0];
            dvVelocityU   = daPropertyDimensional[1];
            dvVelocityV   = daPropertyDimensional[2];
            dvVelocityW   = daPropertyDimensional[3];
            dvTemperature = daPropertyDimensional[4];
            dvQ2          = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            break;
        default:
            error("EOS_Compute_Properties: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    
    // Update the Output
    dpPropertyOut[ 0] = dvRho;
    dpPropertyOut[ 1] = dvPressure;
    dpPropertyOut[ 2] = dvTemperature;
    dpPropertyOut[ 3] = dvTotalEnergy;
    dpPropertyOut[ 4] = dvTotalEnthalpy;
    dpPropertyOut[ 5] = dvSpeedSound;
    dpPropertyOut[ 6] = dvMach;
    dpPropertyOut[ 7] = dvVelocityU;
    dpPropertyOut[ 8] = dvVelocityV;
    dpPropertyOut[ 9] = dvVelocityW;
    dpPropertyOut[10] = dvQ2;
}

