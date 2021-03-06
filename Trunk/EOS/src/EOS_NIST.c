/*******************************************************************************
 * File:        EOS_NIST.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "NISTThermo_Extension.h"
#include "EOS.h"
#include "EOS_NIST.h"
#include "EOS_Internal.h"
#include "Trim_Utils.h"

/* Common Shared Data */
/* NIST Thermodynamics Input-Output Helpers */
EOS_S_NISTHelper SogNISTHelper;

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void EOS_NIST_Init() {
    int i;
    char *NIST_PATH;
    
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
        error("EOS_NIST_Init: Environment Variable NIST_PATH Not Found");
}

//------------------------------------------------------------------------------
//! Finalize EOS NIST
//------------------------------------------------------------------------------
void EOS_NIST_Finalize() {
    
}

//------------------------------------------------------------------------------
//! Setup only for pure fluid
//  Improve to handle mixture fluids
//------------------------------------------------------------------------------
void EOS_NIST_Set(const char* csFluidName) {
    
    // Fluid Models provided by NIST
    SogNISTHelper.ivNumberComponent = 1;
    SogNISTHelper.ivError           = 0;
    strncpy(SogNISTHelper.csFiles, csFluidName, strlen(csFluidName));
    strcat(SogNISTHelper.csFiles,  ".fld");
    strncpy(SogNISTHelper.csFmix,  "hmx.bnc", 7);
    strncpy(SogNISTHelper.csRef,   "DEF",     3);
    strncpy(SogNISTHelper.csError, "Ok",      2);

    // Call NISTThermo Setup
#if defined(WIN32) && !defined(IFORT)
    SETUP(&SogNISTHelper.ivNumberComponent, SogNISTHelper.csFiles, NISTTHERMO_FILE_STR_LEN, SogNISTHelper.csFmix, 
            NISTTHERMO_GEN_STR_LEN, SogNISTHelper.csRef, NISTTHERMO_REF_STR_LEN, &SogNISTHelper.ivError, 
            SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
#else
    SETUP(&SogNISTHelper.ivNumberComponent, SogNISTHelper.csFiles, SogNISTHelper.csFmix, SogNISTHelper.csRef, 
            &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_FILE_STR_LEN, 
            NISTTHERMO_GEN_STR_LEN, NISTTHERMO_REF_STR_LEN, NISTTHERMO_GEN_STR_LEN);
#endif
    if (SogNISTHelper.ivError != 0)
        error("EOS_Set:1: %s", SogNISTHelper.csError);

    // Set the Fluid Name
    strncpy(SogFluidInfo.csFluidName, SogNISTHelper.csFiles , strlen(SogNISTHelper.csFiles));

    // Get the NIST Fluid Information
    EOS_NIST_Set_Fluid_Information(SogNISTHelper.ivNumberComponent);

    // Print the NIST Fluid Information
    EOS_NIST_Print_Fluid_Information(SogNISTHelper.ivNumberComponent);
}

//------------------------------------------------------------------------------
//! Initialize Fluid Information from NIST Thermodynamic Library
//  All Thermodynamic variables are in SI Units
//------------------------------------------------------------------------------
void EOS_NIST_Set_Fluid_Information(int ivComp) {
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
        error("EOS_Set_NIST_Fluid_Information: %s", SogNISTHelper.csError);
    
    // Convert into International Standard Units (SI)
    SogFluidInfo.dvPressure_Crit = EOS_NIST_Get_NIST_To_SI_Pressure(SogFluidInfo.dvPressure_Crit);
    SogFluidInfo.dvPressure_Max  = EOS_NIST_Get_NIST_To_SI_Pressure(SogFluidInfo.dvPressure_Max);
    SogFluidInfo.dvDensity_Crit  = EOS_NIST_Get_NIST_To_SI_Density(SogFluidInfo.dvDensity_Crit);
    SogFluidInfo.dvDensity_Max   = EOS_NIST_Get_NIST_To_SI_Density(SogFluidInfo.dvDensity_Max);
    SogFluidInfo.dvGasConstant  /= 0.001*SogFluidInfo.dvMolecularWeight;
}

//------------------------------------------------------------------------------
//! Print Fluid Information from NIST Thermodynamic Library
//------------------------------------------------------------------------------
void EOS_NIST_Print_Fluid_Information(int ivComp) {
    if (ivComp != 1)
        error("Only Single fluid implemented !");
    printf("=============================================================================\n");
    info("FLUID INFORMATION : %s", SogFluidInfo.csFluidName);
    printf("-----------------------------------------------------------------------------\n");
    info("Molecular Weight --------------------------: %15.6f kg/kmol", SogFluidInfo.dvMolecularWeight);
    info("Gas Constant ------------------------------: %15.6f J/kg-K",  SogFluidInfo.dvGasConstant);
    info("Temperature Critical ----------------------: %15.6f K",       SogFluidInfo.dvTemperature_Crit);
    info("Pressure Critical -------------------------: %15.6f kPa",     EOS_NIST_Get_SI_To_NIST_Pressure(SogFluidInfo.dvPressure_Crit));
    info("Compressibility @ Critical Point ----------: %15.6f ",        SogFluidInfo.dvZ_Crit);
    info("Density Critical --------------------------: %15.6f kg/m^3",  SogFluidInfo.dvDensity_Crit);
    info("Triple Point Temperature ------------------: %15.6f K",       SogFluidInfo.dvTPTemperature);
    info("Normal Boiling Point Temperature ----------: %15.6f K",       SogFluidInfo.dvNBPTemperature);
    info("Acentric Factor ---------------------------: %15.6f ",        SogFluidInfo.dvAcentric_Fact);
    info("Dipole Moment -----------------------------: %15.6f debye",   SogFluidInfo.dvDipole_Mom);
    info("Temperature Minimum -----------------------: %15.6f K",       SogFluidInfo.dvTemperature_Min);
    info("Temperature Maximum -----------------------: %15.6f K",       SogFluidInfo.dvTemperature_Max);
    info("Pressure Maximum --------------------------: %15.6f kPa",     EOS_NIST_Get_SI_To_NIST_Pressure(SogFluidInfo.dvPressure_Max));
    info("Density Maximum ---------------------------: %15.6f kg/m^3",  SogFluidInfo.dvDensity_Max);
}

//------------------------------------------------------------------------------
//! Set the Reference Property for Fluid : Generic (Reference Density is Computed)
//  Pressure (Pa), Temperature (K) and Length (m)
//------------------------------------------------------------------------------
void EOS_NIST_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength) {
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
    dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
    
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
    dvgDensity_Ref           = EOS_NIST_Get_NIST_To_SI_Density(dvDensity);
    dvgSpeedSound_Ref        = dvSpeedSound;
//    dvgSpeedSound_Ref = sqrt(dvSpeedSound*sqrt(dvPressure/dvgDensity_Ref));
    dvgPressure_Ref          = dvgDensity_Ref*dvgSpeedSound_Ref*dvgSpeedSound_Ref; // Pressure in Pa  
/*
    dvgPressure_Ref          = dvPressure;
    dvgTemperature_Ref       = dvTemperature;
    dvgDensity_Ref           = dvgPressure_Ref/dvgTemperature_Ref;
    dvgSpeedSound_Ref        = sqrt(dvgPressure_Ref/dvgDensity_Ref);
*/  
    dvgLength_Ref            = dvLength;
    // Derived Reference Properties
    dvgRatioSpecificHeat_Ref = dvCp/dvCv;
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
//! Compute the NIST to SI compatible Density
//  Density in mol/L to kg/m^3 
//------------------------------------------------------------------------------
double EOS_NIST_Get_NIST_To_SI_Density(double dvDensity) {
    return dvDensity*SogFluidInfo.dvMolecularWeight;
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Pressure
//  Pressure in kPa to Pa
//------------------------------------------------------------------------------
double EOS_NIST_Get_NIST_To_SI_Pressure(double dvPressure) {
    return 1000.0*dvPressure;
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Enthalpy
//  Pressure in J/mol to J/kg
//------------------------------------------------------------------------------
double EOS_NIST_Get_NIST_To_SI_Enthalpy(double dvEnthalpy) {
    return dvEnthalpy*(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Entropy
//  Pressure in J/(mol-K) to J/(kg-K)
//------------------------------------------------------------------------------
double EOS_NIST_Get_NIST_To_SI_Entropy(double dvEntropy) {
    return dvEntropy*(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Entropy
//  Pressure in J/mol to J/kg
//------------------------------------------------------------------------------
double EOS_NIST_Get_NIST_To_SI_InternalEnergy(double dvInternalEnergy) {
    return dvInternalEnergy*(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible HeatCapacityCv
//  Pressure in J/(mol-K) to J/(kg-K)
//------------------------------------------------------------------------------
double EOS_NIST_Get_NIST_To_SI_HeatCapacityCv(double dvHeatCapacityCv) {
    return dvHeatCapacityCv*(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible HeatCapacityCp
//  Pressure in J/(mol-K) to J/(kg-K)
//------------------------------------------------------------------------------
double EOS_NIST_Get_NIST_To_SI_HeatCapacityCp(double dvHeatCapacityCp) {
    return dvHeatCapacityCp*(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Properties
//  
//------------------------------------------------------------------------------
void EOS_NIST_Get_NIST_To_SI_Property(double *dpProperty) {
    double dvTmp1;
    dvTmp1 = 1000.0/SogFluidInfo.dvMolecularWeight;
    
    dpProperty[ 0] = dpProperty[ 0]*SogFluidInfo.dvMolecularWeight; // Density
    dpProperty[ 1] = dpProperty[ 1]*SogFluidInfo.dvMolecularWeight; // Density Liquid
    dpProperty[ 2] = dpProperty[ 2]*SogFluidInfo.dvMolecularWeight; // Density Vapor
    dpProperty[ 3] = dpProperty[ 3]*1000.0;                         // Pressure
    dpProperty[ 4] = dpProperty[ 4];                                // Temperature
    dpProperty[ 5] = dpProperty[ 5];                                // Speed of Sound
    dpProperty[ 6] = dpProperty[ 6]*dvTmp1;                         // Entropy
    dpProperty[ 7] = dpProperty[ 7]*dvTmp1;                         // Enthalpy
    dpProperty[ 8] = dpProperty[ 8]*dvTmp1;                         // Internal Energy
    dpProperty[ 9] = dpProperty[ 9]*dvTmp1;                         // Heat Capacity Cv
    dpProperty[10] = dpProperty[10]*dvTmp1;                         // Heat Capacity Cp
}

//------------------------------------------------------------------------------
//! Compute the NIST to SI compatible Properties
//  
//------------------------------------------------------------------------------
void EOS_NIST_Get_NIST_To_SI_Extended_Property(double *dpProperty) {
    double dvTmp1, dvTmp2, dvTmp3;
    dvTmp1 = 1000.0/SogFluidInfo.dvMolecularWeight;
    dvTmp2 = SogFluidInfo.dvMolecularWeight/1000.0;
    dvTmp3 = 1000.0/(SogFluidInfo.dvMolecularWeight*SogFluidInfo.dvMolecularWeight);
    
    dpProperty[ 0] = dpProperty[ 0]*SogFluidInfo.dvMolecularWeight; // Density
    dpProperty[ 1] = dpProperty[ 1]*SogFluidInfo.dvMolecularWeight; // Density Liquid
    dpProperty[ 2] = dpProperty[ 2]*SogFluidInfo.dvMolecularWeight; // Density Vapor
    dpProperty[ 3] = dpProperty[ 3]*1000.0;                         // Pressure
    dpProperty[ 4] = dpProperty[ 4];                                // Temperature
    dpProperty[ 5] = dpProperty[ 5];                                // Speed of Sound
    dpProperty[ 6] = dpProperty[ 6]*dvTmp1;                         // Entropy
    dpProperty[ 7] = dpProperty[ 7]*dvTmp1;                         // Enthalpy
    dpProperty[ 8] = dpProperty[ 8]*dvTmp1;                         // Internal Energy
    dpProperty[ 9] = dpProperty[ 9]*dvTmp1;                         // Heat Capacity Cv
    dpProperty[10] = dpProperty[10]*dvTmp1;                         // Heat Capacity Cp
    dpProperty[11] = dpProperty[11]*dvTmp1;                         // DPDRho
    dpProperty[12] = dpProperty[12]*1000.0;                         // DPDT
    dpProperty[13] = dpProperty[13]*SogFluidInfo.dvMolecularWeight; // DRhoDT
    dpProperty[14] = dpProperty[14]*dvTmp2;                         // DRhoDP;
    dpProperty[15] = dpProperty[15]*dvTmp3;                         // D2PDRho2
    dpProperty[16] = dpProperty[16]*1000.0;                         // D2PDT2
    dpProperty[17] = dpProperty[17]*dvTmp1;                         // D2PDTRho
    dpProperty[18] = dpProperty[18]*dvTmp1;                         // DHdT_Rho
    dpProperty[19] = dpProperty[19]*dvTmp1;                         // DHdT_P
    dpProperty[20] = dpProperty[20]*dvTmp3;                         // DHDRho_T
    dpProperty[21] = dpProperty[21]*dvTmp3;                         // DHDRho_P
    dpProperty[22] = dpProperty[22]/SogFluidInfo.dvMolecularWeight; // DHDP_T
    dpProperty[23] = dpProperty[23]/SogFluidInfo.dvMolecularWeight; // DHDP_Rho
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Density
//  Density in kg/m^3 to mol/L
//------------------------------------------------------------------------------
double EOS_NIST_Get_SI_To_NIST_Density(double dvDensity) {
    return dvDensity/SogFluidInfo.dvMolecularWeight;
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Pressure
//  Pressure in Pa to kPa
//------------------------------------------------------------------------------
double EOS_NIST_Get_SI_To_NIST_Pressure(double dvPressure) {
    return 0.001*dvPressure;
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Enthalpy
//  Pressure in J/kg to J/mol
//------------------------------------------------------------------------------
double EOS_NIST_Get_SI_To_NIST_Enthalpy(double dvEnthalpy) {
    return dvEnthalpy/(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Entropy
//  Pressure in J/(kg-K) to J/(mol-K)
//------------------------------------------------------------------------------
double EOS_NIST_Get_SI_To_NIST_Entropy(double dvEntropy) {
    return dvEntropy/(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Entropy
//  Pressure in J/kg to J/mol
//------------------------------------------------------------------------------
double EOS_NIST_Get_SI_To_NIST_InternalEnergy(double dvInternalEnergy) {
    return dvInternalEnergy/(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible HeatCapacityCv
//  Pressure in J/(kg-K) to J/(mol-K)
//------------------------------------------------------------------------------
double EOS_NIST_Get_SI_To_NIST_HeatCapacityCv(double dvHeatCapacityCv) {
    return dvHeatCapacityCv/(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible HeatCapacityCp
//  Pressure in J/(kg-K) to J/(mol-K)
//------------------------------------------------------------------------------
double EOS_NIST_Get_SI_To_NIST_HeatCapacityCp(double dvHeatCapacityCp) {
    return dvHeatCapacityCp/(1000.0/SogFluidInfo.dvMolecularWeight);
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Properties
//  
//------------------------------------------------------------------------------
void EOS_NIST_Get_SI_To_NIST_Property(double *dpProperty) {
    double dvTmp1;
    dvTmp1 = 1000.0/SogFluidInfo.dvMolecularWeight;
    
    dpProperty[ 0] = dpProperty[ 0]/SogFluidInfo.dvMolecularWeight; // Density
    dpProperty[ 1] = dpProperty[ 1]/SogFluidInfo.dvMolecularWeight; // Density Liquid
    dpProperty[ 2] = dpProperty[ 2]/SogFluidInfo.dvMolecularWeight; // Density Vapor
    dpProperty[ 3] = dpProperty[ 3]/1000.0;                         // Pressure
    dpProperty[ 4] = dpProperty[ 4];                                // Temperature
    dpProperty[ 5] = dpProperty[ 5];                                // Speed of Sound
    dpProperty[ 6] = dpProperty[ 6]/dvTmp1;                         // Entropy
    dpProperty[ 7] = dpProperty[ 7]/dvTmp1;                         // Enthalpy
    dpProperty[ 8] = dpProperty[ 8]/dvTmp1;                         // Internal Energy
    dpProperty[ 9] = dpProperty[ 9]/dvTmp1;                         // Heat Capacity Cv
    dpProperty[10] = dpProperty[10]/dvTmp1;                         // Heat Capacity Cp
}

//------------------------------------------------------------------------------
//! Compute the SI to NIST compatible Properties
//  
//------------------------------------------------------------------------------
void EOS_NIST_Get_SI_To_NIST_Extended_Property(double *dpProperty) {
    double dvTmp1, dvTmp2, dvTmp3;
    dvTmp1 = 1000.0/SogFluidInfo.dvMolecularWeight;
    dvTmp2 = SogFluidInfo.dvMolecularWeight/1000.0;
    dvTmp3 = 1000.0/(SogFluidInfo.dvMolecularWeight*SogFluidInfo.dvMolecularWeight);
    
    dpProperty[ 0] = dpProperty[ 0]/SogFluidInfo.dvMolecularWeight; // Density
    dpProperty[ 1] = dpProperty[ 1]/SogFluidInfo.dvMolecularWeight; // Density Liquid
    dpProperty[ 2] = dpProperty[ 2]/SogFluidInfo.dvMolecularWeight; // Density Vapor
    dpProperty[ 3] = dpProperty[ 3]/1000.0;                         // Pressure
    dpProperty[ 4] = dpProperty[ 4];                                // Temperature
    dpProperty[ 5] = dpProperty[ 5];                                // Speed of Sound
    dpProperty[ 6] = dpProperty[ 6]/dvTmp1;                         // Entropy
    dpProperty[ 7] = dpProperty[ 7]/dvTmp1;                         // Enthalpy
    dpProperty[ 8] = dpProperty[ 8]/dvTmp1;                         // Internal Energy
    dpProperty[ 9] = dpProperty[ 9]/dvTmp1;                         // Heat Capacity Cv
    dpProperty[10] = dpProperty[10]/dvTmp1;                         // Heat Capacity Cp
    dpProperty[11] = dpProperty[11]/dvTmp1;                         // DPDRho
    dpProperty[12] = dpProperty[12]/1000.0;                         // DPDT
    dpProperty[13] = dpProperty[13]/SogFluidInfo.dvMolecularWeight; // DRhoDT
    dpProperty[14] = dpProperty[14]/dvTmp2;                         // DRhoDP;
    dpProperty[15] = dpProperty[15]/dvTmp3;                         // D2PDRho2
    dpProperty[16] = dpProperty[16]/1000.0;                         // D2PDT2
    dpProperty[17] = dpProperty[17]/dvTmp1;                         // D2PDTRho
    dpProperty[18] = dpProperty[18]/dvTmp1;                         // DHdT_Rho
    dpProperty[19] = dpProperty[19]/dvTmp1;                         // DHdT_P
    dpProperty[20] = dpProperty[20]/dvTmp3;                         // DHDRho_T
    dpProperty[21] = dpProperty[21]/dvTmp3;                         // DHDRho_P
    dpProperty[22] = dpProperty[22]*SogFluidInfo.dvMolecularWeight; // DHDP_T
    dpProperty[23] = dpProperty[23]*SogFluidInfo.dvMolecularWeight; // DHDP_Rho
}

//------------------------------------------------------------------------------
//! Compute the EOS NIST Properties based on Variable Type
//! Input and Output property are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_NIST_Get_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    int i;
    double dvRho, dvRhoL, dvRhoV, dvPressure, dvTemperature;
    double dvVelocityU = 0.0, dvVelocityV = 0.0, dvVelocityW = 0.0, dvQ2 = 0.0, dvSpeedSound, dvMach;
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
            error("EOS_NIST_Get_Properties:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Properties:2: Conservative Variables Not Implemented");
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
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Properties:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Properties:4: EOS_VARIABLE_PUS Variables Not Implemented");
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
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Properties:5: Undefined Variable Type - %d", ivVariableType);
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
    EOS_NIST_Get_NIST_To_SI_Property(dpPropertyOut);
    
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
            error("EOS_NIST_Get_Properties:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
}

//------------------------------------------------------------------------------
//! Compute the EOS NIST Extended Properties based on Variable Type
//! Input and Output property are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_NIST_Get_Extended_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    int i;
    double dvRho, dvRhoL, dvRhoV, dvPressure, dvTemperature;
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST;
    double dvVelocityU = 0.0, dvVelocityV = 0.0, dvVelocityW = 0.0, dvQ2 = 0.0, dvSpeedSound, dvMach;
    double dvEntropy, dvEnthalpy, dvInternalEnergy, dvTotalEnergy, dvTotalEnthalpy;
    double dvCv, dvCp, dvDPDRho, dvDPDT, dvDRhoDT, dvDRhoDP, dvD2PDRho2, dvD2PDT2, dvD2PDTRho;
    double dvDHDT_Rho, dvDHDT_P, dvDHDRho_T, dvDHDRho_P, dvDHDP_T, dvDHDP_Rho;
    double dvDEDT_Rho, dvDEDT_P, dvDEDRho_T, dvDEDRho_P, dvDEDP_T, dvDEDP_Rho;
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
            error("EOS_NIST_Get_Extended_Properties:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Extended_Properties:2: Conservative Variables Not Implemented");
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
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Extended Quantities
            PDEPDFLSH2(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, SogNISTHelper.daXLiq, SogNISTHelper.daXVap,
                    dpPropertyOut, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Extended_Properties:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Extended_Properties:4: EOS_VARIABLE_PUS Variables Not Implemented");
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
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH2(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, SogNISTHelper.daXLiq, SogNISTHelper.daXVap, 
                    dpPropertyOut, &SogNISTHelper.ivError, SogNISTHelper.csError, NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Extended_Properties:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Properties in NIST units
    switch (ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Extended_Properties:6: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvTemperature = dpPropertyOut[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Extended_Properties:7: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Extended_Properties:8: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvPressureNIST = dpPropertyOut[0];
            break;
        default:
            error("EOS_NIST_Get_Extended_Properties:9: Undefined Variable Type - %d", ivVariableType);
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
    dvDHDT_Rho       = dpPropertyOut[23];
    dvDHDT_P         = dpPropertyOut[24];
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
    dpPropertyOut[18] = dvDHDT_Rho;
    dpPropertyOut[19] = dvDHDT_P;
    dpPropertyOut[20] = dvDHDRho_T;
    dpPropertyOut[21] = dvDHDRho_P;
    dpPropertyOut[22] = dvDHDP_T;
    dpPropertyOut[23] = dvDHDP_Rho;
    EOS_NIST_Get_NIST_To_SI_Extended_Property(dpPropertyOut);
    
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
    dvDHDT_Rho       = dpPropertyOut[18];
    dvDHDT_P         = dpPropertyOut[19];
    dvDHDRho_T       = dpPropertyOut[20];
    dvDHDRho_P       = dpPropertyOut[21];
    dvDHDP_T         = dpPropertyOut[22];
    dvDHDP_Rho       = dpPropertyOut[23];
    
    // Compute the Remaining Properties in SI units
    dvMach          = sqrt(dvQ2)/dvSpeedSound;
    dvTotalEnthalpy = dvEnthalpy + 0.5*dvQ2;
    dvTotalEnergy   = dvInternalEnergy + 0.5*dvQ2;
    dvDEDT_Rho      = dvDHDT_Rho - dvDPDT/dvRho;
    dvDEDT_P        = dvDHDT_P + dvPressure*dvDRhoDT/(dvRho*dvRho);
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
    dpPropertyOut[25] = dvDHDT_Rho;
    dpPropertyOut[26] = dvDHDT_P;
    dpPropertyOut[27] = dvDHDRho_T;
    dpPropertyOut[28] = dvDHDRho_P;
    dpPropertyOut[29] = dvDHDP_T;
    dpPropertyOut[30] = dvDHDP_Rho;
    dpPropertyOut[31] = dvDEDT_Rho;
    dpPropertyOut[32] = dvDEDT_P;
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
            error("EOS_NIST_Get_Extended_Properties:10: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Density(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho;
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
            error("EOS_NIST_Get_Density:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvRho = 0.0;
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Density:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Density:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Density:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            break;
        default:
            error("EOS_NIST_Get_Density:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvRho = EOS_Internal_NonDimensionalize_Density(dvRho);
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
            dvRho = EOS_Internal_NonDimensionalize_Density(dvRho);
            break;
        default:
            error("EOS_NIST_Get_Density:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvRho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Density_Liquid(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvRhoL, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Density_Liquid:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Density_Liquid:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Density_Liquid:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Density_Liquid:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Density_Liquid:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Convert Properties From NIST to SI units
    dvRhoL = EOS_NIST_Get_NIST_To_SI_Density(dvRhoLNIST);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvRhoL = EOS_Internal_NonDimensionalize_Density(dvRhoL);
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
            dvRhoL = EOS_Internal_NonDimensionalize_Density(dvRhoL);
            break;
        default:
            error("EOS_NIST_Get_Density_Liquid:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvRhoL;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Density_Vapor(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvRhoV, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Density_Vapor:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Density_Vapor:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Density_Vapor:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_DensityVapor:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Density_Vapor:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Convert Properties From NIST to SI units
    dvRhoV = EOS_NIST_Get_NIST_To_SI_Density(dvRhoVNIST);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvRhoV = EOS_Internal_NonDimensionalize_Density(dvRhoV);
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
            dvRhoV = EOS_Internal_NonDimensionalize_Density(dvRhoV);
            break;
        default:
            error("EOS_NIST_Get_Density_Vapor:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvRhoV;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void EOS_NIST_Get_Density_All(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpDensityOut) {
    int i;
    double dvRho, dvRhoL, dvRhoV, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Density_All:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Density_All:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Density_All:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Density_All:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Density_All:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Convert Properties From NIST to SI units
    dvRho  = EOS_NIST_Get_NIST_To_SI_Density(dvRhoNIST);
    dvRhoL = EOS_NIST_Get_NIST_To_SI_Density(dvRhoLNIST);
    dvRhoV = EOS_NIST_Get_NIST_To_SI_Density(dvRhoVNIST);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvRho  = EOS_Internal_NonDimensionalize_Density(dvRho);
            dvRhoL = EOS_Internal_NonDimensionalize_Density(dvRhoL);
            dvRhoV = EOS_Internal_NonDimensionalize_Density(dvRhoV);
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
            dvRho  = EOS_Internal_NonDimensionalize_Density(dvRho);
            dvRhoL = EOS_Internal_NonDimensionalize_Density(dvRhoL);
            dvRhoV = EOS_Internal_NonDimensionalize_Density(dvRhoV);
            break;
        default:
            error("EOS_NIST_Get_Density_All:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    dpDensityOut[0] = dvRho;
    dpDensityOut[1] = dvRhoL;
    dpDensityOut[2] = dvRhoV;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Pressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure = 0.0, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Pressure:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Pressure:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvPressure  = daVariableDimensional[4];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Pressure:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Pressure:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            // Compute Property in SI Units
            dvPressure = EOS_NIST_Get_NIST_To_SI_Pressure(dvPressureNIST);
            break;
        default:
            error("EOS_NIST_Get_Pressure:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvPressure = EOS_Internal_NonDimensionalize_Pressure(dvPressure);
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
            dvPressure = EOS_Internal_NonDimensionalize_Pressure(dvPressure);
            break;
        default:
            error("EOS_NIST_Get_Pressure:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvPressure;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Temperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Temperature:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Temperature:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Temperature:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Temperature:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvTemperature = daVariableDimensional[4];
            break;
        default:
            error("EOS_NIST_Get_Temperature:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
      
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvTemperature = EOS_Internal_NonDimensionalize_Temperature(dvTemperature);
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
            dvTemperature = EOS_Internal_NonDimensionalize_Temperature(dvTemperature);
            break;
        default:
            error("EOS_NIST_Get_Temperature:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvTemperature;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_SpeedSound(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_SpeedSound:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_SpeedSound:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_SpeedSound:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_SpeedSound:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_SpeedSound:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvSpeedSound = EOS_Internal_NonDimensionalize_SpeedSound(dvSpeedSound);
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
            dvSpeedSound = EOS_Internal_NonDimensionalize_SpeedSound(dvSpeedSound);
            break;
        default:
            error("EOS_NIST_Get_SpeedSound:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvSpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Mach(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU = 0.0, dvVelocityV = 0.0, dvVelocityW = 0.0, dvQ2 = 0.0, dvSpeedSound, dvMach;
    double dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Mach:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Mach:2: Conservative Variables Not Implemented");
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
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Mach:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Mach:4: EOS_VARIABLE_PUS Variables Not Implemented");
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
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Mach:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvMach          = sqrt(dvQ2)/dvSpeedSound;
    
    //-END----------------Dimensional Computations------------------------------
    return dvMach;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Entropy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Entropy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Entropy:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Entropy:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Entropy:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Entropy:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Get the Properties in SI units
    dvEntropy = EOS_NIST_Get_NIST_To_SI_Entropy(dvEntropy);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvEntropy = EOS_Internal_NonDimensionalize_Entropy(dvEntropy);
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
            dvEntropy = EOS_Internal_NonDimensionalize_Entropy(dvEntropy);
            break;
        default:
            error("EOS_NIST_Get_Entropy:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvEntropy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Enthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_Enthalpy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Enthalpy:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Enthalpy:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Enthalpy:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Enthalpy:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Get the Properties in SI units
    dvEnthalpy = EOS_NIST_Get_NIST_To_SI_Enthalpy(dvEnthalpy);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvEnthalpy = EOS_Internal_NonDimensionalize_Enthalpy(dvEnthalpy);
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
            dvEnthalpy = EOS_Internal_NonDimensionalize_Enthalpy(dvEnthalpy);
            break;
        default:
            error("EOS_NIST_Get_Enthalpy:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvEnthalpy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_InternalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_InternalEnergy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_InternalEnergy:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_InternalEnergy:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_InternalEnergy:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_InternalEnergy:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Get the Properties in SI units
    dvInternalEnergy = EOS_NIST_Get_NIST_To_SI_InternalEnergy(dvInternalEnergy);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvInternalEnergy = EOS_Internal_NonDimensionalize_InternalEnergy(dvInternalEnergy);
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
            dvInternalEnergy = EOS_Internal_NonDimensionalize_InternalEnergy(dvInternalEnergy);
            break;
        default:
            error("EOS_NIST_Get_InternalEnergy:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvInternalEnergy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_TotalEnthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU = 0.0, dvVelocityV = 0.0, dvVelocityW = 0.0, dvQ2 = 0.0, dvSpeedSound;
    double dvEntropy, dvEnthalpy, dvInternalEnergy, dvTotalEnthalpy;
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
            error("EOS_NIST_Get_TotalEnthalpy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_TotalEnthalpy:2: Conservative Variables Not Implemented");
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
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_TotalEnthalpy:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_TotalEnthalpy:4: EOS_VARIABLE_PUS Variables Not Implemented");
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
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_TotalEnthalpy:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Get the Properties in SI units
    dvEnthalpy = EOS_NIST_Get_NIST_To_SI_Enthalpy(dvEnthalpy);
    
    // Compute the Remaining Properties in SI units
    dvTotalEnthalpy = dvEnthalpy + 0.5*dvQ2;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvTotalEnthalpy = EOS_Internal_NonDimensionalize_TotalEnthalpy(dvTotalEnthalpy);
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
            dvTotalEnthalpy = EOS_Internal_NonDimensionalize_TotalEnthalpy(dvTotalEnthalpy);
            break;
        default:
            error("EOS_NIST_Get_TotalEnthalpy:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvTotalEnthalpy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_TotalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU = 0.0, dvVelocityV = 0.0, dvVelocityW = 0.0, dvQ2 = 0.0, dvSpeedSound;
    double dvEntropy, dvEnthalpy, dvInternalEnergy, dvTotalEnergy;
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
            error("EOS_NIST_Get_TotalEnergy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_TotalEnergy:2: Conservative Variables Not Implemented");
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
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_TotalEnergy:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_TotalEnergy:4: EOS_VARIABLE_PUS Variables Not Implemented");
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
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_TotalEnergy:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Get the Properties in SI units
    dvInternalEnergy = EOS_NIST_Get_NIST_To_SI_InternalEnergy(dvInternalEnergy);
    
    // Compute the Remaining Properties in SI units
    dvTotalEnergy   = dvInternalEnergy + 0.5*dvQ2;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvTotalEnergy = EOS_Internal_NonDimensionalize_TotalEnergy(dvTotalEnergy);
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
            dvTotalEnergy = EOS_Internal_NonDimensionalize_TotalEnergy(dvTotalEnergy);
            break;
        default:
            error("EOS_NIST_Get_TotalEnergy:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvTotalEnergy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_HeatCapacityCv(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_HeatCapacityCv:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_HeatCapacityCv:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_HeatCapacityCv:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_HeatCapacityCv:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_HeatCapacityCv:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Get the Properties in SI units;
    dvCv = EOS_NIST_Get_NIST_To_SI_HeatCapacityCv(dvCv);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvCv = EOS_Internal_NonDimensionalize_HeatCapacityCv(dvCv);
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
            dvCv = EOS_Internal_NonDimensionalize_HeatCapacityCv(dvCv);
            break;
        default:
            error("EOS_NIST_Get_HeatCapacityCv:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvCv;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_HeatCapacityCp(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound, dvEntropy, dvEnthalpy, dvInternalEnergy;
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
            error("EOS_NIST_Get_HeatCapacityCp:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_HeatCapacityCp:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_HeatCapacityCp:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_HeatCapacityCp:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_HeatCapacityCp:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Get the Properties in SI units
    dvCp = EOS_NIST_Get_NIST_To_SI_HeatCapacityCv(dvCp);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvCp = EOS_Internal_NonDimensionalize_HeatCapacityCp(dvCp);
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
            dvCp = EOS_Internal_NonDimensionalize_HeatCapacityCp(dvCp);
            break;
        default:
            error("EOS_NIST_Get_HeatCapacityCp:6: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvCp;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_Quality(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvSpeedSound;
    double dvEntropy, dvEnthalpy, dvInternalEnergy;
    double dvCv, dvCp;
    double daVariableDimensional[EOS_NEQUATIONS];
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST, dvQualityNIST;
    
    dvQualityNIST = 0.0;
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
            error("EOS_NIST_Get_Quality:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    switch(ivVariableType) {
        // Conservative Variables
        case EOS_VARIABLE_CON:
            error("EOS_NIST_Get_Quality:2: Conservative Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho       = daVariableDimensional[0];
            dvPressure  = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST      = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
            // Compute Thermodynamic Quantities
            PDEPDFLSH(&dvPressureNIST, &dvRhoNIST, SogNISTHelper.daX, &dvTemperature, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            error("EOS_NIST_Get_Quality:3: EOS_VARIABLE_PUT Variables Not Implemented");
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUS:
            error("EOS_NIST_Get_Quality:4: EOS_VARIABLE_PUS Variables Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho         = daVariableDimensional[0];
            dvTemperature = daVariableDimensional[4];
            // Compute NIST compatible values
            dvRhoNIST     = EOS_NIST_Get_SI_To_NIST_Density(dvRho);
            // Compute Thermodynamic Quantities
            TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                    SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                    &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                    NISTTHERMO_GEN_STR_LEN);
            break;
        default:
            error("EOS_NIST_Get_Quality:5: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    return dvQualityNIST;
    //-END----------------Dimensional Computations------------------------------
}

// Compute First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DPressureDDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DPressureDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DDensityDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

// Compute Second Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_D2PressureDDensity2(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_D2PressureDTemperature2(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_D2PressureDTemperatureDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

// Enthalpy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DEnthalpyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DEnthalpyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DEnthalpyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DEnthalpyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DEnthalpyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DEnthalpyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

// Internal Energy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DInternalEnergyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DInternalEnergyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DInternalEnergyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DInternalEnergyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DInternalEnergyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DInternalEnergyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double result = 0.0;
    
    return result;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void EOS_NIST_Get_PT_Density(int ivDimIOType, double dvPressure, double dvTemperature, double *dpDensityOut) {
    double dvRho, dvRhoL, dvRhoV;
    double dvSpeedSound;
    double dvEntropy, dvEnthalpy, dvInternalEnergy;
    double dvCv, dvCp;
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST, dvQualityNIST;
    
    // Dimensionalize the Input Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvPressure    = EOS_Internal_Dimensionalize_Pressure(dvPressure);
            dvTemperature = EOS_Internal_Dimensionalize_Temperature(dvTemperature);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            dvPressure    = EOS_Internal_Dimensionalize_Pressure(dvPressure);
            dvTemperature = EOS_Internal_Dimensionalize_Temperature(dvTemperature);
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        default:
            error("EOS_NIST_Get_PT_Density:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    // Compute NIST compatible values
    dvPressureNIST = EOS_NIST_Get_SI_To_NIST_Pressure(dvPressure);
    // Compute Thermodynamic Quantities
    TPFLSH(&dvTemperature, &dvPressureNIST, SogNISTHelper.daX, &dvRhoNIST, &dvRhoLNIST, &dvRhoVNIST,
            SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
            &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
            NISTTHERMO_GEN_STR_LEN);
    
    // Convert Properties From NIST to SI units
    dvRho  = EOS_NIST_Get_NIST_To_SI_Density(dvRhoNIST);
    dvRhoL = EOS_NIST_Get_NIST_To_SI_Density(dvRhoLNIST);
    dvRhoV = EOS_NIST_Get_NIST_To_SI_Density(dvRhoVNIST);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvRho  = EOS_Internal_NonDimensionalize_Density(dvRho);
            dvRhoL = EOS_Internal_NonDimensionalize_Density(dvRhoL);
            dvRhoV = EOS_Internal_NonDimensionalize_Density(dvRhoV);
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
            dvRho  = EOS_Internal_NonDimensionalize_Density(dvRho);
            dvRhoL = EOS_Internal_NonDimensionalize_Density(dvRhoL);
            dvRhoV = EOS_Internal_NonDimensionalize_Density(dvRhoV);
            break;
        default:
            error("EOS_NIST_Get_PT_Density:2: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    
    // Set the Output
    dpDensityOut[0] = dvRho;
    dpDensityOut[1] = dvRhoL;
    dpDensityOut[2] = dvRhoV;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DH_SpeedSound(int ivDimIOType, double dvDensity, double dvEnthalpy) {
    double dvSpeedSound, dvTemperature;
    double dvEntropy, dvInternalEnergy;
    double dvCv, dvCp;
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST, dvQualityNIST;
    
    // Dimensionalize the Input Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDensity  = EOS_Internal_Dimensionalize_Density(dvDensity);
            dvEnthalpy = EOS_Internal_Dimensionalize_Enthalpy(dvEnthalpy);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            dvDensity  = EOS_Internal_Dimensionalize_Density(dvDensity);
            dvEnthalpy = EOS_Internal_Dimensionalize_Enthalpy(dvEnthalpy);
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        default:
            error("EOS_NIST_Get_DH_SpeedSound:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    // Compute NIST compatible values
    dvRhoNIST  = EOS_NIST_Get_SI_To_NIST_Density(dvDensity);
    dvEnthalpy = EOS_NIST_Get_SI_To_NIST_Enthalpy(dvEnthalpy);
    dvSpeedSound = 0.0;
    // Compute Thermodynamic Quantities
    DHFLSH(&dvRhoNIST, &dvEnthalpy, SogNISTHelper.daX, &dvTemperature, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
            SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy,
            &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
            NISTTHERMO_GEN_STR_LEN);
    
    // Validate Speed of Sound
    if (dvSpeedSound <= 0.0) {
        // Multiphase Region
        // Compute Thermodynamic Quantities
        TDETDFLSH(&dvTemperature, &dvRhoNIST, SogNISTHelper.daX, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
                SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy, &dvEnthalpy,
                &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
                NISTTHERMO_GEN_STR_LEN);
        
        if (dvSpeedSound <= 0.0)
            error("EOS_NIST_Get_DH_SpeedSound:2: Unable to compute speed of sound");
    }
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvSpeedSound = EOS_Internal_NonDimensionalize_SpeedSound(dvSpeedSound);
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
            dvSpeedSound = EOS_Internal_NonDimensionalize_SpeedSound(dvSpeedSound);
            break;
        default:
            error("EOS_NIST_Get_DH_SpeedSound:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvSpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DH_Pressure(int ivDimIOType, double dvDensity, double dvEnthalpy) {
    double dvPressure, dvSpeedSound, dvTemperature;
    double dvEntropy, dvInternalEnergy;
    double dvCv, dvCp;
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST, dvQualityNIST;
    
    // Dimensionalize the Input Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDensity  = EOS_Internal_Dimensionalize_Density(dvDensity);
            dvEnthalpy = EOS_Internal_Dimensionalize_Enthalpy(dvEnthalpy);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            dvDensity  = EOS_Internal_Dimensionalize_Density(dvDensity);
            dvEnthalpy = EOS_Internal_Dimensionalize_Enthalpy(dvEnthalpy);
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        default:
            error("EOS_NIST_Get_DH_Pressure:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    // Compute NIST compatible values
    dvRhoNIST  = EOS_NIST_Get_SI_To_NIST_Density(dvDensity);
    dvEnthalpy = EOS_NIST_Get_SI_To_NIST_Enthalpy(dvEnthalpy);
    // Compute Thermodynamic Quantities
    DHFLSH(&dvRhoNIST, &dvEnthalpy, SogNISTHelper.daX, &dvTemperature, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
            SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy,
            &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
            NISTTHERMO_GEN_STR_LEN);
    
    // Compute Property in SI Units
    dvPressure = EOS_NIST_Get_NIST_To_SI_Pressure(dvPressureNIST);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvPressure = EOS_Internal_NonDimensionalize_Pressure(dvPressure);
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
            dvPressure = EOS_Internal_NonDimensionalize_Pressure(dvPressure);
            break;
        default:
            error("EOS_NIST_Get_DH_Pressure:2: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvPressure;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_NIST_Get_DH_Temperature(int ivDimIOType, double dvDensity, double dvEnthalpy) {
    double dvSpeedSound, dvTemperature;
    double dvEntropy, dvInternalEnergy;
    double dvCv, dvCp;
    double dvRhoNIST, dvRhoLNIST, dvRhoVNIST, dvPressureNIST, dvQualityNIST;
    
    // Dimensionalize the Input Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDensity  = EOS_Internal_Dimensionalize_Density(dvDensity);
            dvEnthalpy = EOS_Internal_Dimensionalize_Enthalpy(dvEnthalpy);
            break;
        // Type Input: Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_D_D:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        // Type Input: Non-Dimensional  Output: Dimensional
        case EOS_DIMENSIONAL_IO_ND_D:
            dvDensity  = EOS_Internal_Dimensionalize_Density(dvDensity);
            dvEnthalpy = EOS_Internal_Dimensionalize_Enthalpy(dvEnthalpy);
            break;
        // Type Input: Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_D_ND:
            // Do Nothing Input is Dimensional (SI Units)
            break;
        default:
            error("EOS_NIST_Get_DH_Temperature:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    // Compute Properties Based on Input Variable Types
    // Compute NIST compatible values
    dvRhoNIST  = EOS_NIST_Get_SI_To_NIST_Density(dvDensity);
    dvEnthalpy = EOS_NIST_Get_SI_To_NIST_Enthalpy(dvEnthalpy);
    // Compute Thermodynamic Quantities
    DHFLSH(&dvRhoNIST, &dvEnthalpy, SogNISTHelper.daX, &dvTemperature, &dvPressureNIST, &dvRhoLNIST, &dvRhoVNIST,
            SogNISTHelper.daXLiq, SogNISTHelper.daXVap, &dvQualityNIST, &dvInternalEnergy,
            &dvEntropy, &dvCv, &dvCp, &dvSpeedSound, &SogNISTHelper.ivError, SogNISTHelper.csError,
            NISTTHERMO_GEN_STR_LEN);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvTemperature = EOS_Internal_NonDimensionalize_Temperature(dvTemperature);
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
            dvTemperature = EOS_Internal_NonDimensionalize_Temperature(dvTemperature);
            break;
        default:
            error("EOS_NIST_Get_DH_Temperature:2: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvTemperature;
}

