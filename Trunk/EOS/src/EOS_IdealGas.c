/*******************************************************************************
 * File:        EOS_IdealGas.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "EOS.h"
#include "EOS_IdealGas.h"
#include "EOS_Internal.h"
#include "Trim_Utils.h"

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
void EOS_IdealGas_Init() {
    
}

//------------------------------------------------------------------------------
//! Finalize EOS Ideal Gas
//------------------------------------------------------------------------------
void EOS_IdealGas_Finalize() {
    
}

//------------------------------------------------------------------------------
//! Setup only for ideal gas
//------------------------------------------------------------------------------
void EOS_IdealGas_Set() {
    
    // Fluid Model provided by Ideal Gas
    // Set the Fluid Name
    strncpy(SogFluidInfo.csFluidName, "Ideal Gas" , 9);
    
    // Get the Fluid Information
    EOS_IdealGas_Set_Fluid_Information();
    
    // Print the Fluid Information
    EOS_IdealGas_Print_Fluid_Information();
}

//------------------------------------------------------------------------------
//! Initialize Fluid Information from NIST Thermodynamic Library
//  All Thermodynamic variables are in SI Units
//------------------------------------------------------------------------------
void EOS_IdealGas_Set_Fluid_Information() {
    // Set the Information about the ideal gas fluid in SI Units
    SogFluidInfo.dvMolecularWeight = 28.966;
    SogFluidInfo.dvGasConstant     = 8314.47215/SogFluidInfo.dvMolecularWeight;
}

//------------------------------------------------------------------------------
//! Print Fluid Information Ideal Gas
//------------------------------------------------------------------------------
void EOS_IdealGas_Print_Fluid_Information() {
    printf("=============================================================================\n");
    info("FLUID INFORMATION : %s", SogFluidInfo.csFluidName);
    printf("=============================================================================\n");
    info("Molecular Weight                  = %15.6f kg/kmol", SogFluidInfo.dvMolecularWeight);
    info("Gas Constant                      = %15.6f J/kg-K",  SogFluidInfo.dvGasConstant);
}

//------------------------------------------------------------------------------
//! Set the Reference Property for Fluid : Generic (Reference Density is Computed)
//  Pressure (Pa), Temperature (K) and Length (m)
//------------------------------------------------------------------------------
void EOS_IdealGas_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength) {
    // Set the Reference Properties in SI Units
    dvgTemperature_Ref       = dvTemperature;
    dvgLength_Ref            = dvLength;
    dvgRatioSpecificHeat_Ref = 1.4;
    dvgDensity_Ref           = dvPressure/(SogFluidInfo.dvGasConstant*dvTemperature);
    dvgSpeedSound_Ref        = sqrt(SogFluidInfo.dvGasConstant*dvgTemperature_Ref);
    
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
//! Compute Equation of State Properties for Ideal Gas
//------------------------------------------------------------------------------
void EOS_IdealGas_Get_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2, dvSpeedSound, dvMach;
    double dvEntropy, dvEnthalpy, dvInternalEnergy, dvTotalEnergy, dvTotalEnthalpy;
    double dvCv, dvCp, dvR, dvGamma;
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
            error("EOS_IdealGas_Get_Properties:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        default:
            error("EOS_IdealGas_Get_Properties:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvInternalEnergy = dvTotalEnergy - 0.5*dvQ2;
    dvTotalEnthalpy  = dvTotalEnergy + dvPressure/dvRho;
    dvEnthalpy       = dvTotalEnthalpy - 0.5*dvQ2;
    dvSpeedSound     = sqrt((dvGamma * dvPressure) / dvRho);
    dvMach           = sqrt(dvQ2)/dvSpeedSound;
    dvEntropy        = dvPressure/pow(dvRho, dvGamma);
    dvCp             = dvGamma*dvR/(dvGamma - 1.0);
    dvCv             = dvCp/dvGamma;
    
    // Update the Output (is Dimensional)
    dpPropertyOut[ 0] = dvRho;
    dpPropertyOut[ 1] = dvRho;
    dpPropertyOut[ 2] = dvRho;
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
            error("EOS_IdealGas_Get_Properties:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
}

//------------------------------------------------------------------------------
//! Compute the EOS Ideal Gas Extended Properties based on Variable Type
//! Input and Output property are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_IdealGas_Get_Extended_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2, dvSpeedSound, dvMach;
    double dvEntropy, dvEnthalpy, dvInternalEnergy, dvTotalEnergy, dvTotalEnthalpy;
    double dvCv, dvCp, dvR, dvGamma;
    double dvDPDRho, dvDPDT, dvDRhoDT, dvDRhoDP, dvD2PDRho2, dvD2PDT2, dvD2PDTRho;
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
            error("EOS_IdealGas_Get_Properties:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        default:
            error("EOS_IdealGas_Get_Properties:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvInternalEnergy = dvTotalEnergy - 0.5*dvQ2;
    dvTotalEnthalpy  = dvTotalEnergy + dvPressure/dvRho;
    dvEnthalpy       = dvTotalEnthalpy - 0.5*dvQ2;
    dvSpeedSound     = sqrt((dvGamma * dvPressure) / dvRho);
    dvMach           = sqrt(dvQ2)/dvSpeedSound;
    dvEntropy        = dvPressure/pow(dvRho, dvGamma);
    dvCp             = dvGamma*dvR/(dvGamma - 1.0);
    dvCv             = dvCp/dvGamma;
    // Compute the derivative
    dvDPDRho         = dvR*dvTemperature;
    dvDPDT           = dvR*dvRho;
    dvDRhoDT         = -dvRho/dvTemperature;
    dvDRhoDP         = 1.0/dvDPDRho;
    dvD2PDRho2       = 0.0;
    dvD2PDT2         = 0.0;
    dvD2PDTRho       = dvR;
    dvDHdT_Rho       = dvCp;
    dvDHdT_P         = dvCp;
    dvDHDRho_T       = 0.0;
    dvDHDRho_P       = -dvCp*dvTemperature/dvRho;
    dvDHDP_T         = 0.0;
    dvDHDP_Rho       = dvGamma/((dvGamma - 1.0)*dvRho);
    dvDEdT_Rho       = dvDHdT_Rho - dvDPDT/dvRho;
    dvDEdT_P         = dvDHdT_P + dvPressure*dvDRhoDT/(dvRho*dvRho);
    dvDEDRho_T       = dvDHDRho_T + dvPressure/(dvRho*dvRho) - dvDPDRho/dvRho;
    dvDEDRho_P       = dvDHDRho_P + dvPressure/(dvRho*dvRho);
    dvDEDP_T         = dvDHDP_T - 1.0/dvRho + dvPressure*dvDRhoDP/(dvRho*dvRho);
    dvDEDP_Rho       = dvDHDP_Rho - 1.0/dvRho;
    
    // Update the Output (is Dimensional)
    dpPropertyOut[ 0] = dvRho;
    dpPropertyOut[ 1] = dvRho;
    dpPropertyOut[ 2] = dvRho;
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
            error("EOS_IdealGas_Get_Properties:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_Density(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature, dvR;
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
            error("EOS_IdealGas_Get_Density:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvR   = SogFluidInfo.dvGasConstant;
    dvRho = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            break;
        default:
            error("EOS_IdealGas_Get_Density:2: Undefined Variable Type - %d", ivVariableType);
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
            error("EOS_IdealGas_Get_Density:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvRho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DensityLiquid(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    return EOS_IdealGas_Get_Density(ivDimIOType, ivVariableType, dpVariableIn);
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DensityVapor(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    return EOS_IdealGas_Get_Density(ivDimIOType, ivVariableType, dpVariableIn);
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_Pressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_Pressure:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma    = dvgRatioSpecificHeat_Ref;
    dvR        = SogFluidInfo.dvGasConstant;
    dvPressure = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvPressure      = daVariableDimensional[4];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            break;
        default:
            error("EOS_IdealGas_Get_Pressure:2: Undefined Variable Type - %d", ivVariableType);
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
            error("EOS_IdealGas_Get_Pressure:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvPressure;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_Temperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_Temperature:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma       = dvgRatioSpecificHeat_Ref;
    dvR           = SogFluidInfo.dvGasConstant;
    dvTemperature = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvTemperature   = daVariableDimensional[4];
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvTemperature   = daVariableDimensional[4];
            break;
        default:
            error("EOS_IdealGas_Get_Temperature:2: Undefined Variable Type - %d", ivVariableType);
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
            error("EOS_IdealGas_Get_Temperature:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvTemperature;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_SpeedSound(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2, dvSpeedSound;
    double dvTotalEnergy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_SpeedSound:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma      = dvgRatioSpecificHeat_Ref;
    dvR          = SogFluidInfo.dvGasConstant;
    dvSpeedSound = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            break;
        default:
            error("EOS_IdealGas_Get_SpeedSound:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvSpeedSound     = sqrt((dvGamma * dvPressure) / dvRho);
    
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
            error("EOS_IdealGas_Get_SpeedSound:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvSpeedSound;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_Mach(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2, dvSpeedSound, dvMach;
    double dvTotalEnergy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_Mach:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    dvMach  = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvPressure      = daVariableDimensional[4];
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            break;
        default:
            error("EOS_IdealGas_Get_Mach:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvSpeedSound     = sqrt((dvGamma * dvPressure) / dvRho);
    dvMach           = sqrt(dvQ2)/dvSpeedSound;
    
    //-END----------------Dimensional Computations------------------------------
    return dvMach;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_Entropy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvEntropy, dvTotalEnergy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_Entropy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma   = dvgRatioSpecificHeat_Ref;
    dvR       = SogFluidInfo.dvGasConstant;
    dvEntropy = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            break;
        default:
            error("EOS_IdealGas_Get_Entropy:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvEntropy        = dvPressure/pow(dvRho, dvGamma);
        
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
            error("EOS_IdealGas_Get_Entropy:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvEntropy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_Enthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvEnthalpy,dvTotalEnergy, dvTotalEnthalpy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_Enthalpy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma    = dvgRatioSpecificHeat_Ref;
    dvR        = SogFluidInfo.dvGasConstant;
    dvEnthalpy = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvPressure      = daVariableDimensional[4];
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        default:
            error("EOS_IdealGas_Get_Enthalpy:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvTotalEnthalpy  = dvTotalEnergy + dvPressure/dvRho;
    dvEnthalpy       = dvTotalEnthalpy - 0.5*dvQ2;
        
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
            error("EOS_IdealGas_Get_Enthalpy:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvEnthalpy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_InternalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvInternalEnergy, dvTotalEnergy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_InternalEnergy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma          = dvgRatioSpecificHeat_Ref;
    dvR              = SogFluidInfo.dvGasConstant;
    dvInternalEnergy = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvPressure      = daVariableDimensional[4];
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        default:
            error("EOS_IdealGas_Get_InternalEnergy:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvInternalEnergy = dvTotalEnergy - 0.5*dvQ2;
        
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
            error("EOS_IdealGas_Get_InternalEnergy:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvInternalEnergy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_TotalEnthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy, dvTotalEnthalpy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_TotalEnthalpy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma         = dvgRatioSpecificHeat_Ref;
    dvR             = SogFluidInfo.dvGasConstant;
    dvTotalEnthalpy = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvPressure      = daVariableDimensional[4];
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        default:
            error("EOS_IdealGas_Get_TotalEnthalpy:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvTotalEnthalpy  = dvTotalEnergy + dvPressure/dvRho;
        
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
            error("EOS_IdealGas_Get_TotalEnthalpy:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvTotalEnthalpy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_TotalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvR, dvGamma;
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
            error("EOS_IdealGas_Get_TotalEnergy:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma       = dvgRatioSpecificHeat_Ref;
    dvR           = SogFluidInfo.dvGasConstant;
    dvTotalEnergy = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvPressure      = daVariableDimensional[4];
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1];
            dvVelocityV     = daVariableDimensional[2];
            dvVelocityW     = daVariableDimensional[3];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvTotalEnergy   = dvPressure/((dvGamma - 1.0)*dvRho) + 0.5*dvQ2;
            break;
        default:
            error("EOS_IdealGas_Get_TotalEnergy:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
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
            error("EOS_IdealGas_Get_TotalEnergy:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvTotalEnergy;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_HeatCapacityCv(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvCv, dvCp, dvR, dvGamma;
      
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    dvCv    = 0.0;
    
    // Compute the Remaining Properties in SI units
    dvCp = dvGamma*dvR/(dvGamma - 1.0);
    dvCv = dvCp/dvGamma;
        
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
            error("EOS_IdealGas_Get_HeatCapacityCv:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvCv;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_HeatCapacityCp(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvCp, dvR, dvGamma;
      
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    dvCp    = 0.0;
    
    // Compute the Remaining Properties in SI units
    dvCp = dvGamma*dvR/(dvGamma - 1.0);
        
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
            error("EOS_IdealGas_Get_HeatCapacityCv:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvCp;
}

// Compute First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DPressureDDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvR, dvGamma;
    double dvDPDRho;
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
            error("EOS_IdealGas_Get_DPressureDDensity:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma  = dvgRatioSpecificHeat_Ref;
    dvR      = SogFluidInfo.dvGasConstant;
    dvDPDRho = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvTemperature   = daVariableDimensional[4];
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvTemperature   = daVariableDimensional[4];
            break;
        default:
            error("EOS_IdealGas_Get_DPressureDDensity:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDPDRho         = dvR*dvTemperature;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDPDRho = EOS_Internal_NonDimensionalize_DPressureDDensity(dvDPDRho);
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
            dvDPDRho = EOS_Internal_NonDimensionalize_DPressureDDensity(dvDPDRho);
            break;
        default:
            error("EOS_IdealGas_Get_DPressureDDensity:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDPDRho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DPressureDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvR;
    double dvDPDT;
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
            error("EOS_IdealGas_Get_DPressureDTemperature:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvR     = SogFluidInfo.dvGasConstant;
    dvDPDT  = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            break;
        default:
            error("EOS_IdealGas_Get_DPressureDTemperature:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDPDT           = dvR*dvRho;
        
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDPDT = EOS_Internal_NonDimensionalize_DPressureDTemperature(dvDPDT);
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
            dvDPDT = EOS_Internal_NonDimensionalize_DPressureDTemperature(dvDPDT);
            break;
        default:
            error("EOS_IdealGas_Get_DPressureDTemperature:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDPDT;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DDensityDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvR, dvGamma;
    double dvDRhoDT;
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
            error("EOS_IdealGas_Get_DDensityDTemperature:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma  = dvgRatioSpecificHeat_Ref;
    dvR      = SogFluidInfo.dvGasConstant;
    dvDRhoDT = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            break;
        default:
            error("EOS_IdealGas_Get_DDensityDTemperature:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDRhoDT         = -dvRho/dvTemperature;
        
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDRhoDT = EOS_Internal_NonDimensionalize_DDensityDTemperature(dvDRhoDT);
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
            dvDRhoDT = EOS_Internal_NonDimensionalize_DDensityDTemperature(dvDRhoDT);
            break;
        default:
            error("EOS_IdealGas_Get_DDensityDTemperature:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDRhoDT;
}

// Compute Second Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_D2PressureDDensity2(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvD2PDRho2;
    
    //-START--------------Dimensional Computations------------------------------
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvD2PDRho2 = 0.0;
        
    //-END----------------Dimensional Computations------------------------------
    return dvD2PDRho2;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_D2PressureDTemperature2(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvD2PDT2;
       
    //-START--------------Dimensional Computations------------------------------
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvD2PDT2         = 0.0;
    
    //-END----------------Dimensional Computations------------------------------
    return dvD2PDT2;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_D2PressureDTemperatureDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvR;
    double dvD2PDTRho;
      
    //-START--------------Dimensional Computations------------------------------
    dvR = SogFluidInfo.dvGasConstant;
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvD2PDTRho = dvR;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvD2PDTRho = EOS_Internal_NonDimensionalize_D2PressureDTemperatureDensity(dvD2PDTRho);
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
            dvD2PDTRho = EOS_Internal_NonDimensionalize_D2PressureDTemperatureDensity(dvD2PDTRho);
            break;
        default:
            error("EOS_IdealGas_Get_D2PressureDTemperatureDensity:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvD2PDTRho;
}

// Enthalpy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DEnthalpyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvCp, dvR, dvGamma;
    double dvDHdT_Rho;
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    dvDHdT_Rho = 0.0;
        
    // Compute the Remaining Properties in SI units
    dvCp       = dvGamma*dvR/(dvGamma - 1.0);
    // Compute the derivative
    dvDHdT_Rho = dvCp;   
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDHdT_Rho = EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CDensity(dvDHdT_Rho);
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
            dvDHdT_Rho = EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CDensity(dvDHdT_Rho);
            break;
        default:
            error("EOS_IdealGas_Get_DEnthalpyDTemperature_CDensity:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDHdT_Rho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DEnthalpyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvCp, dvR, dvGamma;
    double dvDHdT_P;
       
    //-START--------------Dimensional Computations------------------------------
    dvGamma  = dvgRatioSpecificHeat_Ref;
    dvR      = SogFluidInfo.dvGasConstant;
    dvDHdT_P = 0.0;
    
    // Compute the Remaining Properties in SI units
    dvCp     = dvGamma*dvR/(dvGamma - 1.0);
    // Compute the derivative
    dvDHdT_P = dvCp;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDHdT_P = EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CPressure(dvDHdT_P);
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
            dvDHdT_P = EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CPressure(dvDHdT_P);
            break;
        default:
            error("EOS_IdealGas_Get_DEnthalpyDTemperature_CPressure:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDHdT_P;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DEnthalpyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvDHDRho_T;
        
    //-START--------------Dimensional Computations------------------------------   
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDHDRho_T  = 0.0;
           
    //-END----------------Dimensional Computations------------------------------
    return dvDHDRho_T;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DEnthalpyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvCp, dvR, dvGamma;
    double dvDHDRho_P;
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
            error("EOS_IdealGas_Get_DEnthalpyDDensity_CPressure:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            break;
        default:
            error("EOS_IdealGas_Get_DEnthalpyDDensity_CPressure:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvCp       = dvGamma*dvR/(dvGamma - 1.0);
    // Compute the derivative
    dvDHDRho_P = -dvCp*dvTemperature/dvRho;
                
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDHDRho_P = EOS_Internal_NonDimensionalize_DEnthalpyDDensity_CPressure(dvDHDRho_P);
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
            dvDHDRho_P = EOS_Internal_NonDimensionalize_DEnthalpyDDensity_CPressure(dvDHDRho_P);
            break;
        default:
            error("EOS_IdealGas_Get_DEnthalpyDDensity_CPressure:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDHDRho_P;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DEnthalpyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    double dvDHDP_T;
    
    //-START--------------Dimensional Computations------------------------------
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDHDP_T = 0.0;
    
    //-END----------------Dimensional Computations------------------------------
    return dvDHDP_T;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DEnthalpyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvR, dvGamma;
    double dvDHDP_Rho;
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
            error("EOS_IdealGas_Get_DEnthalpyDPressure_CDensity:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma = dvgRatioSpecificHeat_Ref;
    dvR     = SogFluidInfo.dvGasConstant;
    
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            break;
        default:
            error("EOS_IdealGas_Get_DEnthalpyDPressure_CDensity:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDHDP_Rho  = dvGamma/((dvGamma - 1.0)*dvRho);
        
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDHDP_Rho = EOS_Internal_NonDimensionalize_DEnthalpyDPressure_CDensity(dvDHDP_Rho);
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
            dvDHDP_Rho = EOS_Internal_NonDimensionalize_DEnthalpyDPressure_CDensity(dvDHDP_Rho);
            break;
        default:
            error("EOS_IdealGas_Get_DEnthalpyDPressure_CDensity:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDHDP_Rho;
}

// Internal Energy First Derivatives
//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DInternalEnergyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvCp, dvR, dvGamma;
    double dvDPDT;
    double dvDHdT_Rho;
    double dvDEdT_Rho;
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
            error("EOS_IdealGas_Get_DInternalEnergyDTemperature_CDensity:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma    = dvgRatioSpecificHeat_Ref;
    dvR        = SogFluidInfo.dvGasConstant;
    dvDEdT_Rho = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDTemperature_CDensity:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvCp       = dvGamma*dvR/(dvGamma - 1.0);
    // Compute the derivative
    dvDPDT     = dvR*dvRho;
    dvDHdT_Rho = dvCp;
    dvDEdT_Rho = dvDHdT_Rho - dvDPDT/dvRho;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDEdT_Rho = EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CDensity(dvDEdT_Rho);
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
            dvDEdT_Rho = EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CDensity(dvDEdT_Rho);
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDTemperature_CDensity:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDEdT_Rho;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DInternalEnergyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvCp, dvR, dvGamma;
    double dvDRhoDT;
    double dvDHdT_P;
    double dvDEdT_P;
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
            error("EOS_IdealGas_Get_DInternalEnergyDTemperature_CPressure:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma  = dvgRatioSpecificHeat_Ref;
    dvR      = SogFluidInfo.dvGasConstant;
    dvDEdT_P = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDTemperature_CPressure:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvCp     = dvGamma*dvR/(dvGamma - 1.0);
    // Compute the derivative
    dvDRhoDT = -dvRho/dvTemperature;
    dvDHdT_P = dvCp;
    dvDEdT_P = dvDHdT_P + dvPressure*dvDRhoDT/(dvRho*dvRho);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDEdT_P = EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CPressure(dvDEdT_P);
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
            dvDEdT_P = EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CPressure(dvDEdT_P);
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDTemperature_CPressure:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDEdT_P;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DInternalEnergyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvR, dvGamma;
    double dvDPDRho;
    double dvDHDRho_T;
    double dvDEDRho_T;
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
            error("EOS_IdealGas_Get_DInternalEnergyDDensity_CTemperature:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma    = dvgRatioSpecificHeat_Ref;
    dvR        = SogFluidInfo.dvGasConstant;
    dvDEDRho_T = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDDensity_CTemperature:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDPDRho         = dvR*dvTemperature;
    dvDHDRho_T       = 0.0;
    dvDEDRho_T       = dvDHDRho_T + dvPressure/(dvRho*dvRho) - dvDPDRho/dvRho;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDEDRho_T = EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CTemperature(dvDEDRho_T);
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
            dvDEDRho_T = EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CTemperature(dvDEDRho_T);
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDDensity_CTemperature:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDEDRho_T;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DInternalEnergyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvCp, dvR, dvGamma;
    double dvDHDRho_P;
    double dvDEDRho_P;
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
            error("EOS_IdealGas_Get_DInternalEnergyDDensity_CPressure:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma    = dvgRatioSpecificHeat_Ref;
    dvR        = SogFluidInfo.dvGasConstant;
    dvDEDRho_P = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDDensity_CPressure:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    dvCp       = dvGamma*dvR/(dvGamma - 1.0);
    // Compute the derivative
    dvDHDRho_P = -dvCp*dvTemperature/dvRho;
    dvDEDRho_P = dvDHDRho_P + dvPressure/(dvRho*dvRho);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDEDRho_P = EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CPressure(dvDEDRho_P);
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
            dvDEDRho_P = EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CPressure(dvDEDRho_P);
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDDensity_CPressure:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDEDRho_P;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DInternalEnergyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvVelocityU, dvVelocityV, dvVelocityW, dvQ2;
    double dvTotalEnergy;
    double dvR, dvGamma;
    double dvDPDRho, dvDRhoDP;
    double dvDHDP_T;
    double dvDEDP_T;
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
            error("EOS_IdealGas_Get_DInternalEnergyDPressure_CTemperature:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma  = dvgRatioSpecificHeat_Ref;
    dvR      = SogFluidInfo.dvGasConstant;
    dvDEDP_T = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            dvVelocityU     = daVariableDimensional[1]/dvRho;
            dvVelocityV     = daVariableDimensional[2]/dvRho;
            dvVelocityW     = daVariableDimensional[3]/dvRho;
            dvTotalEnergy   = daVariableDimensional[4]/dvRho;
            dvQ2            = dvVelocityU*dvVelocityU + dvVelocityV*dvVelocityV + dvVelocityW*dvVelocityW;
            dvPressure      = (dvGamma - 1.0)*dvRho*(dvTotalEnergy - 0.5*dvQ2);
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            dvPressure      = daVariableDimensional[4];
            dvTemperature   = dvPressure/(dvR*dvRho);
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvPressure      = dvTemperature*dvR*dvRho;
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDPressure_CTemperature:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDPDRho = dvR*dvTemperature;
    dvDRhoDP = 1.0/dvDPDRho;
    dvDHDP_T = 0.0;
    dvDEDP_T = dvDHDP_T - 1.0/dvRho + dvPressure*dvDRhoDP/(dvRho*dvRho);
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDEDP_T = EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CTemperature(dvDEDP_T);
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
            dvDEDP_T = EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CTemperature(dvDEDP_T);
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDPressure_CTemperature:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDEDP_T;
}

//------------------------------------------------------------------------------
//!
//------------------------------------------------------------------------------
double EOS_IdealGas_Get_DInternalEnergyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn) {
    int i;
    double dvRho, dvPressure, dvTemperature;
    double dvR, dvGamma;
    double dvDHDP_Rho;
    double dvDEDP_Rho;
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
            error("EOS_IdealGas_Get_DInternalEnergyDPressure_CDensity:1: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    
    //-START--------------Dimensional Computations------------------------------
    dvGamma    = dvgRatioSpecificHeat_Ref;
    dvR        = SogFluidInfo.dvGasConstant;
    dvDEDP_Rho = 0.0;
    // Compute the EOS Based on Variable Type of dpVariableIn
    switch (ivVariableType) {
        // Conservative Variable Formulation
        case EOS_VARIABLE_CON:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case EOS_VARIABLE_PUT:
            dvPressure      = daVariableDimensional[0];
            dvTemperature   = daVariableDimensional[4];
            dvRho           = dvPressure/(dvR*dvTemperature);
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case EOS_VARIABLE_RUP:
            dvRho           = daVariableDimensional[0];
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case EOS_VARIABLE_RUT:
            dvRho           = daVariableDimensional[0];
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDPressure_CDensity:2: Undefined Variable Type - %d", ivVariableType);
            break;
    }
    
    // Compute the Remaining Properties in SI units
    // Compute the derivative
    dvDHDP_Rho = dvGamma/((dvGamma - 1.0)*dvRho);
    dvDEDP_Rho = dvDHDP_Rho - 1.0/dvRho;
    
    // Dimensionalize the Output Properties based on I/O Type
    switch (ivDimIOType) {
        // Type Input: Non-Dimensional  Output: Non-Dimensional
        case EOS_DIMENSIONAL_IO_ND_ND:
            dvDEDP_Rho = EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CDensity(dvDEDP_Rho);
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
            dvDEDP_Rho = EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CDensity(dvDEDP_Rho);
            break;
        default:
            error("EOS_IdealGas_Get_DInternalEnergyDPressure_CDensity:3: Undefined Dimensional Input/Output -%d", ivDimIOType);
            break;
    }
    //-END----------------Dimensional Computations------------------------------
    return dvDEDP_Rho;
}

