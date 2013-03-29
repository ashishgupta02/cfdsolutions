/*******************************************************************************
 * File:        EOS_Density_Temperature.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_DENSITY_TEMPERATURE_H
#define	_EOS_DENSITY_TEMPERATURE_H

#ifdef	__cplusplus
extern "C" {
#endif

    // -------------- Density and Temperature Formulation ----------------------
    // This functions are valid for All Regions:
    // 1) Super Critical State
    // 2) Subcooled Compressed Liquid
    // 3) Super Heated Vapor
    // 4) Multi-Phase Liquid-Vapor
    // -------------------------------------------------------------------------
    // Compute EOS Properties
    double EOS_DT_Get_Pressure(double dvDensity, double dvTemperature);
    double EOS_DT_Get_SpeedSound(double dvDensity, double dvTemperature);
    double EOS_DT_Get_Entropy(double dvDensity, double dvTemperature);
    double EOS_DT_Get_Enthalpy(double dvDensity, double dvTemperature);
    double EOS_DT_Get_InternalEnergy(double dvDensity, double dvTemperature);
    double EOS_DT_Get_HeatCapacityCv(double dvDensity, double dvTemperature);
    double EOS_DT_Get_HeatCapacityCp(double dvDensity, double dvTemperature);
    double EOS_DT_Get_TotalEnergy(double dvDensity, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W);
    double EOS_DT_Get_Mach(double dvDensity, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W);
    // Compute First Derivatives
    double EOS_DT_Get_DPressureDDensity(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DPressureDTemperature(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DDensityDPressure(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DDensityDTemperature(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DTemperatureDDensity(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DTemperatureDPressure(double dvDensity, double dvTemperature);
    // Enthalpy First Derivatives
    double EOS_DT_Get_DEnthalpyDTemperature_CDensity(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DEnthalpyDTemperature_CPressure(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DEnthalpyDDensity_CTemperature(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DEnthalpyDDensity_CPressure(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DEnthalpyDPressure_CTemperature(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DEnthalpyDPressure_CDensity(double dvDensity, double dvTemperature);
    // Internal Energy First Derivatives
    double EOS_DT_Get_DInternalEnergyDTemperature_CDensity(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DInternalEnergyDTemperature_CPressure(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DInternalEnergyDDensity_CTemperature(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DInternalEnergyDDensity_CPressure(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DInternalEnergyDPressure_CTemperature(double dvDensity, double dvTemperature);
    double EOS_DT_Get_DInternalEnergyDPressure_CDensity(double dvDensity, double dvTemperature);
    // Compute Second Derivatives
    double EOS_DT_Get_D2PressureDDensity2(double dvDensity, double dvTemperature);
    double EOS_DT_Get_D2PressureDTemperature2(double dvDensity, double dvTemperature);
    double EOS_DT_Get_D2PressureDTemperatureDensity(double dvDensity, double dvTemperature);
    
    // Transformation Matrix
    void EOS_DT_Get_Transformation_Matrix_CON_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_RUT_To_CON(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_CON_To_RUT(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_RUT_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_RUP_To_RUT(int ivDimIOType, double *dpVariableIn, double **Matrix);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_DENSITY_TEMPERATURE_H */

