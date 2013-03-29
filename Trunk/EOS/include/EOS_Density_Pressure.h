/*******************************************************************************
 * File:        EOS_Density_Pressure.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_DENSITY_PRESSURE_H
#define	_EOS_DENSITY_PRESSURE_H

#ifdef	__cplusplus
extern "C" {
#endif

    // -------------- Density and Pressure Formulation -------------------------
    // This functions are valid for All Regions:
    // 1) Super Critical State
    // 2) Subcooled Compressed Liquid
    // 3) Super Heated Vapor
    // 4) Multi-Phase Liquid-Vapor
    // -------------------------------------------------------------------------
    // Compute EOS Properties
    double EOS_DP_Get_Pressure(double dvDensity, double dvPressure);
    double EOS_DP_Get_SpeedSound(double dvDensity, double dvPressure);
    double EOS_DP_Get_Entropy(double dvDensity, double dvPressure);
    double EOS_DP_Get_Enthalpy(double dvDensity, double dvPressure);
    double EOS_DP_Get_InternalEnergy(double dvDensity, double dvPressure);
    double EOS_DP_Get_HeatCapacityCv(double dvDensity, double dvPressure);
    double EOS_DP_Get_HeatCapacityCp(double dvDensity, double dvPressure);
    double EOS_DP_Get_TotalEnergy(double dvDensity, double dvPressure, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W);
    double EOS_DP_Get_Mach(double dvDensity, double dvPressure, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W);
    // Compute First Derivatives
    double EOS_DP_Get_DPressureDDensity(double dvDensity, double dvPressure);
    double EOS_DP_Get_DPressureDTemperature(double dvDensity, double dvPressure);
    double EOS_DP_Get_DDensityDPressure(double dvDensity, double dvPressure);
    double EOS_DP_Get_DDensityDTemperature(double dvDensity, double dvPressure);
    double EOS_DP_Get_DTemperatureDDensity(double dvDensity, double dvPressure);
    double EOS_DP_Get_DTemperatureDPressure(double dvDensity, double dvPressure);
    // Enthalpy First Derivatives
    double EOS_DP_Get_DEnthalpyDTemperature_CDensity(double dvDensity, double dvPressure);
    double EOS_DP_Get_DEnthalpyDTemperature_CPressure(double dvDensity, double dvPressure);
    double EOS_DP_Get_DEnthalpyDDensity_CTemperature(double dvDensity, double dvPressure);
    double EOS_DP_Get_DEnthalpyDDensity_CPressure(double dvDensity, double dvPressure);
    double EOS_DP_Get_DEnthalpyDPressure_CTemperature(double dvDensity, double dvPressure);
    double EOS_DP_Get_DEnthalpyDPressure_CDensity(double dvDensity, double dvPressure);
    // Internal Energy First Derivatives
    double EOS_DP_Get_DInternalEnergyDTemperature_CDensity(double dvDensity, double dvPressure);
    double EOS_DP_Get_DInternalEnergyDTemperature_CPressure(double dvDensity, double dvPressure);
    double EOS_DP_Get_DInternalEnergyDDensity_CTemperature(double dvDensity, double dvPressure);
    double EOS_DP_Get_DInternalEnergyDDensity_CPressure(double dvDensity, double dvPressure);
    double EOS_DP_Get_DInternalEnergyDPressure_CTemperature(double dvDensity, double dvPressure);
    double EOS_DP_Get_DInternalEnergyDPressure_CDensity(double dvDensity, double dvPressure);
    // Compute Second Derivatives
    double EOS_DP_Get_D2PressureDDensity2(double dvDensity, double dvPressure);
    double EOS_DP_Get_D2PressureDTemperature2(double dvDensity, double dvPressure);
    double EOS_DP_Get_D2PressureDTemperatureDensity(double dvDensity, double dvPressure);

    // Transformation Matrix
    void EOS_DP_Get_Transformation_Matrix_RUP_To_CON(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DP_Get_Transformation_Matrix_CON_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_DENSITY_PRESSURE_H */

