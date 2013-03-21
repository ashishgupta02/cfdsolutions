/*******************************************************************************
 * File:        EOS.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_H
#define	_EOS_H

#ifdef	__cplusplus
extern "C" {
#endif
    // Initialize the EOS Data Structure
    void EOS_Init(void);
    // Finalize the EOS Data Structure
    void EOS_Finalize(void);
    // Setup the NIST with fluids
    void EOS_Set(void);
    
    // Set the Reference Property for Fluid
    void EOS_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength);
    // Set the Reference Property using Density, Temperature and Length : Generic
    void EOS_Set_Reference_Properties_DensityTemperatureLength(double dvDensity, double dvTemperature, double dvLength);
    // Set the Reference Property using Density, Pressure and Length : LMRoe
    void EOS_Set_Reference_Properties_DensityPressureLength(double dvDensity, double dvTemperature, double dvLength);
    // Set the Reference Property using Density, Temperature, Velocity and Length : BTW
    void EOS_Set_Reference_Properties_DensityTemperatureVelocityLength(double dvDensity, double dvTemperature, double dvVelocity, double dvLength);
    // Get the Reference Property for Fluid
    void EOS_Get_Reference_Properties();
    
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
    
    // -------------- Pressure and Temperature Formulation ----------------------
    // This functions are valid Only in this Regions:
    // 1) Super Critical State
    // 2) Subcooled Compressed Liquid
    // 3) Super Heated Vapor
    // -------------------------------------------------------------------------
    // Compute EOS Properties
    double EOS_PT_Get_Pressure(double dvPressure, double dvTemperature);
    double EOS_PT_Get_SpeedSound(double dvPressure, double dvTemperature);
    double EOS_PT_Get_Entropy(double dvPressure, double dvTemperature);
    double EOS_PT_Get_Enthalpy(double dvPressure, double dvTemperature);
    double EOS_PT_Get_InternalEnergy(double dvPressure, double dvTemperature);
    double EOS_PT_Get_HeatCapacityCv(double dvPressure, double dvTemperature);
    double EOS_PT_Get_HeatCapacityCp(double dvPressure, double dvTemperature);
    double EOS_PT_Get_TotalEnergy(double dvPressure, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W);
    double EOS_PT_Get_Mach(double dvPressure, double dvTemperature, double dvVelocity_U, double dvVelocity_V, double dvVelocity_W);
    // Compute First Derivatives
    double EOS_PT_Get_DPressureDDensity(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DPressureDTemperature(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DDensityDPressure(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DDensityDTemperature(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DTemperatureDDensity(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DTemperatureDPressure(double dvPressure, double dvTemperature);
    // Enthalpy First Derivatives
    double EOS_PT_Get_DEnthalpyDTemperature_CDensity(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DEnthalpyDTemperature_CPressure(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DEnthalpyDDensity_CTemperature(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DEnthalpyDDensity_CPressure(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DEnthalpyDPressure_CTemperature(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DEnthalpyDPressure_CDensity(double dvPressure, double dvTemperature);
    // Internal Energy First Derivatives
    double EOS_PT_Get_DInternalEnergyDTemperature_CDensity(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DInternalEnergyDTemperature_CPressure(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DInternalEnergyDDensity_CTemperature(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DInternalEnergyDDensity_CPressure(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DInternalEnergyDPressure_CTemperature(double dvPressure, double dvTemperature);
    double EOS_PT_Get_DInternalEnergyDPressure_CDensity(double dvPressure, double dvTemperature);
    // Compute Second Derivatives
    double EOS_PT_Get_D2PressureDDensity2(double dvPressure, double dvTemperature);
    double EOS_PT_Get_D2PressureDTemperature2(double dvPressure, double dvTemperature);
    double EOS_PT_Get_D2PressureDTemperatureDensity(double dvPressure, double dvTemperature);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_H */

