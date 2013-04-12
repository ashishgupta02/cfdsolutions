/*******************************************************************************
 * File:        EOS_IdealGas.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_IDEALGAS_H
#define	_EOS_IDEALGAS_H

#ifdef	__cplusplus
extern "C" {
#endif

    /* Functions */
    void EOS_IdealGas_Init();
    void EOS_IdealGas_Finalize();
    void EOS_IdealGas_Set();
    void EOS_IdealGas_Set_Fluid_Information();
    void EOS_IdealGas_Print_Fluid_Information();
    void EOS_IdealGas_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength);
    
    void EOS_IdealGas_Get_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut);
    void EOS_IdealGas_Get_Extended_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut);
    
    double EOS_IdealGas_Get_Density(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DensityLiquid(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DensityVapor(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_Pressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_Temperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_SpeedSound(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_Mach(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_Entropy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_Enthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_InternalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_TotalEnthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_TotalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_HeatCapacityCv(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_HeatCapacityCp(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DPressureDDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DPressureDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DDensityDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_D2PressureDDensity2(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_D2PressureDTemperature2(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_D2PressureDTemperatureDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DEnthalpyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DEnthalpyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DEnthalpyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DEnthalpyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DEnthalpyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DEnthalpyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DInternalEnergyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DInternalEnergyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DInternalEnergyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DInternalEnergyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DInternalEnergyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_IdealGas_Get_DInternalEnergyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);

    // Compute Additional Properties
    void   EOS_IdealGas_Get_PT_Density(int ivDimIOType, double dvPressure, double dvTemperature, double *dpDensityOut);
    double EOS_IdealGas_Get_DH_SpeedSound(int ivDimIOType, double dvDensity, double dvEnthalpy);
    double EOS_IdealGas_Get_DH_Pressure(int ivDimIOType, double dvDensity, double dvEnthalpy);
    double EOS_IdealGas_Get_DH_Temperature(int ivDimIOType, double dvDensity, double dvEnthalpy);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_IDEALGAS_H */

