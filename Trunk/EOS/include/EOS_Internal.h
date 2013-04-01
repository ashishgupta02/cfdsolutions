/*******************************************************************************
 * File:        EOS_Internal.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_INTERNAL_H
#define	_EOS_INTERNAL_H

#ifdef	__cplusplus
extern "C" {
#endif
    
    /*!
    * \brief Fluid Property Structure
    */
    typedef struct {
        int    ivEOSModelType;
        char   csFluidName[255];
        double dvTemperature_Crit;
        double dvPressure_Crit;
        double dvDensity_Crit;
        double dvMolecularWeight;
        double dvTPTemperature;
        double dvNBPTemperature;
        double dvZ_Crit;
        double dvAcentric_Fact;
        double dvDipole_Mom;
        double dvGasConstant;
        double dvTemperature_Min;
        double dvTemperature_Max;
        double dvPressure_Max;
        double dvDensity_Max;
    } EOS_S_FluidInfo;
    
    /* Common Shared Data */
    /* Fluid Information */
    extern EOS_S_FluidInfo SogFluidInfo;
    /* Reference Properties */
    extern int    ivgSet_Ref;
    extern double dvgTemperature_Ref;
    extern double dvgPressure_Ref;
    extern double dvgDensity_Ref;
    extern double dvgVelocity_Ref;
    extern double dvgLength_Ref;
    extern double dvgMach_Ref;
    extern double dvgSpeedSound_Ref;
    extern double dvgTime_Ref;
    extern double dvgEnthalpy_Ref;
    extern double dvgTotalEnthalpy_Ref;
    extern double dvgInternalEnergy_Ref;
    extern double dvgTotalEnergy_Ref;
    extern double dvgEntropy_Ref;
    extern double dvgEntropyConst_Ref;
    extern double dvgGasConstant_Ref;
    extern double dvgHeatCapacityCv_Ref;
    extern double dvgHeatCapacityCp_Ref;
    extern double dvgRatioSpecificHeat_Ref;
    
    /* Functions */
    void EOS_Internal_Init();
    void EOS_Internal_Finalize();
    void EOS_Internal_Init_Fluid_Information();
    void EOS_Internal_Init_Reference_Properties();
    
    void   EOS_Internal_Dimensionalize_Variables(int ivVariableType, double *dpVariableIn, double *dpVariableOut);
    void   EOS_Internal_Dimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut);
    void   EOS_Internal_Dimensionalize_Extended_Properties(double *dpPropertyIn, double *dpPropertyOut);
    double EOS_Internal_Dimensionalize_Density(double dvDensity);
    double EOS_Internal_Dimensionalize_Pressure(double dvPressure);
    double EOS_Internal_Dimensionalize_Temperature(double dvTemperature);
    double EOS_Internal_Dimensionalize_SpeedSound(double dvSpeedSound);
    double EOS_Internal_Dimensionalize_Entropy(double dvEntropy);
    double EOS_Internal_Dimensionalize_Enthalpy(double dvEnthalpy);
    double EOS_Internal_Dimensionalize_InternalEnergy(double dvInternalEnergy);
    double EOS_Internal_Dimensionalize_TotalEnthalpy(double dvTotalEnthalpy);
    double EOS_Internal_Dimensionalize_TotalEnergy(double dvTotalEnergy);
    double EOS_Internal_Dimensionalize_HeatCapacityCv(double dvHeatCapacityCv);
    double EOS_Internal_Dimensionalize_HeatCapacityCp(double dvHeatCapacityCp);
    double EOS_Internal_Dimensionalize_DPressureDDensity(double dvDPressureDDensity);
    double EOS_Internal_Dimensionalize_DPressureDTemperature(double dvDPressureDTemperature);
    double EOS_Internal_Dimensionalize_DDensityDTemperature(double dvDDensityDTemperature);
    double EOS_Internal_Dimensionalize_D2PressureDDensity2(double dvD2PressureDDensity2);
    double EOS_Internal_Dimensionalize_D2PressureDTemperature2(double dvD2PressureDTemperature2);
    double EOS_Internal_Dimensionalize_D2PressureDTemperatureDensity(double dvD2PressureDTemperatureDensity);
    double EOS_Internal_Dimensionalize_DEnthalpyDTemperature_CDensity(double dvDEnthalpyDTemperature_CDensity);
    double EOS_Internal_Dimensionalize_DEnthalpyDTemperature_CPressure(double dvDEnthalpyDTemperature_CPressure);
    double EOS_Internal_Dimensionalize_DEnthalpyDDensity_CTemperature(double dvDEnthalpyDDensity_CTemperature);
    double EOS_Internal_Dimensionalize_DEnthalpyDDensity_CPressure(double dvDEnthalpyDDensity_CPressure);
    double EOS_Internal_Dimensionalize_DEnthalpyDPressure_CTemperature(double dvDEnthalpyDPressure_CTemperature);
    double EOS_Internal_Dimensionalize_DEnthalpyDPressure_CDensity(double dvDEnthalpyDPressure_CDensity);
    double EOS_Internal_Dimensionalize_DInternalEnergyDTemperature_CDensity(double dvDInternalEnergyDTemperature_CDensity);
    double EOS_Internal_Dimensionalize_DInternalEnergyDTemperature_CPressure(double dvDInternalEnergyDTemperature_CPressure);
    double EOS_Internal_Dimensionalize_DInternalEnergyDDensity_CTemperature(double dvDInternalEnergyDDensity_CTemperature);
    double EOS_Internal_Dimensionalize_DInternalEnergyDDensity_CPressure(double dvDInternalEnergyDDensity_CPressure);
    double EOS_Internal_Dimensionalize_DInternalEnergyDPressure_CTemperature(double dvDInternalEnergyDPressure_CTemperature);
    double EOS_Internal_Dimensionalize_DInternalEnergyDPressure_CDensity(double dvDInternalEnergyDPressure_CDensity);
    
    void   EOS_Internal_NonDimensionalize_Variables(int ivVariableType, double *dpVariableIn, double *dpVariableOut);
    void   EOS_Internal_NonDimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut);
    void   EOS_Internal_NonDimensionalize_Extended_Properties(double *dpPropertyIn, double *dpPropertyOut);
    double EOS_Internal_NonDimensionalize_Density(double dvDensity);
    double EOS_Internal_NonDimensionalize_Pressure(double dvPressure);
    double EOS_Internal_NonDimensionalize_Temperature(double dvTemperature);
    double EOS_Internal_NonDimensionalize_SpeedSound(double dvSpeedSound);
    double EOS_Internal_NonDimensionalize_Entropy(double dvEntropy);
    double EOS_Internal_NonDimensionalize_Enthalpy(double dvEnthalpy);
    double EOS_Internal_NonDimensionalize_InternalEnergy(double dvInternalEnergy);
    double EOS_Internal_NonDimensionalize_TotalEnthalpy(double dvTotalEnthalpy);
    double EOS_Internal_NonDimensionalize_TotalEnergy(double dvTotalEnergy);
    double EOS_Internal_NonDimensionalize_HeatCapacityCv(double dvHeatCapacityCv);
    double EOS_Internal_NonDimensionalize_HeatCapacityCp(double dvHeatCapacityCp);
    double EOS_Internal_NonDimensionalize_DPressureDDensity(double dvDPressureDDensity);
    double EOS_Internal_NonDimensionalize_DPressureDTemperature(double dvDPressureDTemperature);
    double EOS_Internal_NonDimensionalize_DDensityDTemperature(double dvDDensityDTemperature);
    double EOS_Internal_NonDimensionalize_D2PressureDDensity2(double dvD2PressureDDensity2);
    double EOS_Internal_NonDimensionalize_D2PressureDTemperature2(double dvD2PressureDTemperature2);
    double EOS_Internal_NonDimensionalize_D2PressureDTemperatureDensity(double dvD2PressureDTemperatureDensity);
    double EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CDensity(double dvDEnthalpyDTemperature_CDensity);
    double EOS_Internal_NonDimensionalize_DEnthalpyDTemperature_CPressure(double dvDEnthalpyDTemperature_CPressure);
    double EOS_Internal_NonDimensionalize_DEnthalpyDDensity_CTemperature(double dvDEnthalpyDDensity_CTemperature);
    double EOS_Internal_NonDimensionalize_DEnthalpyDDensity_CPressure(double dvDEnthalpyDDensity_CPressure);
    double EOS_Internal_NonDimensionalize_DEnthalpyDPressure_CTemperature(double dvDEnthalpyDPressure_CTemperature);
    double EOS_Internal_NonDimensionalize_DEnthalpyDPressure_CDensity(double dvDEnthalpyDPressure_CDensity);
    double EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CDensity(double dvDInternalEnergyDTemperature_CDensity);
    double EOS_Internal_NonDimensionalize_DInternalEnergyDTemperature_CPressure(double dvDInternalEnergyDTemperature_CPressure);
    double EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CTemperature(double dvDInternalEnergyDDensity_CTemperature);
    double EOS_Internal_NonDimensionalize_DInternalEnergyDDensity_CPressure(double dvDInternalEnergyDDensity_CPressure);
    double EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CTemperature(double dvDInternalEnergyDPressure_CTemperature);
    double EOS_Internal_NonDimensionalize_DInternalEnergyDPressure_CDensity(double dvDInternalEnergyDPressure_CDensity);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_INTERNAL_H */

