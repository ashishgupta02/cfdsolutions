/*******************************************************************************
 * File:        EOS_NIST.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_NIST_H
#define	_EOS_NIST_H

#include "NISTThermo.h"

#ifdef	__cplusplus
extern "C" {
#endif
    
    /*!
    * \brief Thermodynamic Region Type
    */
    enum EOS_REG_TYPE {
        EOS_REG_UNKNOWN                      = 0, /*!< \brief Unknown */
        EOS_REG_SUPER_CRITICAL_STATE         = 1, /*!< \brief Super Critical State */
        EOS_REG_SUBCOOOLED_COMPRESSED_LIQUID = 2, /*!< \brief Subcooled Compressed Liquid */
        EOS_REG_SUPER_HEATED_VAPOR           = 3, /*!< \brief Super Heated Vapor */
        EOS_REG_MULTIPHASE                   = 4  /*!< \brief Multi-Phase Liquid-Vapor */
    };

    /*!
    * \brief NIST Helper Structure
    */
    typedef struct {
        int    ivNumberComponent;
        int    ivError;
        char   csFiles[NISTTHERMO_GEN_STR_LEN*NISTTHERMO_MAX_COMPONENT];
        char   csFmix[NISTTHERMO_GEN_STR_LEN];
        char   csRef[NISTTHERMO_REF_STR_LEN];
        char   csError[NISTTHERMO_GEN_STR_LEN];
        char   csPath[NISTTHERMO_GEN_STR_LEN];
        double daX[NISTTHERMO_MAX_COMPONENT];
        double daXLiq[NISTTHERMO_MAX_COMPONENT];
        double daXVap[NISTTHERMO_MAX_COMPONENT];
    } EOS_S_NISTHelper;

    /* NIST Thermodynamics Input-Output Helpers */
    extern EOS_S_NISTHelper SogNISTHelper;
    
    /* Functions */
    void EOS_NIST_Init();
    void EOS_NIST_Finalize();
    void EOS_NIST_Set(const char* csFluidName);
    void EOS_NIST_Set_Fluid_Information(int ivComp);
    void EOS_NIST_Print_Fluid_Information(int ivComp);
    void EOS_NIST_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength);
    
    double EOS_NIST_Get_SI_To_NIST_Density(double dvDensity);
    double EOS_NIST_Get_NIST_To_SI_Density(double dvDensity);
    double EOS_NIST_Get_SI_To_NIST_Pressure(double dvPressure);
    double EOS_NIST_Get_NIST_To_SI_Pressure(double dvPressure);
    
    void EOS_NIST_Get_NIST_To_SI_Property(double *dpProperty);
    void EOS_NIST_Get_NIST_To_SI_Extended_Property(double *dpProperty);
    void EOS_NIST_Get_SI_To_NIST_Property(double *dpProperty);
    void EOS_NIST_Get_SI_To_NIST_Extended_Property(double *dpProperty);
    
    void EOS_NIST_Get_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut);
    void EOS_NIST_Get_Extended_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut);
    
    
    double EOS_NIST_Get_Density(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DensityLiquid(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DensityVapor(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_Pressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_Temperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_SpeedSound(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_Mach(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_Entropy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_Enthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_InternalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_TotalEnthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_TotalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_HeatCapacityCv(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_HeatCapacityCp(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DPressureDDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DPressureDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DDensityDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_D2PressureDDensity2(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_D2PressureDTemperature2(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_D2PressureDTemperatureDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DEnthalpyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DEnthalpyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DEnthalpyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DEnthalpyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DEnthalpyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DEnthalpyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DInternalEnergyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DInternalEnergyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DInternalEnergyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DInternalEnergyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DInternalEnergyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_NIST_Get_DInternalEnergyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);

#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_NIST_H */

