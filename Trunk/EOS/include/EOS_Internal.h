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
    * \brief Fluid Property Structure
    */
    typedef struct {
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
    
    /* Common Shared Data */
    /* NIST Thermodynamics Input-Output Helpers */
    extern EOS_S_NISTHelper SogNISTHelper;
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
    void EOS_Internal_Init_Fluid_Information();
    void EOS_Internal_Set_NIST_Fluid_Information(int ivComp);
    void EOS_Internal_Print_NIST_Fluid_Information(int ivComp);
    void EOS_Internal_Init_Reference_Properties();
    
    double EOS_Internal_Get_SI_To_NIST_Density(double dvDensity);
    double EOS_Internal_Get_NIST_To_SI_Density(double dvDensity);
    double EOS_Internal_Get_SI_To_NIST_Pressure(double dvPressure);
    double EOS_Internal_Get_NIST_To_SI_Pressure(double dvPressure);
    
    void   EOS_Internal_Get_NIST_To_SI_Property(double *dpProperty);
    void   EOS_Internal_Get_NIST_To_SI_Extended_Property(double *dpProperty);
    void   EOS_Internal_Get_SI_To_NIST_Property(double *dpProperty);
    void   EOS_Internal_Get_SI_To_NIST_Extended_Property(double *dpProperty);
    
    void EOS_Internal_Dimensionalize_Variables(int ivVariableType, double *dpVariableIn, double *dpVariableOut);
    void EOS_Internal_NonDimensionalize_Variables(int ivVariableType, double *dpVariableIn, double *dpVariableOut);
    
    void EOS_Internal_Dimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut);
    void EOS_Internal_NonDimensionalize_Properties(double *dpPropertyIn, double *dpPropertyOut);
    void EOS_Internal_Dimensionalize_Extended_Properties(double *dpPropertyIn, double *dpPropertyOut);
    void EOS_Internal_NonDimensionalize_Extended_Properties(double *dpPropertyIn, double *dpPropertyOut);
    
    void EOS_Internal_Dimensionalize_DT(double *dpDensity, double *dpTemperature);
    void EOS_Internal_Dimensionalize_DP(double *dpDensity, double *dpPressure);
    void EOS_Internal_Dimensionalize_PT(double *dpPressure, double *dpTemperature);
    
    void EOS_Internal_NonDimensionalize_Pressure(double *dpPressure);
    void EOS_Internal_NonDimensionalize_Density(double *dpDensity);
    void EOS_Internal_NonDimensionalize_Temperature(double *dpTemperature);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_INTERNAL_H */

