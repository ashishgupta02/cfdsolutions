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
    * \brief Variable Type
    */
    enum EOS_VARIABLE_TYPE {
        EOS_VARIABLE_NONE = -1,    /*!< \brief Variable Type None. */
        EOS_VARIABLE_CON  =  0,    /*!< \brief Variable Type Conservative. */
        EOS_VARIABLE_PUT  =  1,    /*!< \brief Variable Type Presure Velocity and Temperature. */
        EOS_VARIABLE_RUP  =  2,    /*!< \brief Variable Type Density Velocity and Pressure. */
        EOS_VARIABLE_PUS  =  3,    /*!< \brief Variable Type Pressure Velocity and Entropy. */
        EOS_VARIABLE_RUT  =  4     /*!< \brief Variable Type Density Velocity and Temperature. */
    };

    /*!
    * \brief Thermodynamic Region Type
    */
    enum EOS_REG_TYPE {
        EOS_REG_UKN     = 0, /*!< \brief Unknown */
        EOS_REG_SCS     = 1, /*!< \brief Super Critical State */
        EOS_REG_SCL     = 2, /*!< \brief Subcooled Compressed Liquid */
        EOS_REG_SHV     = 3, /*!< \brief Super Heated Vapor */
        EOS_REG_MP      = 4  /*!< \brief Multi-Phase Liquid-Vapor */
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
    extern double dvTemperature_Ref;
    extern double dvPressure_Ref;
    extern double dvDensity_Ref;
    extern double dvVelocity_Ref;
    extern double dvLength_Ref;
    extern double dvMach_Ref;
    extern double dvSpeedSound_Ref;
    extern double dvTime_Ref;
    extern double dvEnthalpy_Ref;
    extern double dvTotalEnthalpy_Ref;
    extern double dvInternalEnergy_Ref;
    extern double dvTotalEnergy_Ref;
    extern double dvEntropy_Ref;
    extern double dvEntropyConst_Ref;
    extern double dvGasConstant_Ref;
    extern double dvRatioSpecificHeat_Ref;
    
    /* Functions */
    void EOS_Internal_Init();
    void EOS_Internal_Init_Fluid_Information();
    void EOS_Internal_Set_NIST_Fluid_Information(int ivComp);
    void EOS_Internal_Print_NIST_Fluid_Information(int ivComp);
    void EOS_Internal_Init_Reference_Properties();
    void EOS_Internal_Dimensionalize_Properties(int ivVariableType, double *dpPropertyIn, double *dpPropertyOut);
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
