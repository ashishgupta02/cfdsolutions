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
    // No of Equations
    #define EOS_NEQUATIONS      5
    #define EOS_THERM_DIM       18
    #define EOS_EX_THERM_DIM    37

    /*!
    * \brief Equation of State Model Type
    */
    enum EOS_MODAL_TYPE {
        EOS_MODEL_NONE     = -1,   /*!< \brief Equation of State Model Type None. */
        EOS_MODEL_IDEALGAS =  0,   /*!< \brief Equation of State Model Type Ideal Gas */
        EOS_MODEL_NIST     =  1    /*!< \brief Equation of State Model Type NIST */
    };
    
    /*!
    * \brief Variable Type
    */
    enum EOS_VARIABLE_TYPE {
        EOS_VARIABLE_NONE = -1,    /*!< \brief Variable Type None. */
        EOS_VARIABLE_CON  =  0,    /*!< \brief Variable Type Conservative. */
        EOS_VARIABLE_RUP  =  1,    /*!< \brief Variable Type Density Velocity and Pressure. */
        EOS_VARIABLE_PUT  =  2,    /*!< \brief Variable Type Presure Velocity and Temperature. */
        EOS_VARIABLE_PUS  =  3,    /*!< \brief Variable Type Pressure Velocity Entropy. */
        EOS_VARIABLE_RUT  =  4     /*!< \brief Variable Type Density Velocity and Temperature. */
    };

    /*!
    * \brief Dimensional Input/Output Type
    */
    enum EOS_DIMENSIONAL_IO_TYPE {
        EOS_DIMENSIONAL_IO_NONE  = -1,   /*!< \brief Dimensional Input/Output Type None. */
        EOS_DIMENSIONAL_IO_ND_ND =  0,   /*!< \brief Dimensional Input/Output Type Input: Non-Dimensional  Output: Non-Dimensional */
        EOS_DIMENSIONAL_IO_D_D   =  1,   /*!< \brief Dimensional Input/Output Type Input:     Dimensional  Output:     Dimensional */
        EOS_DIMENSIONAL_IO_ND_D  =  2,   /*!< \brief Dimensional Input/Output Type Input: Non-Dimensional  Output:     Dimensional */
        EOS_DIMENSIONAL_IO_D_ND  =  3    /*!< \brief Dimensional Input/Output Type Input:     Dimensional  Output: Non-Dimensional */
    };
    
    // Initialize the EOS Data Structure
    void EOS_Init(int ivEOSModelType);
    // Finalize the EOS Data Structure
    void EOS_Finalize(void);
    // Setup the NIST with fluids
    void EOS_Set(const char* csFluidName);
    void EOS_Get_Fluid_Information(double *dpProperty);
    
    // Set the Reference Property for Fluid : Generic (Reference Density is Computed)
    void EOS_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength);
    // Get the Reference Property for Fluid
    void EOS_Get_Reference_Properties(double *dpProperty);
    // Print the Reference Property for Fluid
    void EOS_Print_Reference_Properties();
    
    // Compute EOS Properties
    void EOS_Get_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut);
    void EOS_Get_Extended_Properties(int ivDimIOType, int ivVariableType, double *dpVariableIn, double *dpPropertyOut);
    
    double EOS_Get_Density(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DensityLiquid(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DensityVapor(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_Pressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_Temperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_SpeedSound(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_Mach(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_Entropy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_Enthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_InternalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_TotalEnthalpy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_TotalEnergy(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_HeatCapacityCv(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_HeatCapacityCp(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DPressureDDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DPressureDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DDensityDTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_D2PressureDDensity2(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_D2PressureDTemperature2(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_D2PressureDTemperatureDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DEnthalpyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DEnthalpyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DEnthalpyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DEnthalpyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DEnthalpyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DEnthalpyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DInternalEnergyDTemperature_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DInternalEnergyDTemperature_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DInternalEnergyDDensity_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DInternalEnergyDDensity_CPressure(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DInternalEnergyDPressure_CTemperature(int ivDimIOType, int ivVariableType, double *dpVariableIn);
    double EOS_Get_DInternalEnergyDPressure_CDensity(int ivDimIOType, int ivVariableType, double *dpVariableIn);

    // Variable Transformations
    void EOS_Get_Transformation_Matrix(int ivDimIOType, int ivVarTypeIn, double *dpPropertyIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_H */

