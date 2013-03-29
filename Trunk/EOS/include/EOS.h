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
    void EOS_Init(void);
    // Finalize the EOS Data Structure
    void EOS_Finalize(void);
    // Setup the NIST with fluids
    void EOS_Set(void);
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
    
    // Variable Transformations
    void EOS_Get_Transformation_Matrix(int ivDimIOType, int ivVarTypeIn, double *dpPropertyIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_H */

