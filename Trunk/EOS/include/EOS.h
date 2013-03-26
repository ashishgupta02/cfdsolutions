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
    
    // Set the Reference Property for Fluid : Generic (Reference Density is Computed)
    void EOS_Set_Reference_Properties(double dvPressure, double dvTemperature, double dvLength);
    // Get the Reference Property for Fluid
    void EOS_Get_Reference_Properties(double *dvProperty_Ref);
    // Print the Reference Property for Fluid
    void EOS_Print_Reference_Properties();
    
    // Compute EOS Properties
    void EOS_Get_Properties(int ivVariableType, double *dpVariableIn, double *dpPropertyOut);
    
    // Variable Transformations
    void EOS_Get_Transformation_Matrix(int ivVarTypeIn, double *dpPropertyIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_H */

