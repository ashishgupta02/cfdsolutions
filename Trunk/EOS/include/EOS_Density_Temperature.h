/*******************************************************************************
 * File:        EOS_Density_Temperature.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_DENSITY_TEMPERATURE_H
#define	_EOS_DENSITY_TEMPERATURE_H

#ifdef	__cplusplus
extern "C" {
#endif

    // -------------- Density and Temperature Formulation ----------------------
    // This functions are valid for All Regions:
    // 1) Super Critical State
    // 2) Subcooled Compressed Liquid
    // 3) Super Heated Vapor
    // 4) Multi-Phase Liquid-Vapor
    // -------------------------------------------------------------------------   
    // Transformation Matrix
    void EOS_DT_Get_Transformation_Matrix_CON_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_RUT_To_CON(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_CON_To_RUT(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_RUT_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DT_Get_Transformation_Matrix_RUP_To_RUT(int ivDimIOType, double *dpVariableIn, double **Matrix);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_DENSITY_TEMPERATURE_H */

