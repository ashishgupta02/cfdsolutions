/*******************************************************************************
 * File:        EOS_Density_Pressure.h
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#ifndef _EOS_DENSITY_PRESSURE_H
#define	_EOS_DENSITY_PRESSURE_H

#ifdef	__cplusplus
extern "C" {
#endif

    // -------------- Density and Pressure Formulation -------------------------
    // This functions are valid for All Regions:
    // 1) Super Critical State
    // 2) Subcooled Compressed Liquid
    // 3) Super Heated Vapor
    // 4) Multi-Phase Liquid-Vapor
    // -------------------------------------------------------------------------
    // Transformation Matrix
    void EOS_DP_Get_Transformation_Matrix_RUP_To_CON(int ivDimIOType, double *dpVariableIn, double **Matrix);
    void EOS_DP_Get_Transformation_Matrix_CON_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix);
    
#ifdef	__cplusplus
}
#endif

#endif	/* _EOS_DENSITY_PRESSURE_H */

