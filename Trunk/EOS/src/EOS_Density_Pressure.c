/*******************************************************************************
 * File:        EOS_Density_Pressure.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "NISTThermo_Extension.h"
#include "EOS.h"
#include "EOS_Internal.h"
#include "EOS_Density_Pressure.h"
#include "Trim_Utils.h"

//------------------- Density and Pressure Formulation -------------------------
// This functions are valid for All Regions:
// 1) Super Critical State
// 2) Subcooled Compressed Liquid
// 3) Super Heated Vapor
// 4) Multi-Phase Liquid-Vapor
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//!
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_DP_Get_Transformation_Matrix_RUP_To_CON(int ivDimIOType, double *dpVariableIn, double **Matrix) {
    int    ivVarType = EOS_VARIABLE_RUP;
    double daProperty[EOS_EX_THERM_DIM];
    double dvU, dvV, dvW, dvRho, dvQ2, dvEt, dvDEDRho_P, dvDEDP_Rho;
    double dvX, dvZ;
    
    // Get the extended properties
    EOS_Get_Extended_Properties(ivDimIOType, ivVarType, dpVariableIn, daProperty);
    
    // Get the Required variables
    dvRho      = daProperty[0];
    dvU        = daProperty[5];
    dvV        = daProperty[6];
    dvW        = daProperty[7];
    dvQ2       = daProperty[8];
    dvEt       = daProperty[15];
    dvDEDRho_P = daProperty[34];
    dvDEDP_Rho = daProperty[36];
    dvX        = dvEt + dvRho*dvDEDRho_P;
    dvZ        = dvRho*dvDEDP_Rho;
            
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = -dvU/dvRho;
    Matrix[1][1] = 1.0/dvRho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = -dvV/dvRho;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = 1.0/dvRho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = -dvW/dvRho;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = 1.0/dvRho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = (dvQ2 - dvX)/dvZ;
    Matrix[4][1] = -dvU/dvZ;
    Matrix[4][2] = -dvV/dvZ;
    Matrix[4][3] = -dvW/dvZ;
    Matrix[4][4] = 1.0/dvZ;
}

//------------------------------------------------------------------------------
//!
//  All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_DP_Get_Transformation_Matrix_CON_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix) {
    int    ivVarType = EOS_VARIABLE_RUP;
    double daProperty[EOS_EX_THERM_DIM];
    double dvU, dvV, dvW, dvRho, dvEt, dvDEDRho_P, dvDEDP_Rho;
    
    // Get the Extended properties
    EOS_Get_Extended_Properties(ivDimIOType, ivVarType, dpVariableIn, daProperty);
    
    // Get the Required variables
    dvRho      = daProperty[0];
    dvU        = daProperty[5];
    dvV        = daProperty[6];
    dvW        = daProperty[7];
    dvEt       = daProperty[15];
    dvDEDRho_P = daProperty[34];
    dvDEDP_Rho = daProperty[36];
    
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = dvU;
    Matrix[1][1] = dvRho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = dvV;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = dvRho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = dvW;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = dvRho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = dvEt + dvRho*dvDEDRho_P;
    Matrix[4][1] = dvRho*dvU;
    Matrix[4][2] = dvRho*dvV;
    Matrix[4][3] = dvRho*dvW;
    Matrix[4][4] = dvRho*dvDEDP_Rho;
}

