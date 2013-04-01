/*******************************************************************************
 * File:        EOS_Density_Temperature.c
 * Author:      Ashish Gupta
 * Revision:    1
 ******************************************************************************/

#include <string.h>

#include "NISTThermo.h"
#include "Trim_Utils.h"
#include "EOS.h"
#include "EOS_NIST.h"
#include "EOS_Internal.h"
#include "EOS_Density_Temperature.h"

// ------------------ Density and Temperature Formulation ----------------------
// This functions are valid for All Regions:
// 1) Super Critical State
// 2) Subcooled Compressed Liquid
// 3) Super Heated Vapor
// 4) Multi-Phase Liquid-Vapor
// -----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//!
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_DT_Get_Transformation_Matrix_CON_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix) {
    int    ivVarType = EOS_VARIABLE_RUT;
    double daProperty[EOS_EX_THERM_DIM];
    double dvU, dvV, dvW, dvRho, dvEt, dvDEDRho_P, dvDEDP_Rho;
    
    // Get the extended properties
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

//------------------------------------------------------------------------------
//!
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_DT_Get_Transformation_Matrix_RUT_To_CON(int ivDimIOType, double *dpVariableIn, double **Matrix) {
    int    ivVarType = EOS_VARIABLE_RUT;
    double daProperty[EOS_EX_THERM_DIM];
    double dvU, dvV, dvW, dvRho, dvQ2, dvEt, dvDEDRho_T, dvDEDT_Rho;
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
    dvDEDRho_T = daProperty[33];
    dvDEDT_Rho = daProperty[31];
    dvX        = dvEt + dvRho*dvDEDRho_T;
    dvZ        = dvRho*dvDEDT_Rho;
            
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
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_DT_Get_Transformation_Matrix_CON_To_RUT(int ivDimIOType, double *dpVariableIn, double **Matrix) {
    int    ivVarType = EOS_VARIABLE_RUT;
    double daProperty[EOS_EX_THERM_DIM];
    double dvU, dvV, dvW, dvRho, dvEt, dvDEDRho_T, dvDEDT_Rho;
    
    // Get the extended properties
    EOS_Get_Extended_Properties(ivDimIOType, ivVarType, dpVariableIn, daProperty);
    
    // Get the Required variables
    dvRho      = daProperty[0];
    dvU        = daProperty[5];
    dvV        = daProperty[6];
    dvW        = daProperty[7];
    dvEt       = daProperty[15];
    dvDEDRho_T = daProperty[33];
    dvDEDT_Rho = daProperty[31];
            
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

    Matrix[4][0] = dvEt + dvRho*dvDEDRho_T;
    Matrix[4][1] = dvRho*dvU;
    Matrix[4][2] = dvRho*dvV;
    Matrix[4][3] = dvRho*dvW;
    Matrix[4][4] = dvRho*dvDEDT_Rho;
}

//------------------------------------------------------------------------------
//!
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_DT_Get_Transformation_Matrix_RUT_To_RUP(int ivDimIOType, double *dpVariableIn, double **Matrix) {
    int    ivVarType = EOS_VARIABLE_RUT;
    double daProperty[EOS_EX_THERM_DIM];
    double dvDTDRho, dvDTDP;
    
    // Get the extended properties
    EOS_Get_Extended_Properties(ivDimIOType, ivVarType, dpVariableIn, daProperty);
    
    // Get the Required variables
    dvDTDRho = 1.0/daProperty[20];
    dvDTDP   = 1.0/daProperty[19];
            
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = 0.0;
    Matrix[1][1] = 1.0;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = 0.0;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = 1.0;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = 0.0;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = 1.0;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = dvDTDRho;
    Matrix[4][1] = 0.0;
    Matrix[4][2] = 0.0;
    Matrix[4][3] = 0.0;
    Matrix[4][4] = dvDTDP;
}

//------------------------------------------------------------------------------
//!
//! All Matrix Output are based on Dimensional Input/Output Type
//------------------------------------------------------------------------------
void EOS_DT_Get_Transformation_Matrix_RUP_To_RUT(int ivDimIOType, double *dpVariableIn, double **Matrix) {
    int    ivVarType = EOS_VARIABLE_RUT;
    double daProperty[EOS_EX_THERM_DIM];
    double dvDPDRho, dvDPDT;
    
    // Get the extended properties
    EOS_Get_Extended_Properties(ivDimIOType, ivVarType, dpVariableIn, daProperty);
    
    // Get the Required variables
    dvDPDRho = 1.0/daProperty[21];
    dvDPDT   = daProperty[19];
            
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = 0.0;
    Matrix[1][1] = 1.0;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = 0.0;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = 1.0;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = 0.0;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = 1.0;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = dvDPDRho;
    Matrix[4][1] = 0.0;
    Matrix[4][2] = 0.0;
    Matrix[4][3] = 0.0;
    Matrix[4][4] = dvDPDT;
}

