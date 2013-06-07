/*******************************************************************************
 * File:        Material_Transformations.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Material.h"
#include "SolverParameters.h"
#include "Trim_Utils.h"
#include "EOS.h"
#include "Solver.h"

//------------------------------------------------------------------------------
//! M01
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_CON_To_RUP(int ivNodeID, double **Matrix) {
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    double DPDRho, DPDT, DRhoDT, DEDRho_P; 
    double DEDRho_T, DEDP_Rho, DEDT_Rho;
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[ivNodeID].Get_Extended_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy,
            DPDRho, DPDT, DRhoDT, DEDRho_P,
            DEDRho_T, DEDP_Rho, DEDT_Rho);
    
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = Velocity_U;
    Matrix[1][1] = Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = Velocity_V;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = Velocity_W;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = TotalEnergy + Rho*DEDRho_P;
    Matrix[4][1] = Rho*Velocity_U;
    Matrix[4][2] = Rho*Velocity_V;
    Matrix[4][3] = Rho*Velocity_W;
    Matrix[4][4] = Rho*DEDP_Rho;
}

//------------------------------------------------------------------------------
//! M10
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_RUP_To_CON(int ivNodeID, double **Matrix) {
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    double DPDRho, DPDT, DRhoDT, DEDRho_P; 
    double DEDRho_T, DEDP_Rho, DEDT_Rho;
    double dvX, dvZ;
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[ivNodeID].Get_Extended_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy,
            DPDRho, DPDT, DRhoDT, DEDRho_P,
            DEDRho_T, DEDP_Rho, DEDT_Rho);
    
    dvX = TotalEnergy + Rho*DEDRho_P;
    dvZ = Rho*DEDP_Rho;
    
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = -Velocity_U/Rho;
    Matrix[1][1] = 1.0/Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = -Velocity_V/Rho;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = 1.0/Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = -Velocity_W/Rho;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = 1.0/Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = (Q2 - dvX)/dvZ;
    Matrix[4][1] = -Velocity_U/dvZ;
    Matrix[4][2] = -Velocity_V/dvZ;
    Matrix[4][3] = -Velocity_W/dvZ;
    Matrix[4][4] = 1.0/dvZ;
}

//------------------------------------------------------------------------------
//! M40
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_RUT_To_CON(int ivNodeID, double **Matrix) {
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    double DPDRho, DPDT, DRhoDT, DEDRho_P; 
    double DEDRho_T, DEDP_Rho, DEDT_Rho;
    double dvX, dvZ;
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[ivNodeID].Get_Extended_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy,
            DPDRho, DPDT, DRhoDT, DEDRho_P,
            DEDRho_T, DEDP_Rho, DEDT_Rho);
    
    // Get the Required variables
    dvX = TotalEnergy + Rho*DEDRho_T;
    dvZ = Rho*DEDT_Rho;
    
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = -Velocity_U/Rho;
    Matrix[1][1] = 1.0/Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = -Velocity_V/Rho;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = 1.0/Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = -Velocity_W/Rho;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = 1.0/Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = (Q2 - dvX)/dvZ;
    Matrix[4][1] = -Velocity_U/dvZ;
    Matrix[4][2] = -Velocity_V/dvZ;
    Matrix[4][3] = -Velocity_W/dvZ;
    Matrix[4][4] = 1.0/dvZ;
}

//------------------------------------------------------------------------------
//! M04
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_CON_To_RUT(int ivNodeID, double **Matrix) {
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    double DPDRho, DPDT, DRhoDT, DEDRho_P; 
    double DEDRho_T, DEDP_Rho, DEDT_Rho;
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[ivNodeID].Get_Extended_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy,
            DPDRho, DPDT, DRhoDT, DEDRho_P,
            DEDRho_T, DEDP_Rho, DEDT_Rho);
            
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = Velocity_U;
    Matrix[1][1] = Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = Velocity_V;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = Velocity_W;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = TotalEnergy + Rho*DEDRho_T;
    Matrix[4][1] = Rho*Velocity_U;
    Matrix[4][2] = Rho*Velocity_V;
    Matrix[4][3] = Rho*Velocity_W;
    Matrix[4][4] = Rho*DEDT_Rho;
}

//------------------------------------------------------------------------------
//! M41
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_RUT_To_RUP(int ivNodeID, double **Matrix) {
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    double DPDRho, DPDT, DRhoDT, DEDRho_P; 
    double DEDRho_T, DEDP_Rho, DEDT_Rho;
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[ivNodeID].Get_Extended_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy,
            DPDRho, DPDT, DRhoDT, DEDRho_P,
            DEDRho_T, DEDP_Rho, DEDT_Rho);
    
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

    Matrix[4][0] = 1.0/DRhoDT;
    Matrix[4][1] = 0.0;
    Matrix[4][2] = 0.0;
    Matrix[4][3] = 0.0;
    Matrix[4][4] = 1.0/DPDT;
}

//------------------------------------------------------------------------------
//! M14
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_RUP_To_RUT(int ivNodeID, double **Matrix) {
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    double DPDRho, DPDT, DRhoDT, DEDRho_P; 
    double DEDRho_T, DEDP_Rho, DEDT_Rho;
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[ivNodeID].Get_Extended_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy,
            DPDRho, DPDT, DRhoDT, DEDRho_P,
            DEDRho_T, DEDP_Rho, DEDT_Rho);
    
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

    Matrix[4][0] = DPDRho;
    Matrix[4][1] = 0.0;
    Matrix[4][2] = 0.0;
    Matrix[4][3] = 0.0;
    Matrix[4][4] = DPDT;
}

//------------------------------------------------------------------------------
//! Transformation Matrix: Mpr = dqp/dqr
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix(int ivNodeID, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix) {
    int i, j;
    
    // Check the Variable Type and Return the proper Transformation Matrix
    // Return Identity
    if (ivVarTypeFrom == ivVarTypeTo) {
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                if (i == j)
                    Matrix[i][j] = 1.0;
                else
                    Matrix[i][j] = 0.0;
            }
        }
    } else {  
        switch (ivVarTypeFrom) {
            case VARIABLE_CON:
                switch (ivVarTypeTo) {
                    case VARIABLE_RUP:
                        Material_Get_Transformation_Matrix_CON_To_RUP(ivNodeID, Matrix);
                        break;
                    case VARIABLE_PUT:
                        error("Material_Get_Transformation_Matrix:1: Not Implemented For VARIABLE_PUT");
                        break;
                    case VARIABLE_PUS:
                        error("Material_Get_Transformation_Matrix:2: Not Implemented For VARIABLE_PUS");
                        break;
                    case VARIABLE_RUT:
                        Material_Get_Transformation_Matrix_CON_To_RUT(ivNodeID, Matrix);
                        break;
                    default:
                        error("Material_Get_Transformation_Matrix:3: Undefined Variable Type - %d", ivVarTypeTo);
                }
                break;
            case VARIABLE_RUP:
                switch (ivVarTypeTo) {
                    case VARIABLE_CON:
                        Material_Get_Transformation_Matrix_RUP_To_CON(ivNodeID, Matrix);
                        break;
                    case VARIABLE_PUT:
                        error("Material_Get_Transformation_Matrix:4: Not Implemented For VARIABLE_PUT");
                        break;
                    case VARIABLE_PUS:
                        error("Material_Get_Transformation_Matrix:5: Not Implemented For VARIABLE_PUS");
                        break;
                    case VARIABLE_RUT:
                        Material_Get_Transformation_Matrix_RUP_To_RUT(ivNodeID, Matrix);
                        break;
                    default:
                        error("Material_Get_Transformation_Matrix:6: Undefined Variable Type - %d", ivVarTypeTo);
                }
                break;
            case VARIABLE_PUT:
                error("Material_Get_Transformation_Matrix:7: Not Implemented For VARIABLE_PUT");
                break;
            case VARIABLE_PUS:
                error("Material_Get_Transformation_Matrix:8: Not Implemented For VARIABLE_PUS");
                break;
            case VARIABLE_RUT:
                switch (ivVarTypeTo) {
                    case VARIABLE_CON:
                        Material_Get_Transformation_Matrix_RUT_To_CON(ivNodeID, Matrix);
                        break;
                    case VARIABLE_RUP:
                        Material_Get_Transformation_Matrix_RUT_To_RUP(ivNodeID, Matrix);
                        break;
                    case VARIABLE_PUT:
                        error("Material_Get_Transformation_Matrix:9: Not Implemented For VARIABLE_PUT");
                        break;
                    case VARIABLE_PUS:
                        error("Material_Get_Transformation_Matrix:10: Not Implemented For VARIABLE_PUS");
                        break;
                    default:
                        error("Material_Get_Transformation_Matrix:11: Undefined Variable Type - %d", ivVarTypeTo);
                }
                break;
            default:
                error("Material_Get_Transformation_Matrix:12: Undefined Variable Type - %d", ivVarTypeFrom);
                break;
        }
    }
}

//------------------------------------------------------------------------------
//! Transformation Matrix: Mpr = dqp/dqr
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix(double *dpVariableIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix) {
    double daEOSVariable[NEQUATIONS];
    
    // This is done to preserve the accuracy of dpVariableIn
    // if passed directly to EOS functions will result in loss of accuracy
    for (int i = 0; i < NEQUATIONS; i++)
        daEOSVariable[i] = dpVariableIn[i];
    
    // Compute the EOS Based on Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            // Nothing to Do
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            daEOSVariable[4] += Gauge_Pressure; // Pressure
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            daEOSVariable[0] += Gauge_Pressure; // Pressure
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            // Nothing to Do
            break;
        default:
            error("Material_Get_Transformation_Matrix:1: Undefined Variable Type - %d", VariableType);
            break;
    }
    
    // Compute the Transformation Matrix
    EOS_Get_Transformation_Matrix(MaterialEOS_IO_Type, VariableType, daEOSVariable, ivVarTypeFrom, ivVarTypeTo, Matrix);
}

//------------------------------------------------------------------------------
//! M01
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_Properties_CON_To_RUP(double *dpExPropertyIn, double **Matrix) {
    double Rho, Velocity_U, Velocity_V, Velocity_W, TotalEnergy;
    double DEDRho_P, DEDP_Rho;
    
    // Get the Required Properties values
    Rho         = dpExPropertyIn[ 0];
    Velocity_U  = dpExPropertyIn[ 5];
    Velocity_V  = dpExPropertyIn[ 6];
    Velocity_W  = dpExPropertyIn[ 7];
    TotalEnergy = dpExPropertyIn[15];
    DEDRho_P    = dpExPropertyIn[34];
    DEDP_Rho    = dpExPropertyIn[36];
    
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = Velocity_U;
    Matrix[1][1] = Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = Velocity_V;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = Velocity_W;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = TotalEnergy + Rho*DEDRho_P;
    Matrix[4][1] = Rho*Velocity_U;
    Matrix[4][2] = Rho*Velocity_V;
    Matrix[4][3] = Rho*Velocity_W;
    Matrix[4][4] = Rho*DEDP_Rho;
}

//------------------------------------------------------------------------------
//! M10
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_Properties_RUP_To_CON(double *dpExPropertyIn, double **Matrix) {
    double Rho, Velocity_U, Velocity_V, Velocity_W, Q2, TotalEnergy;
    double DEDRho_P, DEDP_Rho;
    double dvX, dvZ;
    
    // Get the Required Properties values
    Rho           = dpExPropertyIn[ 0];
    Velocity_U    = dpExPropertyIn[ 5];
    Velocity_V    = dpExPropertyIn[ 6];
    Velocity_W    = dpExPropertyIn[ 7];
    Q2            = dpExPropertyIn[ 8];
    TotalEnergy   = dpExPropertyIn[15];
    DEDRho_P      = dpExPropertyIn[34];
    DEDP_Rho      = dpExPropertyIn[36];
    
    dvX = TotalEnergy + Rho*DEDRho_P;
    dvZ = Rho*DEDP_Rho;
    
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = -Velocity_U/Rho;
    Matrix[1][1] = 1.0/Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = -Velocity_V/Rho;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = 1.0/Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = -Velocity_W/Rho;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = 1.0/Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = (Q2 - dvX)/dvZ;
    Matrix[4][1] = -Velocity_U/dvZ;
    Matrix[4][2] = -Velocity_V/dvZ;
    Matrix[4][3] = -Velocity_W/dvZ;
    Matrix[4][4] = 1.0/dvZ;
}

//------------------------------------------------------------------------------
//! M40
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_Properties_RUT_To_CON(double *dpExPropertyIn, double **Matrix) {
    double Rho, Velocity_U, Velocity_V, Velocity_W, Q2, TotalEnergy;
    double DEDRho_T, DEDT_Rho;
    double dvX, dvZ;
    
    // Get the Required Properties values
    Rho         = dpExPropertyIn[ 0];
    Velocity_U  = dpExPropertyIn[ 5];
    Velocity_V  = dpExPropertyIn[ 6];
    Velocity_W  = dpExPropertyIn[ 7];
    Q2          = dpExPropertyIn[ 8];
    TotalEnergy = dpExPropertyIn[15];
    DEDT_Rho    = dpExPropertyIn[31];
    DEDRho_T    = dpExPropertyIn[33];

    // Get the Required variables
    dvX = TotalEnergy + Rho*DEDRho_T;
    dvZ = Rho*DEDT_Rho;
    
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = -Velocity_U/Rho;
    Matrix[1][1] = 1.0/Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = -Velocity_V/Rho;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = 1.0/Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = -Velocity_W/Rho;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = 1.0/Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = (Q2 - dvX)/dvZ;
    Matrix[4][1] = -Velocity_U/dvZ;
    Matrix[4][2] = -Velocity_V/dvZ;
    Matrix[4][3] = -Velocity_W/dvZ;
    Matrix[4][4] = 1.0/dvZ;
}

//------------------------------------------------------------------------------
//! M04
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_Properties_CON_To_RUT(double *dpExPropertyIn, double **Matrix) {
    double Rho, Velocity_U, Velocity_V, Velocity_W, TotalEnergy;
    double DEDRho_T, DEDT_Rho;
    
    // Get the Required Properties values
    Rho         = dpExPropertyIn[ 0];
    Velocity_U  = dpExPropertyIn[ 5];
    Velocity_V  = dpExPropertyIn[ 6];
    Velocity_W  = dpExPropertyIn[ 7];
    TotalEnergy = dpExPropertyIn[15];
    DEDT_Rho    = dpExPropertyIn[31];
    DEDRho_T    = dpExPropertyIn[33];
            
    // Assemble the Matrix
    Matrix[0][0] = 1.0;
    Matrix[0][1] = 0.0;
    Matrix[0][2] = 0.0;
    Matrix[0][3] = 0.0;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = Velocity_U;
    Matrix[1][1] = Rho;
    Matrix[1][2] = 0.0;
    Matrix[1][3] = 0.0;
    Matrix[1][4] = 0.0;

    Matrix[2][0] = Velocity_V;
    Matrix[2][1] = 0.0;
    Matrix[2][2] = Rho;
    Matrix[2][3] = 0.0;
    Matrix[2][4] = 0.0;

    Matrix[3][0] = Velocity_W;
    Matrix[3][1] = 0.0;
    Matrix[3][2] = 0.0;
    Matrix[3][3] = Rho;
    Matrix[3][4] = 0.0;

    Matrix[4][0] = TotalEnergy + Rho*DEDRho_T;
    Matrix[4][1] = Rho*Velocity_U;
    Matrix[4][2] = Rho*Velocity_V;
    Matrix[4][3] = Rho*Velocity_W;
    Matrix[4][4] = Rho*DEDT_Rho;
}

//------------------------------------------------------------------------------
//! M41
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_Properties_RUT_To_RUP(double *dpExPropertyIn, double **Matrix) {
    double DPDT, DRhoDT; 
    
    // Get the Required Properties values
    DPDT   = dpExPropertyIn[19];
    DRhoDT = dpExPropertyIn[20];
    
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

    Matrix[4][0] = 1.0/DRhoDT;
    Matrix[4][1] = 0.0;
    Matrix[4][2] = 0.0;
    Matrix[4][3] = 0.0;
    Matrix[4][4] = 1.0/DPDT;
}

//------------------------------------------------------------------------------
//! M14
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_Properties_RUP_To_RUT(double *dpExPropertyIn, double **Matrix) {
    double DPDRho, DPDT; 
    
    // Get the Required Properties values
    DPDRho = dpExPropertyIn[18];
    DPDT   = dpExPropertyIn[19];
    
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

    Matrix[4][0] = DPDRho;
    Matrix[4][1] = 0.0;
    Matrix[4][2] = 0.0;
    Matrix[4][3] = 0.0;
    Matrix[4][4] = DPDT;
}

//------------------------------------------------------------------------------
//! Transformation Matrix: Mpr = dqp/dqr
//------------------------------------------------------------------------------
void Material_Get_Transformation_Matrix_Properties(double *dpExPropertyIn, int ivVarTypeFrom, int ivVarTypeTo, double **Matrix) {
    int i, j;
    
    // Check the Variable Type and Return the proper Transformation Matrix
    // Return Identity
    if (ivVarTypeFrom == ivVarTypeTo) {
        for (i = 0; i < NEQUATIONS; i++) {
            for (j = 0; j < NEQUATIONS; j++) {
                if (i == j)
                    Matrix[i][j] = 1.0;
                else
                    Matrix[i][j] = 0.0;
            }
        }
    } else {  
        switch (ivVarTypeFrom) {
            case VARIABLE_CON:
                switch (ivVarTypeTo) {
                    case VARIABLE_RUP:
                        Material_Get_Transformation_Matrix_Properties_CON_To_RUP(dpExPropertyIn, Matrix);
                        break;
                    case VARIABLE_PUT:
                        error("Material_Get_Transformation_Matrix_Properties:1: Not Implemented For VARIABLE_PUT");
                        break;
                    case VARIABLE_PUS:
                        error("Material_Get_Transformation_Matrix_Properties:2: Not Implemented For VARIABLE_PUS");
                        break;
                    case VARIABLE_RUT:
                        Material_Get_Transformation_Matrix_Properties_CON_To_RUT(dpExPropertyIn, Matrix);
                        break;
                    default:
                        error("Material_Get_Transformation_Matrix_Properties:3: Undefined Variable Type - %d", ivVarTypeTo);
                }
                break;
            case VARIABLE_RUP:
                switch (ivVarTypeTo) {
                    case VARIABLE_CON:
                        Material_Get_Transformation_Matrix_Properties_RUP_To_CON(dpExPropertyIn, Matrix);
                        break;
                    case VARIABLE_PUT:
                        error("Material_Get_Transformation_Matrix_Properties:4: Not Implemented For VARIABLE_PUT");
                        break;
                    case VARIABLE_PUS:
                        error("Material_Get_Transformation_Matrix_Properties:5: Not Implemented For VARIABLE_PUS");
                        break;
                    case VARIABLE_RUT:
                        Material_Get_Transformation_Matrix_Properties_RUP_To_RUT(dpExPropertyIn, Matrix);
                        break;
                    default:
                        error("Material_Get_Transformation_Matrix_Properties:6: Undefined Variable Type - %d", ivVarTypeTo);
                }
                break;
            case VARIABLE_PUT:
                error("Material_Get_Transformation_Matrix_Properties:7: Not Implemented For VARIABLE_PUT");
                break;
            case VARIABLE_PUS:
                error("Material_Get_Transformation_Matrix_Properties:8: Not Implemented For VARIABLE_PUS");
                break;
            case VARIABLE_RUT:
                switch (ivVarTypeTo) {
                    case VARIABLE_CON:
                        Material_Get_Transformation_Matrix_Properties_RUT_To_CON(dpExPropertyIn, Matrix);
                        break;
                    case VARIABLE_RUP:
                        Material_Get_Transformation_Matrix_Properties_RUT_To_RUP(dpExPropertyIn, Matrix);
                        break;
                    case VARIABLE_PUT:
                        error("Material_Get_Transformation_Matrix_Properties:9: Not Implemented For VARIABLE_PUT");
                        break;
                    case VARIABLE_PUS:
                        error("Material_Get_Transformation_Matrix_Properties:10: Not Implemented For VARIABLE_PUS");
                        break;
                    default:
                        error("Material_Get_Transformation_Matrix_Properties:11: Undefined Variable Type - %d", ivVarTypeTo);
                }
                break;
            default:
                error("Material_Get_Transformation_Matrix_Properties:12: Undefined Variable Type - %d", ivVarTypeFrom);
                break;
        }
    }
}

