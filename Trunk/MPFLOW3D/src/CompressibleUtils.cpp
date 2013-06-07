/*******************************************************************************
 * File:        CompressibleUtils.cpp
 * Author:      Ashish Gupta
 * Revision:    1
 * Package:     Multi-Phase Computational Fluid Dynamics Flow Solver
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Vector3D.h"
#include "Trim_Utils.h"
#include "Material.h"
#include "Solver.h"
#include "SolverParameters.h"
#include "CompressibleUtils.h"
#include "EOS.h"

//------------------------------------------------------------------------------
//! Compute Euler Convective Flux
//  This function should be used if Reevaluation of Properties are required
//------------------------------------------------------------------------------
void Compute_Flux_Euler_Convective(double *Q, Vector3D areavec, double *Flux) {
    double nx, ny, nz;
    double Rho, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Ubar, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    
    // Get area vector
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    // Compute Equation of State
    // Note: Based on VariableType Pressure can be Perturbation of them
    Material_Get_ControlVolume_Properties(Q, Rho, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy);
    Ubar = Velocity_U*nx + Velocity_V*ny + Velocity_W*nz;
    
    // Compute Flux
    Flux[0] = Rho*Ubar;
    Flux[1] = Velocity_U*Flux[0] + Pressure*nx;
    Flux[2] = Velocity_V*Flux[0] + Pressure*ny;
    Flux[3] = Velocity_W*Flux[0] + Pressure*nz;
    Flux[4] = TotalEnthalpy*Flux[0];
}

//------------------------------------------------------------------------------
//! Compute Euler Convective Flux
//! Note: Node's Q is always first order
//------------------------------------------------------------------------------
void Compute_Flux_Euler_Convective(int node_ID, Vector3D areavec, double *Flux) {
    double nx, ny, nz;
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Ubar, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    
    // Get area vector
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[node_ID].Get_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy);
    Ubar = Velocity_U*nx + Velocity_V*ny + Velocity_W*nz;
    
    // Compute Flux
    Flux[0] = Rho*Ubar;
    Flux[1] = Velocity_U*Flux[0] + Pressure*nx;
    Flux[2] = Velocity_V*Flux[0] + Pressure*ny;
    Flux[3] = Velocity_W*Flux[0] + Pressure*nz;
    Flux[4] = TotalEnthalpy*Flux[0];
}

//------------------------------------------------------------------------------
//! Compute Euler Flux Convective Jacobian
//  This function should be used if Reevaluation of Properties are required
//------------------------------------------------------------------------------
void Compute_Flux_Jacobian_Euler_Convective(double *Q, Vector3D areavec, double **Jacobian_Conv) {
    double nx, ny, nz;
    double Ubar, Pressure;
    double Prop[EOS_EX_THERM_DIM];
    
    // Get area vector
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    // Get the Extended Property for the Q's
    // Note: Based on VariableType Pressure can be Perturbation of them
    Material_Get_Extended_Properties(Q, Prop);
    Ubar     = Prop[5]*nx + Prop[6]*ny + Prop[7]*nz;
    Pressure = Prop[3] + Gauge_Pressure;
    
    // Compute the Jacobian Based on Non Dimensionalization and Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            error("Compute_Flux_Jacobian_Euler_Convective:1: Conservative Variable Type Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            // Row:0
            Jacobian_Conv[0][0] = Ubar;
            Jacobian_Conv[0][1] = Q[0]*nx;
            Jacobian_Conv[0][2] = Q[0]*ny;
            Jacobian_Conv[0][3] = Q[0]*nz;
            Jacobian_Conv[0][4] = 0.0;
            // Row:1
            Jacobian_Conv[1][0] = Q[1]*Ubar;
            Jacobian_Conv[1][1] = Q[0]*(Ubar + Q[1]*nx);
            Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
            Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
            Jacobian_Conv[1][4] = nx;
            // Row:2
            Jacobian_Conv[2][0] = Q[2]*Ubar;
            Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
            Jacobian_Conv[2][2] = Q[0]*(Ubar + Q[2]*ny);
            Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
            Jacobian_Conv[2][4] = ny;
            // Row:3
            Jacobian_Conv[3][0] = Q[3]*Ubar;
            Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
            Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
            Jacobian_Conv[3][3] = Q[0]*(Ubar + Q[3]*nz);
            Jacobian_Conv[3][4] = nz;
            // Row:4
            Jacobian_Conv[4][0] = Ubar*(Prop[15] + Q[0]*Prop[34]);
            Jacobian_Conv[4][1] = (Q[0]*Prop[15] + Pressure)*nx + Q[0]*Q[1]*Ubar; // Pressure is not a perturbation
            Jacobian_Conv[4][2] = (Q[0]*Prop[15] + Pressure)*ny + Q[0]*Q[2]*Ubar; // Pressure is not a perturbation
            Jacobian_Conv[4][3] = (Q[0]*Prop[15] + Pressure)*nz + Q[0]*Q[3]*Ubar; // Pressure is not a perturbation
            Jacobian_Conv[4][4] = Ubar*(1.0 + Q[0]*Prop[36]);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            error("Compute_Flux_Jacobian_Euler_Convective:2: PUT Variable Type Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            // Row:0
            Jacobian_Conv[0][0] = Ubar;
            Jacobian_Conv[0][1] = Q[0]*nx;
            Jacobian_Conv[0][2] = Q[0]*ny;
            Jacobian_Conv[0][3] = Q[0]*nz;
            Jacobian_Conv[0][4] = 0.0;
            // Row:1
            Jacobian_Conv[1][0] = Q[1]*Ubar + Prop[18]*nx;
            Jacobian_Conv[1][1] = Q[0]*(Ubar + Q[1]*nx);
            Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
            Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
            Jacobian_Conv[1][4] = Prop[19]*nx;
            // Row:2
            Jacobian_Conv[2][0] = Q[2]*Ubar + Prop[18]*ny;
            Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
            Jacobian_Conv[2][2] = Q[0]*(Ubar + Q[2]*ny);
            Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
            Jacobian_Conv[2][4] = Prop[19]*ny;
            // Row:3
            Jacobian_Conv[3][0] = Q[3]*Ubar + Prop[18]*nz;
            Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
            Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
            Jacobian_Conv[3][3] = Q[0]*(Ubar + Q[3]*nz);
            Jacobian_Conv[3][4] = Prop[19]*nz;
            // Row:4
            Jacobian_Conv[4][0] = Ubar*(Prop[15] + Q[0]*Prop[33] + Prop[18]);
            Jacobian_Conv[4][1] = (Q[0]*Prop[15] + Pressure)*nx + Q[0]*Q[1]*Ubar;
            Jacobian_Conv[4][2] = (Q[0]*Prop[15] + Pressure)*ny + Q[0]*Q[2]*Ubar;
            Jacobian_Conv[4][3] = (Q[0]*Prop[15] + Pressure)*nz + Q[0]*Q[3]*Ubar;
            Jacobian_Conv[4][4] = Ubar*(Q[0]*Prop[31] + Prop[19]);
            break;
        default:
            error("Compute_Flux_Jacobian_Euler_Convective:3: Undefined Variable Type - %d - Error-3", VariableType);
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute Euler Flux Convective Jacobian
//! Note: Node's Q is always first order
//------------------------------------------------------------------------------
void Compute_Flux_Jacobian_Euler_Convective(int node_ID, Vector3D areavec, double **Jacobian_Conv) {
    double nx, ny, nz;
    double Rho, RhoL, RhoV, Pressure, Temperature;
    double Velocity_U, Velocity_V, Velocity_W, Ubar, Q2;
    double SpeedSound, Mach, TotalEnergy, TotalEnthalpy;
    double DPDRho, DPDT, DRhoDT, DEDRho_P; 
    double DEDRho_T, DEDP_Rho, DEDT_Rho;
    
    // Get area vector
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    // Get the Extended Property for the node
    // Note: Based on VariableType Pressure can be Perturbation of them
    CogSolver.CpNodeDB[node_ID].Get_Extended_Properties(Rho, RhoL, RhoV, Pressure, 
            Temperature, Velocity_U, Velocity_V, Velocity_W, Q2,
            SpeedSound, Mach, TotalEnergy, TotalEnthalpy,
            DPDRho, DPDT, DRhoDT, DEDRho_P,
            DEDRho_T, DEDP_Rho, DEDT_Rho);
    Ubar = Velocity_U*nx + Velocity_V*ny + Velocity_W*nz;
    Pressure += Gauge_Pressure;
    
    // Compute the Jacobian Based on Non Dimensionalization and Variable Type of Q
    switch (VariableType) {
        // Conservative Variable Formulation
        case VARIABLE_CON:
            error("Compute_Flux_Jacobian_Euler_Convective:1: Conservative Variable Type Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Pressure
        case VARIABLE_RUP:
            // Row:0
            Jacobian_Conv[0][0] = Ubar;
            Jacobian_Conv[0][1] = Rho*nx;
            Jacobian_Conv[0][2] = Rho*ny;
            Jacobian_Conv[0][3] = Rho*nz;
            Jacobian_Conv[0][4] = 0.0;
            // Row:1
            Jacobian_Conv[1][0] = Velocity_U*Ubar;
            Jacobian_Conv[1][1] = Rho*(Ubar + Velocity_U*nx);
            Jacobian_Conv[1][2] = Rho*Velocity_U*ny;
            Jacobian_Conv[1][3] = Rho*Velocity_U*nz;
            Jacobian_Conv[1][4] = nx;
            // Row:2
            Jacobian_Conv[2][0] = Velocity_V*Ubar;
            Jacobian_Conv[2][1] = Rho*Velocity_V*nx;
            Jacobian_Conv[2][2] = Rho*(Ubar + Velocity_V*ny);
            Jacobian_Conv[2][3] = Rho*Velocity_V*nz;
            Jacobian_Conv[2][4] = ny;
            // Row:3
            Jacobian_Conv[3][0] = Velocity_W*Ubar;
            Jacobian_Conv[3][1] = Rho*Velocity_W*nx;
            Jacobian_Conv[3][2] = Rho*Velocity_W*ny;
            Jacobian_Conv[3][3] = Rho*(Ubar + Velocity_W*nz);
            Jacobian_Conv[3][4] = nz;
            // Row:4
            Jacobian_Conv[4][0] = Ubar*(TotalEnergy + Rho*DEDRho_P);
            Jacobian_Conv[4][1] = (Rho*TotalEnergy + Pressure)*nx + Rho*Velocity_U*Ubar; // Pressure is not a perturbation
            Jacobian_Conv[4][2] = (Rho*TotalEnergy + Pressure)*ny + Rho*Velocity_V*Ubar; // Pressure is not a perturbation
            Jacobian_Conv[4][3] = (Rho*TotalEnergy + Pressure)*nz + Rho*Velocity_W*Ubar; // Pressure is not a perturbation
            Jacobian_Conv[4][4] = Ubar*(1.0 + Rho*DEDP_Rho);
            break;
        // Primitive Variable Formulation Pressure Velocity Temperature
        case VARIABLE_PUT:
            error("Compute_Flux_Jacobian_Euler_Convective:2: PUT Variable Type Not Implemented");
            break;
        // Primitive Variable Formulation Density Velocity Temperature
        case VARIABLE_RUT:
            // Row:0
            Jacobian_Conv[0][0] = Ubar;
            Jacobian_Conv[0][1] = Rho*nx;
            Jacobian_Conv[0][2] = Rho*ny;
            Jacobian_Conv[0][3] = Rho*nz;
            Jacobian_Conv[0][4] = 0.0;
            // Row:1
            Jacobian_Conv[1][0] = Velocity_U*Ubar + DPDRho*nx;
            Jacobian_Conv[1][1] = Rho*(Ubar + Velocity_U*nx);
            Jacobian_Conv[1][2] = Rho*Velocity_U*ny;
            Jacobian_Conv[1][3] = Rho*Velocity_U*nz;
            Jacobian_Conv[1][4] = DPDT*nx;
            // Row:2
            Jacobian_Conv[2][0] = Velocity_V*Ubar + DPDRho*ny;
            Jacobian_Conv[2][1] = Rho*Velocity_V*nx;
            Jacobian_Conv[2][2] = Rho*(Ubar + Velocity_V*ny);
            Jacobian_Conv[2][3] = Rho*Velocity_V*nz;
            Jacobian_Conv[2][4] = DPDT*ny;
            // Row:3
            Jacobian_Conv[3][0] = Velocity_W*Ubar + DPDRho*nz;
            Jacobian_Conv[3][1] = Rho*Velocity_W*nx;
            Jacobian_Conv[3][2] = Rho*Velocity_W*ny;
            Jacobian_Conv[3][3] = Rho*(Ubar + Velocity_W*nz);
            Jacobian_Conv[3][4] = DPDT*nz;
            // Row:4
            Jacobian_Conv[4][0] = Ubar*(TotalEnergy + Rho*DEDRho_T + DPDRho);
            Jacobian_Conv[4][1] = (Rho*TotalEnergy + Pressure)*nx + Rho*Velocity_U*Ubar;
            Jacobian_Conv[4][2] = (Rho*TotalEnergy + Pressure)*ny + Rho*Velocity_V*Ubar;
            Jacobian_Conv[4][3] = (Rho*TotalEnergy + Pressure)*nz + Rho*Velocity_W*Ubar;
            Jacobian_Conv[4][4] = Ubar*(Rho*DEDT_Rho + DPDT);
            break;
        default:
            error("Compute_Flux_Jacobian_Euler_Convective:3: Undefined Variable Type - %d - Error-3", VariableType);
            break;
    }
}

