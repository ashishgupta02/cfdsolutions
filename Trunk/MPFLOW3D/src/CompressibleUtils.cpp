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

//------------------------------------------------------------------------------
//! Compute Euler Convective Flux
//------------------------------------------------------------------------------
void Compute_Flux_Euler_Convective(double *Q, Vector3D areavec, double *Flux) {
    double nx, ny, nz;
    double rho, u, v, w, et, p, c, ht, ubar, q2, T, mach;
    
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    // Compute Equation of State
    // Note: Based on VariableType Pressure and Temperature can be Perturbation of them
    Compute_EOS_Variables_Face(Q, nx, ny, nz, rho, p, T, u, v, w, q2, c, mach, ubar, et, ht);
    
    // Compute Flux
    Flux[0] = rho*ubar;
    Flux[1] =   u*Flux[0] + p*nx;
    Flux[2] =   v*Flux[0] + p*ny;
    Flux[3] =   w*Flux[0] + p*nz;
    Flux[4] =  ht*Flux[0];
}

//------------------------------------------------------------------------------
//! Compute Euler Convective Flux
//! Note: Q is always first order
//------------------------------------------------------------------------------
void Compute_Flux_Euler_Convective(int node_ID, Vector3D areavec, double *Flux) {
    double Q[NEQUATIONS];
    
    // Get the Variables
    Q[0] = Q1[node_ID];
    Q[1] = Q2[node_ID];
    Q[2] = Q3[node_ID];
    Q[3] = Q4[node_ID];
    Q[4] = Q5[node_ID];
    
    // Compute the convective flux
    Compute_Flux_Euler_Convective(Q, areavec, Flux);
}

//------------------------------------------------------------------------------
//! Compute Euler Flux Convective Jacobian
//! Note: Q is always first order
//------------------------------------------------------------------------------
void Compute_Flux_Jacobian_Euler_Convective(double *Q, Vector3D areavec, double **Jacobian_Conv) {
    double nx, ny, nz;
    
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    // Based on variable type update the perturbation variables
    switch (VariableType) {
        case VARIABLE_CONSERVATIVE:
            // Do Nothing
            break;
        case VARIABLE_PRIMITIVE_PUT:
            Q[0] += Gauge_Pressure;
            break;
        case VARIABLE_PRIMITIVE_RUP:
            Q[4] += Gauge_Pressure;
            break;
        case VARIABLE_PRIMITIVE_RUT:
            break;
    }
    
    // Compute the Jacobian Based on Non Dimensionalization and Variable Type of Q
    switch (NonDimensionalMethod) {
        // Default Non Dimensionalization
        case NONDIMENSIONAL_METHOD_GENERIC:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    // Row:0
                    Jacobian_Conv[0][0] = 0.0;
                    Jacobian_Conv[0][1] = nx;
                    Jacobian_Conv[0][2] = ny;
                    Jacobian_Conv[0][3] = nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = (-2.0*Q[1]*(Q[2]*ny + Q[3]*nz) + nx*((Gamma - 3.0)*Q[1]*Q[1] + (Gamma - 1.0)*(Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[1][1] = (-(Gamma - 3.0)*Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[0];
                    Jacobian_Conv[1][2] = (Q[1]*ny - (Gamma - 1.0)*Q[2]*nx)/Q[0];
                    Jacobian_Conv[1][3] = (Q[1]*nz - (Gamma - 1.0)*Q[3]*nx)/Q[0];
                    Jacobian_Conv[1][4] = (Gamma - 1.0)*nx;
                    // Row:2
                    Jacobian_Conv[2][0] = (-2.0*Q[2]*(Q[1]*nx + Q[3]*nz) + ny*((Gamma - 3.0)*Q[2]*Q[2] + (Gamma - 1.0)*(Q[1]*Q[1] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[2][1] = (Q[2]*nx - (Gamma - 1.0)*Q[1]*ny)/Q[0];
                    Jacobian_Conv[2][2] = (Q[1]*nx - (Gamma - 3.0)*Q[2]*ny + Q[3]*nz)/Q[0];
                    Jacobian_Conv[2][3] = (Q[2]*nz - (Gamma - 1.0)*Q[3]*ny)/Q[0];
                    Jacobian_Conv[2][4] = (Gamma - 1.0)*ny;
                    // Row:3
                    Jacobian_Conv[3][0] = (-2.0*Q[3]*(Q[1]*nx + Q[2]*ny) + nz*((Gamma - 3.0)*Q[3]*Q[3] + (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[3][1] = (Q[3]*nx - (Gamma - 1.0)*Q[1]*nz)/Q[0];
                    Jacobian_Conv[3][2] = (Q[3]*ny - (Gamma - 1.0)*Q[2]*nz)/Q[0];
                    Jacobian_Conv[3][3] = (Q[1]*nx + Q[2]*ny - (Gamma - 3.0)*Q[3]*nz)/Q[0];
                    Jacobian_Conv[3][4] = (Gamma - 1.0)*nz;
                    // Row:4
                    Jacobian_Conv[4][0] = (Q[1]*nx + Q[2]*ny + Q[3]*nz)*((Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) - Gamma*Q[0]*Q[4])/(Q[0]*Q[0]*Q[0]);
                    Jacobian_Conv[4][1] = (-2.0*(Gamma - 1.0)*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nx*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][2] = (-2.0*(Gamma - 1.0)*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + ny*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][3] = (-2.0*(Gamma - 1.0)*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nz*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][4] = Gamma*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[0];
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    // Row:0
                    Jacobian_Conv[0][0] = Gamma*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[0][1] = Gamma*Q[0]*nx/Q[4];
                    Jacobian_Conv[0][2] = Gamma*Q[0]*ny/Q[4];
                    Jacobian_Conv[0][3] = Gamma*Q[0]*nz/Q[4];
                    Jacobian_Conv[0][4] = -Gamma*Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:1
                    Jacobian_Conv[1][0] = nx + Gamma*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[1][1] = Gamma*Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[1][2] = Gamma*Q[0]*Q[1]*ny/Q[4];
                    Jacobian_Conv[1][3] = Gamma*Q[0]*Q[1]*nz/Q[4];
                    Jacobian_Conv[1][4] = -Gamma*Q[0]*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:2
                    Jacobian_Conv[2][0] = ny + Gamma*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[2][1] = Gamma*Q[0]*Q[2]*nx/Q[4];
                    Jacobian_Conv[2][2] = Gamma*Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[2][3] = Gamma*Q[0]*Q[2]*nz/Q[4];
                    Jacobian_Conv[2][4] = -Gamma*Q[0]*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:3
                    Jacobian_Conv[3][0] = nz + Gamma*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[3][1] = Gamma*Q[0]*Q[3]*nx/Q[4];
                    Jacobian_Conv[3][2] = Gamma*Q[0]*Q[3]*ny/Q[4];
                    Jacobian_Conv[3][3] = Gamma*Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz)/Q[4];
                    Jacobian_Conv[3][4] = -Gamma*Q[0]*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:4
                    Jacobian_Conv[4][0] = Gamma*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma - 1.0))/(2.0*Q[4]);
                    Jacobian_Conv[4][1] = (Gamma*Q[0]*(Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nx*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma -1.0))))/Q[4];
                    Jacobian_Conv[4][2] = (Gamma*Q[0]*(Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*ny*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma -1.0))))/Q[4];
                    Jacobian_Conv[4][3] = (Gamma*Q[0]*(Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nz*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma -1.0))))/Q[4];
                    Jacobian_Conv[4][4] = -Gamma*Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/(2.0*Q[4]*Q[4]);
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    // Row:0
                    Jacobian_Conv[0][0] = Q[1]*nx + Q[2]*ny + Q[3]*nz;
                    Jacobian_Conv[0][1] = Q[0]*nx;
                    Jacobian_Conv[0][2] = Q[0]*ny;
                    Jacobian_Conv[0][3] = Q[0]*nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][1] = Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
                    Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
                    Jacobian_Conv[1][4] = nx;
                    // Row:2
                    Jacobian_Conv[2][0] = Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
                    Jacobian_Conv[2][2] = Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
                    Jacobian_Conv[2][4] = ny;
                    // Row:3
                    Jacobian_Conv[3][0] = Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);;
                    Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
                    Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
                    Jacobian_Conv[3][3] = Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz);;
                    Jacobian_Conv[3][4] = nz;
                    // Row:4
                    Jacobian_Conv[4][0] = 0.5*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]);
                    Jacobian_Conv[4][1] = (2.0*(Gamma - 1.0)*Q[0]*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nx*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]))/(2.0*(Gamma - 1.0));
                    Jacobian_Conv[4][2] = (2.0*(Gamma - 1.0)*Q[0]*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + ny*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]))/(2.0*(Gamma - 1.0));
                    Jacobian_Conv[4][3] = (2.0*(Gamma - 1.0)*Q[0]*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nz*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]))/(2.0*(Gamma - 1.0));
                    Jacobian_Conv[4][4] = Gamma*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Gamma - 1.0);
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    // Row:0
                    Jacobian_Conv[0][0] = Q[1]*nx + Q[2]*ny + Q[3]*nz;
                    Jacobian_Conv[0][1] = Q[0]*nx;
                    Jacobian_Conv[0][2] = Q[0]*ny;
                    Jacobian_Conv[0][3] = Q[0]*nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*nx/Gamma;
                    Jacobian_Conv[1][1] = Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
                    Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
                    Jacobian_Conv[1][4] = Q[0]*nx/Gamma;
                    // Row:2
                    Jacobian_Conv[2][0] = Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*ny/Gamma;
                    Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
                    Jacobian_Conv[2][2] = Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
                    Jacobian_Conv[2][4] = Q[0]*ny/Gamma;
                    // Row:3
                    Jacobian_Conv[3][0] = Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*nz/Gamma;
                    Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
                    Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
                    Jacobian_Conv[3][3] = Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz);
                    Jacobian_Conv[3][4] = Q[0]*nz/Gamma;
                    // Row:4
                    Jacobian_Conv[4][0] = 0.5*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma - 1.0));
                    Jacobian_Conv[4][1] = Q[0]*(Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nx*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma - 1.0)));
                    Jacobian_Conv[4][2] = Q[0]*(Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*ny*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma - 1.0)));
                    Jacobian_Conv[4][3] = Q[0]*(Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nz*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]/(Gamma - 1.0)));
                    Jacobian_Conv[4][4] = Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Gamma - 1.0);
                    break;
                default:
                    error("Compute_Flux_Jacobian_Euler_Convective: Undefined Variable Type - %d - Error-1", VariableType);
                    break;
            }
            break;
        // Briley Taylor Whitfield Non Dimensionalization
        case NONDIMENSIONAL_METHOD_BTW:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    // Row:0
                    Jacobian_Conv[0][0] = 0.0;
                    Jacobian_Conv[0][1] = nx;
                    Jacobian_Conv[0][2] = ny;
                    Jacobian_Conv[0][3] = nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = (-2.0*Q[1]*(Q[2]*ny + Q[3]*nz) + nx*((Gamma - 3.0)*Q[1]*Q[1] + (Gamma - 1.0)*(Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[1][1] = (-(Gamma - 3.0)*Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[0];
                    Jacobian_Conv[1][2] = (Q[1]*ny - (Gamma - 1.0)*Q[2]*nx)/Q[0];
                    Jacobian_Conv[1][3] = (Q[1]*nz - (Gamma - 1.0)*Q[3]*nx)/Q[0];
                    Jacobian_Conv[1][4] = nx/(Ref_Mach*Ref_Mach);
                    // Row:2
                    Jacobian_Conv[2][0] = (-2.0*Q[2]*(Q[1]*nx + Q[3]*nz) + ny*((Gamma - 3.0)*Q[2]*Q[2] + (Gamma - 1.0)*(Q[1]*Q[1] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[2][1] = (Q[2]*nx - (Gamma - 1.0)*Q[1]*ny)/Q[0];
                    Jacobian_Conv[2][2] = (Q[1]*nx - (Gamma - 3.0)*Q[2]*ny + Q[3]*nz)/Q[0];
                    Jacobian_Conv[2][3] = (Q[2]*nz - (Gamma - 1.0)*Q[3]*ny)/Q[0];
                    Jacobian_Conv[2][4] = ny/(Ref_Mach*Ref_Mach);
                    // Row:3
                    Jacobian_Conv[3][0] = (-2.0*Q[3]*(Q[1]*nx + Q[2]*ny) + nz*((Gamma - 3.0)*Q[3]*Q[3] + (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[3][1] = (Q[3]*nx - (Gamma - 1.0)*Q[1]*nz)/Q[0];
                    Jacobian_Conv[3][2] = (Q[3]*ny - (Gamma - 1.0)*Q[2]*nz)/Q[0];
                    Jacobian_Conv[3][3] = (Q[1]*nx + Q[2]*ny - (Gamma - 3.0)*Q[3]*nz)/Q[0];
                    Jacobian_Conv[3][4] = nz/(Ref_Mach*Ref_Mach);
                    // Row:4
                    Jacobian_Conv[4][0] = (Q[1]*nx + Q[2]*ny + Q[3]*nz)*((Gamma - 1.0)*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) - Gamma*Q[0]*Q[4])/(Q[0]*Q[0]*Q[0]);
                    Jacobian_Conv[4][1] = (-2.0*(Gamma - 1.0)*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nx*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][2] = (-2.0*(Gamma - 1.0)*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + ny*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][3] = (-2.0*(Gamma - 1.0)*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nz*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][4] = Gamma*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[0];
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    // Row:0
                    Jacobian_Conv[0][0] = Gamma*Ref_Mach*Ref_Mach*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[0][1] = Gamma*Ref_Mach*Ref_Mach*Q[0]*nx/Q[4];
                    Jacobian_Conv[0][2] = Gamma*Ref_Mach*Ref_Mach*Q[0]*ny/Q[4];
                    Jacobian_Conv[0][3] = Gamma*Ref_Mach*Ref_Mach*Q[0]*nz/Q[4];
                    Jacobian_Conv[0][4] = -Gamma*Ref_Mach*Ref_Mach*Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:1
                    Jacobian_Conv[1][0] = nx + Gamma*Ref_Mach*Ref_Mach*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[1][1] = Gamma*Ref_Mach*Ref_Mach*Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[1][2] = Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[1]*ny/Q[4];
                    Jacobian_Conv[1][3] = Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[1]*nz/Q[4];
                    Jacobian_Conv[1][4] = -Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:2
                    Jacobian_Conv[2][0] = ny + Gamma*Ref_Mach*Ref_Mach*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[2][1] = Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[2]*nx/Q[4];
                    Jacobian_Conv[2][2] = Gamma*Ref_Mach*Ref_Mach*Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[2][3] = Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[2]*nz/Q[4];
                    Jacobian_Conv[2][4] = -Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:3
                    Jacobian_Conv[3][0] = nz + Gamma*Ref_Mach*Ref_Mach*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[3][1] = Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[3]*nx/Q[4];
                    Jacobian_Conv[3][2] = Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[3]*ny/Q[4];
                    Jacobian_Conv[3][3] = Gamma*Ref_Mach*Ref_Mach*Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz)/Q[4];
                    Jacobian_Conv[3][4] = -Gamma*Ref_Mach*Ref_Mach*Q[0]*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:4
                    Jacobian_Conv[4][0] = Gamma*Ref_Mach*Ref_Mach*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*((Gamma - 1.0)*Ref_Mach*Ref_Mach*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Q[4])/(2.0*Q[4]);
                    Jacobian_Conv[4][1] = Gamma*Ref_Mach*Ref_Mach*Q[0]*((Gamma - 1.0)*Ref_Mach*Ref_Mach*(Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nx*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])) + Q[4]*nx)/Q[4];
                    Jacobian_Conv[4][2] = Gamma*Ref_Mach*Ref_Mach*Q[0]*((Gamma - 1.0)*Ref_Mach*Ref_Mach*(Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*ny*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])) + Q[4]*ny)/Q[4];
                    Jacobian_Conv[4][3] = Gamma*Ref_Mach*Ref_Mach*Q[0]*((Gamma - 1.0)*Ref_Mach*Ref_Mach*(Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nz*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])) + Q[4]*nz)/Q[4];
                    Jacobian_Conv[4][4] = -Gamma*(Gamma - 1.0)*Ref_Mach*Ref_Mach*Ref_Mach*Ref_Mach*Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/(2.0*Q[4]*Q[4]);
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    // Row:0
                    Jacobian_Conv[0][0] = Q[1]*nx + Q[2]*ny + Q[3]*nz;
                    Jacobian_Conv[0][1] = Q[0]*nx;
                    Jacobian_Conv[0][2] = Q[0]*ny;
                    Jacobian_Conv[0][3] = Q[0]*nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][1] = Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
                    Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
                    Jacobian_Conv[1][4] = nx;
                    // Row:2
                    Jacobian_Conv[2][0] = Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
                    Jacobian_Conv[2][2] = Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
                    Jacobian_Conv[2][4] = ny;
                    // Row:3
                    Jacobian_Conv[3][0] = Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);;
                    Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
                    Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
                    Jacobian_Conv[3][3] = Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz);;
                    Jacobian_Conv[3][4] = nz;
                    // Row:4
                    Jacobian_Conv[4][0] = 0.5*(Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]);
                    Jacobian_Conv[4][1] = 0.5*Ref_Mach*Ref_Mach*(2.0*(Gamma - 1.0)*Q[0]*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nx*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]));
                    Jacobian_Conv[4][2] = 0.5*Ref_Mach*Ref_Mach*(2.0*(Gamma - 1.0)*Q[0]*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + ny*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]));
                    Jacobian_Conv[4][3] = 0.5*Ref_Mach*Ref_Mach*(2.0*(Gamma - 1.0)*Q[0]*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nz*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]));
                    Jacobian_Conv[4][4] = Gamma*Ref_Mach*Ref_Mach*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    // Row:0
                    Jacobian_Conv[0][0] = Q[1]*nx + Q[2]*ny + Q[3]*nz;
                    Jacobian_Conv[0][1] = Q[0]*nx;
                    Jacobian_Conv[0][2] = Q[0]*ny;
                    Jacobian_Conv[0][3] = Q[0]*nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*nx/(Gamma*Ref_Mach*Ref_Mach);
                    Jacobian_Conv[1][1] = Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
                    Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
                    Jacobian_Conv[1][4] = Q[0]*nx/(Gamma*Ref_Mach*Ref_Mach);
                    // Row:2
                    Jacobian_Conv[2][0] = Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*ny/(Gamma*Ref_Mach*Ref_Mach);
                    Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
                    Jacobian_Conv[2][2] = Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
                    Jacobian_Conv[2][4] = Q[0]*ny/(Gamma*Ref_Mach*Ref_Mach);
                    // Row:3
                    Jacobian_Conv[3][0] = Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*nz/(Gamma*Ref_Mach*Ref_Mach);
                    Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
                    Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
                    Jacobian_Conv[3][3] = Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz);
                    Jacobian_Conv[3][4] = Q[0]*nz/(Gamma*Ref_Mach*Ref_Mach);
                    // Row:4
                    Jacobian_Conv[4][0] = 0.5*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*((Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Q[4]);
                    Jacobian_Conv[4][1] = Q[0]*((Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nx*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])) + Q[4]*nx);
                    Jacobian_Conv[4][2] = Q[0]*((Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*ny*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])) + Q[4]*ny);
                    Jacobian_Conv[4][3] = Q[0]*((Gamma - 1.0)*(Ref_Mach*Ref_Mach)*(Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nz*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])) + Q[4]*nz);
                    Jacobian_Conv[4][4] = Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    break;
                default:
                    error("Compute_Flux_Jacobian_Euler_Convective: Undefined Variable Type - %d - Error-2", VariableType);
                    break;
            }
            break;
        // LMRoe Non Dimensionalization
        case NONDIMENSIONAL_METHOD_LMROE:
            switch (VariableType) {
                // Conservative Variable Formulation
                case VARIABLE_CONSERVATIVE:
                    // Row:0
                    Jacobian_Conv[0][0] = 0.0;
                    Jacobian_Conv[0][1] = nx;
                    Jacobian_Conv[0][2] = ny;
                    Jacobian_Conv[0][3] = nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = (-2.0*Q[1]*(Q[2]*ny + Q[3]*nz) + nx*((Gamma - 3.0)*Q[1]*Q[1] + (Gamma - 1.0)*(Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[1][1] = (-(Gamma - 3.0)*Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[0];
                    Jacobian_Conv[1][2] = (Q[1]*ny - (Gamma - 1.0)*Q[2]*nx)/Q[0];
                    Jacobian_Conv[1][3] = (Q[1]*nz - (Gamma - 1.0)*Q[3]*nx)/Q[0];
                    Jacobian_Conv[1][4] = (Gamma - 1.0)*nx;
                    // Row:2
                    Jacobian_Conv[2][0] = (-2.0*Q[2]*(Q[1]*nx + Q[3]*nz) + ny*((Gamma - 3.0)*Q[2]*Q[2] + (Gamma - 1.0)*(Q[1]*Q[1] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[2][1] = (Q[2]*nx - (Gamma - 1.0)*Q[1]*ny)/Q[0];
                    Jacobian_Conv[2][2] = (Q[1]*nx - (Gamma - 3.0)*Q[2]*ny + Q[3]*nz)/Q[0];
                    Jacobian_Conv[2][3] = (Q[2]*nz - (Gamma - 1.0)*Q[3]*ny)/Q[0];
                    Jacobian_Conv[2][4] = (Gamma - 1.0)*ny;
                    // Row:3
                    Jacobian_Conv[3][0] = (-2.0*Q[3]*(Q[1]*nx + Q[2]*ny) + nz*((Gamma - 3.0)*Q[3]*Q[3] + (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[3][1] = (Q[3]*nx - (Gamma - 1.0)*Q[1]*nz)/Q[0];
                    Jacobian_Conv[3][2] = (Q[3]*ny - (Gamma - 1.0)*Q[2]*nz)/Q[0];
                    Jacobian_Conv[3][3] = (Q[1]*nx + Q[2]*ny - (Gamma - 3.0)*Q[3]*nz)/Q[0];
                    Jacobian_Conv[3][4] = (Gamma - 1.0)*nz;
                    // Row:4
                    Jacobian_Conv[4][0] = (Q[1]*nx + Q[2]*ny + Q[3]*nz)*((Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) - Gamma*Q[0]*Q[4])/(Q[0]*Q[0]*Q[0]);
                    Jacobian_Conv[4][1] = (-2.0*(Gamma - 1.0)*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nx*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][2] = (-2.0*(Gamma - 1.0)*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + ny*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][3] = (-2.0*(Gamma - 1.0)*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nz*(2.0*Gamma*Q[0]*Q[4] - (Gamma - 1.0)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])))/(2.0*Q[0]*Q[0]);
                    Jacobian_Conv[4][4] = Gamma*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[0];
                    break;
                // Primitive Variable Formulation Pressure Velocity Temperature
                case VARIABLE_PRIMITIVE_PUT:
                    // Row:0
                    Jacobian_Conv[0][0] = (Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[0][1] = Q[0]*nx/Q[4];
                    Jacobian_Conv[0][2] = Q[0]*ny/Q[4];
                    Jacobian_Conv[0][3] = Q[0]*nz/Q[4];
                    Jacobian_Conv[0][4] = -Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:1
                    Jacobian_Conv[1][0] = nx + Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[1][1] = Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[1][2] = Q[0]*Q[1]*ny/Q[4];
                    Jacobian_Conv[1][3] = Q[0]*Q[1]*nz/Q[4];
                    Jacobian_Conv[1][4] = -Q[0]*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:2
                    Jacobian_Conv[2][0] = ny + Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[2][1] = Q[0]*Q[2]*nx/Q[4];
                    Jacobian_Conv[2][2] = Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[2][3] = Q[0]*Q[2]*nz/Q[4];
                    Jacobian_Conv[2][4] = -Q[0]*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:3
                    Jacobian_Conv[3][0] = nz + Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/Q[4];
                    Jacobian_Conv[3][1] = Q[0]*Q[3]*nx/Q[4];
                    Jacobian_Conv[3][2] = Q[0]*Q[3]*ny/Q[4];
                    Jacobian_Conv[3][3] = Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz)/Q[4];
                    Jacobian_Conv[3][4] = -Q[0]*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Q[4]*Q[4]);
                    // Row:4
                    Jacobian_Conv[4][0] = (Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0))/(2.0*Q[4]);
                    Jacobian_Conv[4][1] = (Q[0]*(Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nx*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0))))/Q[4];
                    Jacobian_Conv[4][2] = (Q[0]*(Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*ny*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0))))/Q[4];
                    Jacobian_Conv[4][3] = (Q[0]*(Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nz*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0))))/Q[4];
                    Jacobian_Conv[4][4] = -Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3])/(2.0*Q[4]*Q[4]);
                    break;
                // Primitive Variable Formulation Density Velocity Pressure
                case VARIABLE_PRIMITIVE_RUP:
                    // Row:0
                    Jacobian_Conv[0][0] = Q[1]*nx + Q[2]*ny + Q[3]*nz;
                    Jacobian_Conv[0][1] = Q[0]*nx;
                    Jacobian_Conv[0][2] = Q[0]*ny;
                    Jacobian_Conv[0][3] = Q[0]*nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][1] = Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
                    Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
                    Jacobian_Conv[1][4] = nx;
                    // Row:2
                    Jacobian_Conv[2][0] = Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
                    Jacobian_Conv[2][2] = Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
                    Jacobian_Conv[2][4] = ny;
                    // Row:3
                    Jacobian_Conv[3][0] = Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz);;
                    Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
                    Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
                    Jacobian_Conv[3][3] = Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz);;
                    Jacobian_Conv[3][4] = nz;
                    // Row:4
                    Jacobian_Conv[4][0] = 0.5*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]);
                    Jacobian_Conv[4][1] = (2.0*(Gamma - 1.0)*Q[0]*Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nx*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]))/(2.0*(Gamma - 1.0));
                    Jacobian_Conv[4][2] = (2.0*(Gamma - 1.0)*Q[0]*Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + ny*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]))/(2.0*(Gamma - 1.0));
                    Jacobian_Conv[4][3] = (2.0*(Gamma - 1.0)*Q[0]*Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + nz*((Gamma - 1.0)*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) + 2.0*Gamma*Q[4]))/(2.0*(Gamma - 1.0));
                    Jacobian_Conv[4][4] = Gamma*(Q[1]*nx + Q[2]*ny + Q[3]*nz)/(Gamma - 1.0);
                    break;
                // Primitive Variable Formulation Density Velocity Temperature
                case VARIABLE_PRIMITIVE_RUT:
                    // Row:0
                    Jacobian_Conv[0][0] = Q[1]*nx + Q[2]*ny + Q[3]*nz;
                    Jacobian_Conv[0][1] = Q[0]*nx;
                    Jacobian_Conv[0][2] = Q[0]*ny;
                    Jacobian_Conv[0][3] = Q[0]*nz;
                    Jacobian_Conv[0][4] = 0.0;
                    // Row:1
                    Jacobian_Conv[1][0] = Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*nx;
                    Jacobian_Conv[1][1] = Q[0]*(2.0*Q[1]*nx + Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[1][2] = Q[0]*Q[1]*ny;
                    Jacobian_Conv[1][3] = Q[0]*Q[1]*nz;
                    Jacobian_Conv[1][4] = Q[0]*nx;
                    // Row:2
                    Jacobian_Conv[2][0] = Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*ny;
                    Jacobian_Conv[2][1] = Q[0]*Q[2]*nx;
                    Jacobian_Conv[2][2] = Q[0]*(Q[1]*nx + 2.0*Q[2]*ny + Q[3]*nz);
                    Jacobian_Conv[2][3] = Q[0]*Q[2]*nz;
                    Jacobian_Conv[2][4] = Q[0]*ny;
                    // Row:3
                    Jacobian_Conv[3][0] = Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + Q[4]*nz;
                    Jacobian_Conv[3][1] = Q[0]*Q[3]*nx;
                    Jacobian_Conv[3][2] = Q[0]*Q[3]*ny;
                    Jacobian_Conv[3][3] = Q[0]*(Q[1]*nx + Q[2]*ny + 2.0*Q[3]*nz);
                    Jacobian_Conv[3][4] = Q[0]*nz;
                    // Row:4
                    Jacobian_Conv[4][0] = 0.5*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0));
                    Jacobian_Conv[4][1] = Q[0]*(Q[1]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nx*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0)));
                    Jacobian_Conv[4][2] = Q[0]*(Q[2]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*ny*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0)));
                    Jacobian_Conv[4][3] = Q[0]*(Q[3]*(Q[1]*nx + Q[2]*ny + Q[3]*nz) + 0.5*nz*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + 2.0*Q[4]*Gamma/(Gamma - 1.0)));
                    Jacobian_Conv[4][4] = Q[0]*(Q[1]*nx + Q[2]*ny + Q[3]*nz)*Gamma/(Gamma - 1.0);
                    break;
                default:
                    error("Compute_Flux_Jacobian_Euler_Convective: Undefined Variable Type - %d - Error-3", VariableType);
                    break;
            }
            break;
        default:
            error("Compute_Flux_Jacobian_Euler_Convective: Undefined Non Dimensional Method - %d - Error-4", NonDimensionalMethod);
            break;
    }
    
    // Based on variable type update the perturbation variables
    switch (VariableType) {
        case VARIABLE_CONSERVATIVE:
            // Do Nothing
            break;
        case VARIABLE_PRIMITIVE_PUT:
            Q[0] -= Gauge_Pressure;
            break;
        case VARIABLE_PRIMITIVE_RUP:
            Q[4] -= Gauge_Pressure;
            break;
        case VARIABLE_PRIMITIVE_RUT:
            break;
    }
}

//------------------------------------------------------------------------------
//! Compute Euler Convective Jacobian
//! Note: Q is always first order
//------------------------------------------------------------------------------
void Compute_Flux_Jacobian_Euler_Convective(int node_ID, Vector3D areavec, double **Jacobian_Conv) {
    double Q[NEQUATIONS];
    
    // Get the Variables
    Q[0] = Q1[node_ID];
    Q[1] = Q2[node_ID];
    Q[2] = Q3[node_ID];
    Q[3] = Q4[node_ID];
    Q[4] = Q5[node_ID];
    
    // Compute the Flux Jacobian of Convective Term
    Compute_Flux_Jacobian_Euler_Convective(Q, areavec, Jacobian_Conv);
}
