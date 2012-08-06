/*******************************************************************************
 * File:        CompressibleUtils.cpp
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

// Custom header files
#include "Vector3D.h"
#include "CompressibleUtils.h"

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
double ConservativeGetPressure(double *Qc, double gamma, double gauge_pressure) {
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    return (gamma - 1.0)*(Qc[4] - (0.5/Qc[0])*(Qc[1]*Qc[1] + Qc[2]*Qc[2] + Qc[3]*Qc[3])) - gauge_pressure;
}

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
double ConservativeGetTemperature(double *Qc, double gamma, double gauge_temperature) {
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    // Temperature = gamma*Pressure/Rho
    return gamma*(gamma - 1.0)*(Qc[4] - (0.5/Qc[0])*(Qc[1]*Qc[1] + Qc[2]*Qc[2] + Qc[3]*Qc[3]))/Qc[0] - gauge_temperature;
}

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Primitive Variable [Rho, U, V, W, P]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void ConservativeToRhoVelocityPressure(double *Qc, double *Qp, double gamma, double gauge_pressure) {
    // Rho
    Qp[0] = Qc[0];
    // U - Velocity
    Qp[1] = Qc[1]/Qc[0];
    // V - Velocity
    Qp[2] = Qc[2]/Qc[0];
    // W - Velocity
    Qp[3] = Qc[3]/Qc[0];
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    Qp[4] = (gamma - 1.0)*(Qc[4] - 0.5*Qp[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3])) - gauge_pressure;
}

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Primitive Variable [P, U, V, W, T]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void ConservativeToPressureVelocityTemperature(double *Qc, double *Qp, double gamma, double gauge_pressure, double gauge_temperature) {
    // U - Velocity
    Qp[1] = Qc[1]/Qc[0];
    // V - Velocity
    Qp[2] = Qc[2]/Qc[0];
    // W - Velocity
    Qp[3] = Qc[3]/Qc[0];
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    Qp[0] = (gamma - 1.0)*(Qc[4] - 0.5*Qc[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]));
    // Temperature = gamma*Pressure/Rho
    Qp[4] = gamma*Qp[0]/Qc[0] - gauge_temperature;
    Qp[0] = Qp[0] - gauge_pressure;
}

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Primitive Variable [Rho, U, V, W, T]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void ConservativeToRhoVelocityTemperature(double *Qc, double *Qp, double gamma, double gauge_temperature) {
    // Rho
    Qp[0] = Qc[0];
    // U - Velocity
    Qp[1] = Qc[1]/Qc[0];
    // V - Velocity
    Qp[2] = Qc[2]/Qc[0];
    // W - Velocity
    Qp[3] = Qc[3]/Qc[0];
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    Qp[4] = (gamma - 1.0)*(Qc[4] - 0.5*Qp[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]));
    // Temperature = gamma*Pressure/Rho
    Qp[4] = gamma*Qp[4]/Qc[0] - gauge_temperature;
}

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Primitive Variable [Rho, U, V, W, P]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void RhoVelocityPressureToConservative(double *Qp, double *Qc, double gamma, double gauge_pressure) {
    // Rho
    Qc[0] = Qp[0];
    // Rho*U
    Qc[1] = Qp[0]*Qp[1];
    // Rho*V
    Qc[2] = Qp[0]*Qp[2];
    // Rho*W
    Qc[3] = Qp[0]*Qp[3];
    // Rho*Et = Pressure/(gamma - 1) + 0.5*Rho*(U*U + V*V + W*W)
    Qc[4] = (Qp[4] + gauge_pressure)/(gamma - 1.0) + 0.5*Qp[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]);
}

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Primitive Variable [P, U, V, W, T]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void PressureVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma, double gauge_pressure, double gauge_temperature) {
    // Rho = gamma*Pressure/Temperature
    Qc[0] = gamma*(Qp[0] + gauge_pressure)/(Qp[4] + gauge_temperature);
    // Rho*U
    Qc[1] = Qc[0]*Qp[1];
    // Rho*V
    Qc[2] = Qc[0]*Qp[2];
    // Rho*W
    Qc[3] = Qc[0]*Qp[3];
    // Rho*Et = Pressure/(gamma - 1) + 0.5*Rho*(U*U + V*V + W*W)
    Qc[4] = (Qp[0] + gauge_pressure)/(gamma - 1.0) + 0.5*Qc[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]);
}

//------------------------------------------------------------------------------
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Primitive Variable [Rho, U, V, W, T]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void RhoVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma, double gauge_temperature){
    // Rho
    Qc[0] = Qp[0];
    // Rho*U
    Qc[1] = Qp[0]*Qp[1];
    // Rho*V
    Qc[2] = Qp[0]*Qp[2];
    // Rho*W
    Qc[3] = Qp[0]*Qp[3];
    // Pressure = Rho*Temperature/gamma
    Qc[4] = Qp[0]*(Qp[4] + gauge_temperature)/gamma;
    // Rho*Et = Pressure/(gamma - 1) + 0.5*Rho*(U*U + V*V + W*W)
    Qc[4] = Qc[4]/(gamma - 1.0) + 0.5*Qp[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]);
}

//------------------------------------------------------------------------------
//! Compute Euler Flux using conservative variable
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
void ConservativeEulerFlux(double *Qc, Vector3D areavec, double *Flux, double gamma) {
    double nx, ny, nz;
    double rho, u, v, w, rhoet, p, ht, ubar;
    
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    rho   = Qc[0];
    u     = Qc[1] / rho;
    v     = Qc[2] / rho;
    w     = Qc[3] / rho;
    rhoet = Qc[4];
    p     = (gamma - 1.0)*(rhoet - 0.5*rho*(u*u + v*v + w*w));
    ht    = (rhoet + p)/rho;
    ubar  = u*nx + v*ny + w*nz;
    
    // Compute Flux
    Flux[0] = rho*ubar;
    Flux[1] =   u*Flux[0] + p*nx;
    Flux[2] =   v*Flux[0] + p*ny;
    Flux[3] =   w*Flux[0] + p*nz;
    Flux[4] =  ht*Flux[0];
}

//------------------------------------------------------------------------------
//! Compute Euler Flux Jacobian using conservative variable
//! Qc : Conservative Variable [Rho, RhoU, RhoV, RhoEt]
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//! AIAA 2001-2609
//------------------------------------------------------------------------------
void ConservativeEulerFluxJacobian(double *Qc, Vector3D areavec, double **Jacobian, double gamma) {
    double nx, ny, nz;
    double rho, u, v, w, rhoet, p, ht, ubar;
    double ek;
    
    areavec.normalize();
    nx = areavec.vec[0];
    ny = areavec.vec[1];
    nz = areavec.vec[2];
    
    rho   = Qc[0];
    u     = Qc[1] / rho;
    v     = Qc[2] / rho;
    w     = Qc[3] / rho;
    rhoet = Qc[4];
    ek    = 0.5*(u*u + v*v + w*w);
    p     = (gamma - 1.0)*(rhoet - rho*ek);
    ht    = (rhoet + p)/rho;
    ubar  = u*nx + v*ny + w*nz;
    
    // Compute Jacobian
    Jacobian[0][0] = 0.0;
    Jacobian[1][0] = (gamma - 1.0)*ek*nx - u*ubar;
    Jacobian[2][0] = (gamma - 1.0)*ek*ny - v*ubar;
    Jacobian[3][0] = (gamma - 1.0)*ek*nz - w*ubar;
    Jacobian[4][0] = ((gamma - 1.0)*ek - ht)*ubar;
    
    Jacobian[0][1] = nx;
    Jacobian[1][1] = ubar  - (gamma - 2.0)*u*nx;
    Jacobian[2][1] = v*nx  - (gamma - 1.0)*u*ny;
    Jacobian[3][1] = w*nx  - (gamma - 1.0)*u*nz;
    Jacobian[4][1] = ht*nx - (gamma - 1.0)*u*ubar;
    
    Jacobian[0][2] = ny;
    Jacobian[1][2] = u*ny  - (gamma - 1.0)*v*nx;
    Jacobian[2][2] = ubar  - (gamma - 2.0)*v*ny;
    Jacobian[3][2] = w*ny  - (gamma - 1.0)*v*nz;
    Jacobian[4][2] = ht*ny - (gamma - 1.0)*v*ubar;
    
    Jacobian[0][3] = nz;
    Jacobian[1][3] = u*nz  - (gamma - 1.0)*w*nx;
    Jacobian[2][3] = v*nz  - (gamma - 1.0)*w*ny;
    Jacobian[3][3] = ubar  - (gamma - 2.0)*w*nz;
    Jacobian[4][3] = ht*nz - (gamma - 1.0)*w*ubar;
    
    Jacobian[0][4] = 0.0;
    Jacobian[1][4] = (gamma - 1.0)*nx;
    Jacobian[2][4] = (gamma - 1.0)*ny;
    Jacobian[3][4] = (gamma - 1.0)*nz;
    Jacobian[4][4] = gamma*ubar;
}

