/*******************************************************************************
 * File:        CompressibleUtils.cpp
 * Author:      Ashish Gupta
 * Revision:    3
 ******************************************************************************/

// Custom header files
#include "CompressibleUtils.h"

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//------------------------------------------------------------------------------
double ConservativeGetPressure(double *Qc, double gamma) {
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    return (gamma - 1.0)*(Qc[4] - (0.5/Qc[0])*(Qc[1]*Qc[1] + Qc[2]*Qc[2] + Qc[3]*Qc[3]));
}

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//------------------------------------------------------------------------------
double ConservativeGetTemperature(double *Qc, double gamma) {
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    // Temperature = gamma*Pressure/Rho
    return gamma*((gamma - 1.0)*(Qc[4] - (0.5/Qc[0])*(Qc[1]*Qc[1] + Qc[2]*Qc[2] + Qc[3]*Qc[3])))/Qc[0];
}

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Premitive Variable [Rho, U, V, W, P]
//------------------------------------------------------------------------------
void ConservativeToRhoVelocityPressure(double *Qc, double *Qp, double gamma) {
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
}

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Premitive Variable [P, U, V, W, T]
//------------------------------------------------------------------------------
void ConservativeToPressureVelocityTemperature(double *Qc, double *Qp, double gamma) {
    // U - Velocity
    Qp[1] = Qc[1]/Qc[0];
    // V - Velocity
    Qp[2] = Qc[2]/Qc[0];
    // W - Velocity
    Qp[3] = Qc[3]/Qc[0];
    // Pressure = (gamma - 1)*(Rho*Et - 0.5*Rho*(U*U + V*V + W*W))
    Qp[0] = (gamma - 1.0)*(Qc[4] - 0.5*Qc[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]));
    // Temperature = gamma*Pressure/Rho
    Qp[4] = gamma*Qp[0]/Qc[0];
}

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Premitive Variable [Rho, U, V, W, T]
//------------------------------------------------------------------------------
void ConservativeToRhoVelocityTemperature(double *Qc, double *Qp, double gamma) {
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
    Qp[4] = gamma*Qp[4]/Qc[0];
}

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Premitive Variable [Rho, U, V, W, P]
//------------------------------------------------------------------------------
void RhoVelocityPressureToConservative(double *Qp, double *Qc, double gamma) {
    // Rho
    Qc[0] = Qp[0];
    // Rho*U
    Qc[1] = Qp[0]*Qp[1];
    // Rho*V
    Qc[2] = Qp[0]*Qp[2];
    // Rho*W
    Qc[3] = Qp[0]*Qp[3];
    // Rho*Et = Pressure/(gamma - 1) + 0.5*Rho*(U*U + V*V + W*W)
    Qc[4] = Qp[4]/(gamma - 1.0) + 0.5*Qp[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]);
}

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Premitive Variable [P, U, V, W, T]
//------------------------------------------------------------------------------
void PressureVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma) {
    // Rho = gamma*Pressure/Temperature
    Qc[0] = gamma*Qp[0]/Qp[4];
    // Rho*U
    Qc[1] = Qc[0]*Qp[1];
    // Rho*V
    Qc[2] = Qc[0]*Qp[2];
    // Rho*W
    Qc[3] = Qc[0]*Qp[3];
    // Rho*Et = Pressure/(gamma - 1) + 0.5*Rho*(U*U + V*V + W*W)
    Qc[4] = Qp[0]/(gamma - 1.0) + 0.5*Qc[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]);
}

//------------------------------------------------------------------------------
//! Qc : Conservatiove Variable [Rho, RhoU, RhoV, RhoEt]
//! Qp : Premitive Variable [Rho, U, V, W, T]
//------------------------------------------------------------------------------
void RhoVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma){
    // Rho
    Qc[0] = Qp[0];
    // Rho*U
    Qc[1] = Qp[0]*Qp[1];
    // Rho*V
    Qc[2] = Qp[0]*Qp[2];
    // Rho*W
    Qc[3] = Qp[0]*Qp[3];
    // Pressure = Rho*Temperature/gamma
    Qc[4] = Qp[0]*Qp[4]/gamma;
    // Rho*Et = Pressure/(gamma - 1) + 0.5*Rho*(U*U + V*V + W*W)
    Qc[4] = Qc[4]/(gamma - 1.0) + 0.5*Qp[0]*(Qp[1]*Qp[1] + Qp[2]*Qp[2] + Qp[3]*Qp[3]);
}

