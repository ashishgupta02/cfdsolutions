/*******************************************************************************
 * File:        CompressibleUtils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _COMPRESSIBLEUTILS_H
#define	_COMPRESSIBLEUTILS_H

//------------------------------------------------------------------------------
//! Note: Only Primitive Variable are Gauge Pressure and Temperature Adjusted 
//!       not Conservative Variables
//------------------------------------------------------------------------------
double ConservativeGetPressure(double *Qc, double gamma, double gauge_pressure);
double ConservativeGetTemperature(double *Qc, double gamma, double gauge_temperature);

void ConservativeToRhoVelocityPressure(double *Qc, double *Qp, double gamma, double gauge_pressure);
void ConservativeToPressureVelocityTemperature(double *Qc, double *Qp, double gamma, double gauge_pressure, double gauge_temperature);
void ConservativeToRhoVelocityTemperature(double *Qc, double *Qp, double gamma, double gauge_temperature);

void RhoVelocityPressureToConservative(double *Qp, double *Qc, double gamma, double gauge_pressure);
void PressureVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma, double gauge_pressure, double guage_temperature);
void RhoVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma, double gauge_temperature);

void ConservativeEulerFlux(double *Qc, Vector3D areavec, double *Flux);
void ConservativeEulerFluxJacobian(double *Qc, Vector3D areavec, double **Jacobian, double gamma);

#endif	/* _COMPRESSIBLEUTILS_H */

