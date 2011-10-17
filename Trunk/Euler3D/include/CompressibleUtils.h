/*******************************************************************************
 * File:        CompressibleUtils.h
 * Author:      Ashish Gupta
 * Revision:    4
 ******************************************************************************/

#include "License.h"

#ifndef _COMPRESSIBLEUTILS_H
#define	_COMPRESSIBLEUTILS_H

double ConservativeGetPressure(double *Qc, double gamma, double gauge_pressure);
double ConservativeGetTemperature(double *Qc, double gamma, double gauge_pressure);

void ConservativeToRhoVelocityPressure(double *Qc, double *Qp, double gamma, double gauge_pressure);
void ConservativeToPressureVelocityTemperature(double *Qc, double *Qp, double gamma, double gauge_pressure);
void ConservativeToRhoVelocityTemperature(double *Qc, double *Qp, double gamma, double gauge_pressure);

void RhoVelocityPressureToConservative(double *Qp, double *Qc, double gamma, double gauge_pressure);
void PressureVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma, double gauge_pressure);
void RhoVelocityTemperatureToConservative(double *Qp, double *Qc, double gamma, double gauge_pressure);

void ConservativeEulerFlux(double *Qc, Vector3D areavec, double *Flux, double gauge_pressure);
void ConservativeEulerFluxJacobian(double *Qc, Vector3D areavec, double **Jacobian, double gamma, double gauge_pressure);

#endif	/* _COMPRESSIBLEUTILS_H */

